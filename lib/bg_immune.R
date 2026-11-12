#!/usr/bin/env Rscript
# =============================================================================
# dsh-geobulk :: 免疫浸润（IOBR 优先）R 脚本
#
# 由 lib/rbridge.py 调用：
#   Rscript --vanilla lib/bg_immune.R --in <bundle> --out <r_out> --params '<json>'
# bundle（rbridge 写）：
#   expr.tsv      样本×基因（首列样本名，首行基因名）
#   info.csv      样本信息（可选）
#   params.json   params 的 JSON
#   extra/signatures.tsv  ssGSEA 用签名表（列 cell / gene）
# r_out（本脚本写）：
#   result.json   {"summary":{"title","bullets","metrics"},"warnings":[],"notes":[],"methods":[...]}
#   tables/immune_<method>.csv  每方法宽表（首列 sample，其余细胞/分数列）
#   tables/ips.csv
#
# 方法：cibersort / epic / mcpcounter / xcell / quantiseq / estimate / timer（IOBR::deconvo_tme）
#       ssgsea（GSVA::gsva ssgseaParam + extra/signatures.tsv 的签名）
#       ips（IOBR::deconvo_ips，单独表）
# 每个方法独立 tryCatch：失败跳过并记 warning，不中断其它方法。
# 退出码：0 成功（即使部分方法失败）；3 数据不满足前置（如 bundle 缺失/样本不足）。
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(key) {
  i <- match(key, args)
  if (is.na(i) || i >= length(args)) return(NULL)
  args[i + 1L]
}

in_dir <- getarg("--in")
out_dir <- getarg("--out")
params_arg <- getarg("--params")

fail_prereq <- function(msg) {
  cat(msg, "\n", file = stderr())
  quit(status = 3L, save = "no")
}

if (is.null(in_dir) || is.null(out_dir)) fail_prereq("R 脚本缺少 --in / --out 参数")
if (!dir.exists(in_dir)) fail_prereq(paste0("R bundle 目录不存在：", in_dir))
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# IOBR 在 0.99.0 里以自由变量引用包内数据集（lm22/mcp_genes/...），且 EPIC 需要
# IOBR 处于 search path 上（as.environment("package:IOBR")）。Rscript --vanilla 下
# 必须显式 library() + data() 加载，否则各方法会报 "object 'lm22' not found"。
IOBR_OK <- suppressWarnings(suppressPackageStartupMessages(
  requireNamespace("IOBR", quietly = TRUE)
))
if (IOBR_OK) {
  suppressWarnings(suppressPackageStartupMessages(try(library(IOBR), silent = TRUE)))
  for (.iobr_data in c("lm22", "mcp_genes", "mcp_probesets", "xCell.data",
                       "quantiseq_data", "common_genes", "SI_geneset",
                       "ips_gene_set", "TRef", "BRef", "mRNA_cell_default")) {
    if (!exists(.iobr_data, envir = globalenv(), inherits = FALSE)) {
      suppressWarnings(try(
        utils::data(list = .iobr_data, package = "IOBR", envir = globalenv()),
        silent = TRUE
      ))
    }
  }
}

read_params <- function() {
  p <- list()
  if (!is.null(params_arg) && nzchar(params_arg)) {
    p <- tryCatch(jsonlite::fromJSON(params_arg), error = function(e) list())
  }
  pfile <- file.path(in_dir, "params.json")
  if (file.exists(pfile)) {
    p2 <- tryCatch(jsonlite::fromJSON(pfile), error = function(e) list())
    for (nm in names(p2)) if (is.null(p[[nm]])) p[[nm]] <- p2[[nm]]
  }
  p
}

params <- read_params()
as_chr <- function(x, default = NULL) {  if (is.null(x) || length(x) == 0) return(default)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) default else x
}
as_int <- function(x, default) {
  if (is.null(x) || length(x) == 0) return(as.integer(default))
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) as.integer(default) else v
}

# 避免 IOBR 把 ESTIMATE 中间文件写到调用方 cwd（引擎工作目录）
expr_path <- file.path(in_dir, "expr.tsv")
if (!file.exists(expr_path)) fail_prereq("R bundle 缺少 expr.tsv，无法做免疫浸润")
expr_df <- read.delim(expr_path, row.names = 1, check.names = FALSE, sep = "\t")
eset <- t(as.matrix(expr_df))
storage.mode(eset) <- "numeric"
if (ncol(eset) < 3) fail_prereq(sprintf("样本数不足（%d < 3），无法做免疫浸润", ncol(eset)))
if (anyDuplicated(rownames(eset))) eset <- eset[!duplicated(rownames(eset)), , drop = FALSE]
if (nrow(eset) < 20) fail_prereq(sprintf("基因数不足（%d < 20），无法做免疫浸润", nrow(eset)))
genes <- rownames(eset)
samples <- colnames(eset)
setwd(tempdir())

ALL_METHODS <- c("cibersort", "epic", "mcpcounter", "xcell",
                 "quantiseq", "estimate", "ssgsea", "timer")
methods <- as_chr(params$methods, ALL_METHODS)
methods <- unique(methods[methods %in% ALL_METHODS])
if (length(methods) == 0) methods <- ALL_METHODS
do_ips <- is.null(params$ips) || isTRUE(as.logical(params$ips)[1])
perm <- as_int(params$perm, 100L)
perm <- max(0L, min(perm, 1000L))
indications <- as_chr(params$indications, "LUAD")[1]
platform <- as_chr(params$platform, "affymetrix")[1]

stderr_note <- function(...) cat(paste0("[bg_immune] ", ...), "\n", file = stderr())

signature_path <- file.path(in_dir, "extra", "signatures.tsv")

# ---------------------------------------------------------------------------
# ssGSEA：GSVA 2.x 的 ssgseaParam；签名由 Python 侧传入（来自 bg_genesets IMMUNE_28 或兜底）
# ---------------------------------------------------------------------------
run_ssgsea <- function() {
  if (!file.exists(signature_path)) stop("ssGSEA 缺少签名文件 extra/signatures.tsv")
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("GSVA 包不可用，无法做 ssGSEA")
  sdf <- read.delim(signature_path, row.names = 1, check.names = FALSE, sep = "\t")
  if (ncol(sdf) < 2) stop("ssGSEA 签名表格式不正确（需要 cell / gene 两列）")
  cell_col <- if ("cell" %in% colnames(sdf)) "cell" else colnames(sdf)[1]
  gene_col <- if ("gene" %in% colnames(sdf)) "gene" else colnames(sdf)[2]
  genes_by_cell <- split(as.character(sdf[[gene_col]]), as.character(sdf[[cell_col]]))
  upper_genes <- toupper(genes)
  genes_by_cell <- lapply(genes_by_cell, function(g) {
    idx <- match(toupper(unique(g)), upper_genes)
    idx <- idx[!is.na(idx)]
    unique(genes[idx])
  })
  genes_by_cell <- genes_by_cell[lengths(genes_by_cell) >= 3]
  if (length(genes_by_cell) < 3) stop("ssGSEA 签名基因与表达矩阵重叠不足（<3 个细胞类型可算）")
  par <- GSVA::ssgseaParam(eset, genes_by_cell, normalize = TRUE)
  mat <- GSVA::gsva(par, verbose = FALSE)
  data.frame(ID = colnames(mat), t(mat), check.names = FALSE)
}

# ---------------------------------------------------------------------------
# IOBR 方法分发
# ---------------------------------------------------------------------------
run_iobr <- function(method) {
  if (!isTRUE(IOBR_OK)) stop("IOBR 包不可用")
  if (method == "cibersort") {
    return(IOBR::deconvo_tme(eset, method = "cibersort", arrays = FALSE, perm = perm))
  }
  if (method == "quantiseq") {
    # IOBR 0.99.0 的 deconvo_tme 未给 scale_mrna 默认值，必须显式传入
    return(IOBR::deconvo_tme(eset, method = "quantiseq", arrays = FALSE, scale_mrna = TRUE))
  }
  if (method == "timer") {
    # deconvo_timer 要求 indications 是「每样本一个癌种」的向量（长度 = 样本数）
    return(IOBR::deconvo_tme(eset, method = "timer", arrays = FALSE,
                             group_list = rep(indications, ncol(eset))))
  }
  if (method == "estimate") {
    return(IOBR::deconvo_tme(eset, method = "estimate", arrays = FALSE, platform = platform))
  }
  IOBR::deconvo_tme(eset, method = method, arrays = FALSE)
}

# 统一整理方法输出：首列 → sample；按 samples 顺序对齐；其余列转 numeric
tidy_method_output <- function(obj) {
  df <- as.data.frame(obj, check.names = FALSE)
  if (ncol(df) < 2) stop("方法输出不足两列")
  id_col <- if ("ID" %in% names(df)) "ID" else names(df)[1]
  names(df)[names(df) == id_col] <- "sample"
  df$sample <- as.character(df$sample)
  for (cn in setdiff(names(df), "sample")) {
    df[[cn]] <- suppressWarnings(as.numeric(df[[cn]]))
  }
  keep <- c("sample", setdiff(names(df), "sample"))
  df <- df[, keep, drop = FALSE]
  if (!all(samples %in% df$sample)) stop("方法输出的样本名与输入不一致")
  df <- df[match(samples, df$sample), , drop = FALSE]
  rownames(df) <- NULL
  df
}

write_method <- function(df, method) {
  path <- file.path(out_dir, "tables", paste0("immune_", method, ".csv"))
  write.csv(df, path, row.names = FALSE)
  list(
    method = method,
    ok = TRUE,
    n_cells = as.integer(ncol(df) - 1L),
    cells = setdiff(names(df), "sample"),
    n_samples = as.integer(nrow(df)),
    file = basename(path)
  )
}

t0_all <- Sys.time()
status <- list()
warnings <- list()
notes <- list()

for (method in methods) {
  t0 <- Sys.time()
  run <- tryCatch({
    obj <- if (method == "ssgsea") run_ssgsea() else run_iobr(method)
    df <- tidy_method_output(obj)
    list(ok = TRUE, meta = write_method(df, method))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  if (isTRUE(run$ok)) {
    meta <- run$meta
    meta$elapsed_s <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    status[[method]] <- meta
    stderr_note(sprintf("OK  %-12s %2d cells  %.1fs", method, meta$n_cells, meta$elapsed_s))
  } else {
    status[[method]] <- list(method = method, ok = FALSE, error = as.character(run$error))
    warnings[[length(warnings) + 1L]] <- sprintf("%s 计算失败已跳过：%s", method, as.character(run$error))
    stderr_note(sprintf("FAIL %-12s %s", method, as.character(run$error)))
  }
}

# ---------------------------------------------------------------------------
# IPS（Immunophenoscore，IOBR::deconvo_ips）
# ---------------------------------------------------------------------------
ips_meta <- list(ok = FALSE, method = "ips", error = "not requested")
if (do_ips) {
  t0 <- Sys.time()
  ips_meta <- tryCatch({
    if (!isTRUE(IOBR_OK)) stop("IOBR 包不可用")
    obj <- IOBR::deconvo_ips(eset, plot = FALSE)
    df <- tidy_method_output(obj)
    write.csv(df, file.path(out_dir, "tables", "ips.csv"), row.names = FALSE)
    list(method = "ips", ok = TRUE, n_cells = as.integer(ncol(df) - 1L),
         cells = setdiff(names(df), "sample"), n_samples = as.integer(nrow(df)),
         elapsed_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2))
  }, error = function(e) {
    warnings[[length(warnings) + 1L]] <<- paste0("IPS 计算失败已跳过：", conditionMessage(e))
    list(method = "ips", ok = FALSE, error = conditionMessage(e))
  })
  if (isTRUE(ips_meta$ok)) stderr_note("OK  ips")
}

# ---------------------------------------------------------------------------
# result.json
# ---------------------------------------------------------------------------
ok_methods <- names(status)[vapply(status, function(x) isTRUE(x$ok), logical(1))]
failed <- names(status)[!vapply(status, function(x) isTRUE(x$ok), logical(1))]
n_ok <- length(ok_methods)
elapsed_all <- round(as.numeric(difftime(Sys.time(), t0_all, units = "secs")), 2)

bullets <- c(
  sprintf("IOBR/GSVA 成功 %d 个免疫浸润方法：%s", n_ok,
          if (n_ok > 0) paste(ok_methods, collapse = "、") else "无"),
  sprintf("样本 %d × 基因 %d；单方法耗时合计 %.1fs", length(samples), length(genes), elapsed_all)
)
if (length(failed) > 0) {
  bullets <- c(bullets, sprintf("失败/跳过：%s", paste(failed, collapse = "、")))
}
if (isTRUE(ips_meta$ok)) bullets <- c(bullets, "IPS（Immunophenoscore）已计算")

metrics <- list(
  n_methods = n_ok,
  methods = ok_methods,
  failed = failed,
  n_samples = length(samples),
  n_genes = length(genes),
  ips = isTRUE(ips_meta$ok),
  engine = "IOBR/GSVA (R)"
)

notes <- c(
  "cibersort=IOBR::deconvo_tme cibersort（LM22，相对比例）；epic/quantiseq 为比例类方法；",
  "mcpcounter/timer 为 marker 打分（非比例）；xcell/estimate 为打分类方法；ssgsea 为 GSVA ssgseaParam 签名得分；ips=IOBR::deconvo_ips。",
  "estimate 平台参数默认 affymetrix；timer 需要癌种 indications（默认 LUAD）。"
)

out <- list(
  summary = list(
    title = "免疫浸润（IOBR 优先）",
    bullets = bullets,
    metrics = metrics
  ),
  warnings = unlist(warnings, use.names = FALSE),
  notes = notes,
  methods = unname(status),
  ips = ips_meta,
  engine = "IOBR/GSVA (R)"
)
if (length(out$warnings) == 0) out$warnings <- character(0)

jsonlite::write_json(out, file.path(out_dir, "result.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("[bg_immune] done: %d methods ok, %d failed\n", n_ok, length(failed)), file = stderr())
quit(status = 0L, save = "no")
