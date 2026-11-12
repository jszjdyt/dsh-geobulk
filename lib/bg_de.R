#!/usr/bin/env Rscript
# dsh-geobulk — 差异表达 R 脚本（limma / DESeq2 官方实现）
#
# 调用：Rscript --vanilla lib/bg_de.R --in <bundle> --out <out> --params '<json>'
# bundle: expr.tsv（样本×基因）、info.csv（样本注释）、params.json
# 输出：<out>/result.json、<out>/tables/de.csv
# 退出码：0 成功；3 缺前置（stderr 给中文原因）；其它失败
#
# 约定见 CONTRACT.md §3.2。图由 Python 侧统一绘制，本脚本只产出差异结果表。

suppressPackageStartupMessages({
  ok_json <- requireNamespace("jsonlite", quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(key) {
  idx <- match(key, args)
  if (is.na(idx) || idx >= length(args)) return(NULL)
  args[idx + 1]
}

fail_prereq <- function(message) {
  cat(message, "\n", file = stderr())
  quit(save = "no", status = 3)
}
fail_error <- function(message) {
  cat(message, "\n", file = stderr())
  quit(save = "no", status = 1)
}

in_dir <- getarg("--in")
out_dir <- getarg("--out")
params_raw <- getarg("--params")
if (is.null(in_dir) || is.null(out_dir)) fail_prereq("R 脚本缺少 --in / --out 参数")
if (!dir.exists(in_dir)) fail_prereq(paste0("bundle 目录不存在：", in_dir))
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

params <- list()
if (!is.null(params_raw)) {
  params <- tryCatch({
    if (requireNamespace("jsonlite", quietly = TRUE)) jsonlite::fromJSON(params_raw, simplifyVector = TRUE)
    else list()
  }, error = function(e) list())
}

read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

expr <- read_tsv(file.path(in_dir, "expr.tsv"))
if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) fail_prereq("bundle 里没有可用的 expr.tsv")
info <- read_tsv(file.path(in_dir, "info.csv"))
if (is.null(info)) {
  info <- data.frame(row.names = rownames(expr))
}

method <- if (!is.null(params$method)) tolower(as.character(params$method)) else "limma"
group_col <- if (!is.null(params$group_col)) as.character(params$group_col) else ""
case_label <- if (!is.null(params$case)) as.character(params$case) else ""
control_label <- if (!is.null(params$control)) as.character(params$control) else ""
covars <- if (!is.null(params$covars)) as.character(unlist(params$covars)) else character(0)
paired <- isTRUE(params$paired)
pair_col <- if (!is.null(params$pair_col)) as.character(params$pair_col) else ""

if (group_col == "" || !(group_col %in% colnames(info))) {
  fail_prereq(paste0("临床/样本表里没有分组列：", group_col))
}
groups <- as.character(info[[group_col]])
sample_order <- rownames(expr)
if (!all(sample_order %in% rownames(info))) {
  missing <- setdiff(sample_order, rownames(info))
  fail_prereq(paste0("样本注释缺少这些样本的分组：", paste(head(missing, 5), collapse = ", ")))
}
groups <- groups[match(sample_order, rownames(info))]

if (case_label == "" || !(case_label %in% groups)) {
  fail_prereq(paste0("实验组不存在或为空：", case_label))
}
if (control_label == "" || !(control_label %in% groups)) {
  fail_prereq(paste0("对照组不存在或为空：", control_label))
}
if (sum(groups == case_label) < 2 || sum(groups == control_label) < 2) {
  fail_prereq("实验组或对照组样本数不足 2 个")
}

keep_samples <- groups %in% c(case_label, control_label, if (paired) unique(groups) else character(0))
if (all(groups %in% c(case_label, control_label))) {
  keep_samples <- rep(TRUE, length(groups))
}
expr <- expr[keep_samples, , drop = FALSE]
info <- info[match(rownames(expr), rownames(info)), , drop = FALSE]
groups <- groups[keep_samples]
if (nrow(expr) < 4) fail_prereq("可用样本不足 4 个")

result <- list()

if (method == "deseq2") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) fail_prereq("本机没有安装 DESeq2")
  suppressPackageStartupMessages(library(DESeq2))
  mat <- t(round(as.matrix(expr)))  # DESeq2 要 genes × samples 的整数矩阵
  storage.mode(mat) <- "integer"
  col_data <- data.frame(
    group = factor(groups, levels = c(control_label, case_label)),
    row.names = rownames(expr)
  )
  for (covariate in covars) {
    if (covariate %in% colnames(info)) {
      col_data[[covariate]] <- factor(as.character(info[[covariate]]))
    }
  }
  formula_text <- paste(c("group", covars[covars %in% colnames(col_data)]), collapse = " + ")
  design <- as.formula(paste("~", formula_text))
  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = col_data, design = design),
    error = function(e) fail_prereq(paste0("DESeq2 构建对象失败：", conditionMessage(e)))
  )
  dds <- tryCatch(
    DESeq2::DESeq(dds, quiet = TRUE),
    error = function(e) fail_prereq(paste0("DESeq2 拟合失败：", conditionMessage(e)))
  )
  res <- DESeq2::results(dds, contrast = c("group", case_label, control_label))
  frame <- as.data.frame(res)
  frame$gene <- rownames(frame)
  frame <- frame[, c("gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
  frame$B <- 0
  method_used <- "deseq2"
} else {
  if (!requireNamespace("limma", quietly = TRUE)) fail_prereq("本机没有安装 limma")
  suppressPackageStartupMessages(library(limma))
  mat <- t(as.matrix(expr))  # genes × samples
  group_factor <- factor(groups)
  group_factor <- relevel(group_factor, ref = control_label)
  design_terms <- "group_factor"
  covar_frames <- list()
  if (length(covars) > 0) {
    for (covariate in covars) {
      if (covariate %in% colnames(info)) {
        values <- as.character(info[[covariate]])
        if (length(unique(values)) < 2) next
        col_name <- paste0("cov_", make.names(covariate))
        covar_frames[[col_name]] <- factor(values)
      }
    }
    if (length(covar_frames) > 0) {
      design_terms <- paste(c(design_terms, names(covar_frames)), collapse = " + ")
    }
  }
  if (paired) {
    if (pair_col == "" || !(pair_col %in% colnames(info))) {
      fail_prereq("配对分析需要提供 pair_col 且该列存在于样本表")
    }
    block <- factor(as.character(info[[pair_col]]))
    design <- model.matrix(as.formula(paste("~", paste(c("0", "group_factor", names(covar_frames)), collapse = " + "))))
    colnames(design) <- make.names(colnames(design))
    fit <- tryCatch(
      limma::lmFit(mat, design),
      error = function(e) fail_prereq(paste0("limma::lmFit 失败：", conditionMessage(e)))
    )
    contrast_col <- grep("^group_factor", colnames(design), value = TRUE)
    case_col <- contrast_col[grepl(make.names(case_label), contrast_col)]
    ctrl_col <- contrast_col[grepl(make.names(control_label), contrast_col)]
    if (length(case_col) == 0 || length(ctrl_col) == 0) fail_prereq("无法构造实验组/对照组对比")
    contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
    rownames(contrast_matrix) <- colnames(design)
    contrast_matrix[case_col[1], 1] <- 1
    contrast_matrix[ctrl_col[1], 1] <- -1
    colnames(contrast_matrix) <- paste0(case_label, "-", control_label)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    table <- limma::topTable(fit2, number = Inf, sort.by = "P")
    frame <- data.frame(
      gene = rownames(table),
      logFC = table$logFC,
      AveExpr = table$AveExpr,
      t = table$t,
      P.Value = table$P.Value,
      adj.P.Val = table$adj.P.Val,
      B = table$B,
      stringsAsFactors = FALSE
    )
  } else {
    data_list <- c(list(group_factor = group_factor), covar_frames)
    design <- model.matrix(as.formula(paste("~", paste(c("0", names(data_list)), collapse = " + "))))
    colnames(design) <- make.names(colnames(design))
    fit <- tryCatch(
      limma::lmFit(mat, design),
      error = function(e) fail_prereq(paste0("limma::lmFit 失败：", conditionMessage(e)))
    )
    case_col <- grep(paste0("^group_factor", make.names(case_label)), colnames(design), value = TRUE)
    ctrl_col <- grep(paste0("^group_factor", make.names(control_label)), colnames(design), value = TRUE)
    if (length(case_col) == 0 || length(ctrl_col) == 0) fail_prereq("无法构造实验组/对照组对比")
    contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
    rownames(contrast_matrix) <- colnames(design)
    contrast_matrix[case_col[1], 1] <- 1
    contrast_matrix[ctrl_col[1], 1] <- -1
    colnames(contrast_matrix) <- "case_vs_control"
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    table <- limma::topTable(fit2, number = Inf, sort.by = "P")
    frame <- data.frame(
      gene = rownames(table),
      logFC = table$logFC,
      AveExpr = table$AveExpr,
      t = table$t,
      P.Value = table$P.Value,
      adj.P.Val = table$adj.P.Val,
      B = table$B,
      stringsAsFactors = FALSE
    )
  }
  method_used <- "limma"
}

if (is.null(frame) || nrow(frame) == 0) fail_error("R 差异分析没有产生结果行")
frame <- frame[!is.na(frame$P.Value), , drop = FALSE]
frame$regulation <- ifelse(
  frame$adj.P.Val <= 0.05 & frame$logFC > 0, "up",
  ifelse(frame$adj.P.Val <= 0.05 & frame$logFC < 0, "down", "ns")
)

write.csv(frame, file.path(out_dir, "tables", "de.csv"), row.names = FALSE, quote = FALSE)

n_up <- sum(frame$regulation == "up")
n_down <- sum(frame$regulation == "down")
result <- list(
  summary = list(
    title = if (method_used == "deseq2") "差异表达（DESeq2）" else "差异表达（limma）",
    bullets = c(
      sprintf("%s：%s vs %s（%d vs %d 样本）", method_used, case_label, control_label,
              sum(groups == case_label), sum(groups == control_label)),
      sprintf("差异基因：上调 %d、下调 %d（FDR<0.05）", n_up, n_down)
    ),
    metrics = list(
      method = method_used,
      n_genes = nrow(frame),
      n_up = n_up,
      n_down = n_down,
      n_case = sum(groups == case_label),
      n_ctrl = sum(groups == control_label)
    )
  ),
  tables = list(list(
    id = "de",
    title = "差异基因",
    columns = lapply(c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "regulation"),
                     function(name) list(key = name, label = name,
                                         type = if (name %in% c("gene", "regulation")) "string" else "number")),
    rows = lapply(seq_len(min(nrow(frame), 500)), function(i) as.list(frame[i, ])),
    total = nrow(frame)
  )),
  warnings = list(),
  notes = list(sprintf("limma/DESeq2 官方 R 实现；design = %s", if (method == "deseq2") "~group(+covars)" else "~0+group(+covars)"))
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(result, file.path(out_dir, "result.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
} else {
  # 没有 jsonlite 时手写最小 JSON
  writeLines(
    sprintf('{"summary":{"title":"%s","bullets":["R 分析完成"],"metrics":{"method":"%s","n_genes":%d}}}',
            result$summary$title, method_used, nrow(frame)),
    file.path(out_dir, "result.json")
  )
}

quit(save = "no", status = 0)
