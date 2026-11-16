#!/usr/bin/env Rscript
# dsh-geobulk — GSEA（fgsea 官方实现）
#
# 调用：Rscript --vanilla lib/bg_gsea.R --in <bundle> --out <out> --params '<json>'
# bundle/extra：ranked.tsv（gene, logFC）、sets.tsv（term, gene）
# 输出：<out>/tables/gsea.csv、<out>/result.json
# 退出码：0 成功；3 缺前置（没有 fgsea / 数据不满足）；4 运行期程序错误
#
# 修复记录（GEO-08）：
#   插件自实现的 "NES" 是 ES/std(所有置换 ES)，与官方 GSEA 定义（按符号分别用置换 ES
#   的均值归一化）不同。这里接官方 fgsea，由调用方记录 actual_method / 版本 /
#   基因集来源 / 覆盖率 / 置换策略。

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(key) {
  idx <- match(key, args)
  if (is.na(idx) || idx >= length(args)) return(NULL)
  args[idx + 1]
}
fail_prereq <- function(message) { cat(message, "\n", file = stderr()); quit(save = "no", status = 3) }
fail_runtime <- function(message) { cat(message, "\n", file = stderr()); quit(save = "no", status = 4) }

in_dir <- getarg("--in")
out_dir <- getarg("--out")
params_raw <- getarg("--params")
if (is.null(in_dir) || is.null(out_dir)) fail_prereq("R 脚本缺少 --in / --out 参数")
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

params <- list()
if (!is.null(params_raw) && requireNamespace("jsonlite", quietly = TRUE)) {
  params <- tryCatch(jsonlite::fromJSON(params_raw, simplifyVector = TRUE), error = function(e) list())
}

if (!requireNamespace("fgsea", quietly = TRUE)) fail_prereq("本机没有安装 fgsea")
suppressPackageStartupMessages(library(fgsea))

read_tsv <- function(path) {
  if (!file.exists(path)) fail_prereq(paste0("bundle 里缺少 ", basename(path)))
  tryCatch(read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
           error = function(e) fail_prereq(paste0("读取 ", basename(path), " 失败：", conditionMessage(e))))
}

ranked <- read_tsv(file.path(in_dir, "extra", "ranked.tsv"))
sets_raw <- read_tsv(file.path(in_dir, "extra", "sets.tsv"))
if (is.null(ranked$gene) || is.null(ranked$logFC)) fail_prereq("ranked.tsv 需要 gene / logFC 两列")
if (is.null(sets_raw$term) || is.null(sets_raw$gene)) fail_prereq("sets.tsv 需要 term / gene 两列")

stats <- as.numeric(ranked$logFC)
names(stats) <- as.character(ranked$gene)
stats <- stats[!is.na(stats) & !duplicated(names(stats))]
if (length(stats) < 50) fail_prereq(sprintf("排序基因只有 %d 个，太少，无法做 GSEA", length(stats)))

paths <- split(as.character(sets_raw$gene), as.character(sets_raw$term))
paths <- lapply(paths, function(x) unique(x[!is.na(x)]))
min_size <- if (!is.null(params$min_size)) as.integer(params$min_size) else 15L
max_size <- if (!is.null(params$max_size)) as.integer(params$max_size) else 500L
seed <- if (!is.null(params$seed)) as.integer(params$seed) else 0L
nperm <- if (!is.null(params$nperm)) as.integer(params$nperm) else 1000L

set.seed(seed)
res <- tryCatch(
  suppressWarnings(fgsea::fgseaMultilevel(
    pathways = paths,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    nPermSimple = max(1000L, nperm),
    nproc = 1
  )),
  error = function(e) fail_runtime(paste0("fgsea::fgseaMultilevel 失败：", conditionMessage(e)))
)
if (is.null(res) || nrow(res) == 0) fail_prereq("fgsea 没有返回任何基因集结果（检查基因名与集合大小）")

out <- data.frame(
  pathway = as.character(res$pathway),
  Description = as.character(res$pathway),
  SetSize = as.integer(res$size),
  ES = as.numeric(res$ES),
  NES = as.numeric(res$NES),
  pvalue = as.numeric(res$pval),
  p.adjust = as.numeric(res$padj),
  log2err = as.numeric(res$log2err),
  geneID = vapply(res$leadingEdge, function(x) paste(head(as.character(x), 200), collapse = "/"), character(1)),
  stringsAsFactors = FALSE
)
out <- out[order(out$p.adjust, out$pvalue), , drop = FALSE]
write.csv(out, file.path(out_dir, "tables", "gsea.csv"), row.names = FALSE, quote = FALSE)

version <- tryCatch(as.character(utils::packageVersion("fgsea")), error = function(e) "unknown")
result <- list(
  summary = list(
    title = "GSEA（fgsea 官方实现）",
    bullets = c(
      sprintf("fgsea::fgseaMultilevel：%d 个基因集通过大小筛选（%d–%d）", nrow(out), min_size, max_size),
      sprintf("NES 为官方定义（按符号分别用置换 ES 均值归一化）；显著 %d 条（padj≤0.05）", sum(out$p.adjust <= 0.05, na.rm = TRUE))
    ),
    metrics = list(
      actual_method = "fgsea::fgseaMultilevel",
      software = "fgsea",
      software_version = version,
      permutation = "fgseaMultilevel 自适应多层置换（nPermSimple 预设 + 多层抽样）",
      n_genes = length(stats),
      n_sets = nrow(out),
      min_size = min_size,
      max_size = max_size
    )
  ),
  warnings = list(),
  notes = list(sprintf("GSEA 由 R fgsea %s 执行；数值不可与自实现 ES/std 版本互换比较", version))
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(result, file.path(out_dir, "result.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
}
quit(save = "no", status = 0)
