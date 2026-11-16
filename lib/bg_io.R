#!/usr/bin/env Rscript
# dsh-geobulk — R 对象读取桥（仅被 bg_io.py 调用；rds/RData → MTX + 文本）
#
# 用法：
#   Rscript --vanilla bg_io.R --in <file.rds|.RData> --out <dir> [--object NAME] [--inspect]
#   Rscript --vanilla bg_io.R --annotate --genes <genes.txt> --out <dir> --species human|mouse \
#           --keytype SYMBOL|ENSEMBL|ENTREZID|REFSEQ
#
# --annotate（基因标识符注释，任务单第七节建议 2）：用 AnnotationDbi +
#   org.Hs.eg.db / org.Mm.eg.db 把 ID 映射到 SYMBOL，输出：
#   annotation.tsv  input_id\tsymbol（映射不到的 symbol 留空，不猜测）
#   meta.tsv        status ok|unavailable|error / database / database_version / keytype /
#                   n_input / n_mapped / one_to_many / many_to_one
#   注释包不存在时 status=unavailable（退出码 0），由 Python 侧退化成纯报告，不硬造映射。
#
# 输出（都在 --out 目录）：
#   meta.tsv        key\tvalue（status/class/n_genes/n_samples/has_counts/notes/warnings/error）
#   matrix.mtx      基因 × 样本 稀疏矩阵（Matrix::writeMM）
#   genes.tsv       基因名（每行一个）
#   samples.tsv     样本名（每行一个）
#   obs.csv         样本表（pData/colData/DGEList$samples，行名 = 样本名）
#   candidates.tsv  多对象时的候选列表（object/class/n_genes/n_samples/label）
#
# 支持对象：ExpressionSet / SummarizedExperiment / Seurat / DGEList /
#           matrix / Matrix / data.frame / 常见 list 容器（counts/data/expr/expression）
# 退出码：0 成功；1 失败（meta.tsv 里 status=error）；2 参数错误
# 不依赖 jsonlite，元信息统一用 TSV 传回，避免转义问题。

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) return(default)
  args[[i + 1]]
}
has_flag <- function(flag) flag %in% args

in_path <- get_arg("--in")
out_dir <- get_arg("--out")
object_name <- get_arg("--object")
inspect_only <- has_flag("--inspect")
annotate_mode <- has_flag("--annotate")

usage <- function() {
  cat("用法：\n",
      "  Rscript --vanilla bg_io.R --in <file> --out <dir> [--object NAME] [--inspect]\n",
      "  Rscript --vanilla bg_io.R --annotate --genes <genes.txt> --out <dir>",
      " --species human|mouse --keytype SYMBOL|ENSEMBL|ENTREZID|REFSEQ\n",
      file = stderr())
}

if (is.null(out_dir)) {
  usage()
  quit(save = "no", status = 2)
}
if (!annotate_mode) {
  if (is.null(in_path)) {
    usage()
    quit(save = "no", status = 2)
  }
  if (!file.exists(in_path)) {
    cat(paste0("文件不存在：", in_path, "\n"), file = stderr())
    quit(save = "no", status = 1)
  }
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

NOTES <- character(0)
WARNINGS <- character(0)

add_note <- function(...) {
  NOTES <<- c(NOTES, paste0(...))
}
add_warning <- function(...) {
  WARNINGS <<- c(WARNINGS, paste0(...))
}

write_meta <- function(meta) {
  keys <- names(meta)
  lines <- vapply(keys, function(k) {
    v <- meta[[k]]
    if (is.null(v) || length(v) == 0) v <- ""
    if (is.logical(v)) v <- ifelse(is.na(v), "", ifelse(v, "TRUE", "FALSE"))
    text <- gsub("[\r\n]+", " ", paste(v, collapse = " | "))
    paste0(k, "\t", text)
  }, character(1))
  writeLines(lines, file.path(out_dir, "meta.tsv"), useBytes = TRUE)
}

fail <- function(msg) {
  write_meta(list(status = "error", error = msg))
  cat(paste0(msg, "\n"), file = stderr())
  quit(save = "no", status = 1)
}

# ---------------------------------------------------------------------------
# --annotate：基因标识符 ID → SYMBOL 注释报告（本机 AnnotationDbi / org.*.eg.db）
# ---------------------------------------------------------------------------

annotation_unavailable <- function(reason, database = "", version = "", keytype = "",
                                   n_input = 0) {
  write_meta(list(status = "unavailable", database = database,
                  database_version = version, keytype = keytype,
                  n_input = n_input, n_mapped = 0, one_to_many = 0, many_to_one = 0,
                  error = reason))
  cat(paste0(reason, "\n"), file = stderr())
  quit(save = "no", status = 0)
}

run_annotation <- function() {
  genes_path <- get_arg("--genes")
  species <- tolower(get_arg("--species", "human"))
  keytype <- toupper(get_arg("--keytype", "SYMBOL"))
  if (is.null(genes_path) || !file.exists(genes_path)) {
    annotation_unavailable("--annotate 模式缺少 --genes 文件", keytype = keytype)
  }
  genes <- tryCatch(readLines(genes_path, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) character(0))
  genes <- unique(genes[nzchar(trimws(genes))])
  if (length(genes) == 0) {
    annotation_unavailable("基因列表为空，没有可映射的 ID", keytype = keytype)
  }
  pkg <- switch(species,
                human = "org.Hs.eg.db",
                mouse = "org.Mm.eg.db",
                NULL)
  if (is.null(pkg)) {
    annotation_unavailable(paste0("不支持的物种（只支持 human / mouse）：", species),
                           keytype = keytype, n_input = length(genes))
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    annotation_unavailable(paste0("缺少 R 注释包：", pkg), database = pkg,
                           keytype = keytype, n_input = length(genes))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    annotation_unavailable("缺少 R 注释包：AnnotationDbi", database = pkg,
                           keytype = keytype, n_input = length(genes))
  }
  db <- tryCatch(getExportedValue(pkg, pkg), error = function(e) NULL)
  if (is.null(db)) {
    annotation_unavailable(paste0("无法加载注释库对象：", pkg), database = pkg,
                           keytype = keytype, n_input = length(genes))
  }
  valid_types <- tryCatch(AnnotationDbi::keytypes(db), error = function(e) character(0))
  if (!keytype %in% valid_types) {
    annotation_unavailable(paste0(pkg, " 不支持 keytype：", keytype),
                           database = pkg, keytype = keytype, n_input = length(genes))
  }
  version <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "")
  one_to_many <- 0L
  if (keytype == "SYMBOL") {
    # 注意：mapIds(column="SYMBOL", keytype="SYMBOL") 会把未知 key 原样返回（看起来
    # 100% 覆盖）。必须用 keys(db, "SYMBOL") 校验，未知 symbol 留空，不伪造映射。
    valid_keys <- tryCatch(AnnotationDbi::keys(db, keytype = "SYMBOL"),
                           error = function(e) character(0))
    if (length(valid_keys) == 0) {
      annotation_unavailable(paste0("无法读取 ", pkg, " 的 SYMBOL keys"),
                             database = pkg, version = version, keytype = keytype,
                             n_input = length(genes))
    }
    idx <- match(toupper(genes), toupper(valid_keys))
    symbols <- ifelse(is.na(idx), "", valid_keys[idx])
  } else {
    mapped <- tryCatch(
      suppressMessages(AnnotationDbi::mapIds(db, keys = genes, column = "SYMBOL",
                                             keytype = keytype, multiVals = "first")),
      error = function(e) NULL)
    if (is.null(mapped)) {
      annotation_unavailable(paste0("AnnotationDbi::mapIds 失败（", keytype, " → SYMBOL）"),
                             database = pkg, version = version, keytype = keytype,
                             n_input = length(genes))
    }
    multi <- tryCatch(
      suppressMessages(AnnotationDbi::mapIds(db, keys = genes, column = "SYMBOL",
                                             keytype = keytype, multiVals = "list")),
      error = function(e) NULL)
    one_to_many <- if (is.null(multi)) 0L else
      sum(vapply(multi, function(values) length(unique(values)) > 1, logical(1)))
    symbols <- as.character(mapped)
    symbols[is.na(symbols)] <- ""
  }
  keep <- nzchar(symbols)
  n_mapped <- sum(keep)
  many_to_one <- if (n_mapped > 0) sum(duplicated(symbols[keep])) else 0L
  frame <- data.frame(input_id = genes,
                      symbol = symbols,
                      stringsAsFactors = FALSE)
  write.table(frame, file = file.path(out_dir, "annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  write_meta(list(status = "ok", database = pkg, database_version = version,
                  keytype = keytype, n_input = length(genes), n_mapped = n_mapped,
                  one_to_many = one_to_many, many_to_one = many_to_one,
                  notes = "AnnotationDbi::mapIds(multiVals='first')；无 GPL 平台表时不做 probe 猜测",
                  warnings = ""))
  quit(save = "no", status = 0)
}

if (annotate_mode) {
  run_annotation()
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# 读取对象
# ---------------------------------------------------------------------------

read_objects <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.rds$", lower)) {
    obj <- tryCatch(readRDS(path), error = function(e) fail(paste0("readRDS 失败：", conditionMessage(e))))
    return(list(values = list(obj), names = c(basename(path))))
  }
  env <- new.env(parent = emptyenv())
  ok <- tryCatch({
    load(path, envir = env)
    TRUE
  }, error = function(e) {
    fail(paste0("RData 加载失败：", conditionMessage(e)))
    FALSE
  })
  if (!ok) return(NULL)
  objs <- ls(env)
  if (length(objs) == 0) fail("RData 里没有任何对象")
  values <- lapply(objs, function(n) get(n, envir = env))
  names(values) <- objs
  list(values = values, names = objs)
}

# ---------------------------------------------------------------------------
# 统一成 基因 × 样本 数值矩阵
# ---------------------------------------------------------------------------

.need <- function(pkg) requireNamespace(pkg, quietly = TRUE)

.extract <- function(obj) {
  # 返回 list(mat=基因×样本 matrix, pdata=data.frame|NULL, class=..., label=...)
  cls <- paste(class(obj), collapse = "/")
  mat <- NULL
  pdata <- NULL
  label <- NULL
  try({
    if (inherits(obj, "ExpressionSet")) {
      if (!.need("Biobase")) fail("缺少 Biobase 包，无法读取 ExpressionSet")
      mat <- Biobase::exprs(obj)
      pdata <- as.data.frame(Biobase::pData(obj))
      fd <- tryCatch(Biobase::fData(obj), error = function(e) NULL)
      if (!is.null(fd) && nrow(fd) == nrow(mat)) {
        rn <- rownames(fd)
        if (!is.null(rn) && all(nzchar(rn))) rownames(mat) <- rn
      }
      if (is.null(rownames(pdata)) || !all(nzchar(rownames(pdata)))) {
        rownames(pdata) <- colnames(mat)
      }
    } else if (methods::is(obj, "SummarizedExperiment")) {
      if (!.need("SummarizedExperiment")) fail("缺少 SummarizedExperiment 包，无法读取该对象")
      assays <- SummarizedExperiment::assays(obj)
      if (length(assays) == 0) fail("SummarizedExperiment 没有 assay")
      mat <- as.matrix(assays[[1]])
      cd <- as.data.frame(SummarizedExperiment::colData(obj))
      if (nrow(cd) == ncol(mat)) {
        pdata <- cd
        if (is.null(rownames(pdata)) || !all(nzchar(rownames(pdata)))) rownames(pdata) <- colnames(mat)
      }
      rd <- rownames(obj)
      if (!is.null(rd) && length(rd) == nrow(mat)) rownames(mat) <- rd
    } else if (inherits(obj, "Seurat")) {
      if (!.need("SeuratObject")) fail("缺少 SeuratObject 包，无法读取 Seurat 对象")
      get_layer <- function() {
        tryCatch(
          SeuratObject::GetAssayData(obj, assay = "RNA", slot = "data"),
          error = function(e) {
            if (exists("LayerData", envir = asNamespace("SeuratObject"))) {
              SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
            } else {
              stop(e)
            }
          }
        )
      }
      mat <- get_layer()
      meta <- tryCatch(obj@meta.data, error = function(e) NULL)
      if (!is.null(meta)) {
        pdata <- as.data.frame(meta)
        if (is.null(rownames(pdata)) || !all(nzchar(rownames(pdata)))) rownames(pdata) <- colnames(mat)
      }
    } else if (inherits(obj, "DGEList")) {
      mat <- obj$counts
      sm <- tryCatch(as.data.frame(obj$samples), error = function(e) NULL)
      if (!is.null(sm)) pdata <- sm
    } else if (is.matrix(obj) || (requireNamespace("Matrix", quietly = TRUE) && inherits(obj, "Matrix"))) {
      mat <- as.matrix(obj)
    } else if (is.data.frame(obj)) {
      frame <- obj
      numeric_cols <- vapply(frame, function(col) is.numeric(col) || is.integer(col), logical(1))
      if (sum(numeric_cols) < 2) {
        return(NULL)
      }
      first <- frame[[1]]
      if (!is.numeric(first) && !is.integer(first) && !is.logical(first)) {
        rn <- as.character(first)
        frame <- frame[, numeric_cols, drop = FALSE]
        rownames(frame) <- rn
        mat <- as.matrix(frame)
      } else {
        mat <- as.matrix(frame[, numeric_cols, drop = FALSE])
      }
    } else if (is.list(obj) && !is.data.frame(obj)) {
      candidates <- c("counts", "data", "expr", "expression", "mat", "matrix", "assay", "tpm", "fpkm")
      hit <- intersect(candidates, names(obj))
      if (length(hit) > 0) {
        inner <- obj[[hit[1]]]
        if (is.matrix(inner) || is.data.frame(inner) ||
            (requireNamespace("Matrix", quietly = TRUE) && inherits(inner, "Matrix"))) {
          mat <- as.matrix(inner)
          for (info_name in c("samples", "sample_info", "pdata", "phenoData", "metadata")) {
            if (info_name %in% names(obj)) {
              info <- obj[[info_name]]
              if (is.data.frame(info)) {
                pdata <- info
                break
              } else if (is.list(info) && !is.null(info$samples)) {
                pdata <- as.data.frame(info$samples)
                break
              }
            }
          }
        }
      }
    }
  }, silent = TRUE)
  if (is.null(mat)) return(NULL)
  mat <- as.matrix(mat)
  if (nrow(mat) < 1 || ncol(mat) < 1) return(NULL)
  # 非数值矩阵尝试数值化
  if (!is.numeric(mat)) {
    converted <- suppressWarnings(matrix(as.numeric(as.character(mat)),
                                         nrow = nrow(mat), ncol = ncol(mat)))
    if (all(is.na(converted))) {
      return(NULL)
    }
    dimnames(converted) <- dimnames(mat)
    mat <- converted
  }
  storage.mode(mat) <- "double"
  # 方向：基因行数通常 >= 样本列数；若行数 < 列数，认定为样本×基因，转置
  if (nrow(mat) < ncol(mat)) {
    add_note("矩阵原始方向为 样本×基因，已转置")
    mat <- t(mat)
  }
  if (is.null(rownames(mat)) || !all(nzchar(rownames(mat)))) {
    rownames(mat) <- paste0("GENE", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat)) || !all(nzchar(colnames(mat)))) {
    colnames(mat) <- paste0("S", seq_len(ncol(mat)))
  }
  rownames(mat) <- make.unique(rownames(mat))
  colnames(mat) <- make.unique(colnames(mat))
  if (!is.null(pdata)) {
    rownames(pdata) <- as.character(rownames(pdata))
    common <- intersect(colnames(mat), rownames(pdata))
    if (length(common) >= 0.5 * ncol(mat)) {
      pdata <- pdata[match(colnames(mat), rownames(pdata)), , drop = FALSE]
      rownames(pdata) <- colnames(mat)
    } else {
      add_warning("pData 的行名与矩阵列名匹配不上，样本表已丢弃")
      pdata <- NULL
    }
  }
  list(mat = mat, pdata = pdata, class = cls)
}

# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

loaded <- read_objects(in_path)
values <- loaded$values
obj_names <- loaded$names

extracted <- list()
candidate_rows <- list()
for (i in seq_along(values)) {
  result <- tryCatch(.extract(values[[i]]), error = function(e) NULL)
  if (is.null(result)) next
  extracted[[obj_names[i]]] <- result
  candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
    object = obj_names[i],
    class = result$class,
    n_genes = nrow(result$mat),
    n_samples = ncol(result$mat),
    label = paste0(obj_names[i], "（", result$class, "，",
                   nrow(result$mat), " 基因 × ", ncol(result$mat), " 样本）"),
    stringsAsFactors = FALSE
  )
}

if (length(extracted) == 0) {
  fail(paste0("文件里没有可识别的表达矩阵对象（支持 ExpressionSet / SummarizedExperiment / Seurat /",
              " DGEList / matrix / data.frame）"))
}

chosen_name <- NULL
if (!is.null(object_name) && nzchar(object_name)) {
  if (!object_name %in% names(extracted)) {
    fail(paste0("找不到对象「", object_name, "」；可用对象：", paste(names(extracted), collapse = "、")))
  }
  chosen_name <- object_name
} else if (length(extracted) == 1) {
  chosen_name <- names(extracted)[1]
} else {
  # 多个候选：优先唯一一个 ExpressionSet/SummarizedExperiment/Seurat
  preferred <- names(extracted)[vapply(extracted, function(item) {
    grepl("ExpressionSet|SummarizedExperiment|Seurat", item$class)
  }, logical(1))]
  if (length(preferred) == 1) {
    chosen_name <- preferred[1]
    add_note(paste0("文件里有多个对象，自动选择唯一的高阶对象：", chosen_name))
  } else {
    cand <- do.call(rbind, candidate_rows)
    write.table(cand, file = file.path(out_dir, "candidates.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write_meta(list(status = "ambiguous",
                    error = paste0("文件里有多个可用对象：", paste(names(extracted), collapse = "、")),
                    n_candidates = nrow(cand)))
    quit(save = "no", status = 0)
  }
}

chosen <- extracted[[chosen_name]]
mat <- chosen$mat
pdata <- chosen$pdata
add_note(paste0("对象名：", chosen_name, "；R class：", chosen$class))
has_counts <- all(is.finite(mat)) && all(mat == round(mat)) && min(mat) >= 0

if (inspect_only) {
  write_meta(list(status = "ok", class = chosen$class, object = chosen_name,
                  orientation = "genes_by_samples",
                  n_genes = nrow(mat), n_samples = ncol(mat),
                  has_counts = has_counts,
                  notes = NOTES, warnings = WARNINGS))
  quit(save = "no", status = 0)
}

# 丢弃全 NA 行/列，缺失按 0 填充
if (any(!is.finite(mat))) {
  bad_rows <- apply(mat, 1, function(row) all(!is.finite(row)))
  bad_cols <- apply(mat, 2, function(col) all(!is.finite(col)))
  if (any(bad_rows) || any(bad_cols)) {
    add_warning(paste0("丢弃 ", sum(bad_rows), " 个全缺失基因行 / ", sum(bad_cols), " 个全缺失样本列"))
    mat <- mat[!bad_rows, !bad_cols, drop = FALSE]
    if (!is.null(pdata)) {
      pdata <- pdata[intersect(colnames(mat), rownames(pdata)), , drop = FALSE]
    }
  }
  n_missing <- sum(!is.finite(mat))
  if (n_missing > 0) {
    add_warning(paste0(n_missing, " 个缺失值按 0 填充"))
    mat[!is.finite(mat)] <- 0
  }
}
if (nrow(mat) < 1 || ncol(mat) < 1) {
  fail("矩阵过滤缺失后为空")
}

if (!.need("Matrix")) {
  # 极端兜底：没有 Matrix 包时写稠密文本（Matrix 是 R 推荐包，正常都有）
  write.table(mat, file = file.path(out_dir, "matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  write_meta(list(status = "ok", class = chosen$class, object = chosen_name,
                  orientation = "genes_by_samples", n_genes = nrow(mat), n_samples = ncol(mat),
                  has_counts = has_counts, dense = TRUE, notes = NOTES, warnings = WARNINGS))
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(Matrix))
sparse <- Matrix::Matrix(mat, sparse = TRUE)
Matrix::writeMM(sparse, file = file.path(out_dir, "matrix.mtx"))
writeLines(rownames(mat), file.path(out_dir, "genes.tsv"), useBytes = TRUE)
writeLines(colnames(mat), file.path(out_dir, "samples.tsv"), useBytes = TRUE)
if (!is.null(pdata)) {
  write.csv(pdata, file = file.path(out_dir, "obs.csv"), quote = TRUE)
}

write_meta(list(status = "ok", class = chosen$class, object = chosen_name,
                orientation = "genes_by_samples",
                n_genes = nrow(mat), n_samples = ncol(mat),
                has_counts = has_counts,
                notes = NOTES, warnings = WARNINGS))
quit(save = "no", status = 0)
