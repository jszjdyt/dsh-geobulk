#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — Bulk 转录组全格式导入层（CONTRACT §7）

把本机各种 bulk 表达数据统一读成：

    X            np.ndarray float32 (samples × genes)
    genes        list[str]（列名 = 基因）
    samples      list[str]（行名 = 样本）
    sample_info  pandas.DataFrame(index=samples)
    notes/warnings list[str]
    meta         {"n_samples","n_genes","species","has_counts","value_kind", ...}

支持的格式（`SUPPORTED_FORMATS`）：
    matrix_text / xlsx / geo_series / archive / rds / rdata
    + 非编码模块数据 circ / mirna / lncrna / fusion / splicing

对外 API（Host / registry 按它接线）：
    inspect_source(source, workdir=None, options=None) -> dict      # 不加载整表
    load_source(source, workdir=None, options=None, progress=None)  -> dict   # §7.1 payload
    load_module(source, kind=None, workdir=None, options=None)      -> dict   # attach 用

只依赖标准库 + numpy/pandas + 本机 Rscript（仅 rds/RData 用）。
"""

from __future__ import annotations

import bz2
import csv
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from bg_common import EngineError


# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

SUPPORTED_FORMATS = [
    "matrix_text", "xlsx", "geo_series", "archive", "rds", "rdata",
    "circ", "mirna", "lncrna", "fusion", "splicing", "auto",
]

MODULE_KINDS = ["circ", "mirna", "lncrna", "fusion", "splicing"]

FORMAT_INFO = {
    "matrix_text": {
        "label": "表达矩阵（CSV/TSV/TXT，可 gz）",
        "extensions": [".csv", ".tsv", ".txt", ".tab", ".matrix", ".counts", ".gz", ".bz2"],
        "description": "首行/首列自动识别，支持逗号/制表符/分号/空白分隔，自动判断样本×基因方向",
    },
    "xlsx": {
        "label": "Excel 表达矩阵",
        "extensions": [".xlsx", ".xlsm", ".xls"],
        "description": "读取第一个可用 sheet，可用 sheet 参数指定；多 sheet 会列出候选",
    },
    "geo_series": {
        "label": "GEO series matrix",
        "extensions": [".txt", ".txt.gz"],
        "description": "含 !series_matrix_table_begin 段的 GEO 系列矩阵；!Sample_* 元数据转成样本表",
    },
    "archive": {
        "label": "压缩包（zip/tar/tar.gz）",
        "extensions": [".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tar.xz"],
        "description": "解压后自动寻找系列矩阵/表达矩阵；多个候选时用 member 指定",
    },
    "rds": {
        "label": "R 对象（.rds）",
        "extensions": [".rds", ".rds.gz"],
        "description": "ExpressionSet / SummarizedExperiment / Seurat / DGEList / matrix / data.frame",
    },
    "rdata": {
        "label": "R 工作空间（.RData）",
        "extensions": [".RData", ".rdata", ".rda"],
        "description": "加载后从多个对象里选目标对象，可用 object 参数指定",
    },
    "circ": {
        "label": "circRNA 表达 / 鉴定表",
        "extensions": [".tsv", ".csv", ".txt"],
        "description": "识别 chr1:12345-12789:+ / hsa_circ_0001234 等 ID 与样本列",
    },
    "mirna": {
        "label": "miRNA 表达矩阵",
        "extensions": [".tsv", ".csv", ".txt"],
        "description": "miRBase / mature 名（hsa-miR-21-5p），内置常见 precursor→mature 映射",
    },
    "lncrna": {
        "label": "lncRNA 表达矩阵",
        "extensions": [".tsv", ".csv", ".txt"],
        "description": "ENSG/GENCODE ID 或 LINC… / …-AS1 符号，启发式标注 lncRNA 类型",
    },
    "fusion": {
        "label": "融合基因结果表",
        "extensions": [".tsv", ".csv", ".txt"],
        "description": "STAR-Fusion / Arriba / SOAPfuse 输出，识别融合对、junction/spanning reads",
    },
    "splicing": {
        "label": "可变剪接结果表",
        "extensions": [".tsv", ".txt", ".csv"],
        "description": "rMATS（SE/A5SS/A3SS/MXE/RI，IJC_SAMPLE_1/IncLevel1）或 SUPPA/通用 PSI 表",
    },
    "auto": {"label": "自动识别", "extensions": [], "description": "按扩展名与文件内容自动判断格式"},
}

TEXT_SUFFIXES = (".csv", ".tsv", ".txt", ".tab", ".matrix", ".counts", ".tsv2")
ARCHIVE_SUFFIXES = (".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz", ".tbz2", ".tar.xz", ".txz")
XLSX_SUFFIXES = (".xlsx", ".xlsm", ".xls")
R_SUFFIXES = (".rds", ".rdata", ".rda")

MAX_UNPACK_BYTES = 20 * 1024 * 1024 * 1024   # 解压总大小上限 20GB
MAX_DOWNLOAD_BYTES = 8 * 1024 * 1024 * 1024  # 单文件下载上限 8GB
MAX_TEXT_BYTES_FOR_FULL = 1024 ** 3          # >1GB 文本矩阵直接拒绝
R_TIMEOUT = 300
DOWNLOAD_TIMEOUT = 120
_UA = "dsh-geobulk/1.0 (bulk RNA-seq workbench)"

# 强基因名特征（大写符号 / Ensembl / 芯片探针 / miRBase）
_STRONG_GENE_RE = re.compile(r"^(?:MT-[A-Z0-9]+|[A-Z][A-Z0-9_.\-]{2,}|ENS[A-Z]*G\d{6,}|"
                             r"(?:hsa|mmu|rno|ptr|gga)-[A-Za-z0-9.\-]+|\d+_[a-z]+_at)$")
# 样本名特征（GEO / SRA / TCGA / 常见分组命名）
_SAMPLE_RE = re.compile(
    r"^(?:GSM\d+|GSE\d+|GS[EM]\d+|TCGA[-_][A-Za-z0-9\-]+|SRR\d+|ERR\d+|DRR\d+|"
    r"(?i:ctrl|control|case|treat|treated|normal|tumor|tumour|sample|patient|pt|rep|"
    r"replicate|mock|blank|wt|ko|sh\d*|oe|disease|healthy)[-_]?\d*)$"
)
_SAFE_MEMBER_RE = re.compile(r"^[^/\\]*$")


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------


def _err(code: str, message: str, **details) -> EngineError:
    return EngineError(code, message, **details)


def _emit(progress, message, pct=None):
    """进度回调容错：同时兼容 progress(message,pct) / progress(pct)。"""
    if not progress:
        return
    try:
        progress(message, pct)
        return
    except TypeError:
        pass
    except Exception:  # noqa: BLE001
        return
    try:
        progress(pct if pct is not None else message)
    except Exception:  # noqa: BLE001
        pass


def _workdirs(workdir) -> dict:
    base = Path(workdir).expanduser() if workdir else Path(tempfile.mkdtemp(prefix="bg-io-"))
    try:
        base.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise _err("bad_request", f"无法创建工作目录：{base}（{error}）")
    dirs = {
        "base": base,
        "unpack": base / "unpack",
        "tmp": base / "tmp",
        "downloads": base / "downloads",
    }
    for key in ("unpack", "tmp", "downloads"):
        try:
            dirs[key].mkdir(parents=True, exist_ok=True)
        except OSError:
            dirs[key] = Path(tempfile.mkdtemp(prefix=f"{key}-"))
    return dirs


def _looks_like_url(source) -> bool:
    return bool(re.match(r"^(?:https?|ftp)://", str(source or "").strip(), re.I))


def _is_plain_filename(name: str) -> bool:
    text = str(name or "")
    return bool(text) and bool(_SAFE_MEMBER_RE.match(text)) and ".." not in text


def _dedupe(names, prefix="gene") -> list[str]:
    """去空、去引号、同名去重（保持顺序）。"""
    out = []
    seen = {}
    for i, raw in enumerate(names):
        text = str(raw).strip().strip('"').strip("'").strip()
        if not text:
            text = f"{prefix}_{i + 1}"
        if text in seen:
            seen[text] += 1
            text = f"{text}_{seen[text]}"
        else:
            seen[text] = 1
        out.append(text)
    return out


def _open_text(path):
    """文本读取句柄（gz/bz2 透明解压）。"""
    path = Path(path)
    name = path.name.lower()
    if name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if name.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _read_head(path, size: int = 4 * 1024 * 1024) -> str:
    with _open_text(path) as handle:
        return handle.read(size)


def _is_numeric_text(value) -> bool:
    try:
        text = str(value).strip().replace(",", "")
        if text == "" or text.lower() in ("na", "nan", "null", "none", "-"):
            return False
        float(text)
        return True
    except Exception:  # noqa: BLE001
        return False


def _split_line(line: str, delimiter) -> list[str]:
    if delimiter == "whitespace":
        return [item for item in re.split(r"\s+", line.strip()) if item != ""]
    try:
        return next(csv.reader([line], delimiter=delimiter, quotechar='"'))
    except Exception:  # noqa: BLE001
        return line.split(delimiter)


def _sniff_delimiter(lines) -> str:
    """分隔符嗅探：逗号 / 制表符 / 分号 / 空白。返回 ',' '\t' ';' 或 'whitespace'。"""
    best, best_score = "\t", -1.0
    for delimiter in (",", "\t", ";", "whitespace"):
        counts = []
        for line in lines[:25]:
            if not line.strip():
                continue
            counts.append(len(_split_line(line, delimiter)))
        if not counts or max(counts) < 2:
            continue
        mode = max(set(counts), key=counts.count)
        consistency = counts.count(mode) / len(counts)
        # 字段数×一致度；字段数是硬信号，避免单列被误判
        score = mode * consistency
        # 首行（通常是表头）应与数据行字段数一致：rMATS 里 IncLevel1 的值形如
        # "0.8,0.85,0.82"，逗号会让数据行字段数暴增，用首行一致性把误判压下去
        if counts[0] != mode:
            score *= 0.3
        if score > best_score:
            best_score, best = score, delimiter
    return best


def _detect_header(cells) -> bool:
    """首行除首列外是否全为数值 → 无表头（首行即数据）。"""
    rest = [cell for cell in list(cells)[1:] if str(cell).strip() != ""]
    if not rest:
        return True
    return not all(_is_numeric_text(cell) for cell in rest)


def _strong_gene(text: str) -> bool:
    return bool(_STRONG_GENE_RE.match(text)) and not bool(_SAMPLE_RE.match(text))


def _sample_like(text: str) -> bool:
    return bool(_SAMPLE_RE.match(text))


def _score(labels, predicate) -> float:
    hits = checked = 0
    for raw in list(labels)[:200]:
        text = str(raw).strip()
        if not text:
            continue
        checked += 1
        if predicate(text):
            hits += 1
    return hits / checked if checked else 0.0


def _infer_rows_are_genes(row_labels, col_labels, shape, notes) -> bool:
    """判断矩阵的行是否代表基因（True ⇒ 需要转置成 samples×genes）。"""
    row_gene = _score(row_labels, _strong_gene)
    col_gene = _score(col_labels, _strong_gene)
    row_sample = _score(row_labels, _sample_like)
    col_sample = _score(col_labels, _sample_like)

    if row_sample >= 0.5 and col_sample < 0.2:
        return False
    if col_sample >= 0.5 and row_sample < 0.2:
        return True
    if row_gene >= 0.75 and col_gene < 0.30:
        return True
    if col_gene >= 0.75 and row_gene < 0.30:
        return False
    # 形状启发式：bulk 的基因数通常远多于样本数
    if shape[0] > shape[1]:
        return True
    if shape[1] > shape[0]:
        return False
    notes.append("矩阵为方阵，方向无法从数据判断，默认按「行 = 基因」处理；如方向反了请指定 transpose 参数")
    return True


def _infer_species(genes) -> str:
    sample = [str(g) for g in list(genes)[:2000]]
    if not sample:
        return "unknown"
    mouse_ens = sum(1 for g in sample if g.upper().startswith("ENSMUSG"))
    human_ens = sum(1 for g in sample if g.upper().startswith("ENSG"))
    if mouse_ens and mouse_ens >= human_ens:
        return "mouse"
    if human_ens:
        return "human"
    mouse_case = sum(1 for g in sample if re.fullmatch(r"[A-Z][a-z][A-Za-z0-9\-]*", g))
    upper_case = sum(1 for g in sample if _strong_gene(g))
    if mouse_case > 0.35 * len(sample) and mouse_case > upper_case:
        return "mouse"
    if upper_case > 0.35 * len(sample):
        return "human"
    return "unknown"


def _detect_value_kind(X) -> tuple[str, bool, str]:
    """表达值类型判定 → (value_kind, has_counts, 说明)。"""
    arr = np.asarray(X)
    if arr.ndim != 2 or arr.size == 0:
        return "unknown", False, ""
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return "unknown", False, ""
    maximum = float(finite.max())
    minimum = float(finite.min())
    if minimum < 0:
        return "lognorm", False, "存在负值，按标准化后的 log/ratio 数据处理"
    integer_like = bool(np.all(np.abs(finite - np.round(finite)) < 1e-6))
    if integer_like and maximum >= 1:
        return "counts", True, "整数且非负，判定为原始 counts"
    if maximum < 50:
        return "lognorm", False, "最大值 < 50，判定为已 log 标准化数据"
    if maximum > 100:
        return "tpm", False, "最大值 > 100 且全正，判定为 TPM/FPKM（两者无法从数值区分）"
    return "fpkm", False, "数值范围像 FPKM/TPM"


def _align_sample_info(sample_info, samples) -> pd.DataFrame:
    samples = [str(s) for s in samples]
    index = pd.Index(samples, name="sample")
    if sample_info is None:
        return pd.DataFrame(index=index)
    frame = sample_info.copy()
    frame.index = [str(item) for item in frame.index]
    if list(frame.index) == samples:
        frame.index = index
        return frame
    lower_map = {}
    for label in frame.index:
        lower_map.setdefault(str(label).strip().lower(), label)
    rows = []
    for sample in samples:
        key = lower_map.get(sample.strip().lower())
        rows.append(frame.loc[key] if key is not None else pd.Series(dtype=object, name=sample))
    out = pd.DataFrame(rows, index=index)
    for column in frame.columns:
        if column not in out.columns:
            out[column] = np.nan
    out = out.reindex(columns=list(frame.columns))
    return out


def _payload(fmt: str, source, X, genes, samples, sample_info=None, notes=None, warnings=None, meta=None) -> dict:
    """构造 §7.1 的统一导入结果。"""
    if X is None:
        matrix = None
        n_samples, n_genes = 0, 0
    else:
        matrix = np.asarray(X, dtype=np.float32)
        if matrix.ndim != 2:
            raise _err("internal", "导入结果矩阵维度异常（应为 样本×基因）")
        n_samples, n_genes = matrix.shape
    samples = _dedupe(samples if samples is not None else [f"S{i + 1}" for i in range(n_samples or 0)], "sample")
    genes = _dedupe(genes if genes is not None else [f"GENE{i + 1}" for i in range(n_genes or 0)], "gene")
    if matrix is not None and (len(samples), len(genes)) != matrix.shape:
        raise _err("internal", f"导入结果名称与矩阵不匹配：{len(samples)}×{len(genes)} vs {matrix.shape}")
    value_kind, has_counts, kind_note = _detect_value_kind(matrix if matrix is not None else [])
    species = _infer_species(genes)
    out_meta = {
        "n_samples": int(n_samples),
        "n_genes": int(n_genes),
        "species": species,
        "has_counts": bool(has_counts),
        "value_kind": value_kind,
        "source": str(source),
    }
    if kind_note:
        out_meta["value_kind_note"] = kind_note
    for key, value in (meta or {}).items():
        out_meta[key] = value
    out_notes = list(notes or [])
    if kind_note:
        out_notes.append(kind_note)
    return {
        "format": str(fmt),
        "name": Path(urllib.parse.urlparse(str(source)).path or str(source)).name or str(source),
        "X": matrix,
        "genes": genes,
        "samples": samples,
        "sample_info": _align_sample_info(sample_info, samples) if matrix is not None else (
            sample_info if isinstance(sample_info, pd.DataFrame) else pd.DataFrame()),
        "notes": [str(item) for item in out_notes],
        "warnings": [str(item) for item in (warnings or [])],
        "meta": out_meta,
    }


# ---------------------------------------------------------------------------
# 样本信息表 / 临床表 join
# ---------------------------------------------------------------------------


def _read_any_table(path, options=None) -> pd.DataFrame:
    """读取 CSV/TSV/TXT/XLSX 为 DataFrame（不强制数值化）。"""
    path = Path(path)
    options = options or {}
    name = path.name.lower()
    if name.endswith(XLSX_SUFFIXES):
        frame, _ = _xlsx_frame(path, options)
        return frame
    head = _read_head(path, 1024 * 1024)
    delimiter = options.get("delimiter") or _sniff_delimiter(head.splitlines())
    try:
        if delimiter == "whitespace":
            # engine="python" 不支持 low_memory，必须去掉（否则 pandas 直接 ValueError）
            with _open_text(path) as handle:
                return pd.read_csv(handle, sep=r"\s+", engine="python")
        return pd.read_csv(path, sep=delimiter, low_memory=False)
    except (UnicodeDecodeError, LookupError):
        raise _err("unsupported", "文件不是可读文本（可能是二进制或非 UTF-8 编码），"
                                  "请确认格式或用 data_format 指定")
    except pd.errors.ParserError as error:
        raise _err("bad_request", f"无法解析表格（{Path(path).name}）：{error}")
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", f"无法读取样本信息表（{Path(path).name}）：{error}")


def _pick_meta_key(frame: pd.DataFrame, samples) -> tuple[str | None, pd.DataFrame]:
    """找样本名匹配列：优先行名，其次 sample/ID 类列，最后第一列。"""
    samples_lower = {str(s).strip().lower() for s in samples}

    def _match_ratio(values) -> float:
        values = [str(v).strip().lower() for v in values]
        if not values:
            return 0.0
        hits = sum(1 for v in values if v in samples_lower)
        return hits / len(values)

    if _match_ratio(list(frame.index)) >= 0.5:
        return None, frame
    candidates = ["sample", "sample_id", "sampleid", "samples", "id", "geo_accession", "gsm",
                  "name", "sample_name", "title"]
    lower_cols = {str(c).strip().lower(): c for c in frame.columns}
    for key in candidates:
        if key in lower_cols and _match_ratio(list(frame[lower_cols[key]])) >= 0.3:
            return lower_cols[key], frame
    # 兜底：逐列试匹配
    best_col, best_ratio = None, 0.0
    for column in frame.columns:
        ratio = _match_ratio(list(frame[column]))
        if ratio > best_ratio:
            best_col, best_ratio = column, ratio
    if best_col is not None and best_ratio >= 0.3:
        return best_col, frame
    first = frame.columns[0] if len(frame.columns) else None
    return first, frame


def _join_metadata(sample_info: pd.DataFrame, samples, path, label, warnings) -> pd.DataFrame:
    """把样本信息表/临床表按样本名 join 进 sample_info（大小写不敏感）。"""
    try:
        source = _read_any_table(path)
    except EngineError:
        raise
    if source.empty:
        warnings.append(f"{label}表为空：{Path(path).name}")
        return sample_info
    key_column, source = _pick_meta_key(source, samples)
    if key_column is None:
        mapping = {str(index).strip().lower(): index for index in source.index}
    else:
        mapping = {}
        for _, row in source.iterrows():
            mapping[str(row[key_column]).strip().lower()] = row
    out = sample_info.copy()
    matched = 0
    columns = [c for c in source.columns if c != key_column]
    for column in columns:
        if column in out.columns:
            target = f"{column}_{label}"
            suffix = 2
            while target in out.columns:
                target = f"{column}_{label}{suffix}"
                suffix += 1
        else:
            target = column
        values = []
        for sample in samples:
            key = str(sample).strip().lower()
            if key not in mapping:
                values.append(np.nan)
                continue
            matched += 1
            if key_column is None:
                values.append(source.loc[mapping[key], column] if column in source.columns else np.nan)
            else:
                values.append(mapping[key][column])
        out[target] = values
    unmatched = sorted({str(s) for s in samples if str(s).strip().lower() not in mapping})
    if unmatched:
        warnings.append(f"{label}表有 {len(unmatched)} 个样本未匹配：{', '.join(unmatched[:5])}"
                        + ("…" if len(unmatched) > 5 else ""))
    if matched == 0:
        warnings.append(f"{label}表与样本名完全匹配不上（已忽略该表）")
    return out


# ---------------------------------------------------------------------------
# 文本矩阵 / GEO series matrix
# ---------------------------------------------------------------------------


def _frame_from_text(path, options, notes, warnings):
    """读取任意分隔文本矩阵 → (DataFrame, row_labels, col_labels)（尚未判断方向）。"""
    options = options or {}
    head = _read_head(path, 4 * 1024 * 1024)
    if not head.strip():
        raise _err("bad_request", f"文件是空的：{Path(path).name}")
    delimiter = options.get("delimiter") or _sniff_delimiter(head.splitlines())
    if options.get("delimiter") and delimiter != "whitespace" and len(delimiter) > 1:
        raise _err("bad_request", f"delimiter 参数不合法：{delimiter!r}（应为逗号/制表符/分号或 whitespace）")
    lines = head.splitlines()
    first_cells = _split_line(lines[0], delimiter) if lines else [""]
    header = _detect_header(first_cells)
    data_lines = [line for line in lines[1:41] if line.strip()]
    has_index = False
    if data_lines:
        first_cells_data = [_split_line(line, delimiter)[0] if _split_line(line, delimiter) else "" for line in data_lines]
        non_numeric = sum(1 for cell in first_cells_data if not _is_numeric_text(cell))
        has_index = non_numeric >= max(1, int(0.6 * len(first_cells_data)))
    first_header = str(first_cells[0]).strip().lower() if first_cells else ""
    if not has_index and header and first_header in (
            "", "gene", "gene_id", "geneid", "id", "probe", "probe_id", "probeid",
            "transcript", "transcript_id", "symbol", "ensembl", "ensembl_id", "features", "feature"):
        has_index = True  # 数值型探针 ID / 基因 ID 列（列名是 gene/id 等）

    read_kwargs = dict(header=0 if header else None, low_memory=False)
    read_kwargs["index_col"] = 0 if has_index else None
    sep = r"\s+" if delimiter == "whitespace" else delimiter
    try:
        if delimiter == "whitespace":
            # engine="python" 不支持 low_memory，必须去掉（否则 pandas 直接 ValueError）
            read_kwargs.update(sep=sep, engine="python")
            read_kwargs.pop("low_memory", None)
            with _open_text(path) as handle:
                frame = pd.read_csv(handle, **read_kwargs)
        else:
            read_kwargs["sep"] = delimiter
            frame = pd.read_csv(path, **read_kwargs)
    except EngineError:
        raise
    except (UnicodeDecodeError, LookupError):
        raise _err("unsupported", "文件不是可读文本（可能是二进制或非 UTF-8 编码），"
                                  "请确认格式或用 data_format 指定")
    except pd.errors.ParserError as error:
        raise _err("bad_request", f"无法解析矩阵（{Path(path).name}）：{error}")
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", f"读取矩阵失败（{Path(path).name}）：{error}")

    if options.get("gene_col") is not None and header:
        gene_col = options["gene_col"]
        if isinstance(gene_col, int) or (isinstance(gene_col, str) and gene_col.isdigit()):
            gene_col = frame.columns[int(gene_col)]
        if gene_col not in frame.columns:
            raise _err("bad_request", f"gene_col={gene_col!r} 不在表头里；可用列：{list(frame.columns)[:10]}")
        gene_values = frame[gene_col].astype(str).tolist()
        frame = frame.drop(columns=[gene_col])
        frame.index = pd.Index(gene_values)
        has_index = True  # 行名 = 基因
        notes.append(f"按 gene_col={gene_col} 指定基因列")
    notes.append("分隔符：%s" % ("TAB" if delimiter == "\t" else ("空白" if delimiter == "whitespace" else delimiter)))
    if not header:
        notes.append("未检测到表头，自动生成列名")
    if not has_index:
        notes.append("未检测到行名，自动生成行名")
    # 清理：去空列/空行
    frame = frame.loc[:, [str(c).strip() != "" for c in frame.columns]]
    frame.columns = [str(c) for c in frame.columns]
    frame.index = [str(i) for i in frame.index]
    return frame


def _finalize_matrix(frame: pd.DataFrame, fmt: str, source, options, notes, warnings, sample_info=None) -> dict:
    """数值化 + 方向判定 + 组装 payload。frame 允许非数值字符串。"""
    options = options or {}
    try:
        numeric = frame.apply(pd.to_numeric, errors="coerce")
    except Exception:
        numeric = frame
    total_cells = max(1, numeric.shape[0] * numeric.shape[1])
    nan_before = int(numeric.isna().sum().sum())
    drop_rows = int(numeric.isna().all(axis=1).sum()) if numeric.shape[1] else numeric.shape[0]
    drop_cols = int(numeric.isna().all(axis=0).sum()) if numeric.shape[0] else numeric.shape[1]
    if drop_cols and numeric.shape[1] > 1:
        numeric = numeric.loc[:, ~numeric.isna().all(axis=0)]
    if drop_rows and numeric.shape[0] > 1:
        numeric = numeric.loc[~numeric.isna().all(axis=1)]
    if numeric.empty or numeric.shape[0] == 0 or numeric.shape[1] == 0:
        raise _err("bad_request", f"矩阵里没有可解析的数值：{Path(str(source)).name}")
    if nan_before > 0:
        warnings.append(f"有 {nan_before} 个非数值/缺失单元格（{nan_before / total_cells:.2%}）已按 0 填充")
    if drop_rows:
        warnings.append(f"丢弃 {drop_rows} 个全空/全非数值行")
    if drop_cols:
        warnings.append(f"丢弃 {drop_cols} 个全空/全非数值列")
    # 少量 NaN 用 0 填充（GEO 常见）
    numeric = numeric.fillna(0.0)

    row_labels = [str(item) for item in numeric.index]
    col_labels = [str(item) for item in numeric.columns]
    transpose = options.get("transpose")
    if transpose is None:
        transpose = _infer_rows_are_genes(row_labels, col_labels, numeric.shape, notes)
    else:
        transpose = bool(transpose)
        notes.append("按用户指定 transpose=%s 处理" % transpose)
    values = numeric.to_numpy(dtype=np.float64)
    if transpose:
        values = values.T
        genes, samples = row_labels, col_labels
    else:
        genes, samples = col_labels, row_labels
    return _payload(fmt, source, values, genes, samples, sample_info=sample_info,
                    notes=notes, warnings=warnings,
                    meta={"orientation": "genes_by_samples" if transpose else "samples_by_genes"})


def _parse_geo_series(path, source) -> dict:
    """解析 GEO series matrix：table 段 → X，!Sample_* → sample_info。"""
    sample_lines: list[str] = []
    header_cells: list[str] | None = None
    rows: list[list[str]] = []
    in_table = False
    notes = ["从 !series_matrix_table_begin 段解析；!Sample_* 行已转成样本表"]
    warnings: list[str] = []
    size = Path(path).stat().st_size if Path(path).exists() else 0
    if size > MAX_TEXT_BYTES_FOR_FULL:
        raise _err("bad_request", f"GEO series matrix 过大（>{MAX_TEXT_BYTES_FOR_FULL // 1024 ** 3}GB），请先裁剪基因/样本再导入")
    try:
        with _open_text(path) as handle:
            for line in handle:
                stripped = line.strip()
                if stripped.startswith("!series_matrix_table_begin"):
                    in_table = True
                    continue
                if stripped.startswith("!series_matrix_table_end"):
                    break
                if in_table:
                    cells = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
                    if not any(cell.strip() for cell in cells):
                        continue
                    if header_cells is None:
                        header_cells = cells
                    else:
                        rows.append(cells)
                elif stripped.startswith("!Sample_"):
                    sample_lines.append(line.rstrip("\n"))
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", f"GEO series matrix 读取失败（{Path(path).name}）：{error}")
    if header_cells is None or not rows:
        raise _err("bad_request", f"GEO series matrix 里没有可用的表格段：{Path(path).name}")
    table_samples = _dedupe([str(cell).strip() for cell in header_cells[1:]], "sample")
    n_table_samples = len(table_samples)
    values = np.full((len(rows), n_table_samples), np.nan, dtype=np.float64)
    genes_raw = []
    for i, row in enumerate(rows):
        if not row or not str(row[0]).strip():
            continue
        genes_raw.append(str(row[0]).strip())
        cells = row[1:1 + n_table_samples]
        for j, cell in enumerate(cells):
            text = str(cell).strip().strip('"')
            if text in ("", "NA", "null", "NaN", "nan", "-"):
                continue
            try:
                values[i, j] = float(text)
            except Exception:  # noqa: BLE001
                values[i, j] = np.nan
    keep = np.array([bool(str(row[0]).strip()) for row in rows], dtype=bool)
    if keep.sum() != len(rows):
        values = values[keep]
        warnings.append(f"丢弃 {int((~keep).sum())} 个 ID 为空的行")
    if np.isnan(values).all():
        raise _err("bad_request", f"GEO series matrix 表格里没有可解析的数值：{Path(path).name}")
    if np.isnan(values).any():
        warnings.append("表格里有缺失/非数值（NA 等）按 0 填充")
    values = np.nan_to_num(values, nan=0.0)

    # ---- 解析 !Sample_* 元数据 ----
    entries: list[tuple[str, list[str]]] = []
    for line in sample_lines:
        try:
            parts = next(csv.reader([line], delimiter="\t", quotechar='"'))
        except Exception:  # noqa: BLE001
            continue
        if len(parts) < 2:
            continue
        key = parts[0].strip()
        if key.startswith("!Sample_"):
            key = key[len("!Sample_"):]
        vals = [str(v).strip().strip('"') for v in parts[1:]]
        if len(vals) < n_table_samples:
            vals = vals + [""] * (n_table_samples - len(vals))
        entries.append((key, vals[:n_table_samples]))

    accessions: list[str] = []
    titles: list[str] = []
    info_data: dict[str, list] = {}
    used = set()

    def _unique_name(base: str) -> str:
        name = base
        counter = 1
        while name in used:
            counter += 1
            name = f"{base}_{counter}"
        used.add(name)
        return name

    for key, vals in entries:
        if key == "geo_accession":
            accessions = vals
            continue
        if key == "title":
            titles = vals
            continue
        if key == "characteristics_ch1":
            prefix = None
            if vals and ":" in vals[0]:
                prefix = vals[0].split(":", 1)[0].strip()
            if prefix and all((":" in v) or v == "" for v in vals):
                column = _unique_name(prefix)
                info_data[column] = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]
            else:
                column = _unique_name(key)
                info_data[column] = vals
            continue
        if key in ("supplementary_file", "relation", "sample_attributes"):
            info_data[_unique_name(key)] = vals
            continue
        info_data[_unique_name(key)] = vals

    # 样本名：优先表格表头（与矩阵列严格对齐），geo_accession 保留为列
    samples = table_samples
    info_index = list(samples)
    if accessions and len(set(accessions)) == len(accessions):
        if titles and len(set(titles)) == len(titles) and set(accessions) == set(table_samples):
            # 表头是 GSM 号时保持 accession，title 作为列
            pass
        info_data.setdefault("geo_accession", accessions)
        if titles:
            info_data.setdefault("title", titles)
    elif titles:
        if len(set(titles)) == len(titles):
            # 表头与 title 不一致时，仍用表头（矩阵列）做样本名，title 入列
            info_data.setdefault("title", titles)
        else:
            warnings.append("!Sample_title 有重复，样本名改用矩阵表头")
    sample_info = pd.DataFrame(index=pd.Index(info_index, name="sample"))
    for column, vals in info_data.items():
        if len(vals) == len(samples):
            sample_info[column] = vals
        elif len(vals) > len(samples):
            sample_info[column] = vals[:len(samples)]
        else:
            sample_info[column] = vals + [""] * (len(samples) - len(vals))
    notes.append("样本名取自 series matrix 表头（通常 = GEO accession），!Sample_title 保留在 title 列")
    payload = _payload("geo_series", source, values.T, genes_raw, samples, sample_info=sample_info,
                       notes=notes, warnings=warnings,
                       meta={"series": {"n_samples": len(samples), "n_genes": len(genes_raw)}})
    return payload


def _load_text(path, options) -> dict:
    options = options or {}
    head = _read_head(path, 2 * 1024 * 1024)
    if not head.strip():
        raise _err("bad_request", f"文件是空的：{Path(path).name}")
    if "!series_matrix_table_begin" in head:
        return _parse_geo_series(path, path)
    notes: list[str] = []
    warnings: list[str] = []
    frame = _frame_from_text(path, options, notes, warnings)
    return _finalize_matrix(frame, "matrix_text", path, options, notes, warnings)


# ---------------------------------------------------------------------------
# xlsx
# ---------------------------------------------------------------------------


def _open_workbook(path):
    try:
        import openpyxl
    except Exception:  # noqa: BLE001
        raise _err("unsupported", "缺少 openpyxl，无法读取 Excel 文件；请另存为 CSV/TSV 后导入")
    try:
        return openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", f"无法读取 Excel 文件（{Path(path).name}）：{error}")


def _xlsx_frame(path, options) -> tuple[pd.DataFrame, list[str]]:
    options = options or {}
    workbook = _open_workbook(path)
    try:
        sheet_names = list(workbook.sheetnames)
        requested = options.get("sheet")
        candidates = []
        if requested is not None:
            if isinstance(requested, int) or (isinstance(requested, str) and requested.isdigit()):
                index = int(requested)
                if not (0 <= index < len(sheet_names)):
                    raise _err("not_found", f"sheet 序号越界：{index}（共 {len(sheet_names)} 个 sheet）")
                sheet_names = [sheet_names[index]]
            elif requested not in sheet_names:
                raise _err("not_found", f"工作表不存在：{requested}（可用：{', '.join(sheet_names[:10])}）")
            else:
                sheet_names = [requested]
        for name in sheet_names:
            worksheet = workbook[name]
            rows = []
            for row in worksheet.iter_rows(values_only=True):
                if row is None:
                    continue
                if all(cell is None or str(cell).strip() == "" for cell in row):
                    continue
                rows.append(["" if cell is None else cell for cell in row])
            if len(rows) < 2 or len(rows[0]) < 2:
                continue
            candidates.append((name, rows))
        if not candidates:
            raise _err("bad_request", f"Excel 文件里没有可用的二维表：{Path(path).name}")
        if options.get("sheet") is None and len(candidates) > 1:
            # 选择第一个「数值列 ≥2」的 sheet
            chosen = candidates[0]
            for name, rows in candidates:
                numeric_cols = 0
                for column in range(1, min(len(rows[0]), 6)):
                    numbers = 0
                    for row in rows[1:6]:
                        if column < len(row) and _is_numeric_text(row[column]):
                            numbers += 1
                    if numbers >= 2:
                        numeric_cols += 1
                if numeric_cols >= 2:
                    chosen = (name, rows)
                    break
            chosen_name = chosen[0]
        else:
            chosen = candidates[0]
            chosen_name = chosen[0]
        rows = chosen[1]
        header = _detect_header(rows[0])
        first = [str(cell) for cell in rows[0]]
        data_rows = rows[1:] if header else rows
        has_index = False
        if data_rows:
            non_numeric = sum(1 for row in data_rows[:20]
                              if row and not _is_numeric_text(row[0]))
            has_index = non_numeric >= max(1, int(0.6 * min(len(data_rows), 20)))
        first_header = str(first[0]).strip().lower() if first else ""
        if not has_index and header and first_header in (
                "", "gene", "gene_id", "geneid", "id", "probe", "probe_id", "probeid",
                "transcript", "transcript_id", "symbol", "ensembl", "ensembl_id", "features", "feature"):
            has_index = True
        width = max(len(row) for row in rows)
        columns = _dedupe(first[:width], "col") if header else [f"col_{i + 1}" for i in range(width)]
        body = data_rows
        padded = []
        for row in body:
            values = list(row) + [""] * (width - len(row))
            padded.append(values[:width])
        frame = pd.DataFrame(padded, columns=columns)
        if has_index:
            frame = frame.set_index(frame.columns[0])
        notes = [f"Excel：{Path(path).name}，读取 sheet「{chosen_name}」"]
        if options.get("sheet") is None and len(candidates) > 1:
            notes.append(f"文件共 {len(candidates)} 个可用 sheet（{', '.join(n for n, _ in candidates[:5])}），已选第一个；可用 sheet 参数指定")
        return frame, notes
    finally:
        try:
            workbook.close()
        except Exception:  # noqa: BLE001
            pass


def _load_xlsx(path, options) -> dict:
    options = options or {}
    frame, notes = _xlsx_frame(path, options)
    warnings: list[str] = []
    return _finalize_matrix(frame, "xlsx", path, options, notes, warnings)


# ---------------------------------------------------------------------------
# 压缩包
# ---------------------------------------------------------------------------


def _classify_member(name: str) -> str | None:
    lower = str(name).lower()
    if any(lower.endswith(suffix) for suffix in ARCHIVE_SUFFIXES):
        return None
    while lower.endswith((".gz", ".bz2")):
        lower = lower[: -3]
    if lower.endswith(XLSX_SUFFIXES):
        return "xlsx"
    if lower.endswith(R_SUFFIXES):
        return "rdata" if lower.endswith((".rdata", ".rda")) else "rds"
    if lower.endswith((".csv", ".tsv", ".txt", ".tab", ".matrix", ".counts")):
        if "series_matrix" in lower or lower.endswith(("series_matrix.txt", "series_matrix.txt")):
            return "geo_series"
        return "matrix_text"
    return None


def _extract_archive(path, dirs, progress=None) -> Path:
    path = Path(path)
    name = path.name.lower()
    root = Path(tempfile.mkdtemp(prefix="unpack-", dir=str(dirs["unpack"])))
    total = 0
    _emit(progress, f"正在解压 {path.name}", 10)
    try:
        if name.endswith(".zip"):
            with zipfile.ZipFile(path) as archive:
                for info in archive.infolist():
                    if info.is_dir():
                        continue
                    member = info.filename
                    if not _is_safe_member(member):
                        continue
                    total += int(info.file_size)
                    if total > MAX_UNPACK_BYTES:
                        raise _err("bad_request", "压缩包解压后超过 20GB 上限，请先裁剪后重新打包")
                    target = _safe_member_path(root, member)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(info) as source, open(target, "wb") as sink:
                        shutil.copyfileobj(source, sink)
        elif name.endswith(ARCHIVE_SUFFIXES) or tarfile.is_tarfile(path):
            with tarfile.open(path, "r:*") as archive:
                for member in archive.getmembers():
                    if not member.isfile():
                        continue
                    if not _is_safe_member(member.name):
                        continue
                    total += int(member.size)
                    if total > MAX_UNPACK_BYTES:
                        raise _err("bad_request", "压缩包解压后超过 20GB 上限，请先裁剪后重新打包")
                    target = _safe_member_path(root, member.name)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    handle = archive.extractfile(member)
                    if handle is None:
                        continue
                    with handle as source, open(target, "wb") as sink:
                        shutil.copyfileobj(source, sink)
        else:
            raise _err("unsupported", f"不是支持的压缩格式：{path.name}")
    except EngineError:
        raise
    except (zipfile.BadZipFile, tarfile.TarError) as error:
        raise _err("bad_request", f"压缩包损坏或格式不识别（{path.name}）：{error}")
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", f"解压失败（{path.name}）：{error}")
    _emit(progress, f"解压完成：{total / 1024 ** 2:.1f} MB", 60)
    return root


def _is_safe_member(name: str) -> bool:
    text = str(name or "").replace("\\", "/")
    if not text or text.startswith("/"):
        return False
    parts = text.split("/")
    return ".." not in parts


def _safe_member_path(root, member) -> Path:
    root_resolved = Path(root).resolve()
    target = (root_resolved / str(member)).resolve()
    if target != root_resolved and root_resolved not in target.parents:
        raise _err("bad_request", f"压缩包含越权路径：{member}")
    return target


def _archive_candidates(root: Path) -> list[dict]:
    candidates = []
    for file_path in sorted(Path(root).rglob("*")):
        if not file_path.is_file() or file_path.name.startswith("."):
            continue
        kind = _classify_member(file_path.name)
        if kind is None:
            continue
        size = int(file_path.stat().st_size)
        if size < 8:
            continue
        relative = str(file_path.relative_to(root))
        priority = 0
        if kind == "geo_series" or "series_matrix" in file_path.name.lower():
            priority = 3
        elif any(token in file_path.name.lower() for token in ("matrix", "counts", "tpm", "fpkm",
                                                               "expression", "expr")):
            priority = 2
        candidates.append({"name": relative, "member": relative, "size": size,
                           "kind": kind, "format": kind, "priority": priority})
    candidates.sort(key=lambda item: (-item["priority"], item["name"]))
    for item in candidates:
        item.pop("priority", None)
    # 文本候选再看一眼文件头：可能是 GEO series matrix（扩展名是 .txt）
    for item in candidates:
        if item["kind"] != "matrix_text":
            continue
        try:
            head = _read_head(Path(root) / item["name"], 256 * 1024)
            if "!series_matrix_table_begin" in head:
                item["kind"] = item["format"] = "geo_series"
        except Exception:  # noqa: BLE001
            continue
    return candidates


def _member_matches(candidate_name: str, member: str) -> bool:
    if candidate_name == member:
        return True
    if Path(candidate_name).name == member:
        return True
    return candidate_name.replace("\\", "/").endswith(member.replace("\\", "/"))


def _load_archive(path, dirs, options, progress=None) -> dict:
    options = dict(options or {})
    root = _extract_archive(path, dirs, progress=progress)
    candidates = _archive_candidates(root)
    if not candidates:
        raise _err("bad_request", f"压缩包里没有找到可识别的表达矩阵：{Path(path).name}")
    member = options.get("member")
    if member:
        target = next((item for item in candidates if _member_matches(item["name"], str(member))), None)
        if target is None:
            raise _err("not_found", f"压缩包里没有 {member}；可用文件：{', '.join(item['name'] for item in candidates[:10])}")
        selected = target
    elif len(candidates) == 1:
        selected = candidates[0]
    else:
        raise _err("ambiguous", f"压缩包里找到 {len(candidates)} 个候选数据文件，请用 member 指定",
                   candidates=[{"name": item["name"], "size": item["size"], "kind": item["kind"]}
                               for item in candidates])
    _emit(progress, f"正在解析压缩包内文件：{selected['name']}", 75)
    inner_path = root / selected["name"]
    result = _load_path(inner_path, dirs, {**options, "data_format": selected["kind"]}, progress)
    result["notes"] = list(result.get("notes") or []) + [
        f"来自压缩包 {Path(path).name}（已解压 {len(candidates)} 个候选文件）"
    ]
    result["format"] = selected["kind"]
    return result


# ---------------------------------------------------------------------------
# R 对象（rds / RData）—— bg_io.R
# ---------------------------------------------------------------------------


def _rscript() -> str | None:
    candidates = [
        os.environ.get("DS_BG_RSCRIPT") or "",
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
        "/usr/local/bin/Rscript",
        "/opt/homebrew/bin/Rscript",
        "/usr/bin/Rscript",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    found = shutil.which("Rscript")
    return found or None


def _bg_io_r_script() -> Path:
    return Path(__file__).resolve().parent / "bg_io.R"


def _decompress_r(path, dirs) -> Path:
    name = Path(path).name.lower()
    if not name.endswith(".gz"):
        return Path(path)
    target = Path(tempfile.mkdtemp(prefix="rds-", dir=str(dirs["tmp"]))) / Path(path).name[:-3]
    with gzip.open(path, "rb") as source, open(target, "wb") as sink:
        shutil.copyfileobj(source, sink)
    return target


def _load_r(path, dirs, options) -> dict:
    options = dict(options or {})
    rscript = _rscript()
    if not rscript:
        raise _err("unsupported", "本机没有找到 Rscript，无法读取 .rds/.RData；请用 R 另存为 CSV/TSV 后再导入，"
                                  "或安装 R（需要 Biobase/SummarizedExperiment）")
    script = _bg_io_r_script()
    if not script.exists():
        raise _err("internal", f"缺少 R 桥脚本：{script}")
    work = Path(tempfile.mkdtemp(prefix="r-read-", dir=str(dirs["tmp"])))
    command = [rscript, "--vanilla", str(script), "--in", str(Path(path).resolve()), "--out", str(work)]
    if options.get("object"):
        command += ["--object", str(options["object"])]
    if options.get("inspect"):
        command += ["--inspect"]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=R_TIMEOUT)
    except subprocess.TimeoutExpired:
        raise _err("bad_request", f"R 读取超时（>{R_TIMEOUT} 秒）：{Path(path).name}")
    except Exception as error:  # noqa: BLE001
        raise _err("internal", f"无法调用 Rscript：{error}")
    meta_path = work / "meta.tsv"
    meta = _read_meta_tsv(meta_path)
    status = meta.get("status", "error")
    if completed.returncode != 0 or status == "error":
        detail = (meta.get("error") or "").strip() or (completed.stderr or "")[-1200:] or "未知错误"
        raise _err("bad_request", f"R 读取失败（{Path(path).name}）：{detail}")
    if status == "ambiguous":
        candidates = _read_candidates_tsv(work / "candidates.tsv")
        names = "、".join(str(item.get("object")) for item in candidates[:8])
        raise _err("ambiguous", f"R 文件里有多个可用对象，请用 object 参数指定：{names}", candidates=candidates)
    n_genes = int(float(meta.get("n_genes") or 0))
    n_samples = int(float(meta.get("n_samples") or 0))
    r_class = meta.get("class") or "R object"
    notes = [f"R 对象：{r_class}"]
    if meta.get("assay"):
        notes.append(f"assay：{meta['assay']}")
    if meta.get("notes"):
        notes.extend(str(meta["notes"]).split(" | "))
    warnings = [str(item) for item in (meta.get("warnings") or "").split(" | ") if str(item).strip()]
    if options.get("inspect"):
        return {"format": "rdata" if Path(path).name.lower().endswith((".rdata", ".rda")) else "rds",
                "source": str(path), "class": r_class, "n_samples": n_samples, "n_genes": n_genes,
                "orientation": meta.get("orientation"), "notes": notes, "warnings": warnings,
                "candidates": None}
    matrix_path = work / "matrix.mtx"
    genes_path = work / "genes.tsv"
    samples_path = work / "samples.tsv"
    if not matrix_path.exists() or not genes_path.exists():
        raise _err("internal", f"R 桥没有写出矩阵（{Path(path).name}）")
    try:
        from scipy.io import mmread

        matrix = mmread(str(matrix_path))
        matrix = np.asarray(matrix.todense() if hasattr(matrix, "todense") else matrix, dtype=np.float64)
    except Exception as error:  # noqa: BLE001
        raise _err("internal", f"矩阵文件读取失败：{error}")
    genes = _read_lines(genes_path)
    samples = _read_lines(samples_path) if samples_path.exists() else None
    orientation = str(meta.get("orientation") or "genes_by_samples")
    if orientation == "genes_by_samples":
        matrix = matrix.T
    if samples is None:
        samples = [f"S{i + 1}" for i in range(matrix.shape[0])]
    sample_info = None
    obs_path = work / "obs.csv"
    if obs_path.exists():
        try:
            sample_info = pd.read_csv(obs_path, index_col=0, low_memory=False)
        except Exception:  # noqa: BLE001
            sample_info = None
    fmt = "rdata" if Path(path).name.lower().endswith((".rdata", ".rda")) else "rds"
    payload = _payload(fmt, path, matrix, genes, samples, sample_info=sample_info,
                       notes=notes, warnings=warnings,
                       meta={"r_class": r_class, "orientation": orientation,
                             "image": str(matrix_path)})
    if meta.get("has_counts") and str(meta["has_counts"]).upper() == "TRUE":
        payload["meta"]["has_counts"] = True
    return payload


def _read_meta_tsv(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if not line or "\t" not in line:
                    continue
                key, value = line.split("\t", 1)
                out[key.strip()] = value.strip()
    except Exception:  # noqa: BLE001
        return out
    return out


def _read_candidates_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        frame = pd.read_csv(path, sep="\t")
    except Exception:  # noqa: BLE001
        return []
    candidates = []
    for _, row in frame.iterrows():
        candidates.append({
            "object": str(row.get("object", "")),
            "name": str(row.get("object", "")),
            "class": str(row.get("class", "")),
            "label": str(row.get("label", row.get("object", ""))),
            "n_obs": _to_int(row.get("n_samples")),
            "n_vars": _to_int(row.get("n_genes")),
        })
    return candidates


def _to_int(value):
    try:
        return int(float(value))
    except Exception:  # noqa: BLE001
        return None


def _read_lines(path: Path) -> list[str]:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return [line.rstrip("\n").rstrip("\r") for line in handle if line.strip()]
    except Exception:  # noqa: BLE001
        return []


# ---------------------------------------------------------------------------
# 非编码模块数据（circ / mirna / lncrna / fusion / splicing）
# ---------------------------------------------------------------------------

_MIRNA_PRECURSOR_MAP = {
    "hsa-mir-21": "hsa-miR-21-5p", "hsa-mir-155": "hsa-miR-155-5p",
    "hsa-mir-let-7a-1": "hsa-let-7a-5p", "hsa-mir-let-7b": "hsa-let-7b-5p",
    "hsa-mir-10b": "hsa-miR-10b-5p", "hsa-mir-17": "hsa-miR-17-5p",
    "hsa-mir-18a": "hsa-miR-18a-5p", "hsa-mir-19a": "hsa-miR-19a-3p",
    "hsa-mir-20a": "hsa-miR-20a-5p", "hsa-mir-25": "hsa-miR-25-3p",
    "hsa-mir-26a-1": "hsa-miR-26a-5p", "hsa-mir-27a": "hsa-miR-27a-3p",
    "hsa-mir-29a": "hsa-miR-29a-3p", "hsa-mir-29b-1": "hsa-miR-29b-3p",
    "hsa-mir-30a": "hsa-miR-30a-5p", "hsa-mir-31": "hsa-miR-31-5p",
    "hsa-mir-34a": "hsa-miR-34a-5p", "hsa-mir-92a-1": "hsa-miR-92a-3p",
    "hsa-mir-93": "hsa-miR-93-5p", "hsa-mir-101-1": "hsa-miR-101-3p",
    "hsa-mir-125a": "hsa-miR-125a-5p", "hsa-mir-125b-1": "hsa-miR-125b-5p",
    "hsa-mir-126": "hsa-miR-126-3p", "hsa-mir-133a-1": "hsa-miR-133a-3p",
    "hsa-mir-141": "hsa-miR-141-3p", "hsa-mir-143": "hsa-miR-143-3p",
    "hsa-mir-145": "hsa-miR-145-5p", "hsa-mir-146a": "hsa-miR-146a-5p",
    "hsa-mir-148a": "hsa-miR-148a-3p", "hsa-mir-181a": "hsa-miR-181a-5p",
    "hsa-mir-191": "hsa-miR-191-5p", "hsa-mir-200b": "hsa-miR-200b-3p",
    "hsa-mir-200c": "hsa-miR-200c-3p", "hsa-mir-203": "hsa-miR-203a-3p",
    "hsa-mir-205": "hsa-miR-205-5p", "hsa-mir-210": "hsa-miR-210-3p",
    "hsa-mir-214": "hsa-miR-214-3p", "hsa-mir-221": "hsa-miR-221-3p",
    "hsa-mir-222": "hsa-miR-222-3p", "hsa-mir-375": "hsa-miR-375-3p",
    "hsa-mir-451a": "hsa-miR-451a", "hsa-mir-486-1": "hsa-miR-486-5p",
    "hsa-mir-21-5p": "hsa-miR-21-5p", "hsa-mir-155-5p": "hsa-miR-155-5p",
}

_CIRC_ID_RE = re.compile(r"^(?:chr[0-9XYM]+[:：]\d+[-|]\d+[:：][+-]|hsa_circ_\d+|circ_\d+|"
                         r"chr[0-9XYM]+_\d+_\d+_[+-])")
_LNCRNA_ID_RE = re.compile(r"(?:^LINC\d+|^.*-AS\d*$|^.*-IT\d*$|^ENS[A-Z]*G\d{6,})", re.I)


def _read_table_any(path, options=None) -> pd.DataFrame:
    """通用表格（保留字符串列，不做数值化）。"""
    path = Path(path)
    options = options or {}
    head = _read_head(path, 2 * 1024 * 1024)
    if not head.strip():
        raise _err("bad_request", f"文件是空的：{Path(path).name}")
    if "!series_matrix_table_begin" in head:
        raise _err("bad_request", f"这是 GEO series matrix（{Path(path).name}），请用 data_format=geo_series 导入")
    delimiter = options.get("delimiter") or _sniff_delimiter(head.splitlines())
    try:
        if delimiter == "whitespace":
            with _open_text(path) as handle:
                return pd.read_csv(handle, sep=r"\s+", engine="python", low_memory=False)
        return pd.read_csv(path, sep=delimiter, low_memory=False)
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", f"表格读取失败（{Path(path).name}）：{error}")


def _first_non_numeric_column(frame: pd.DataFrame):
    for column in frame.columns:
        series = frame[column].astype(str)
        numeric_ratio = sum(1 for value in series if _is_numeric_text(value)) / max(1, len(series))
        if numeric_ratio < 0.3:
            return column
    return None


def _feature_matrix(frame: pd.DataFrame, id_candidates) -> tuple[np.ndarray, list[str], list[str], pd.DataFrame]:
    """从宽表提取 samples×features 矩阵与特征注释。"""
    lower_map = {str(column).strip().lower(): column for column in frame.columns}
    id_column = None
    for candidate in id_candidates:
        if candidate.lower() in lower_map:
            id_column = lower_map[candidate.lower()]
            break
    if id_column is None:
        first = frame.columns[0]
        numeric_ratio = sum(1 for value in frame[first].astype(str) if _is_numeric_text(value)) / max(1, len(frame))
        if numeric_ratio < 0.5:
            id_column = first
    if id_column is None:
        raise _err("bad_request", "表格里没有找到特征 ID 列（circRNA/miRNA/lncRNA 名称列）")
    numeric_columns = []
    for column in frame.columns:
        if column == id_column:
            continue
        series = frame[column]
        numeric_ratio = sum(1 for value in series.astype(str) if _is_numeric_text(value)) / max(1, len(series))
        if numeric_ratio >= 0.8:
            numeric_columns.append(column)
    if not numeric_columns:
        raise _err("bad_request", "表格里没有找到样本表达列（数值列）")
    features = _dedupe([str(value) for value in frame[id_column]], "feature")
    info_columns = [column for column in frame.columns if column != id_column and column not in numeric_columns]
    info = frame[info_columns].copy()
    info.index = pd.Index(features, name="feature")
    values = frame[numeric_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
    samples = _dedupe([str(column) for column in numeric_columns], "sample")
    return values, features, samples, info


def _pivot_long_table(frame: pd.DataFrame) -> pd.DataFrame | None:
    """长表（feature, sample, value）→ 宽表。"""
    lower_map = {str(column).strip().lower(): column for column in frame.columns}
    feature_col = next((lower_map[key] for key in
                        ("feature", "circ_id", "circrna", "mirna", "lncrna", "gene", "id", "junction")
                        if key in lower_map), None)
    sample_col = next((lower_map[key] for key in ("sample", "sample_id", "samples", "gsm")
                       if key in lower_map), None)
    value_col = next((lower_map[key] for key in ("value", "count", "counts", "expression", "expr", "tpm", "rpkm")
                      if key in lower_map), None)
    if not (feature_col and sample_col and value_col):
        return None
    try:
        return frame.pivot_table(index=feature_col, columns=sample_col, values=value_col, aggfunc="sum")
    except Exception:  # noqa: BLE001
        return None


def _load_feature_kind(path, kind, options) -> dict:
    options = options or {}
    frame = _read_table_any(path, options)
    if frame.empty:
        raise _err("bad_request", f"{kind} 表格为空：{Path(path).name}")
    pivoted = _pivot_long_table(frame)
    if pivoted is not None:
        frame = pivoted.reset_index()
    id_candidates = {
        "circ": ["circ_id", "circRNA", "circrna", "circ", "junction", "id", "feature", "gene"],
        "mirna": ["mirna", "mirna_id", "mature_id", "name", "id", "feature", "gene"],
        "lncrna": ["lncrna", "lncrna_id", "gene", "gene_id", "id", "feature", "name"],
    }[kind]
    values, features, samples, info = _feature_matrix(frame, id_candidates)
    notes = [f"{kind} 表达矩阵：{values.shape[1]} 个样本 × {len(features)} 个特征"]
    warnings: list[str] = []
    if kind == "circ":
        recognized = sum(1 for feature in features if _CIRC_ID_RE.match(feature))
        notes.append(f"circRNA ID 识别：{recognized}/{len(features)} 形如 chr:start-end:strand 或 hsa_circ_0001234")
        if recognized < len(features) * 0.5:
            warnings.append("超过一半的特征名不像标准 circRNA ID，请确认表格是否选对")
    elif kind == "mirna":
        mapped = 0
        renamed = []
        for feature in features:
            key = feature.lower().replace("_", "-")
            if key in _MIRNA_PRECURSOR_MAP:
                renamed.append(_MIRNA_PRECURSOR_MAP[key])
                mapped += 1
            else:
                renamed.append(feature)
        features = _dedupe(renamed, "mirna")
        notes.append(f"miRNA precursor→mature 映射：{mapped} 个（内置 {len(_MIRNA_PRECURSOR_MAP)} 条）")
        recognized = sum(1 for feature in features if "mir" in feature.lower())
        if recognized < len(features) * 0.5:
            warnings.append("超过一半的特征名不像 miRNA 名称（hsa-miR-…）")
    else:  # lncrna
        recognized = sum(1 for feature in features if _LNCRNA_ID_RE.search(feature))
        notes.append(f"lncRNA 识别（LINC… / …-AS1 / …-IT1 / ENSG）：{recognized}/{len(features)}")
        info = info.copy()
        info["lncrna_class"] = [
            ("lincRNA" if str(feature).upper().startswith("LINC") else
             "antisense" if str(feature).upper().endswith(("-AS1", "-AS2", "-AS")) else
             "intronic" if str(feature).upper().endswith(("-IT1", "-IT2")) else
             "ensembl" if str(feature).upper().startswith("ENS") else "other")
            for feature in features
        ]
    matrix = values.T  # samples × features
    sample_info = pd.DataFrame(index=pd.Index(samples, name="sample"))
    payload = _payload(kind, path, matrix, features, samples, sample_info=sample_info,
                       notes=notes, warnings=warnings,
                       meta={"kind": kind, "feature_kind": kind, "has_counts": np.all(np.abs(values - np.round(values)) < 1e-6)})
    payload["info"] = info
    return payload


def _load_fusion(path, options) -> dict:
    options = options or {}
    frame = _read_table_any(path, options)
    lower_map = {str(column).strip().lower(): column for column in frame.columns}

    def pick(*names):
        for name in names:
            if name.lower() in lower_map:
                return lower_map[name.lower()]
        return None

    fusion_col = pick("#fusionname", "fusionname", "fusion_name", "fusion", "name", "fusions")
    gene1_col = pick("leftgene", "gene1", "genea", "left_gene", "gene_1", "5gene", "five_gene")
    gene2_col = pick("rightgene", "gene2", "geneb", "right_gene", "gene_2", "3gene", "three_gene")
    junction_col = pick("junctionreadcount", "junctionreadscount", "junction_reads", "junctionreads")
    spanning_col = pick("spanningfragcount", "spanningfragments", "spanning_reads", "spanningreads")
    breakpoint_col = pick("breakpoint", "junction", "chromosome", "chr")
    frame_name = Path(path).name
    if fusion_col is None and gene1_col is None:
        raise _err("bad_request", f"没识别出融合结果表列（{frame_name}）；需要 FusionName 或 gene1/gene2 列")
    def _split_fusion(text: str) -> tuple[str, str]:
        raw = str(text).strip()
        for token in ("^", "--", "-"):
            if token in raw:
                left, right = raw.split(token, 1)
                return left.strip(), right.strip()
        return raw, ""

    pairs = []
    for _, row in frame.iterrows():
        if fusion_col is not None:
            genes = _split_fusion(row[fusion_col])
        else:
            genes = (str(row[gene1_col]).strip(), str(row[gene2_col]).strip())
        pairs.append({
            "fusion": f"{genes[0]}--{genes[1]}",
            "gene1": genes[0], "gene2": genes[1],
            "junction_reads": _to_int(row[junction_col]) if junction_col else None,
            "spanning_reads": _to_int(row[spanning_col]) if spanning_col else None,
            "breakpoint": str(row[breakpoint_col]) if breakpoint_col else "",
        })
    info = pd.DataFrame(pairs)
    info.index = pd.Index(info["fusion"], name="fusion")
    # 样本：优先找样本名列，否则按单样本处理
    sample_col = pick("sample", "sample_id", "samples")
    samples = _dedupe([str(value) for value in frame[sample_col]], "sample") if sample_col else [Path(path).stem]
    notes = [f"融合结果表：{len(info)} 个融合事件；识别列 "
             f"{'FusionName' if fusion_col else 'gene1/gene2'}"
             f"{'，junction reads' if junction_col else ''}{'，spanning reads' if spanning_col else ''}"]
    warnings = []
    if not sample_col:
        warnings.append("表里没有样本列，按单样本文本处理；如需多样本，请提供 wide 矩阵或加 sample 列")
    payload = _payload("fusion", path, None, [], samples, sample_info=None,
                       notes=notes, warnings=warnings,
                       meta={"kind": "fusion", "n_events": int(len(info)), "n_samples": len(samples), "n_genes": 0})
    payload["samples"] = samples
    payload["info"] = info
    payload["fusion_table"] = frame
    return payload


_RMATS_COLUMNS = ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2"]


def _parse_comma_numbers(text) -> list[float]:
    values = []
    for token in str(text).replace(";", ",").split(","):
        token = token.strip()
        if token == "":
            continue
        try:
            values.append(float(token))
        except Exception:  # noqa: BLE001
            continue
    return values


def _load_splicing(path, options) -> dict:
    options = options or {}
    frame = _read_table_any(path, options)
    columns = [str(column) for column in frame.columns]
    upper = {column.upper(): column for column in columns}
    # rMATS 列名可能带下划线（IJC_SAMPLE_1）或写成驼峰/空格，统一去掉非字母数字后比较
    normalized = {re.sub(r"[^A-Z0-9]", "", column.upper()): column for column in columns}
    rmats_keys = ("IJCSAMPLE1", "SJCSAMPLE1", "INCLEVEL1", "IJCSAMPLE2", "INCLEVEL2", "INCLEVELDIFFERENCE")
    is_rmats = any(key in normalized for key in rmats_keys)
    event_type = "SE"
    match = re.search(r"(SE|A5SS|A3SS|MXE|RI)", Path(path).name.upper())
    if match:
        event_type = match.group(1)
    notes = []
    warnings = []
    if is_rmats:
        def col(name):
            return normalized.get(name)

        inc1 = col("INCLEVEL1")
        inc2 = col("INCLEVEL2")
        inc_diff = col("INCLEVELDIFFERENCE")
        if inc1 is None and inc_diff is None:
            raise _err("bad_request", f"rMATS 表缺少 IncLevel1 / IncLevelDifference 列：{Path(path).name}")
        diff_values = None
        if inc_diff is not None:
            parsed_diff = [_parse_comma_numbers(value) for value in frame[inc_diff]]
            diff_values = np.array([np.mean(item) if item else np.nan for item in parsed_diff],
                                   dtype=np.float64)
        if inc1 is None:
            # 只有 IncLevelDifference 的差异剪接表：矩阵退化为单行 dPSI，仍可挂载与分析
            features = []
            for i in range(len(frame)):
                if "ID" in upper:
                    features.append(str(frame[upper["ID"]].iloc[i]))
                elif "GENE" in upper:
                    features.append(str(frame[upper["GENE"]].iloc[i]))
                else:
                    features.append(f"{event_type}_event_{i + 1}")
            samples = ["dPSI"]
            matrix = np.nan_to_num((diff_values or np.zeros(len(frame)))[None, :], nan=0.0)
            info = pd.DataFrame({"event_type": [event_type] * len(features),
                                 "ID": features,
                                 "dPSI": diff_values if diff_values is not None else np.nan})
            info.index = pd.Index(_dedupe(features, "event"), name="event")
            notes.append(f"rMATS 差异剪接表（{event_type}）：{len(frame)} 个事件（仅 IncLevelDifference）")
            payload = _payload("splicing", path, matrix, info.index.tolist(), samples,
                               sample_info=pd.DataFrame(index=pd.Index(samples, name="sample")),
                               notes=notes, warnings=warnings,
                               meta={"kind": "splicing", "event_type": event_type, "n_events": len(frame)})
            payload["info"] = info
            return payload
        levels1 = [_parse_comma_numbers(value) for value in frame[inc1]]
        levels2 = [_parse_comma_numbers(value) for value in frame[inc2]] if inc2 is not None else []
        n1 = max((len(item) for item in levels1), default=1) or 1
        n2 = max((len(item) for item in levels2), default=0)
        samples = [f"group1_{i + 1}" for i in range(n1)] + [f"group2_{i + 1}" for i in range(n2)]
        matrix = np.full((len(samples), len(frame)), np.nan, dtype=np.float64)
        for i, item in enumerate(levels1):
            for j, value in enumerate(item[:n1]):
                matrix[j, i] = value
        for i, item in enumerate(levels2):
            for j, value in enumerate(item[:n2]):
                matrix[n1 + j, i] = value
        matrix = np.nan_to_num(matrix, nan=0.0)
        sample_info = pd.DataFrame(index=pd.Index(samples, name="sample"))
        if n2:
            sample_info["group"] = ["group1"] * n1 + ["group2"] * n2
        features = []
        for i in range(len(frame)):
            if "ID" in upper:
                features.append(str(frame[upper["ID"]].iloc[i]))
            elif "GENE" in upper:
                features.append(str(frame[upper["GENE"]].iloc[i]))
            else:
                features.append(f"{event_type}_event_{i + 1}")
        if n2:
            dpsi = matrix[:n1].mean(axis=0) - matrix[n1:].mean(axis=0)
        elif diff_values is not None:
            dpsi = diff_values
        else:
            dpsi = np.full(len(features), np.nan)
        info = pd.DataFrame({"event_type": [event_type] * len(features),
                             "ID": features,
                             "dPSI": dpsi})
        info.index = pd.Index(_dedupe(features, "event"), name="event")
        notes.append(f"rMATS 表（{event_type}）：{len(frame)} 个事件，group1 {n1} 例 / group2 {n2} 例")
        if inc_diff is not None and not n2:
            notes.append("dPSI 取自表内 IncLevelDifference 列")
        notes.append("样本名按 rMATS 列内逗号顺序自动编号（group1_1…）；如需真实样本名请提供 sample_sheet")
        payload = _payload("splicing", path, matrix, info.index.tolist(), samples, sample_info=sample_info,
                           notes=notes, warnings=warnings,
                           meta={"kind": "splicing", "event_type": event_type, "n_events": len(frame)})
        payload["info"] = info
        return payload
    # SUPPA / 通用 PSI 宽表
    payload = _load_feature_kind_like_matrix(frame, path, "splicing")
    payload["meta"]["kind"] = "splicing"
    payload["notes"].append("按通用 PSI 矩阵解析（特征 ID × 样本）")
    return payload


def _load_feature_kind_like_matrix(frame: pd.DataFrame, path, kind) -> dict:
    lower_map = {str(column).strip().lower(): column for column in frame.columns}
    id_col = next((lower_map[key] for key in ("id", "event_id", "feature", "name", "gene") if key in lower_map),
                  frame.columns[0])
    numeric_columns = []
    for column in frame.columns:
        if column == id_col:
            continue
        series = frame[column].astype(str)
        ratio = sum(1 for value in series if _is_numeric_text(value)) / max(1, len(series))
        if ratio >= 0.8:
            numeric_columns.append(column)
    if not numeric_columns:
        raise _err("bad_request", f"{kind} 表里没有数值样本列：{Path(path).name}")
    features = _dedupe([str(value) for value in frame[id_col]], "feature")
    values = frame[numeric_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
    samples = _dedupe([str(column) for column in numeric_columns], "sample")
    info = frame[[c for c in frame.columns if c != id_col and c not in numeric_columns]].copy()
    info.index = pd.Index(features, name="feature")
    payload = _payload(kind, path, values.T, features, samples, sample_info=None,
                       notes=[f"{kind} 宽表：{len(samples)} 个样本 × {len(features)} 个特征"],
                       warnings=[], meta={"kind": kind})
    payload["info"] = info
    return payload


def load_module(source, kind=None, workdir=None, options=None) -> dict:
    """挂载非编码/融合/剪接模块数据（attach 用）。

    返回结构（与 bg_registry.attach_module 对齐）：
      {"kind","format","X":ndarray|None,"genes":[str],"samples":[str],
       "table":DataFrame|None,"sample_info":DataFrame|None,"notes":"warnings"}
    矩阵型模态（circ/mirna/lncrna/splicing）给 X（样本×特征）+ genes（特征名）；
    表格型（fusion，以及可选的 splicing 事件注释）放 table。
    """
    options = dict(options or {})
    dirs = _workdirs(workdir)
    path, resolve_notes = _resolve_source(source, dirs, options)
    if kind is None:
        kind = options.get("kind") or options.get("data_format")
    if kind in (None, "auto"):
        kind = _guess_module_kind(path, options)
    if kind not in MODULE_KINDS:
        raise _err("bad_request", f"不支持的非编码模块类型：{kind}（可选：{', '.join(MODULE_KINDS)}）")
    if kind in ("circ", "mirna", "lncrna"):
        payload = _load_feature_kind(path, kind, options)
    elif kind == "fusion":
        payload = _load_fusion(path, options)
    else:
        payload = _load_splicing(path, options)
    X = payload.get("X")
    if X is not None:
        X = np.asarray(X, dtype=np.float32)
    genes = payload.get("genes") or []
    samples = payload.get("samples") or []
    sample_info = payload.get("sample_info")
    if not isinstance(sample_info, pd.DataFrame):
        sample_info = None
    table = None
    if kind == "fusion":
        table = payload.get("fusion_table")
    elif kind == "splicing":
        table = payload.get("info") if isinstance(payload.get("info"), pd.DataFrame) else None
    if not isinstance(table, pd.DataFrame):
        table = None
    return {
        "kind": kind,
        "format": payload.get("format", kind),
        "X": X,
        "genes": [str(g) for g in genes],
        "samples": [str(s) for s in samples],
        "table": table,
        "sample_info": sample_info,
        "notes": list(resolve_notes) + list(payload.get("notes") or []),
        "warnings": list(payload.get("warnings") or []),
    }


def _guess_module_kind(path, options) -> str:
    options = options or {}
    forced = options.get("data_format") or options.get("kind")
    if forced and forced != "auto":
        return str(forced)
    head = _read_head(path, 1024 * 1024)
    flat = " ".join(head.splitlines()[:5])
    if re.search(r"FusionName|LeftGene|RightGene|SpanningFragCount", head, re.I):
        return "fusion"
    if re.search(r"IJC_SAMPLE_1|SJC_SAMPLE_1|IncLevel1|IncLevel2", head):
        return "splicing"
    if re.search(r"hsa[-_]mir[-_]", head, re.I):
        return "mirna"
    if re.search(r"hsa_circ_|chr[0-9XYM]+:\d+-\d+:[+-]", head):
        return "circ"
    if re.search(r"\bLINC\d+|\b\S+-AS1\b", head):
        return "lncrna"
    raise _err("ambiguous", "无法自动判断非编码数据类型，请显式指定 kind（circ/mirna/lncrna/fusion/splicing）",
               candidates=[{"name": kind, "size": None, "kind": kind} for kind in MODULE_KINDS])


# ---------------------------------------------------------------------------
# 格式识别 / 源解析
# ---------------------------------------------------------------------------


def _strip_compression(name: str) -> str:
    lower = str(name).lower()
    while lower.endswith((".gz", ".bz2")):
        lower = lower[:-3]
    return lower


def _detect_format(path: Path, options=None) -> str:
    options = options or {}
    forced = options.get("data_format") or options.get("format")
    if forced and str(forced) not in ("auto", ""):
        return str(forced)
    name = path.name
    lower = name.lower()
    if any(lower.endswith(suffix) for suffix in ARCHIVE_SUFFIXES):
        return "archive"
    inner = _strip_compression(lower)
    if inner.endswith(XLSX_SUFFIXES):
        return "xlsx"
    if inner.endswith((".rdata", ".rda")):
        return "rdata"
    if inner.endswith(".rds"):
        return "rds"
    if inner.endswith((".csv", ".tsv", ".txt", ".tab", ".matrix", ".counts")):
        head = _read_head(path, 2 * 1024 * 1024)
        if "!series_matrix_table_begin" in head:
            return "geo_series"
        return "matrix_text"
    # 无扩展名：看内容
    if path.is_file():
        head = _read_head(path, 512 * 1024)
        if "!series_matrix_table_begin" in head:
            return "geo_series"
        if head.startswith(b"\x89HDF") or head.startswith("PK\x03\x04"):
            raise _err("unsupported", f"文件看起来是 HDF5/zip 容器，请显式指定 data_format：{path.name}")
        if "\t" in head or "," in head or ";" in head:
            return "matrix_text"
    raise _err("unsupported", f"无法识别文件格式：{path.name}；请用 data_format 参数显式指定")


def _download_to(url: str, dest_dir: Path, progress=None, filename: str | None = None) -> Path:
    """流式下载到本地（urllib，单文件 8GB 上限）。"""
    dest_dir.mkdir(parents=True, exist_ok=True)
    parsed = urllib.parse.urlparse(url)
    name = filename or Path(urllib.parse.unquote(parsed.path or "")).name or "download.bin"
    if not name:
        name = "download.bin"
    target = dest_dir / name
    request = urllib.request.Request(url, headers={"User-Agent": _UA})
    try:
        with urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT) as response:
            disposition = response.headers.get("Content-Disposition") or ""
            match = re.search(r'filename="?([^";]+)"?', disposition)
            if match and filename is None:
                target = dest_dir / Path(match.group(1)).name
            total = _to_int(response.headers.get("Content-Length")) or 0
            written = 0
            with open(target, "wb") as sink:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    written += len(chunk)
                    if written > MAX_DOWNLOAD_BYTES:
                        sink.close()
                        target.unlink(missing_ok=True)
                        raise _err("bad_request", f"下载文件超过 8GB 上限：{name}")
                    sink.write(chunk)
                    if total:
                        _emit(progress, f"正在下载 {name}", min(99, int(written * 100 / total)))
            _emit(progress, f"下载完成：{name}（{written / 1024 ** 2:.1f} MB）", 100)
            return target
    except EngineError:
        raise
    except urllib.error.HTTPError as error:
        raise _err("unsupported", f"下载失败（HTTP {error.code}）：{url}")
    except urllib.error.URLError as error:
        raise _err("unsupported", f"下载失败（网络不可达）：{getattr(error, 'reason', error)}")
    except Exception as error:  # noqa: BLE001
        raise _err("unsupported", f"下载失败：{error}")


def _resolve_source(source, dirs, options, inspect=False) -> tuple[Path, list[str]]:
    """把 path/url 解析成本地文件路径；目录则找唯一数据文件。"""
    options = options or {}
    if source is None or str(source).strip() == "":
        raise _err("bad_request", "缺少数据源：请提供 path 或 url")
    text = str(source).strip()
    notes = []
    if _looks_like_url(text):
        target = _download_to(text, dirs["downloads"])
        notes.append(f"已从 URL 下载：{target.name}")
        return target, notes
    path = Path(text).expanduser()
    try:
        path = path.resolve()
    except Exception:  # noqa: BLE001
        pass
    if not path.exists():
        hint = ""
        if text.startswith("~"):
            hint = "（~ 展开后的路径也不存在）"
        raise _err("not_found", f"路径不存在：{text}{hint}")
    if path.is_dir():
        entries = sorted(item for item in path.iterdir()
                         if item.is_file() and not item.name.startswith(".") and item.stat().st_size >= 8)
        recognizable = [item for item in entries if _classify_member(item.name) or
                        item.name.lower().endswith(ARCHIVE_SUFFIXES)]
        if len(recognizable) == 1:
            notes.append(f"目录里只有 {recognizable[0].name}，自动选中")
            return recognizable[0], notes
        if recognizable:
            raise _err("ambiguous", f"目录里有 {len(recognizable)} 个数据文件，请指定具体文件或 member",
                       candidates=[{"name": item.name, "size": int(item.stat().st_size),
                                    "kind": _classify_member(item.name) or "archive"}
                                   for item in recognizable])
        raise _err("bad_request", f"目录里没有可识别的数据文件：{path}")
    return path, notes


# ---------------------------------------------------------------------------
# inspect（不加载整表）
# ---------------------------------------------------------------------------


def _inspect_text(path: Path, options) -> dict:
    options = options or {}
    notes = []
    warnings = []
    head = _read_head(path, 2 * 1024 * 1024)
    if not head.strip():
        raise _err("bad_request", f"文件是空的：{path.name}")
    size = int(path.stat().st_size)
    if "!series_matrix_table_begin" in head:
        n_samples = n_genes = 0
        header_seen = False
        in_table = False
        try:
            with _open_text(path) as handle:
                for line in handle:
                    stripped = line.strip()
                    if stripped.startswith("!series_matrix_table_begin"):
                        in_table = True
                        continue
                    if stripped.startswith("!series_matrix_table_end"):
                        break
                    if in_table and stripped:
                        cells = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
                        if not header_seen:
                            n_samples = max(0, len(cells) - 1)
                            header_seen = True
                        elif any(cell.strip() for cell in cells):
                            n_genes += 1
        except Exception:  # noqa: BLE001
            pass
        notes.append("GEO series matrix（含 !series_matrix_table_begin）")
        notes.append("inspect 仅流式计数行列，未构建矩阵")
        return {"format": "geo_series", "n_samples": n_samples, "n_genes": n_genes,
                "size": size, "notes": notes, "warnings": warnings, "candidates": None}
    delimiter = options.get("delimiter") or _sniff_delimiter(head.splitlines())
    lines = head.splitlines()
    first_cells = _split_line(lines[0], delimiter) if lines else [""]
    header = _detect_header(first_cells)
    data_lines = [line for line in lines[1:41] if line.strip()]
    has_index = False
    if data_lines:
        first_cells_data = [_split_line(line, delimiter)[0] if _split_line(line, delimiter) else ""
                            for line in data_lines]
        non_numeric = sum(1 for cell in first_cells_data if not _is_numeric_text(cell))
        has_index = non_numeric >= max(1, int(0.6 * len(first_cells_data)))
    n_rows = 0
    try:
        with _open_text(path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                if n_rows == 0 and header:
                    n_rows += 1  # 表头行
                    continue
                n_rows += 1
    except Exception:  # noqa: BLE001
        n_rows = 0
    data_rows = max(0, n_rows - (1 if header else 0))
    n_cols = len(first_cells)
    if has_index:
        n_cols = max(0, n_cols - 1)
    row_labels = []
    if has_index:
        for line in data_lines[:100]:
            cells = _split_line(line, delimiter)
            if cells:
                row_labels.append(cells[0])
    col_labels = []
    if header:
        col_labels = [str(cell) for cell in first_cells[(1 if has_index else 0):]]
    elif first_cells:
        col_labels = []
    # 方向
    n_samples, n_genes = data_rows, n_cols
    if has_index and row_labels:
        rows_are_genes = _infer_rows_are_genes(row_labels, col_labels, (data_rows, n_cols), [])
        if rows_are_genes:
            n_samples, n_genes = n_cols, data_rows
    elif col_labels:
        # 没有行名时用列名判断（列名像基因 → 行是样本）
        if _score(col_labels, _strong_gene) >= 0.75 and _score(col_labels, _sample_like) < 0.2:
            n_samples, n_genes = data_rows, n_cols
    # 值域粗判（前若干行）
    values = []
    for line in data_lines[:100]:
        cells = _split_line(line, delimiter)
        cells = cells[1:] if has_index else cells
        for cell in list(cells)[:50]:
            if _is_numeric_text(cell):
                try:
                    values.append(float(str(cell).replace(",", "")))
                except Exception:  # noqa: BLE001
                    pass
    value_kind = "unknown"
    if values:
        arr = np.asarray(values, dtype=np.float64)
        value_kind, _, _ = _detect_value_kind(arr.reshape(1, -1))
    notes.append("分隔符：%s" % ("TAB" if delimiter == "\t" else ("空白" if delimiter == "whitespace" else delimiter)))
    notes.append("inspect 仅读文件头并流式计数行数，未加载整表")
    return {"format": "matrix_text", "n_samples": int(n_samples), "n_genes": int(n_genes),
            "value_kind": value_kind, "size": size, "delimiter": delimiter,
            "notes": notes, "warnings": warnings, "candidates": None}


def _inspect_xlsx(path: Path, options) -> dict:
    workbook = _open_workbook(path)
    try:
        sheets = []
        for name in workbook.sheetnames:
            worksheet = workbook[name]
            rows = int(worksheet.max_row or 0)
            columns = int(worksheet.max_column or 0)
            sheets.append({"name": name, "rows": rows, "cols": columns, "size": None, "kind": "sheet"})
        chosen = None
        wanted = options.get("sheet") if options else None
        for sheet in sheets:
            if wanted is not None and str(wanted) not in (sheet["name"], ""):
                continue
            if sheet["rows"] >= 2 and sheet["cols"] >= 2:
                chosen = sheet
                break
        if chosen is None:
            raise _err("bad_request", f"Excel 里没有可用的二维表：{path.name}")
        worksheet = workbook[chosen["name"]]
        sample_rows = []
        for i, row in enumerate(worksheet.iter_rows(values_only=True)):
            if i >= 20:
                break
            sample_rows.append(["" if cell is None else cell for cell in row])
        header = _detect_header(sample_rows[0]) if sample_rows else True
        data_rows = max(0, chosen["rows"] - (1 if header else 0))
        has_index = False
        if len(sample_rows) > 1:
            non_numeric = sum(1 for row in sample_rows[1:] if row and not _is_numeric_text(row[0]))
            has_index = non_numeric >= max(1, int(0.6 * (len(sample_rows) - 1)))
        first_header = str(sample_rows[0][0]).strip().lower() if sample_rows and sample_rows[0] else ""
        if not has_index and header and first_header in (
                "", "gene", "gene_id", "geneid", "id", "probe", "probe_id", "probeid",
                "transcript", "transcript_id", "symbol", "ensembl", "ensembl_id", "features", "feature"):
            has_index = True
        data_cols = max(0, chosen["cols"] - (1 if has_index else 0))
        n_samples, n_genes = data_rows, data_cols
        if len(sample_rows) > 1:
            row_labels = [str(row[0]) for row in sample_rows[1:] if row and str(row[0]).strip()]
            col_labels = [str(cell) for cell in (sample_rows[0][(1 if has_index else 0):] if sample_rows else [])]
            if row_labels:
                rows_are_genes = _infer_rows_are_genes(row_labels, col_labels,
                                                        (data_rows, data_cols), [])
                if rows_are_genes:
                    n_samples, n_genes = data_cols, data_rows
                elif _score(col_labels, _strong_gene) >= 0.75:
                    n_samples, n_genes = data_rows, data_cols
        notes = [f"Excel：{path.name}，sheet「{chosen['name']}」",
                 f"共 {len(sheets)} 个 sheet；已用 openpyxl read_only 只读维度与前 20 行"]
        return {"format": "xlsx", "n_samples": int(n_samples), "n_genes": int(max(0, n_genes)),
                "size": int(path.stat().st_size), "sheet": chosen["name"],
                "sheets": [sheet["name"] for sheet in sheets],
                "notes": notes, "warnings": [],
                "candidates": sheets if len(sheets) > 1 else None}
    finally:
        try:
            workbook.close()
        except Exception:  # noqa: BLE001
            pass


def _inspect_archive(path: Path) -> dict:
    name = path.name.lower()
    entries = []
    total = 0
    try:
        if name.endswith(".zip"):
            with zipfile.ZipFile(path) as archive:
                for info in archive.infolist():
                    if info.is_dir():
                        continue
                    kind = _classify_member(info.filename)
                    if kind is None:
                        continue
                    total += int(info.file_size)
                    entries.append({"name": info.filename, "size": int(info.file_size), "kind": kind})
        elif name.endswith(ARCHIVE_SUFFIXES) or tarfile.is_tarfile(path):
            with tarfile.open(path, "r:*") as archive:
                for member in archive.getmembers():
                    if not member.isfile():
                        continue
                    kind = _classify_member(member.name)
                    if kind is None:
                        continue
                    total += int(member.size)
                    entries.append({"name": member.name, "size": int(member.size), "kind": kind})
        else:
            raise _err("unsupported", f"不是支持的压缩格式：{path.name}")
    except EngineError:
        raise
    except (zipfile.BadZipFile, tarfile.TarError) as error:
        raise _err("bad_request", f"压缩包损坏或格式不识别（{path.name}）：{error}")
    entries.sort(key=lambda item: (-(3 if item["kind"] == "geo_series" else 2 if "matrix" in item["name"].lower() else 1),
                                   item["name"]))
    if not entries:
        raise _err("bad_request", f"压缩包里没有可识别的数据文件：{path.name}")
    notes = [f"压缩包：{path.name}，共 {len(entries)} 个候选数据文件",
             "inspect 只列目录，未解压"]
    return {"format": "archive", "n_samples": None, "n_genes": None,
            "size": int(path.stat().st_size), "candidates": entries[:50],
            "notes": notes, "warnings": []}


def _http_head_size(url: str) -> int:
    """HEAD 拿 Content-Length（失败返回 0）。"""
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": _UA})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return int(response.headers.get("Content-Length") or 0)
    except Exception:  # noqa: BLE001
        return 0


def _http_range(url: str, max_bytes: int) -> bytes:
    """Range 读取文件头（服务器不支持 Range 时最多也只读 max_bytes）。"""
    request = urllib.request.Request(url, headers={"User-Agent": _UA,
                                                   "Range": f"bytes=0-{max(1, max_bytes) - 1}"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read(max_bytes)


def _inspect_url(url: str, dirs, options) -> dict:
    """URL 输入只读头部预检（CONTRACT §2.2 inspect 支持 url）。"""
    parsed = urllib.parse.urlparse(url)
    name = Path(urllib.parse.unquote(parsed.path or "")).name or "download"
    lower = name.lower()
    inner = _strip_compression(lower)
    if any(lower.endswith(suffix) for suffix in ARCHIVE_SUFFIXES):
        kind = "archive"
    elif inner.endswith(XLSX_SUFFIXES):
        kind = "xlsx"
    elif inner.endswith((".rdata", ".rda")):
        kind = "rdata"
    elif inner.endswith(".rds"):
        kind = "rds"
    else:
        kind = "matrix_text"
    size = _http_head_size(url)
    limit = 64 * 1024 * 1024
    if kind in ("archive", "xlsx", "rds", "rdata"):
        if size and size > limit:
            return {"format": kind, "source": url, "size": size, "n_samples": None, "n_genes": None,
                    "candidates": None, "warnings": [],
                    "notes": [f"URL 文件约 {size / 1024 ** 2:.1f} MB，仅按扩展名识别；下载到本地后可精确 inspect"]}
        target = _download_to(url, dirs["downloads"], filename=name)
        info = inspect_source(str(target), workdir=str(dirs["base"]), options=options)
        info["source"] = url
        info["notes"] = [f"URL 文件已下载到临时目录（{name}）后 inspect"] + list(info.get("notes") or [])
        return info
    # 文本：优先 Range 只读头部
    partial = b""
    try:
        partial = _http_range(url, 2 * 1024 * 1024)
    except Exception:  # noqa: BLE001
        partial = b""
    if partial:
        suffix = Path(name).suffix or ".txt"
        digest = hashlib.md5(url.encode("utf-8")).hexdigest()[:8]
        target = dirs["tmp"] / f"url-head-{digest}{suffix}"
        try:
            target.write_bytes(partial)
            info = _inspect_text(target, options)
            info["source"] = url
            info["size"] = int(size or len(partial))
            info["notes"] = [f"URL 仅读取文件头 {len(partial) / 1024 ** 1:.0f} KB 预检（行数统计可能不完整）"] + \
                            list(info.get("notes") or [])
            return info
        except EngineError:
            pass
    if size and size > limit:
        return {"format": "matrix_text", "source": url, "size": size, "n_samples": None, "n_genes": None,
                "candidates": None, "warnings": [],
                "notes": [f"URL 文件约 {size / 1024 ** 2:.1f} MB，仅按扩展名识别"]}
    target = _download_to(url, dirs["downloads"], filename=name)
    info = inspect_source(str(target), workdir=str(dirs["base"]), options=options)
    info["source"] = url
    info["notes"] = [f"URL 文件已下载到临时目录（{name}）后 inspect"] + list(info.get("notes") or [])
    return info


def inspect_source(source, workdir=None, options=None, **opts) -> dict:
    """只探测格式与维度，不构建矩阵（CONTRACT §2.2 inspect）。"""
    merged = dict(options or {})
    merged.update(opts)
    dirs = _workdirs(workdir)
    if _looks_like_url(source):
        return _inspect_url(str(source).strip(), dirs, merged)
    path, resolve_notes = _resolve_source(source, dirs, merged, inspect=True)
    fmt = _detect_format(path, merged)
    if fmt in ("rds", "rdata"):
        target = _decompress_r(path, dirs) if path.name.lower().endswith(".gz") else path
        info = _load_r(target, dirs, {**merged, "inspect": True})
        info["source"] = str(path)
        info["size"] = int(path.stat().st_size)
        info["notes"] = list(resolve_notes) + list(info.get("notes") or [])
        return info
    if fmt == "archive":
        info = _inspect_archive(path)
    elif fmt == "xlsx":
        info = _inspect_xlsx(path, merged)
    elif fmt in ("geo_series", "matrix_text"):
        info = _inspect_text(path, merged)
    else:
        info = {"format": fmt, "notes": [], "warnings": [], "candidates": None}
    info["format"] = info.get("format") or fmt
    info["source"] = str(path)
    info["size"] = int(info.get("size") or (path.stat().st_size if path.is_file() else 0))
    info["notes"] = list(resolve_notes) + list(info.get("notes") or [])
    info.setdefault("candidates", None)
    info["suggested_import"] = info["format"] if info["format"] not in ("archive",) else "archive"
    return info


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------


def _load_path(path: Path, dirs, options, progress=None) -> dict:
    fmt = _detect_format(path, options)
    if fmt == "archive":
        return _load_archive(path, dirs, options, progress=progress)
    if fmt == "xlsx":
        _emit(progress, f"正在读取 Excel：{path.name}", 30)
        return _load_xlsx(path, options)
    if fmt in ("rds", "rdata"):
        _emit(progress, f"正在用 R 读取对象：{path.name}", 30)
        target = _decompress_r(path, dirs) if path.name.lower().endswith(".gz") else path
        result = _load_r(target, dirs, options)
        result["format"] = fmt
        return result
    if fmt in ("circ", "mirna", "lncrna"):
        return _load_feature_kind(path, fmt, options)
    if fmt == "fusion":
        return _load_fusion(path, options)
    if fmt == "splicing":
        return _load_splicing(path, options)
    if fmt in ("geo_series", "matrix_text"):
        _emit(progress, f"正在解析文本矩阵：{path.name}", 30)
        return _load_text(path, options)
    raise _err("unsupported", f"不支持的导入格式：{fmt}（支持：{', '.join(SUPPORTED_FORMATS)}）")


def load_source(source, workdir=None, options=None, progress=None, **opts) -> dict:
    """全格式导入主入口（CONTRACT §7.1）。"""
    merged = dict(options or {})
    merged.update(opts)
    dirs = _workdirs(workdir)
    _emit(progress, "正在识别数据源", 5)
    path, resolve_notes = _resolve_source(source, dirs, merged)
    _emit(progress, f"已定位文件：{path.name}", 10)
    result = _load_path(path, dirs, merged, progress=progress)
    result["name"] = str(merged.get("name") or (Path(result.get("name") or path.name).stem if path.name else path.name))
    result["notes"] = list(resolve_notes) + list(result.get("notes") or [])
    warnings = list(result.get("warnings") or [])
    if merged.get("species"):
        species = str(merged["species"])
        if species in ("human", "mouse", "unknown"):
            result["meta"]["species"] = species
    # 附加样本信息表 / 临床表 join
    for key, label in (("sample_sheet", "样本信息"), ("clinical", "临床")):
        extra_path = merged.get(key)
        if not extra_path:
            continue
        try:
            result["sample_info"] = _join_metadata(result["sample_info"], result["samples"],
                                                   extra_path, label, warnings)
            result["notes"].append(f"已 join {label}表：{Path(str(extra_path)).name}")
        except EngineError as error:
            warnings.append(f"{label}表读取失败（{error.message}）")
    result["warnings"] = warnings
    n_samples = result["meta"].get("n_samples") or 0
    n_genes = result["meta"].get("n_genes") or 0
    _emit(progress, f"解析完成：{n_samples} 个样本 × {n_genes} 个基因", 100)
    return result
