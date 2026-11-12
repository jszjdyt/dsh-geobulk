"""dsh-geobulk — 免疫浸润（bulk-clin 负责，CONTRACT §2.6 / §8.3）

步骤 `immune`：把归一化后的 bulk 表达矩阵转成免疫细胞比例/打分。

实现优先级
----------
1. **R 官方实现（默认优先）**：`lib/bg_immune.R`
   * `IOBR::deconvo_tme`：cibersort / epic / mcpcounter / xcell / quantiseq / estimate / timer
   * `GSVA::ssgseaParam`：ssgsea（签名来自 `bg_genesets.IMMUNE_28`，Charoentong 28 类）
   * `IOBR::deconvo_ips`：ips（Immunophenoscore）
   每个方法在 R 侧独立 tryCatch，失败跳过但保留 warning。
2. **Python 回退**（R 不可用 / IOBR 缺失 / R 侧全部方法失败）：
   * `ssgsea`：自实现 ssGSEA（rank^α 加权游走积分，等价 GSVA ssgseaParam 逻辑）+ IMMUNE_28
   * `estimate`：ESTIMATE-lite（免疫/基质基因集 z-score 均值 → ImmuneScore/StromalScore/ESTIMATEScore）
   * `mcpcounter`：MCP-counter-lite（marker 基因均值，10 个群体）

输出
----
* 样本 × 细胞比例表（`normalize=relative` 时每行和≈1；`absolute` 时保留原始比例/得分）
* 方法间相关矩阵（每个方法取特征矩阵第一主成分后做 Pearson，避免相对比例"行和=1"的虚假负相关）
* 组间差异表（2 组 Mann-Whitney U；>2 组 Kruskal-Wallis；BH FDR）
* 图：`immune_stack`（堆叠柱）、`immune_corr`（方法相关热图）

约定：只依赖标准库 + numpy/pandas/scipy/matplotlib + bg_common/rbridge/bg_genesets（纯函数）。
"""

from __future__ import annotations

import os
import re
import tempfile

import numpy as np
import pandas as pd

import bg_common as bc
import rbridge
from bg_common import EngineError, export_csv, result, save_fig, table

try:  # 内置签名（Charoentong 28 类）；缺失时退化为 MCP 面板
    from bg_genesets import IMMUNE_28  # type: ignore
except Exception:  # noqa: BLE001
    IMMUNE_28 = {}

__all__ = ["STEPS", "ANALYSES", "run_immune"]

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

ALL_METHODS = [
    "cibersort", "epic", "mcpcounter", "xcell", "quantiseq",
    "estimate", "ssgsea", "timer", "ips",
]
R_METHODS = [m for m in ALL_METHODS if m != "ips"]
PROPORTION_METHODS = ("cibersort", "epic", "quantiseq")  # 天然输出比例（行和≈1）
SCORE_METHODS = ("mcpcounter", "xcell", "estimate", "ssgsea", "timer", "ips")

METHOD_LABELS = {
    "cibersort": "CIBERSORT (IOBR, LM22)",
    "epic": "EPIC (IOBR)",
    "mcpcounter": "MCP-counter (IOBR)",
    "xcell": "xCell (IOBR)",
    "quantiseq": "quanTIseq (IOBR)",
    "estimate": "ESTIMATE (IOBR)",
    "ssgsea": "ssGSEA (GSVA + Charoentong 28)",
    "timer": "TIMER (IOBR)",
    "ips": "IPS / Immunophenoscore (IOBR)",
}

# CIBERSORT 输出里混着 P-value/Correlation/RMSE：不是细胞比例，不参与堆叠/相关
_META_COL_RE = re.compile(r"(p[\s_.-]?value|correlation|rmse|pval|\brho\b)", re.I)

# MCP-counter 式 marker 面板（Python 回退用；官方实现走 R IOBR）
MCP_PANELS = {
    "T cells": ["CD3D", "CD3E", "CD3G", "CD2", "CD6", "SH2D1A", "TRAT1"],
    "CD8 T cells": ["CD8A", "CD8B"],
    "Cytotoxic lymphocytes": ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "KLRD1", "NKG7"],
    "NK cells": ["KIR2DL3", "KIR3DL1", "KLRD1", "KLRF1", "NCAM1"],
    "B lineage": ["CD19", "CD79A", "CD79B", "MS4A1", "TNFRSF17"],
    "Monocytic lineage": ["CD14", "CD68", "CSF1R", "FCGR1A", "LYZ", "FCN1"],
    "Myeloid dendritic cells": ["CD1C", "CD1E", "CLEC10A", "FCER1A"],
    "Neutrophils": ["FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9"],
    "Endothelial cells": ["PECAM1", "VWF", "CDH5", "CLDN5", "ENG"],
    "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRB", "THY1", "POSTN"],
}

STROMAL_MARKERS = [
    "COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "ACTA2", "TAGLN", "FAP",
    "PDGFRB", "THY1", "POSTN", "MMP2", "MMP9", "SPARC", "LUM", "DCN",
]

# ESTIMATE-lite 的免疫分数字典：IMMUNE_28 并集 + MCP 免疫面板 + 常见免疫效应/检查点
_EXTRA_IMMUNE_MARKERS = [
    "PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD4", "MS4A1", "CD19", "CD79A",
    "LYZ", "CD14", "FCGR3A", "ITGAM", "NKG7", "GNLY", "PRF1", "GZMB",
    "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "B2M", "TAP1", "TAP2",
    "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "PDCD1",
]


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------

def _outdir(ctx, outdir) -> str:
    path = outdir or (ctx or {}).get("outdir")
    if not path:
        path = tempfile.mkdtemp(prefix="bg-immune-")
    os.makedirs(path, exist_ok=True)
    return path


def _as_str_list(value) -> list[str]:
    """把 list / 逗号或空格分隔字符串统一成去重小写列表。"""
    if value is None:
        return []
    if isinstance(value, str):
        raw = [x for x in re.split(r"[,\s;|]+", value) if x]
    else:
        try:
            raw = [str(x) for x in value if str(x).strip()]
        except TypeError:
            raw = [str(value)]
    out: list[str] = []
    seen: set[str] = set()
    for item in raw:
        key = item.strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def _methods_from_params(params) -> list[str]:
    requested = _as_str_list((params or {}).get("methods"))
    if not requested:
        return list(R_METHODS)  # 默认 8 个 R 方法（ips 需要显式请求）
    valid = [m for m in requested if m in ALL_METHODS]
    unknown = [m for m in requested if m not in ALL_METHODS]
    if unknown and not valid:
        raise EngineError("bad_request", f"未知免疫浸润方法：{', '.join(unknown)}；可选：{', '.join(ALL_METHODS)}")
    return valid


def _ascii(text, fallback: str = "group") -> str:
    text = str(text)
    return text if text.isascii() else fallback


def _safe_float(value):
    try:
        out = float(value)
        if not np.isfinite(out):
            return None
        return out
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# ctx → 表达矩阵
# ---------------------------------------------------------------------------

def _ctx_data(ctx) -> tuple[np.ndarray, list[str], list[str], pd.DataFrame]:
    """校验并取出 X / genes / samples / sample_info（info 按 samples 重排）。"""
    if not isinstance(ctx, dict):
        raise EngineError("internal", "ctx 不是 dict，无法运行免疫浸润")
    X = np.asarray(ctx.get("X"), dtype=np.float64)
    genes = [str(g) for g in (ctx.get("genes") or [])]
    samples = [str(s) for s in (ctx.get("samples") or [])]
    if X.ndim != 2:
        raise EngineError("bad_request", "表达矩阵不是二维数组")
    if X.shape[0] != len(samples):
        raise EngineError("internal", f"表达矩阵行数（{X.shape[0]}）与样本数（{len(samples)}）不一致")
    if X.shape[1] != len(genes):
        raise EngineError("internal", f"表达矩阵列数（{X.shape[1]}）与基因数（{len(genes)}）不一致")
    info = ctx.get("sample_info")
    if not isinstance(info, pd.DataFrame):
        info = pd.DataFrame(index=samples)
    info = info.copy()
    info.index = [str(i) for i in info.index]
    missing = [s for s in samples if s not in info.index]
    if missing:
        info = pd.concat([info, pd.DataFrame(index=missing)], axis=0)
    info = info.reindex(samples)
    return X, genes, samples, info


def _prepare_expr(ctx) -> tuple[np.ndarray, list[str], list[str], pd.DataFrame, list[str], list[str]]:
    """免疫分析输入：补齐 NaN、负值截断、必要时把线性尺度转 log2。

    返回 (expr_log, genes, samples, info, notes, warnings)。expr_log 为 log 尺度（样本×基因）。
    """
    X, genes, samples, info = _ctx_data(ctx)
    notes: list[str] = []
    warnings: list[str] = []
    if len(samples) < 3:
        raise EngineError("prerequisite", f"免疫浸润至少需要 3 个样本（当前 {len(samples)} 个）")
    if len(genes) < 20:
        raise EngineError("prerequisite", f"免疫浸润至少需要 20 个基因（当前 {len(genes)} 个）")

    expr = np.array(X, dtype=np.float64, copy=True)
    # 全 NaN 基因列 → 0；其余 NaN 用基因中位数填补
    if not np.isfinite(expr).all():
        col_median = np.nanmedian(np.where(np.isfinite(expr), expr, np.nan), axis=0)
        col_median = np.where(np.isfinite(col_median), col_median, 0.0)
        bad = ~np.isfinite(expr)
        if bad.any():
            expr[bad] = np.take(col_median, np.where(bad)[1])
            warnings.append(f"表达矩阵含 {int(bad.sum())} 个非有限值，已用基因中位数填补")
    if float(np.nanmin(expr)) < 0:
        expr = np.clip(expr, 0.0, None)
        notes.append("表达矩阵含负值，已截断为 0（免疫去卷积要求非负）")

    kind = str(ctx.get("X_kind") or "unknown").lower()
    linear = kind in ("counts", "tpm", "fpkm") or float(np.nanmax(expr)) > 30.0
    if linear:
        expr = np.log2(expr + 1.0)
        notes.append("输入为非 log 尺度，已按 log2(x+1) 转换后做免疫浸润")
    return expr, genes, samples, info, notes, warnings


# ---------------------------------------------------------------------------
# R 路径
# ---------------------------------------------------------------------------

def _signature_table() -> pd.DataFrame:
    rows = [(str(cell), str(gene).upper()) for cell, genes in (IMMUNE_28 or {}).items() for gene in genes]
    if not rows:
        return pd.DataFrame(columns=["cell", "gene"])
    return pd.DataFrame(rows, columns=["cell", "gene"])


def _run_r_immune(expr, genes, samples, info, params, outdir):
    methods = _methods_from_params(params)
    r_methods = [m for m in methods if m != "ips"] or ["estimate"]
    payload = {
        "expr": expr,
        "genes": genes,
        "samples": samples,
        "info": info,
        "params": {
            "methods": r_methods,
            "ips": "ips" in methods,
            "perm": int(params.get("perm") or 100),
            "seed": int(params.get("seed") or 0),
            "indications": params.get("indications") or "LUAD",
            "platform": params.get("platform") or "affymetrix",
        },
        "extra": {"signatures": _signature_table()},
    }
    return rbridge.run_r("bg_immune.R", payload, outdir, timeout=int(params.get("timeout") or 1200))


def _coerce_sample_frame(df, samples) -> pd.DataFrame | None:
    """把 R 输出的「sample + 数值列」表转成 index=samples 的数值 DataFrame。"""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return None
    frame = df.copy()
    sample_set = set(samples)
    sample_col = None
    for col in frame.columns:
        try:
            values = [str(v) for v in frame[col].tolist()]
        except Exception:  # noqa: BLE001
            continue
        if len(set(values) & sample_set) >= max(1, len(samples) // 2):
            sample_col = col
            break
    if sample_col is not None:
        frame = frame.set_index(frame[sample_col].astype(str))
        frame = frame.drop(columns=[sample_col])
    else:
        frame.index = [str(i) for i in frame.index]
        if not (set(frame.index) & sample_set):
            return None
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.loc[:, ~numeric.isna().all(axis=0)]
    # 去掉 CIBERSORT 的 P-value/Correlation/RMSE 等非比例列
    keep = [c for c in numeric.columns if not _META_COL_RE.search(str(c))]
    if keep:
        numeric = numeric[keep]
    if numeric.shape[1] == 0:
        return None
    numeric = numeric.reindex([s for s in samples if s in numeric.index])
    if numeric.empty:
        return None
    return numeric


def _parse_r_output(run: dict, samples) -> tuple[dict, pd.DataFrame | None, list[dict], list[str], list[str]]:
    """解析 run_r 的返回：{method: 样本×特征 DataFrame}, ips 表, 方法状态, warnings, notes。"""
    tables = run.get("tables") or {}
    payload = run.get("result") or {}
    per_method: dict[str, pd.DataFrame] = {}
    ips_frame = None
    for name, frame in tables.items():
        key = str(name)
        if key.startswith("immune_"):
            method = key[len("immune_"):]
            parsed = _coerce_sample_frame(frame, samples)
            if parsed is not None and parsed.shape[1] >= 1:
                per_method[method] = parsed
        elif key == "ips":
            parsed = _coerce_sample_frame(frame, samples)
            if parsed is not None:
                ips_frame = parsed
    if ips_frame is not None:
        per_method["ips"] = ips_frame

    statuses: list[dict] = []
    raw_methods = payload.get("methods")
    if isinstance(raw_methods, list):
        for item in raw_methods:
            if isinstance(item, dict) and item.get("method"):
                statuses.append(item)
    elif isinstance(raw_methods, dict):
        for item in raw_methods.values():
            if isinstance(item, dict) and item.get("method"):
                statuses.append(item)
    warnings = [str(w) for w in (payload.get("warnings") or [])]
    notes = [str(n) for n in (payload.get("notes") or [])]
    return per_method, ips_frame, statuses, warnings, notes


# ---------------------------------------------------------------------------
# Python 回退：ssGSEA / ESTIMATE-lite / MCP-counter-lite
# ---------------------------------------------------------------------------

def _gene_index(genes: list[str]) -> dict[str, int]:
    idx: dict[str, int] = {}
    for i, gene in enumerate(genes):
        key = str(gene).strip().upper()
        if key and key not in idx:
            idx[key] = i
    return idx


def _ssgsea(expr: np.ndarray, signatures: dict[str, list[str]], genes: list[str], alpha: float = 0.25) -> pd.DataFrame:
    """自实现 ssGSEA（GSVA ssgseaParam 逻辑的向量化版本）。

    ES = Σ_j [ P_in(j) - P_out(j) ]，其中 j 为按表达升序的基因位置，
    P_in 为 |rank|^alpha 加权累计占比。闭式化简后
        ES_i = Σ_{g∈set} w_g·(n - rank_g + 1)/W_total - (n_out+1)/2
    再按签名跨样本 min-max 归一化（GSVA normalize=TRUE）。
    """
    samples = expr.shape[0]
    n_genes = expr.shape[1]
    idx_map = _gene_index(genes)
    genes_by_cell = {}
    for cell, gene_list in signatures.items():
        idxs = sorted({idx_map[str(g).upper()] for g in gene_list if str(g).upper() in idx_map})
        if len(idxs) >= 3:
            genes_by_cell[str(cell)] = np.array(idxs, dtype=int)
    if len(genes_by_cell) < 2:
        return pd.DataFrame(index=[f"S{i}" for i in range(samples)])

    # 升序秩（1..n），并列取平均
    try:
        from scipy.stats import rankdata
    except Exception:  # noqa: BLE001
        rankdata = None
    if rankdata is not None:
        ranks = np.vstack([rankdata(row) for row in expr]).astype(np.float32)
    else:  # pragma: no cover - 极简回退
        ranks = np.argsort(np.argsort(expr, axis=1), axis=1).astype(np.float32) + 1.0
    weights = np.power(ranks, float(alpha))
    w_total = weights.sum(axis=1)
    w_total = np.where(w_total > 0, w_total, 1.0)

    cells = list(genes_by_cell.keys())
    indicator = np.zeros((n_genes, len(cells)), dtype=np.float32)
    for k, cell in enumerate(cells):
        indicator[genes_by_cell[cell], k] = 1.0
    const = np.array([(n_genes - len(genes_by_cell[c]) + 1) / 2.0 for c in cells], dtype=np.float32)

    es = np.zeros((samples, len(cells)), dtype=np.float64)
    chunk = max(1, int(4_000_000 // max(n_genes, 1)))  # 控制单块内存 < 32MB(float32)
    for start in range(0, samples, chunk):
        stop = min(samples, start + chunk)
        part = (weights[start:stop] * (n_genes - ranks[start:stop] + 1.0)) @ indicator
        part = part / w_total[start:stop, None] - const[None, :]
        es[start:stop] = part

    # 跨样本 min-max 归一化 → [0,1]
    lo = es.min(axis=0, keepdims=True)
    hi = es.max(axis=0, keepdims=True)
    span = np.where(hi - lo > 1e-12, hi - lo, 1.0)
    norm = (es - lo) / span
    return pd.DataFrame(norm, columns=cells)


def _estimate_lite(expr: np.ndarray, genes: list[str]) -> pd.DataFrame:
    """ESTIMATE-lite：免疫/基质基因集 z-score 均值 → 三个分数 + 近似肿瘤纯度。"""
    idx_map = _gene_index(genes)
    immune_pool = set()
    for gene_list in (IMMUNE_28 or {}).values():
        immune_pool.update(str(g).upper() for g in gene_list)
    for panel in MCP_PANELS.values():
        immune_pool.update(str(g).upper() for g in panel)
    immune_pool.update(g.upper() for g in _EXTRA_IMMUNE_MARKERS)
    immune_idx = sorted({idx_map[g] for g in immune_pool if g in idx_map})
    stromal_idx = sorted({idx_map[g] for g in STROMAL_MARKERS if g in idx_map})
    if len(immune_idx) < 5 or len(stromal_idx) < 3:
        return pd.DataFrame()
    mean = np.nanmean(expr, axis=0)
    std = np.nanstd(expr, axis=0, ddof=1)
    std = np.where(std > 1e-9, std, 1.0)
    z = (expr - mean) / std
    immune_score = z[:, immune_idx].mean(axis=1)
    stromal_score = z[:, stromal_idx].mean(axis=1)
    estimate_score = immune_score + stromal_score
    # 简化纯度：ESTIMATEScore 越高 → 纯度越低（sigmoid 映射，非官方 ESTIMATE 绝对值）
    purity = 1.0 - 1.0 / (1.0 + np.exp(-(estimate_score - 0.5)))
    return pd.DataFrame(
        {
            "StromalScore": stromal_score,
            "ImmuneScore": immune_score,
            "ESTIMATEScore": estimate_score,
            "TumorPurity": np.clip(purity, 0.0, 1.0),
        }
    )


def _mcp_lite(expr: np.ndarray, genes: list[str]) -> pd.DataFrame:
    """MCP-counter-lite：每个群体的 marker 均值（log 表达）。"""
    idx_map = _gene_index(genes)
    out: dict[str, np.ndarray] = {}
    for label, markers in MCP_PANELS.items():
        idxs = [idx_map[str(g).upper()] for g in markers if str(g).upper() in idx_map]
        if len(idxs) < 2:
            continue
        out[label] = expr[:, idxs].mean(axis=1)
    if not out:
        return pd.DataFrame()
    return pd.DataFrame(out)


def _fallback_immune(expr: np.ndarray, genes: list[str], warnings: list[str], notes: list[str]) -> dict:
    """R 不可用时的 Python 回退（三种方法，方法名诚实标注 lite/自实现）。"""
    per_method: dict[str, pd.DataFrame] = {}
    signatures = {str(k): [str(g) for g in v] for k, v in (IMMUNE_28 or {}).items()}
    if signatures:
        ssgsea = _ssgsea(expr, signatures, genes)
        if ssgsea.shape[1] >= 3:
            per_method["ssgsea"] = ssgsea
            notes.append("ssGSEA：Python 自实现（rank^0.25 加权游走积分 + 跨样本 min-max），签名 = IMMUNE_28")
        else:
            warnings.append(
                f"ssGSEA 签名与表达矩阵重叠不足（可算细胞类型 {ssgsea.shape[1]} < 3），已跳过"
            )
    estimate = _estimate_lite(expr, genes)
    if not estimate.empty:
        per_method["estimate"] = estimate
        notes.append("ESTIMATE：Python 回退（ESTIMATE-lite，免疫/基质 marker z-score 均值；TumorPurity 为近似值）")
    mcp = _mcp_lite(expr, genes)
    if not mcp.empty:
        per_method["mcpcounter"] = mcp
        notes.append("MCP-counter：Python 回退（marker 均值，非官方 rank 法）")
    if not per_method:
        raise EngineError(
            "prerequisite",
            "免疫浸润失败：R/IOBR 不可用且 Python 回退的签名基因与表达矩阵重叠不足（请确认是人的 bulk 表达数据）",
        )
    return per_method


# ---------------------------------------------------------------------------
# 结果整理：比例 / 相关 / 组间差异
# ---------------------------------------------------------------------------

def _relative_rows(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.astype(np.float64).clip(lower=0.0)
    row_sum = values.sum(axis=1)
    row_sum = row_sum.where(row_sum > 0, 1.0)
    return values.div(row_sum, axis=0)


def _build_proportions(per_method: dict, normalize: str, warnings: list[str], notes: list[str]):
    """选出主比例表：优先天然比例方法，否则用 ssGSEA/MCP 相对化。

    返回 (frame, meta)。
    """
    normalize = "absolute" if str(normalize).lower() == "absolute" else "relative"
    primary = next((m for m in PROPORTION_METHODS if m in per_method), None)
    derived = False
    if primary is not None:
        frame = per_method[primary].astype(np.float64).clip(lower=0.0)
        method_label = primary
        if normalize == "relative":
            frame = _relative_rows(frame)
            notes.append(f"细胞比例表来自 {primary}（相对比例，行和归一化为 1）")
        else:
            notes.append(f"细胞比例表来自 {primary}（absolute：保留方法原始输出，未做行归一化）")
        n_cells = frame.shape[1]
        return frame, {
            "primary_method": primary,
            "primary_label": METHOD_LABELS.get(primary, primary),
            "derived": False,
            "normalize": normalize,
            "n_cells": int(n_cells),
        }

    # 没有天然比例方法 → 用 ssGSEA（优先，28 类且最细）或 MCP-counter 的 marker 得分相对化
    fallback_order = ["ssgsea", "mcpcounter", "xcell", "timer"]
    source = next((m for m in fallback_order if m in per_method), None)
    if source is None:
        # 只剩 estimate/ips：用免疫分数一列，作为相对值的退化表示
        source = next(iter(per_method.keys()))
        frame = per_method[source]
        if source == "estimate" and "ImmuneScore" in frame.columns:
            frame = frame[["ImmuneScore"]]
        frame = frame.astype(np.float64).clip(lower=0.0)
        if normalize == "relative":
            base = frame.iloc[:, 0].values
            lo, hi = np.nanmin(base), np.nanmax(base)
            scaled = (base - lo) / (hi - lo) if hi - lo > 1e-12 else np.full_like(base, 1.0)
            frame = pd.DataFrame({"Immune_relative": scaled}, index=frame.index)
        warnings.append(f"没有天然比例类方法（cibersort/epic/quantiseq），比例表退化为 {source} 的免疫分数")
        return frame, {
            "primary_method": source,
            "primary_label": METHOD_LABELS.get(source, source),
            "derived": True,
            "normalize": normalize,
            "n_cells": int(frame.shape[1]),
        }

    frame = per_method[source].astype(np.float64).clip(lower=0.0)
    if normalize == "relative":
        # 每个特征先跨样本 min-max → 相对得分，再行归一化
        lo = frame.min(axis=0)
        hi = frame.max(axis=0)
        span = (hi - lo).where(hi - lo > 1e-12, 1.0)
        frame = (frame - lo).div(span)
        frame = _relative_rows(frame)
    derived = True
    notes.append(
        f"无天然比例方法，细胞比例表由 {source} 得分相对化得到（{normalize}；方法名已标注为 derived）"
    )
    return frame, {
        "primary_method": source,
        "primary_label": METHOD_LABELS.get(source, source),
        "derived": derived,
        "normalize": normalize,
        "n_cells": int(frame.shape[1]),
    }


def _method_score(per_method: dict) -> pd.DataFrame:
    """每个方法取特征矩阵第一主成分（符号对齐行均值），得到方法级得分向量。"""
    out: dict[str, np.ndarray] = {}
    for method, frame in per_method.items():
        values = frame.astype(np.float64).values
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
        if values.shape[1] == 0:
            continue
        mean = np.nanmean(values, axis=0)
        std = np.nanstd(values, axis=0, ddof=1)
        std = np.where(std > 1e-12, std, 1.0)
        z = (values - mean) / std
        if z.shape[1] == 1:
            score = z[:, 0]
        else:
            try:
                u, s, _vt = np.linalg.svd(z - z.mean(axis=0, keepdims=True), full_matrices=False)
                score = u[:, 0] * s[0]
            except np.linalg.LinAlgError:  # pragma: no cover
                score = z.mean(axis=1)
            if np.corrcoef(score, z.mean(axis=1))[0, 1] < 0:
                score = -score
        out[method] = score
    return pd.DataFrame(out)


def _method_corr(per_method: dict) -> pd.DataFrame:
    scores = _method_score(per_method)
    if scores.empty:
        return pd.DataFrame()
    scores.index = [str(i) for i in scores.index]
    corr = scores.corr(method="pearson")
    return corr.round(4)


def _concat_features(per_method: dict, limit: int = 200) -> pd.DataFrame:
    frames = []
    for method, frame in per_method.items():
        prefixed = frame.copy()
        prefixed.columns = [f"{method}::{str(c)}" for c in prefixed.columns]
        frames.append(prefixed)
    if not frames:
        return pd.DataFrame()
    wide = pd.concat(frames, axis=1)
    if wide.shape[1] > limit:
        wide = wide.iloc[:, :limit]
    return wide


def _group_context(ctx, info, params) -> tuple[str | None, list, list[str]]:
    """确定分组列与分组名。返回 (group_col, labels(n_samples), groups)。"""
    group_col = params.get("group_col") or (ctx or {}).get("group_key")
    samples = [str(s) for s in ((ctx or {}).get("samples") or [])]
    if not group_col or group_col not in getattr(info, "columns", []):
        return None, [], []
    labels = info[group_col].astype(str).values
    hint = [str(g) for g in ((ctx or {}).get("groups") or [])]
    uniq = [str(v) for v in pd.unique(labels) if str(v) not in ("", "nan", "None")]
    groups = [g for g in hint if g in uniq] if hint else []
    for g in uniq:
        if g not in groups:
            groups.append(g)
    if len(groups) < 2:
        return None, [], []
    return str(group_col), [str(v) for v in labels], groups


def _group_diff(prop: pd.DataFrame, labels: list[str], groups: list[str]) -> pd.DataFrame:
    """组间差异：2 组 Mann-Whitney U；>2 组 Kruskal-Wallis；BH FDR。"""
    if not labels or len(labels) != prop.shape[0] or len(groups) < 2:
        return pd.DataFrame()
    labels_arr = np.array(labels)
    rows = []
    for feature in prop.columns:
        values = pd.to_numeric(prop[feature], errors="coerce").values
        masks = [np.isfinite(values) & (labels_arr == group) for group in groups]
        counts = [int(m.sum()) for m in masks]
        if min(counts) < 2:
            continue
        samples_per_group = [values[m] for m in masks]
        try:
            if len(groups) == 2:
                from scipy import stats

                stat, pval = stats.mannwhitneyu(samples_per_group[0], samples_per_group[1], alternative="two-sided")
                test_name = "Mann-Whitney U"
            else:
                from scipy import stats

                stat, pval = stats.kruskal(*samples_per_group)
                test_name = "Kruskal-Wallis"
        except Exception:  # noqa: BLE001
            continue
        means = [float(np.mean(arr)) for arr in samples_per_group]
        row = {
            "feature": str(feature),
            "test": test_name,
            "n_groups": len(groups),
            "statistic": _safe_float(stat),
            "p": _safe_float(pval),
        }
        for group, n, mean in zip(groups, counts, means):
            row[f"n_{group}"] = n
            row[f"mean_{group}"] = mean
        if len(groups) == 2:
            row["delta"] = means[0] - means[1]
            row["direction"] = f"higher_in_{groups[0]}" if means[0] >= means[1] else f"higher_in_{groups[1]}"
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    frame["fdr"] = bc.bh_fdr(frame["p"].values)
    frame = frame.sort_values(["fdr", "p"], na_position="last").reset_index(drop=True)
    return frame


# ---------------------------------------------------------------------------
# 图
# ---------------------------------------------------------------------------

def _cell_color_map(columns) -> dict[str, str]:
    return {str(col): bc.color_for(str(col), i) for i, col in enumerate(columns)}


def _plot_stack(prop: pd.DataFrame, labels: list[str], groups: list[str], outdir: str) -> dict:
    """样本 × 细胞比例堆叠柱（按分组排序，组间留白）。"""
    frame = prop.copy()
    frame.index = [str(i) for i in frame.index]
    if frame.shape[1] > 24:  # 细胞类型过多时合并低频
        means = frame.mean(axis=0).sort_values(ascending=False)
        keep = list(means.index[:20])
        other = [c for c in frame.columns if c not in keep]
        if other:
            frame = pd.concat([frame[keep], frame[other].sum(axis=1).rename("Other")], axis=1)
    order = np.arange(len(frame))
    if labels and len(labels) == len(frame) and groups:
        labels_arr = np.array(labels)
        order = np.concatenate([np.where(labels_arr == group)[0] for group in groups if (labels_arr == group).any()])
    data = frame.iloc[order]
    fig, ax = bc.new_fig(11.0, 5.0)
    bottom = np.zeros(len(data))
    colors = _cell_color_map(data.columns)
    for col in data.columns:
        values = np.nan_to_num(data[col].astype(float).values, nan=0.0)
        ax.bar(np.arange(len(data)), values, bottom=bottom, width=1.0, color=colors[str(col)], label=str(col), linewidth=0)
        bottom += values
    ax.set_xlim(-0.5, len(data) - 0.5)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Estimated proportion / score")
    ax.set_title("Immune cell composition")
    ax.grid(False)
    if len(data) <= 48:
        ax.set_xticks(np.arange(len(data)))
        ax.set_xticklabels([str(v) for v in data.index], rotation=90, fontsize=6)
    else:
        ax.set_xticks([])
    # 分组分界线
    if labels and groups and len(labels) == len(frame):
        start = 0
        for group in groups:
            n = int((np.array(labels) == group).sum())
            if n == 0:
                continue
            if start > 0:
                ax.axvline(start - 0.5, color="#8a919c", linewidth=0.8, linestyle="--")
            ax.text(start + n / 2.0, 1.01, _ascii(group, "group"), transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=8, color="#5a6068")
            start += n
    if data.shape[1] <= 24:
        ax.legend(fontsize=6.5, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return save_fig(fig, outdir, "immune_stack", "Immune cell composition")


def _plot_corr(corr: pd.DataFrame, outdir: str) -> dict:
    """方法间相关热图（Pearson，第一主成分得分）。"""
    fig, ax = bc.new_fig(7.2, 6.2)
    values = corr.values.astype(float)
    im = ax.imshow(values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.index)))
    ax.set_xticklabels([str(c) for c in corr.columns], rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels([str(i) for i in corr.index], fontsize=8)
    ax.grid(False)
    ax.set_title("Method correlation (PC1 scores)")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(value) > 0.6 else "#2b2f36")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="Pearson r")
    return save_fig(fig, outdir, "immune_corr", "Method correlation")


def _plot_group_box(prop: pd.DataFrame, group_frame: pd.DataFrame, labels: list[str], groups: list[str],
                    outdir: str) -> dict | None:
    """组间差异显著的细胞类型箱线图（最多 8 个）。"""
    if group_frame is None or group_frame.empty or not labels or len(groups) != 2:
        return None
    try:
        from scipy import stats  # noqa: F401
    except Exception:  # noqa: BLE001
        return None
    feats = [f for f in group_frame["feature"].tolist() if f in prop.columns][:8]
    if not feats:
        return None
    groups_arr = np.array(labels)
    ncol = min(4, len(feats))
    nrow = int(np.ceil(len(feats) / ncol))
    fig, axes = bc.new_fig(3.2 * ncol, 2.9 * nrow, nrows=nrow, ncols=ncol, squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for k, feature in enumerate(feats):
        ax = axes[k // ncol][k % ncol]
        ax.set_visible(True)
        data = [prop[feature].values[groups_arr == g] for g in groups]
        data = [np.array([v for v in arr if np.isfinite(v)], dtype=float) for arr in data]
        data = [arr for arr in data if arr.size]
        if len(data) < 2:
            continue
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False)
        for patch, group in zip(bp["boxes"], groups):
            patch.set_facecolor(bc.color_for(str(group), groups.index(group)))
            patch.set_alpha(0.75)
            patch.set_edgecolor("#5a6068")
        for median in bp["medians"]:
            median.set_color("#2b2f36")
        row = group_frame[group_frame["feature"] == feature]
        pval = _safe_float(row["fdr"].iloc[0]) if not row.empty else None
        ax.set_title(f"{str(feature)[:30]}" + (f"  FDR={pval:.2g}" if pval is not None else ""), fontsize=8.5)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([_ascii(g, f"G{i + 1}") for i, g in enumerate(groups)], fontsize=7.5)
        ax.grid(axis="y", color="#e8eaee", linewidth=0.6)
    fig.suptitle("Immune features by group", fontsize=10.5)
    return save_fig(fig, outdir, "immune_group_box", "Immune group comparison")


# ---------------------------------------------------------------------------
# 步骤主入口
# ---------------------------------------------------------------------------

def run_immune(ctx, params=None, outdir=None) -> dict:
    """步骤 `immune`：免疫浸润（IOBR 优先，Python 回退）。"""
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    methods = _methods_from_params(params)
    normalize = "absolute" if str(params.get("normalize") or "relative").lower() == "absolute" else "relative"

    expr, genes, samples, info, notes, warnings = _prepare_expr(ctx)
    per_method: dict[str, pd.DataFrame] = {}
    statuses: list[dict] = []
    engine_label = "IOBR/GSVA (R 官方实现)"

    # ---- 1) R 优先 ----
    r_error = None
    try:
        if not rbridge.r_available():
            raise EngineError("no_r", "本机没有可用的 Rscript")
        run = _run_r_immune(expr, genes, samples, info, params, outdir)
        per_method, _ips, statuses, r_warnings, r_notes = _parse_r_output(run, samples)
        warnings.extend(r_warnings)
        notes.extend(r_notes)
        if not per_method:
            raise EngineError("internal", "R 侧未返回任何可用的免疫浸润方法结果")
    except EngineError as error:
        if error.code == "prerequisite":
            raise
        r_error = error.message
        per_method = {}
        statuses = []
        engine_label = "Python 回退（无 R/IOBR）"
        warnings.append(f"R/IOBR 路径不可用：{r_error[:200]}")
    except Exception as error:  # noqa: BLE001 - R 侧任何异常都回退
        r_error = f"{type(error).__name__}: {error}"
        per_method = {}
        statuses = []
        engine_label = "Python 回退（无 R/IOBR）"
        warnings.append(f"R/IOBR 路径异常：{r_error[:200]}")

    r_used = engine_label.startswith("IOBR")
    if not r_used:
        per_method = _fallback_immune(expr, genes, warnings, notes)
        notes.append("引擎：Python 回退（ssGSEA 自实现 + ESTIMATE-lite + MCP-counter-lite），方法名均为 Python 实现")

    # 只保留请求的方法（R 侧可能默认多跑了若干；回退则保留其算出的方法）
    if r_used:
        wanted = set(methods)
        per_method = {m: f for m, f in per_method.items() if m in wanted or (m == "ips" and "ips" in wanted)}
    if not per_method:
        raise EngineError("internal", "免疫浸润没有得到任何可用方法结果")

    # 方法状态（成功/失败原因）
    failed = {}
    for item in statuses:
        method = str(item.get("method"))
        if method in per_method:
            continue
        reason = str(item.get("error") or "计算失败")
        if method in methods or method == "ips":
            failed[method] = reason
    for item in warnings:
        match = re.match(r"^([a-z0-9_]+) 计算失败已跳过[：:](.*)$", item)
        if match and match.group(1) not in per_method:
            failed.setdefault(match.group(1), match.group(2).strip())
    requested_missing = [m for m in methods if m not in per_method and m not in failed]
    for method in requested_missing:
        failed[method] = "R 侧未返回该方法结果"
    if r_used and "ips" in methods and "ips" not in per_method:
        failed.setdefault("ips", "IPS 未成功计算")

    # ---- 2) 比例 / 相关 / 组间差异 ----
    prop, prop_meta = _build_proportions(per_method, normalize, warnings, notes)
    method_corr = _method_corr(per_method)
    wide = _concat_features(per_method)
    group_col, labels, groups = _group_context(ctx, info, params)
    group_frame = pd.DataFrame()
    if group_col:
        group_frame = _group_diff(prop, labels, groups)
        if group_frame.empty:
            warnings.append(f"分组列「{group_col}」每组有效样本不足 2 个，未做组间差异检验")
    else:
        notes.append("未指定分组列（group_col）或分组不足 2 组，跳过组间差异检验")

    # ---- 3) 图 ----
    figures: list[dict] = []
    try:
        figures.append(_plot_stack(prop, labels, groups, outdir))
    except Exception as error:  # noqa: BLE001 - 单图失败不影响结果
        warnings.append(f"immune_stack 绘图失败：{error}")
    if not method_corr.empty and method_corr.shape[0] >= 2:
        try:
            figures.append(_plot_corr(method_corr, outdir))
        except Exception as error:  # noqa: BLE001
            warnings.append(f"immune_corr 绘图失败：{error}")
    if not group_frame.empty:
        try:
            box = _plot_group_box(prop, group_frame, labels, groups, outdir)
            if box:
                figures.append(box)
        except Exception as error:  # noqa: BLE001
            warnings.append(f"immune_group_box 绘图失败：{error}")

    # ---- 4) 表 ----
    prop_records = _frame_records(prop)
    wide_records = _frame_records(wide)
    corr_records = _frame_records(method_corr, index_name="method")
    group_records = group_frame.to_dict(orient="records") if not group_frame.empty else []
    tables = [
        table(
            "immune_proportions",
            f"免疫细胞比例（{prop_meta['primary_method']}）",
            _columns_for(prop),
            prop_records,
            note=("相对比例，每行和≈1" if normalize == "relative" else "原始输出（absolute）")
            + ("；由得分相对化 derived" if prop_meta.get("derived") else ""),
            limit=200,
        ),
    ]
    if not group_frame.empty:
        tables.append(table("immune_group_diff", "免疫浸润组间差异", _plain_columns(group_frame),
                            group_records, note="Mann-Whitney U（2 组）/ Kruskal-Wallis（多组）+ BH FDR", limit=300))
    if not method_corr.empty:
        tables.append(table("immune_method_corr", "方法间相关矩阵", _columns_for(method_corr, index_name="method"),
                            corr_records, note="每个方法取特征矩阵第一主成分后计算 Pearson r", limit=100))
    if not wide.empty:
        tables.append(table("immune_scores", "各方法免疫打分（宽表）", _columns_for(wide),
                            wide_records, note="列名 = 方法::特征", limit=200))

    # ---- 5) 落盘 ----
    data_out = []
    try:
        path = export_csv(prop, outdir, "immune_proportions", index=True)
        data_out.append(bc.data_file("immune_csv", path, "免疫细胞比例表"))
        if not wide.empty:
            path = export_csv(wide, outdir, "immune_scores", index=True)
            data_out.append(bc.data_file("immune_csv", path, "各方法免疫打分宽表"))
        if not group_frame.empty:
            path = export_csv(group_frame, outdir, "immune_group_diff", index=False)
            data_out.append(bc.data_file("immune_csv", path, "免疫组间差异表"))
        if not method_corr.empty:
            path = export_csv(method_corr, outdir, "immune_method_corr", index=True)
            data_out.append(bc.data_file("immune_csv", path, "方法间相关矩阵"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"免疫结果 CSV 导出失败：{error}")

    # ---- 6) 摘要 ----
    methods_ok = list(per_method.keys())
    n_sig = 0
    if not group_frame.empty:
        n_sig = int((group_frame["fdr"] < 0.05).sum())
    bullets = [
        f"{'R/IOBR' if r_used else 'Python 回退'}成功 {len(methods_ok)} 个免疫浸润方法："
        + "、".join(methods_ok),
        f"细胞比例表：{prop.shape[0]} 样本 × {prop.shape[1]} 类（{prop_meta['primary_method']}，"
        + ("行和≈1" if normalize == "relative" else "absolute") + "）",
    ]
    if failed:
        bullets.append("失败/跳过：" + "、".join(sorted(failed.keys())))
    if not group_frame.empty:
        bullets.append(
            f"组间差异（{group_col}）{len(group_frame)} 个特征，FDR<0.05 共 {n_sig} 个"
        )
    metrics = {
        "n_methods": len(methods_ok),
        "methods_ok": methods_ok,
        "methods_failed": sorted(failed.keys()),
        "n_cells": int(prop.shape[1]),
        "n_samples": int(prop.shape[0]),
        "n_genes": int(len(genes)),
        "primary_method": prop_meta["primary_method"],
        "derived": bool(prop_meta.get("derived")),
        "normalize": normalize,
        "group_tests": int(group_frame.shape[0]) if not group_frame.empty else 0,
        "group_significant": n_sig,
        "engine": engine_label,
    }

    store = {
        "engine": engine_label,
        "requested": methods,
        "methods_ok": methods_ok,
        "methods_failed": failed,
        "normalize": normalize,
        "primary_method": prop_meta["primary_method"],
        "primary_label": prop_meta["primary_label"],
        "primary_derived": bool(prop_meta.get("derived")),
        "n_cells": int(prop.shape[1]),
        "n_samples": int(prop.shape[0]),
        "table": prop_records,
        "proportions": prop_records,
        "scores": wide_records,
        "method_corr": corr_records,
        "group_test": group_records,
        "group_col": group_col,
        "groups": groups,
        "notes": list(notes),
    }

    out = result(
        title="免疫浸润",
        bullets=bullets,
        metrics=metrics,
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=data_out,
    )
    out["store"] = {"immune": store}
    return out


# ---------------------------------------------------------------------------
# 表记录辅助
# ---------------------------------------------------------------------------

def _frame_records(frame: pd.DataFrame, index_name: str = "sample") -> list[dict]:
    if frame is None or not isinstance(frame, pd.DataFrame) or frame.empty:
        return []
    work = frame.copy()
    try:
        work.index = [str(i) for i in work.index]
    except Exception:  # noqa: BLE001
        pass
    work = work.reset_index().rename(columns={"index": index_name})
    return bc.jsonable(work.to_dict(orient="records"))


def _plain_columns(frame: pd.DataFrame) -> list[dict]:
    columns = []
    for col in frame.columns:
        numeric = pd.api.types.is_numeric_dtype(frame[col])
        columns.append({"key": str(col), "label": str(col), "type": "number" if numeric else "string"})
    return columns


def _columns_for(frame: pd.DataFrame, index_name: str | None = "sample") -> list[dict]:
    if not index_name:
        return _plain_columns(frame)
    columns: list[dict] = []
    if str(index_name) not in {str(c) for c in frame.columns}:
        columns.append({"key": index_name, "label": "样本" if index_name == "sample" else index_name,
                        "type": "string"})
    return columns + _plain_columns(frame)


def _stored_immune(ctx) -> dict | None:
    stored = ((ctx or {}).get("results") or {}).get("immune")
    return stored if isinstance(stored, dict) else None


def _store_frame(store: dict, key: str = "proportions") -> pd.DataFrame:
    value = store.get(key)
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, list) and value:
        return pd.DataFrame(value)
    # 兼容旧结构
    value = store.get("table")
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, list) and value:
        return pd.DataFrame(value)
    return pd.DataFrame()


def _store_scores(store: dict) -> pd.DataFrame:
    value = store.get("scores")
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, list) and value:
        return pd.DataFrame(value)
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# 分析实现
# ---------------------------------------------------------------------------

def _run_method_analysis(ctx, params, outdir, method: str, title: str) -> dict:
    payload = dict(params or {})
    payload["methods"] = [method]
    payload.setdefault("normalize", "relative")
    out = run_immune(ctx, payload, outdir)
    summary = out.get("summary") or {}
    summary["title"] = title
    out["summary"] = summary
    metrics = summary.get("metrics") or {}
    metrics["analysis_method"] = method
    summary["metrics"] = metrics
    return out


def _analysis_correlation(ctx, params, outdir) -> dict:
    store = _stored_immune(ctx)
    outdir = _outdir(ctx, outdir)
    if not store:
        raise EngineError("prerequisite", "还没有免疫浸润结果，请先运行「免疫浸润」步骤")
    scores = _store_scores(store)
    corr = pd.DataFrame()
    if not scores.empty:
        # 每方法 PC1 → 方法间相关
        by_method: dict[str, pd.DataFrame] = {}
        for col in scores.columns:
            key = str(col)
            if "::" not in key:
                continue
            method = key.split("::", 1)[0]
            by_method.setdefault(method, {})[key] = scores[col]
        if by_method:
            frames = {m: pd.DataFrame(cols) for m, cols in by_method.items()}
            corr = _method_corr(frames)
    if corr.empty:
        table_frame = _store_frame(store)
        if table_frame.empty or table_frame.shape[1] < 2:
            raise EngineError("prerequisite", "免疫结果没有足够特征计算相关矩阵")
        corr = table_frame.select_dtypes(include=[np.number]).corr().round(4)
        corr.index = [str(i) for i in corr.index]
    figure = _plot_corr(corr, outdir) if corr.shape[0] >= 2 else None
    records = _frame_records(corr, index_name="method")
    tables = [table("immune_method_corr", "方法间相关矩阵", _columns_for(corr, index_name="method"),
                    records, note="每个方法取特征矩阵第一主成分后计算 Pearson r", limit=100)]
    return result(
        title="免疫方法间相关性",
        bullets=[f"方法间相关矩阵：{corr.shape[0]} 个方法，r 范围 "
                 f"{np.nanmin(corr.values):.2f} ~ {np.nanmax(corr.values):.2f}"],
        metrics={"n_methods": int(corr.shape[0]), "min_r": _safe_float(np.nanmin(corr.values)),
                 "max_r": _safe_float(np.nanmax(corr.values))},
        tables=tables,
        figures=[figure] if figure else [],
        notes=["相关矩阵基于各方法 PC1 得分（避免相对比例行和恒定的伪相关）"],
        data_out=[bc.data_file("immune_csv", export_csv(corr, outdir, "immune_method_corr", index=True), "方法间相关矩阵")],
    ) | {"store": {"immune_corr": {"method_corr": records}}}


def _analysis_group_diff(ctx, params, outdir) -> dict:
    store = _stored_immune(ctx)
    outdir = _outdir(ctx, outdir)
    if not store:
        raise EngineError("prerequisite", "还没有免疫浸润结果，请先运行「免疫浸润」步骤")
    prop = _store_frame(store)
    if prop.empty:
        raise EngineError("prerequisite", "免疫结果为空，无法做组间差异")
    sample_col = "sample" if "sample" in prop.columns else None
    if sample_col:
        prop = prop.set_index(sample_col)
    prop = prop.select_dtypes(include=[np.number])
    info = (ctx or {}).get("sample_info")
    samples = [str(s) for s in ((ctx or {}).get("samples") or [])]
    if not isinstance(info, pd.DataFrame):
        raise EngineError("prerequisite", "缺少样本信息，无法做组间差异")
    info = info.reindex(samples)
    params = dict(params or {})
    if not params.get("group_col"):
        params["group_col"] = (ctx or {}).get("group_key")
    group_col, labels, groups = _group_context(ctx, info, params)
    if not group_col:
        raise EngineError("prerequisite", "请指定至少含 2 组的分组列（group_col）")
    prop = prop.reindex([s for s in samples if s in prop.index])
    labels = [lab for s, lab in zip(samples, labels) if s in prop.index]
    frame = _group_diff(prop, labels, groups)
    if frame.empty:
        raise EngineError("prerequisite", f"分组列「{group_col}」每组有效样本不足 2 个，无法检验")
    figure = None
    try:
        figure = _plot_group_box(prop, frame, labels, groups, outdir)
    except Exception:  # noqa: BLE001
        figure = None
    n_sig = int((frame["fdr"] < 0.05).sum())
    top = frame.head(5)
    bullets = [
        f"组间差异（{group_col}，{len(groups)} 组）：{len(frame)} 个免疫特征，FDR<0.05 共 {n_sig} 个",
        "最显著：" + "、".join(f"{r['feature']}(FDR={r['fdr']:.2g})" for _, r in top.iterrows()) if n_sig else "无显著差异",
    ]
    payload = result(
        title="免疫浸润组间差异",
        bullets=bullets,
        metrics={"n_features": int(len(frame)), "n_significant": n_sig, "group_col": group_col,
                 "groups": groups},
        tables=[table("immune_group_diff", "免疫浸润组间差异", _plain_columns(frame),
                      frame.to_dict(orient="records"), note="Mann-Whitney U / Kruskal-Wallis + BH FDR", limit=300)],
        figures=[figure] if figure else [],
        notes=["检验为双侧；多组时用 Kruskal-Wallis。免疫比例应结合临床分组解释，相关≠因果。"],
        data_out=[bc.data_file("immune_csv", export_csv(frame, outdir, "immune_group_diff", index=False), "组间差异表")],
    )
    return payload | {"store": {"immune_group": {"group_test": frame.to_dict(orient="records"),
                                                 "group_col": group_col, "groups": groups}}}


def _analysis_tme_cluster(ctx, params, outdir) -> dict:
    """基于免疫比例 KMeans 的 TME 分型（探索性）。"""
    store = _stored_immune(ctx)
    outdir = _outdir(ctx, outdir)
    if not store:
        raise EngineError("prerequisite", "还没有免疫浸润结果，请先运行「免疫浸润」步骤")
    prop = _store_frame(store)
    if prop.empty:
        raise EngineError("prerequisite", "免疫结果为空，无法做 TME 聚类")
    if "sample" in prop.columns:
        prop = prop.set_index("sample")
    prop = prop.select_dtypes(include=[np.number]).dropna(axis=1, how="all")
    if prop.shape[1] < 2:
        raise EngineError("prerequisite", "免疫特征不足 2 列，无法聚类")
    samples = [str(s) for s in prop.index]
    n = len(samples)
    if n < 8:
        raise EngineError("prerequisite", f"样本量不足（{n} < 8），TME 聚类不稳定")
    k = int((params or {}).get("k") or 0)
    values = prop.values.astype(float)
    values = np.nan_to_num(values, nan=float(np.nanmedian(values)))
    mean = values.mean(axis=0)
    std = values.std(axis=0, ddof=1)
    std = np.where(std > 1e-12, std, 1.0)
    z = (values - mean) / std

    k_range = [k] if k >= 2 else list(range(2, min(5, n // 4 if n >= 16 else 3) + 1) or [2])
    best = None
    try:
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score

        for kk in k_range:
            model = KMeans(n_clusters=kk, n_init=20, random_state=0)
            labels = model.fit_predict(z)
            score = float(silhouette_score(z, labels)) if len(set(labels)) > 1 else -1.0
            if best is None or score > best[0]:
                best = (score, kk, labels)
        if best is None:  # pragma: no cover
            raise RuntimeError("KMeans 无可用结果")
    except Exception as error:  # noqa: BLE001 - sklearn 缺失/失败 → 层次聚类回退
        try:
            from scipy.cluster.hierarchy import fcluster, linkage

            kk = k if k >= 2 else 2
            labels = fcluster(linkage(z, method="ward"), t=kk, criterion="maxclust")
            best = (-1.0, kk, labels)
        except Exception:  # noqa: BLE001
            raise EngineError("internal", f"TME 聚类失败：{error}")
    score, kk, labels = best
    labels = [f"TME-{int(lab) + 1}" for lab in labels]
    counts = {name: labels.count(name) for name in dict.fromkeys(labels)}
    # 各簇平均比例
    cluster_frame = pd.DataFrame(prop.values, index=samples, columns=prop.columns)
    cluster_frame["tme_cluster"] = labels
    means = cluster_frame.groupby("tme_cluster").mean(numeric_only=True).round(4)
    means.insert(0, "n", [counts.get(idx, 0) for idx in means.index])

    # 热图：免疫比例 z-score，按簇排序
    order = np.argsort([int(lab.split("-")[1]) for lab in labels], kind="stable")
    fig, ax = bc.new_fig(9.6, 5.2)
    im = ax.imshow(z[order].T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(prop.shape[1]))
    ax.set_yticklabels([str(c)[:32] for c in prop.columns], fontsize=6.5)
    ax.set_xticks([])
    ax.set_xlabel("Samples (sorted by TME cluster)")
    ax.set_title(f"TME clusters (k={kk}, silhouette={score:.2f})")
    ax.grid(False)
    boundaries = []
    running = 0
    for idx in means.index:
        running += counts.get(idx, 0)
        boundaries.append(running)
    for boundary in boundaries[:-1]:
        ax.axvline(boundary - 0.5, color="#2b2f36", linewidth=0.8, linestyle="--")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="z-score")
    figure = save_fig(fig, outdir, "tme_cluster_heat", "TME clusters")

    cluster_table = means.reset_index()
    bullets = [
        f"TME 聚类：{kk} 个簇（{', '.join(f'{k}:{v}' for k, v in counts.items())}）",
        "簇内平均免疫比例见「TME 簇均值」表；分型为探索性结果，不替代已发表的分子分型",
    ]
    payload = result(
        title="TME 免疫分型聚类",
        bullets=bullets,
        metrics={"k": int(kk), "silhouette": _safe_float(score), "n_samples": len(samples),
                 "cluster_sizes": counts},
        tables=[
            table("tme_cluster", "TME 簇均值", _columns_for(cluster_table, index_name="tme_cluster"),
                  cluster_table.to_dict(orient="records"), note="各 TME 簇的平均免疫细胞比例", limit=50),
            table("tme_samples", "样本 TME 分型", _columns_for(cluster_frame, index_name="sample"),
                  _frame_records(cluster_frame, index_name="sample"), note="每个样本的 TME 簇分配", limit=300),
        ],
        figures=[figure],
        notes=["KMeans（z-score 免疫比例，固定 seed=0）或层次聚类回退；簇数选择用 silhouette。"],
        sample_columns=[{"key": "tme_cluster", "kind": "categorical"}],
        sample_values={"tme_cluster": labels},
        data_out=[bc.data_file("immune_csv", export_csv(cluster_frame, outdir, "tme_clusters", index=True), "TME 分型表")],
    )
    return payload | {"store": {"tme_cluster": {"k": int(kk), "labels": dict(zip(samples, labels)),
                                                "silhouette": _safe_float(score)}}}


# ---------------------------------------------------------------------------
# STEPS / ANALYSES
# ---------------------------------------------------------------------------

_METHOD_OPTIONS = [
    {"value": "cibersort", "label": "CIBERSORT（LM22 相对比例）"},
    {"value": "epic", "label": "EPIC（比例）"},
    {"value": "mcpcounter", "label": "MCP-counter（marker 打分）"},
    {"value": "xcell", "label": "xCell（ssGSEA 类打分）"},
    {"value": "quantiseq", "label": "quanTIseq（比例）"},
    {"value": "estimate", "label": "ESTIMATE（免疫/基质/纯度）"},
    {"value": "ssgsea", "label": "ssGSEA（Charoentong 28 类）"},
    {"value": "timer", "label": "TIMER（marker 打分）"},
    {"value": "ips", "label": "IPS（Immunophenoscore）"},
]

_COMMON_PARAMS = [
    {"key": "group_col", "label": "分组列", "type": "column", "default": None},
    {"key": "normalize", "label": "数值归一化", "type": "select", "default": "relative",
     "options": [
         {"value": "relative", "label": "相对比例（每行和≈1）"},
         {"value": "absolute", "label": "原始输出（absolute）"},
     ]},
]

STEPS = [
    {
        "id": "immune",
        "title": "免疫浸润",
        "icon": "immune",
        "requires": ["normalize"],
        "params": [
            {"key": "methods", "label": "方法", "type": "multiselect",
             "default": ["cibersort", "epic", "mcpcounter", "xcell", "quantiseq", "estimate", "ssgsea", "timer"],
             "options": _METHOD_OPTIONS},
            {"key": "group_col", "label": "分组列", "type": "column", "default": "group"},
            {"key": "normalize", "label": "数值归一化", "type": "select", "default": "relative",
             "options": [
                 {"value": "relative", "label": "相对比例（每行和≈1）"},
                 {"value": "absolute", "label": "原始输出（absolute）"},
             ]},
            {"key": "perm", "label": "CIBERSORT 置换次数", "type": "number", "default": 100, "min": 0, "max": 1000, "step": 50},
        ],
        "run": run_immune,
    }
]


def _single_method_analysis(method: str, title: str) -> dict:
    def _run(ctx, params=None, outdir=None):
        return _run_method_analysis(ctx, params, outdir, method, title)

    return {
        "id": f"immune_{method}",
        "title": title,
        "group": "免疫浸润",
        "method": f"{METHOD_LABELS.get(method, method)}（R/IOBR 官方实现优先，Python 回退见 notes）",
        "needs": ["normalize"],
        "params": _COMMON_PARAMS + [
            {"key": "methods", "label": "方法", "type": "multiselect", "default": [method],
             "options": _METHOD_OPTIONS},
        ],
        "run": _run,
    }


ANALYSES = [
    _single_method_analysis("cibersort", "CIBERSORT 免疫细胞比例"),
    _single_method_analysis("epic", "EPIC 免疫/基质比例"),
    _single_method_analysis("mcpcounter", "MCP-counter 免疫群体打分"),
    _single_method_analysis("xcell", "xCell 细胞类型富集打分"),
    _single_method_analysis("quantiseq", "quanTIseq 免疫细胞比例"),
    _single_method_analysis("estimate", "ESTIMATE 免疫/基质评分"),
    _single_method_analysis("ssgsea", "ssGSEA（28 类免疫细胞）"),
    _single_method_analysis("timer", "TIMER 免疫浸润打分"),
    {
        "id": "immune_all",
        "title": "免疫浸润（全方法汇总）",
        "group": "免疫浸润",
        "method": "8 种 IOBR/GSVA 方法 + Python 回退，方法间相关 + 组间差异",
        "needs": ["normalize"],
        "params": [
            {"key": "methods", "label": "方法", "type": "multiselect",
             "default": ["cibersort", "epic", "mcpcounter", "xcell", "quantiseq", "estimate", "ssgsea", "timer"],
             "options": _METHOD_OPTIONS},
            {"key": "group_col", "label": "分组列", "type": "column", "default": None},
            {"key": "normalize", "label": "数值归一化", "type": "select", "default": "relative",
             "options": [
                 {"value": "relative", "label": "相对比例"},
                 {"value": "absolute", "label": "原始输出"},
             ]},
        ],
        "run": run_immune,
    },
    {
        "id": "immune_correlation",
        "title": "免疫方法间相关性",
        "group": "免疫浸润",
        "method": "各方法特征矩阵第一主成分的 Pearson 相关（避免相对比例的伪相关）",
        "needs": ["immune"],
        "params": [],
        "run": _analysis_correlation,
    },
    {
        "id": "immune_group_diff",
        "title": "免疫浸润组间差异",
        "group": "免疫浸润",
        "method": "Mann-Whitney U（2 组）/ Kruskal-Wallis（多组）+ BH FDR",
        "needs": ["immune"],
        "params": [
            {"key": "group_col", "label": "分组列", "type": "column", "default": None},
        ],
        "run": _analysis_group_diff,
    },
    {
        "id": "tme_cluster",
        "title": "TME 免疫分型聚类",
        "group": "免疫浸润",
        "method": "免疫比例 z-score + KMeans（轮廓系数选 k，层次聚类回退）",
        "needs": ["immune"],
        "params": [
            {"key": "k", "label": "簇数（0=自动）", "type": "number", "default": 0, "min": 0, "max": 8, "step": 1},
        ],
        "run": _analysis_tme_cluster,
    },
]
