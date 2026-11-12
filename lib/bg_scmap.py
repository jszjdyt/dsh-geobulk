"""dsh-geobulk — 单细胞映射（bulk-ncrna 负责）

步骤 `scmap`：把 bulk 的差异基因/自定义签名映射到单细胞参考（h5ad）的细胞亚群上。
参考来源优先级：`params.ref_path` → `ctx["sc_ref"]`（Host 由 `/geobulk/sc/attach` 写入）。

打分方法：
  * mean      —— 每个基因在参考细胞内 z-score 后，亚群内取签名基因均值（再跨亚群归一化到 [0,1]）
  * aucell    —— rank-AUC 简化版，AUC = 1 - Σrank(签名基因)/(k·N)，值域 [0,1]
  * addmodule —— Seurat AddModuleScore 式：签名均值 - 50 组表达量匹配随机对照基因集的背景均值

输出：①细胞类型×得分/细胞数 ②目标分子在各亚群的平均表达+表达比例 ③细胞类型×bulk 样本相关矩阵
图：sc_gene_violin（签名得分小提琴）、sc_target_violin（目标基因表达）、sc_corr_heat（相关热图）

约定：本文只依赖标准库 + numpy/pandas/scipy/matplotlib/anndata + bg_common（契约 §3.3）。
只读参考数据，不改写参考的细胞类型注释。
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pandas as pd

import bg_common as bc
from bg_common import EngineError, export_csv, palette, result, save_fig, table

__all__ = ["STEPS", "ANALYSES", "run_scmap"]

SCORE_METHODS = ("mean", "aucell", "addmodule")

# obs 中细胞类型列的优先顺序（小写比较）
CELL_TYPE_COLUMNS = [
    "cell_type",
    "celltype",
    "cell_types",
    "cell_type_annotation",
    "annotation",
    "majorcelltype",
    "major_cell_type",
    "clusters",
    "cluster",
    "leiden",
    "seurat_clusters",
    "louvain",
    "subcluster",
]

# bulk 差异表里可能的列名（大小写合并后匹配）
_DE_GENE_COLS = ("gene", "gene_name", "symbol", "genes", "gene_id", "feature")
_DE_LFC_COLS = ("logfc", "log2fc", "log2_fold_change", "lfc", "foldchange")
_DE_PADJ_COLS = ("adj.p.val", "padj", "fdr", "qvalue", "q_value", "p.adj", "p_adjusted", "adjp")
_DE_PVAL_COLS = ("p.value", "pvalue", "p_val", "p", "pval")


def _outdir(ctx, outdir) -> str:
    path = outdir or (ctx or {}).get("outdir")
    if not path:
        path = tempfile.mkdtemp(prefix="bg-scmap-")
    os.makedirs(path, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# 参数 / 参考路径
# ---------------------------------------------------------------------------

def _ctx_sc_ref(ctx) -> str:
    """从 ctx 里取 Host 写入的 sc_ref（可能是字符串或 dict）。"""
    sc = (ctx or {}).get("sc_ref")
    if isinstance(sc, str) and sc.strip():
        return sc.strip()
    if isinstance(sc, dict):
        for key in ("path", "ref_path", "file", "h5ad", "h5"):
            val = sc.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
    return ""


def _resolve_ref_path(ctx, params) -> tuple[str, str]:
    """返回 (绝对路径, 来源说明)；缺参考 → prerequisite，文件不存在 → not_found。"""
    p = str((params or {}).get("ref_path") or "").strip()
    source = "params.ref_path"
    if not p:
        p = _ctx_sc_ref(ctx)
        source = "dataset.sc_ref"
    if not p:
        raise EngineError(
            "prerequisite",
            "请先指定单细胞参考 h5ad（可在「单细胞」tab 选择本机文件或从单细胞工作台拉取）",
        )
    p = os.path.abspath(os.path.expanduser(p))
    if not os.path.exists(p):
        raise EngineError("not_found", f"单细胞参考文件不存在：{p}")
    if os.path.isdir(p):
        raise EngineError("bad_request", f"单细胞参考应为 .h5ad 文件，而不是目录：{p}")
    if not p.lower().endswith((".h5ad", ".h5")):
        raise EngineError("bad_request", f"单细胞参考只支持 .h5ad（anndata）文件：{p}")
    return p, source


def _as_gene_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        parts = [x for x in value.replace(",", " ").replace(";", " ").split() if x]
        return parts
    try:
        return [str(x).strip() for x in value if str(x).strip()]
    except TypeError:
        return [str(value)]


# ---------------------------------------------------------------------------
# 读取 h5ad 参考
# ---------------------------------------------------------------------------

def _read_h5ad(path: str):
    """读 h5ad；大文件用 backed 模式。返回 (adata, backed)。"""
    try:
        import anndata as ad
    except Exception as exc:  # noqa: BLE001
        raise EngineError("unsupported", f"当前 Python 缺少 anndata，无法读取单细胞参考：{exc}")
    size = os.path.getsize(path)
    use_backed = size > 350 * 1024 * 1024
    try:
        adata = ad.read_h5ad(path, backed="r") if use_backed else ad.read_h5ad(path)
    except MemoryError:
        adata = ad.read_h5ad(path, backed="r")
        use_backed = True
    except Exception as exc:  # noqa: BLE001
        raise EngineError(
            "internal", f"无法读取单细胞参考 {os.path.basename(path)}：{exc}"
        ) from exc
    return adata, use_backed


def _celltype_column(obs: pd.DataFrame, requested) -> str:
    if obs is None or len(obs.columns) == 0:
        raise EngineError("prerequisite", "单细胞参考的 obs 为空，无法确定细胞类型列")
    if requested:
        req = str(requested).strip()
        for col in obs.columns:
            if str(col) == req:
                return str(col)
        for col in obs.columns:
            if str(col).lower() == req.lower():
                return str(col)
        raise EngineError(
            "bad_request",
            f"参考数据的 obs 里没有列「{req}」；可用列：{', '.join(map(str, obs.columns))[:300]}",
        )
    lower = {str(c).lower(): str(c) for c in obs.columns}
    for cand in CELL_TYPE_COLUMNS:
        if cand in lower:
            return lower[cand]
    raise EngineError(
        "prerequisite",
        "参考 h5ad 的 obs 中找不到细胞类型列（cell_type/clusters/leiden/seurat_clusters 等）；"
        "请用 cell_type_col 参数指定",
    )


def _celltype_values(series: pd.Series) -> np.ndarray:
    vals = []
    for val in series.astype(object):
        if val is None:
            vals.append("unassigned")
            continue
        try:
            if isinstance(val, float) and np.isnan(val):
                vals.append("unassigned")
                continue
        except (TypeError, ValueError):
            pass
        text = str(val).strip()
        vals.append(text if text and text.lower() not in ("nan", "none", "na") else "unassigned")
    return np.asarray(vals, dtype=object)


def _first_attr(value) -> str:
    """anndata 的字符串列可能是 categorical / 数组 / str。"""
    try:
        if isinstance(value, str):
            return value
        if hasattr(value, "iloc"):
            return str(value.iloc[0] if len(value) else "")
        return str(value[0]) if len(value) else ""
    except Exception:  # noqa: BLE001
        return ""


def _subset_genes(adata, wanted_idx: list[int], ref_genes: list[str], path: str):
    """从 backed 参考里只读指定基因列；返回 (X float32, genes)。"""
    idx = np.asarray(sorted(set(int(i) for i in wanted_idx)), dtype=int)
    if idx.size == 0:
        raise EngineError("internal", "没有可读取的基因列")
    # 1) anndata 切片
    try:
        sub = adata[:, idx]
        X = sub.X
        if hasattr(X, "to_memory"):
            X = X.to_memory().X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        if X.ndim == 2 and X.shape[1] == idx.size:
            return X, [ref_genes[i] for i in idx]
    except Exception:  # noqa: BLE001
        pass
    # 2) h5py 直接读 dense X
    try:
        import h5py

        with h5py.File(path, "r") as handle:
            node = handle["X"]
            if isinstance(node, h5py.Dataset):
                X = np.asarray(node[:, idx], dtype=np.float32)
                return X, [ref_genes[i] for i in idx]
    except Exception as exc:  # noqa: BLE001
        raise EngineError(
            "unsupported",
            f"参考 h5ad 的 X 为稀疏/特殊存储，backed 模式下无法按基因切片：{exc}",
        ) from exc
    raise EngineError("unsupported", "参考 h5ad 的 X 存储格式暂不支持（既非 dense 也非可切片稀疏）")


def _read_expression(adata, backed: bool, ref_genes: list[str], wanted_idx: list[int] | None, path: str):
    """返回 (X cells×genes float32, genes_subset, approximate 提示或空串)。"""
    n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
    if not backed:
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        if X.ndim != 2:
            raise EngineError("internal", "参考 h5ad 的 X 不是二维矩阵")
        return X, list(ref_genes), ""
    # backed：能整块读就整块读（便于 AUCell 与 marker 分析）
    est = float(n_obs) * float(n_vars) * 4.0
    if est <= 300 * 1024 * 1024:
        try:
            X = adata.X[:]
            if hasattr(X, "toarray"):
                X = X.toarray()
            if hasattr(X, "to_memory"):
                X = X.to_memory().X
            X = np.asarray(X, dtype=np.float32)
            if X.ndim == 2 and X.shape == (n_obs, n_vars):
                return X, list(ref_genes), ""
        except Exception:  # noqa: BLE001
            pass
    if wanted_idx is None or len(wanted_idx) == 0:
        # 没有签名基因时取前 min(2000, n_vars) 列做描述性分析
        wanted_idx = list(range(min(2000, n_vars)))
        note = f"参考以 backed 模式读取，仅载入前 {len(wanted_idx)} 个基因列"
    else:
        note = "参考以 backed 模式读取，仅载入签名基因所在列（AUCell/addModule 为近似值）"
    X, genes = _subset_genes(adata, wanted_idx, ref_genes, path)
    return X, genes, note


def _maybe_log(X: np.ndarray, warnings: list[str]) -> np.ndarray:
    """X 看起来像原始 counts（max>50 且非负）时做 log1p，避免均值打分被大值主导。"""
    if X.size == 0:
        return X
    finite = X[np.isfinite(X)]
    if finite.size == 0:
        return X
    vmax = float(np.nanmax(finite))
    if vmax > 50.0 and float(np.nanmin(finite)) >= 0.0:
        warnings.append("参考表达矩阵数值范围较大（max>50），已按 log1p 处理后再打分")
        return np.log1p(X)
    return X


# ---------------------------------------------------------------------------
# 基因匹配
# ---------------------------------------------------------------------------

def _build_lookup(genes: list[str]) -> dict:
    lookup: dict[str, int] = {}
    for i, gene in enumerate(genes):
        g = str(gene)
        lookup.setdefault(g, i)
        lookup.setdefault(g.upper(), i)
        lookup.setdefault(g.split(".")[0].upper(), i)
    return lookup


def _match_signature(sig_genes: list[str], ref_genes: list[str]) -> tuple[list[tuple[str, int]], list[str]]:
    """把 bulk 基因匹配到参考行下标（symbol / 大小写 / 去 Ensembl 版本号）。返回 [(bulk_gene, ref_idx)], 未命中基因。"""
    lookup = _build_lookup(ref_genes)
    matched: list[tuple[str, int]] = []
    missing: list[str] = []
    seen: set[int] = set()
    for gene in sig_genes:
        g = str(gene).strip()
        if not g:
            continue
        idx = lookup.get(g)
        if idx is None:
            idx = lookup.get(g.upper())
        if idx is None:
            idx = lookup.get(g.split(".")[0].upper())
        if idx is None:
            missing.append(g)
            continue
        if idx in seen:
            continue
        seen.add(idx)
        matched.append((g, int(idx)))
    return matched, missing


# ---------------------------------------------------------------------------
# 差异结果解析（ctx["de"]）
# ---------------------------------------------------------------------------

def _de_table(ctx) -> pd.DataFrame | None:
    de = (ctx or {}).get("de")
    if de is None:
        return None
    if isinstance(de, pd.DataFrame):
        return de
    if isinstance(de, dict):
        frame = de.get("table")
        if isinstance(frame, pd.DataFrame):
            return frame
        for key in ("de_all", "de", "results"):
            frame = de.get(key)
            if isinstance(frame, pd.DataFrame):
                return frame
        tables = de.get("tables")
        if isinstance(tables, list):
            for item in tables:
                if not isinstance(item, dict):
                    continue
                if str(item.get("id")) in ("de", "de_all") and item.get("rows"):
                    try:
                        return pd.DataFrame(item["rows"])
                    except Exception:  # noqa: BLE001
                        continue
        if de.get("rows"):
            try:
                return pd.DataFrame(de["rows"])
            except Exception:  # noqa: BLE001
                return None
        for value in de.values():
            if isinstance(value, pd.DataFrame):
                return value
    return None


def _find_column(frame: pd.DataFrame, candidates) -> str | None:
    lower = {str(c).lower().replace(" ", ""): str(c) for c in frame.columns}
    for cand in candidates:
        key = cand.lower().replace(" ", "")
        if key in lower:
            return lower[key]
    return None


def _de_gene_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """标准化差异表：gene / logFC / padj。"""
    gene_col = _find_column(frame, _DE_GENE_COLS)
    if gene_col is None:
        # 第一列当基因名
        gene_col = str(frame.columns[0])
    lfc_col = _find_column(frame, _DE_LFC_COLS)
    padj_col = _find_column(frame, _DE_PADJ_COLS)
    pval_col = _find_column(frame, _DE_PVAL_COLS)
    out = pd.DataFrame({"gene": frame[gene_col].astype(str)})
    out["logFC"] = pd.to_numeric(frame[lfc_col], errors="coerce") if lfc_col else 0.0
    if padj_col:
        out["padj"] = pd.to_numeric(frame[padj_col], errors="coerce")
    elif pval_col:
        out["padj"] = bc.bh_fdr(pd.to_numeric(frame[pval_col], errors="coerce").fillna(1.0).to_numpy())
    else:
        out["padj"] = 1.0
    out["padj"] = out["padj"].fillna(1.0)
    out["logFC"] = out["logFC"].fillna(0.0)
    return out


def _signature_genes(ctx, params, warnings: list[str]) -> tuple[list[str], str, pd.DataFrame | None]:
    """按 gene_set 取签名基因。返回 (基因列表, 中文说明, 标准化差异表)。"""
    gene_set = str((params or {}).get("gene_set") or "de_up").strip().lower()
    max_genes = int((params or {}).get("max_genes") or 50)
    max_genes = max(3, min(max_genes, 500))
    de_frame = _de_table(ctx)
    de_std = _de_gene_frame(de_frame) if de_frame is not None and len(de_frame) else None

    if gene_set == "custom":
        genes = _as_gene_list((params or {}).get("genes"))
        if not genes:
            raise EngineError("bad_request", "gene_set=custom 时请用 genes 参数给出基因列表（list 或逗号分隔字符串）")
        return genes[:max_genes], f"自定义基因列表（{len(genes[:max_genes])} 个）", de_std

    if gene_set not in ("de_up", "de_down", "de_both"):
        raise EngineError("bad_request", f"不支持的 gene_set：{gene_set}（可选 de_up/de_down/de_both/custom）")

    if de_std is None:
        raise EngineError(
            "prerequisite",
            "尚未运行差异表达（de）步骤，无法用 DE 基因做单细胞签名；请先运行 de，或把 gene_set 设为 custom 并给出 genes",
        )

    sig = de_std[(de_std["padj"] <= 0.05) & (de_std["logFC"].abs() > 0)]
    if sig.empty:
        warnings.append("差异表中没有 FDR<0.05 的显著基因，退化为按 |logFC| 排序取 top 基因作为签名")
        sig = de_std.assign(_abs=de_std["logFC"].abs()).sort_values("_abs", ascending=False).head(max_genes)

    up = sig[sig["logFC"] > 0].sort_values(["padj", "logFC"], ascending=[True, False])
    down = sig[sig["logFC"] < 0].sort_values(["padj", "logFC"], ascending=[True, True])
    if gene_set == "de_up":
        genes = up["gene"].head(max_genes).tolist()
        label = f"DE 上调基因（FDR<0.05，取前 {len(genes)} 个）"
    elif gene_set == "de_down":
        genes = down["gene"].head(max_genes).tolist()
        label = f"DE 下调基因（FDR<0.05，取前 {len(genes)} 个）"
    else:
        half = max(3, max_genes // 2)
        genes = up["gene"].head(half).tolist() + down["gene"].head(half).tolist()
        label = f"DE 上调+下调基因（各取前 {half} 个，共 {len(genes)} 个）"
    if not genes:
        raise EngineError(
            "prerequisite",
            f"差异表里没有符合条件的基因（{label}）；请检查 de 结果或改用 gene_set=custom",
        )
    return genes, label, de_std


# ---------------------------------------------------------------------------
# 打分
# ---------------------------------------------------------------------------

def _zscore_columns(X: np.ndarray) -> np.ndarray:
    mean = np.nanmean(X, axis=0, keepdims=True)
    std = np.nanstd(X, axis=0, keepdims=True)
    std = np.where(np.isfinite(std) & (std > 0), std, 1.0)
    out = (X - mean) / std
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def _aucell(X: np.ndarray, sig_idx: np.ndarray, warnings: list[str]) -> np.ndarray:
    """rank-AUC 简化：AUC = 1 - Σposition/(k·N)，position 为 1-based 表达降序名次。"""
    n_cells, n_genes = X.shape
    k = int(len(sig_idx))
    if k == 0 or n_genes == 0:
        return np.zeros(n_cells)
    sig_idx = np.asarray(sig_idx, dtype=int)
    if float(n_cells) * float(n_genes) <= 4.0e7:
        mat = X
        n = n_genes
        cols = sig_idx
    else:
        # 太大：在签名基因 + 4000 个随机背景基因上近似排名（note 由调用方写入 warnings）
        rng = np.random.default_rng(0)
        bg = rng.choice(n_genes, size=min(n_genes, 4000), replace=False)
        keep = np.unique(np.concatenate([sig_idx, bg]))
        mat = X[:, keep]
        n = int(len(keep))
        cols = np.where(np.isin(keep, sig_idx))[0]
        warnings.append("参考矩阵较大，AUCell 在签名基因+随机背景基因上近似计算")
    if cols.size == 0 or n == 0:
        return np.zeros(n_cells)
    order = np.argsort(-mat, axis=1, kind="stable")
    ranks = np.empty_like(order)
    np.put_along_axis(ranks, order, np.broadcast_to(np.arange(n), order.shape), axis=1)
    rank_sum = ranks[:, cols].sum(axis=1) + cols.size  # 转 1-based 名次
    auc = 1.0 - rank_sum / (cols.size * n)
    return np.clip(auc, 0.0, 1.0)


def _addmodule(X: np.ndarray, sig_idx: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Seurat AddModuleScore 式：签名均值 − 50 组表达量匹配随机对照集的背景均值。"""
    n_cells, n_genes = X.shape
    k = int(len(sig_idx))
    if k == 0:
        return np.zeros(n_cells)
    n_ctrl = 50
    if n_genes < 2 * k + 4:
        return bc.mean_top_score(X, sig_idx) - float(np.mean(X))
    sig_idx = np.asarray(sig_idx, dtype=int)
    mu = np.nanmean(X, axis=0)
    nbins = int(min(24, max(2, n_genes // max(1, 2 * k))))
    order = np.argsort(mu)
    bin_of = np.empty(n_genes, dtype=int)
    for b, chunk in enumerate(np.array_split(order, nbins)):
        bin_of[chunk] = b
    sig_set = np.zeros(n_genes, dtype=bool)
    sig_set[sig_idx] = True
    ctrl_means = np.zeros((n_cells, n_ctrl), dtype=np.float64)
    for c in range(n_ctrl):
        picks = np.empty(k, dtype=int)
        for j, gene in enumerate(sig_idx):
            cand = np.where((bin_of == bin_of[gene]) & (~sig_set))[0]
            if cand.size == 0:
                cand = np.where(~sig_set)[0]
            picks[j] = int(cand[rng.integers(0, cand.size)] if cand.size else gene)
        ctrl_means[:, c] = X[:, picks].mean(axis=1)
    sig_mean = X[:, sig_idx].mean(axis=1)
    return sig_mean - ctrl_means.mean(axis=1)


def _score_cells(prep: dict, method: str, warnings: list[str]) -> np.ndarray:
    X = prep["X"]
    sig_idx = np.asarray([idx for _, idx in prep["matched"]], dtype=int)
    if method == "mean":
        return bc.mean_top_score(_zscore_columns(X), sig_idx)
    if method == "aucell":
        return _aucell(X, sig_idx, warnings)
    if method == "addmodule":
        rng = np.random.default_rng(0)
        return _addmodule(X, sig_idx, rng)
    raise EngineError("bad_request", f"不支持的打分方法：{method}")


def _aggregate_celltypes(scores: np.ndarray, celltypes: np.ndarray) -> pd.DataFrame:
    types = list(pd.unique(pd.Series(celltypes, dtype=object)))
    rows = []
    for ct in types:
        mask = celltypes == ct
        n = int(mask.sum())
        raw = float(np.nanmean(scores[mask])) if n else float("nan")
        rows.append({"cell_type": str(ct), "n_cells": n, "raw_score": raw,
                     "median_score": float(np.nanmedian(scores[mask])) if n else float("nan")})
    frame = pd.DataFrame(rows)
    total = max(1, len(celltypes))
    frame["pct_cells"] = (frame["n_cells"] / total * 100.0).round(2)
    if len(frame) > 1:
        frame["score"] = bc.min_max_scale(frame["raw_score"].to_numpy())
    else:
        frame["score"] = 1.0
    frame = frame.sort_values(["score", "n_cells"], ascending=[False, False]).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    return frame


def _mean_expression_table(prep: dict, genes_limit: int = 20) -> tuple[pd.DataFrame, pd.DataFrame]:
    """目标分子 × 亚群 平均表达 / 表达比例（完整表 + 展示用截断表）。"""
    X = prep["X"]
    celltypes = prep["celltypes"]
    types_order = [str(t) for t in prep["score_frame"]["cell_type"].tolist()]
    matched = prep["matched"][:max(1, genes_limit)]
    rows = []
    for bulk_gene, ref_idx in matched:
        col = X[:, ref_idx]
        for ct in types_order:
            mask = celltypes == ct
            if not mask.any():
                continue
            vals = col[mask]
            rows.append({
                "gene": str(bulk_gene),
                "cell_type": ct,
                "mean_expr": round(float(np.nanmean(vals)), 4),
                "pct_expr": round(float(np.mean(vals > 0) * 100.0), 2),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 相关性
# ---------------------------------------------------------------------------

def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    aa, bb = a[mask], b[mask]
    aa = aa - aa.mean()
    bb = bb - bb.mean()
    denom = float(np.sqrt((aa * aa).sum() * (bb * bb).sum()))
    if denom <= 0:
        return float("nan")
    return float((aa * bb).sum() / denom)


def _celltype_bulk_corr(prep: dict, ctx, warnings: list[str]):
    """细胞类型 × bulk 样本 相关矩阵（签名基因上的表达向量相关）。"""
    bulk_genes = [str(g) for g in (ctx or {}).get("genes") or []]
    X_bulk = (ctx or {}).get("X")
    if not bulk_genes or X_bulk is None:
        warnings.append("ctx 缺少 bulk 表达矩阵，跳过细胞类型×bulk 样本相关分析")
        return None, 0
    X_bulk = np.asarray(X_bulk, dtype=np.float64)
    bulk_lookup = {g: i for i, g in enumerate(bulk_genes)}
    common: list[tuple[int, int]] = []  # (matched 位置, bulk 列)
    for pos, (bulk_gene, _ridx) in enumerate(prep["matched"]):
        col = bulk_lookup.get(str(bulk_gene))
        if col is None:
            col = bulk_lookup.get(str(bulk_gene).upper())
        if col is not None:
            common.append((pos, int(col)))
    if len(common) < 3:
        warnings.append(f"签名基因在 bulk 中只命中 {len(common)} 个（需 ≥3），跳过相关性分析")
        return None, len(common)
    X_ref = prep["X"]
    celltypes = prep["celltypes"]
    samples = [str(s) for s in (ctx or {}).get("samples") or [f"S{i + 1}" for i in range(X_bulk.shape[0])]]
    ref_pos = [p for p, _ in common]
    bulk_cols = [c for _, c in common]
    rows = []
    for ct in [str(t) for t in prep["score_frame"]["cell_type"].tolist()]:
        mask = celltypes == ct
        if not mask.any():
            continue
        ref_vec = X_ref[np.ix_(np.where(mask)[0], [idx for idx in [prep["matched"][p][1] for p in ref_pos]])].mean(axis=0)
        row = {"cell_type": ct}
        for si, sample in enumerate(samples):
            row[sample] = round(_pearson(ref_vec, X_bulk[si, bulk_cols]), 4)
        rows.append(row)
    frame = pd.DataFrame(rows)
    return frame, len(common)


# ---------------------------------------------------------------------------
# 图
# ---------------------------------------------------------------------------

def _violin_by_group(ax, values: np.ndarray, labels: np.ndarray, order: list[str], colors: list[str]):
    data = [values[labels == ct] for ct in order]
    data = [d[np.isfinite(d)] for d in data]
    if any(len(d) == 0 for d in data):
        data = [d if len(d) else np.array([0.0]) for d in data]
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(colors[i % len(colors)])
        body.set_edgecolor("#5a6068")
        body.set_alpha(0.75)
    ax.set_xticks(np.arange(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_ylabel("Signature score")


def _fig_score_violin(prep, scores, outdir, method) -> dict:
    types = [str(t) for t in prep["score_frame"]["cell_type"].tolist()]
    colors = palette(len(types))
    fig, ax = bc.new_fig(8.6, 4.6)
    _violin_by_group(ax, scores, prep["celltypes"], types, colors)
    ax.set_title(f"Signature score by cell type ({method}, n={len(scores)} cells)")
    return save_fig(fig, outdir, "sc_gene_violin", title="Signature score by cell type")


def _fig_target_violin(prep, outdir, top_n: int) -> dict | None:
    matched = prep["matched"][:max(1, top_n)]
    if not matched:
        return None
    X = prep["X"]
    celltypes = prep["celltypes"]
    types = [str(t) for t in prep["score_frame"]["cell_type"].tolist()]
    fig, ax = bc.new_fig(9.2, 4.8)
    positions, data, tick_labels = [], [], []
    offset = 0
    colors = palette(len(types))
    for bulk_gene, ref_idx in matched:
        vals = X[:, ref_idx]
        for j, ct in enumerate(types):
            mask = celltypes == ct
            data.append(vals[mask])
            positions.append(offset + j + 1)
            tick_labels.append(ct)
        offset += len(types) + 1
    parts = ax.violinplot(data, positions=positions, showextrema=False, widths=0.85)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(colors[i % len(colors)])
        body.set_edgecolor("#5a6068")
        body.set_alpha(0.7)
    centers = [offset + (len(types) + 1) / 2.0 for offset in
               [k * (len(types) + 1) for k in range(len(matched))]]
    ax.set_xticks(centers)
    ax.set_xticklabels([g for g, _ in matched], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Expression (log scale)")
    ax.set_title(f"Target genes across cell types (n={len(matched)} genes)")
    return save_fig(fig, outdir, "sc_target_violin", title="Target gene expression by cell type")


def _fig_corr_heat(corr: pd.DataFrame | None, outdir) -> dict | None:
    if corr is None or corr.empty:
        return None
    mat = corr.set_index("cell_type").apply(pd.to_numeric, errors="coerce")
    if mat.shape[1] < 2:
        return None
    fig, ax = bc.new_fig(max(6.0, 0.32 * mat.shape[1] + 2.4), max(3.2, 0.55 * mat.shape[0] + 1.6))
    arr = mat.to_numpy(dtype=float)
    im = ax.imshow(arr, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels([str(c) for c in mat.columns], rotation=90, fontsize=7)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels([str(i) for i in mat.index], fontsize=8)
    ax.set_title("Cell type vs bulk sample correlation (signature genes)")
    ax.grid(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Pearson r")
    return save_fig(fig, outdir, "sc_corr_heat", title="Cell type vs bulk sample correlation")


def _fig_prop_bar(prop: pd.DataFrame, outdir) -> dict:
    frame = prop.sort_values("n_cells", ascending=True)
    fig, ax = bc.new_fig(7.2, max(3.2, 0.38 * len(frame) + 1.8))
    colors = palette(len(frame))
    ax.barh([str(x) for x in frame["cell_type"]], frame["pct_cells"].to_numpy(), color=colors)
    for i, val in enumerate(frame["pct_cells"].to_numpy()):
        ax.text(float(val) + 0.3, i, f"{val:.1f}%", va="center", fontsize=8)
    ax.set_xlabel("Cells (%)")
    ax.set_title("Cell type composition in the single-cell reference")
    return save_fig(fig, outdir, "sc_celltype_prop", title="Cell type composition")


# ---------------------------------------------------------------------------
# 参考准备
# ---------------------------------------------------------------------------

def _prepare(ctx, params, outdir, *, require_sig: bool = True, genes_override: list[str] | None = None) -> dict:
    params = dict(params or {})
    warnings: list[str] = []
    notes: list[str] = []
    ref_path, ref_source = _resolve_ref_path(ctx, params)
    raw_de = _de_table(ctx)
    de_std = _de_gene_frame(raw_de) if raw_de is not None and len(raw_de) else None
    if genes_override is not None:
        sig_genes = list(genes_override)
        sig_label = f"自定义基因列表（{len(sig_genes)} 个）"
    elif require_sig:
        sig_genes, sig_label, de_std = _signature_genes(ctx, params, warnings)
    else:
        sig_genes, sig_label = [], "（未使用签名基因）"

    method = str(params.get("score_method") or "mean").strip().lower()
    if method not in SCORE_METHODS:
        raise EngineError("bad_request", f"不支持的 score_method：{method}（可选 mean/aucell/addmodule）")

    adata, backed = _read_h5ad(ref_path)
    obs = adata.obs.copy()
    ref_genes = [str(g) for g in adata.var_names]
    n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
    if n_obs < 10:
        raise EngineError("prerequisite", f"单细胞参考细胞数过少（{n_obs}），无法做亚群映射")
    ct_col = _celltype_column(obs, params.get("cell_type_col"))
    celltypes = _celltype_values(obs[ct_col])
    matched, missing = _match_signature(sig_genes, ref_genes)
    if missing:
        warnings.append(
            f"签名基因里有 {len(missing)} 个在参考中未找到（例如：{', '.join(missing[:8])}）"
        )
    if require_sig and len(matched) == 0:
        raise EngineError(
            "prerequisite",
            "签名基因在单细胞参考里一个都没找到：请确认参考与 bulk 同物种、基因命名一致"
            "（symbol vs Ensembl ID），或改用 gene_set=custom 指定参考中存在的基因",
        )
    if require_sig and len(matched) < 3:
        warnings.append(f"签名基因仅命中 {len(matched)} 个，打分稳定性有限")
    wanted_idx = [idx for _, idx in matched]
    X, x_genes, approx_note = _read_expression(adata, backed, ref_genes, wanted_idx, ref_path)
    if approx_note:
        notes.append(approx_note)
    X = _maybe_log(X, warnings)
    # 重新定位 matched 到 x_genes（backed 子集时列顺序变化）
    if x_genes != ref_genes:
        x_lookup = {str(g): i for i, g in enumerate(x_genes)}
        rematched = []
        for bulk_gene, _old_idx in matched:
            pos = x_lookup.get(str(ref_genes[_old_idx]))
            if pos is not None:
                rematched.append((bulk_gene, int(pos)))
        matched = rematched
    if require_sig and len(matched) == 0:
        raise EngineError("internal", "签名基因定位失败（参考子集读取后无法匹配）")
    prep = {
        "ref_path": ref_path,
        "ref_source": ref_source,
        "ref_name": os.path.basename(ref_path),
        "n_obs": n_obs,
        "n_vars": n_vars,
        "backed": backed,
        "obs": obs,
        "ct_col": ct_col,
        "celltypes": celltypes,
        "ref_genes": ref_genes,
        "x_genes": x_genes,
        "X": X,
        "matched": matched,
        "missing": missing,
        "sig_genes": sig_genes,
        "sig_label": sig_label,
        "method": method,
        "de_frame": de_std,
        "warnings": warnings,
        "notes": notes,
        "approx": bool(approx_note),
        "outdir": _outdir(ctx, outdir),
    }
    return prep


def _score_frame(prep) -> pd.DataFrame:
    scores = _score_cells(prep, prep["method"], prep["warnings"])
    frame = _aggregate_celltypes(scores, prep["celltypes"])
    prep["scores"] = scores
    prep["score_frame"] = frame
    return frame


# ---------------------------------------------------------------------------
# 步骤 scmap
# ---------------------------------------------------------------------------

def run_scmap(ctx, params=None, outdir=None):
    params = dict(params or {})
    prep = _prepare(ctx, params, outdir)
    outdir = prep["outdir"]
    frame = _score_frame(prep)

    warnings = list(prep["warnings"])
    notes = list(prep["notes"])
    de_label = prep["sig_label"]
    top = frame.iloc[0]

    # ① 细胞类型 × 得分 / 细胞数
    t_score = table(
        "sc_celltype_score",
        "细胞类型 × 签名得分",
        [
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "n_cells", "label": "细胞数", "type": "number"},
            {"key": "pct_cells", "label": "占比(%)", "type": "number"},
            {"key": "score", "label": "归一化得分", "type": "number"},
            {"key": "raw_score", "label": "原始得分", "type": "number"},
            {"key": "rank", "label": "排名", "type": "number"},
        ],
        frame,
        note=f"score 为跨亚群 min-max 归一化（0–1）；打分方法 {prep['method']}",
    )

    # ② 目标分子 × 亚群 平均表达 / 表达比例
    expr_full = _mean_expression_table(prep, genes_limit=int(params.get("max_genes") or 50))
    t_expr = table(
        "sc_gene_expr",
        "目标分子亚群平均表达",
        [
            {"key": "gene", "label": "分子", "type": "string"},
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "mean_expr", "label": "平均表达", "type": "number"},
            {"key": "pct_expr", "label": "表达比例(%)", "type": "number"},
        ],
        expr_full,
        note="mean_expr 为参考内该亚群的平均表达（必要时 log1p）；pct_expr 为表达>0 的细胞比例",
    )

    tables = [t_score, t_expr]
    data_out = []
    if not expr_full.empty:
        path = export_csv(expr_full, outdir, "sc_gene_expr_by_celltype")
        data_out.append(bc.data_file("table_csv", path, "目标分子亚群表达表"))

    # ③ 细胞类型 × bulk 样本 相关
    corr, n_common = _celltype_bulk_corr(prep, ctx, warnings)
    if corr is not None:
        columns = [{"key": "cell_type", "label": "细胞类型", "type": "string"}] + [
            {"key": str(s), "label": str(s), "type": "number"} for s in corr.columns if s != "cell_type"
        ]
        tables.append(table(
            "sc_celltype_bulk_corr",
            "细胞类型 × bulk 样本相关",
            columns,
            corr,
            note=f"用 {n_common} 个共同签名基因的亚群平均表达向量与每个 bulk 样本表达向量做 Pearson 相关",
        ))
        cpath = export_csv(corr, outdir, "sc_celltype_bulk_corr")
        data_out.append(bc.data_file("table_csv", cpath, "细胞类型×bulk 样本相关矩阵"))

    figures = []
    try:
        figures.append(_fig_score_violin(prep, prep["scores"], outdir, prep["method"]))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"签名得分小提琴图绘制失败：{exc}")
    top_violin = int(params.get("top_violin") if params.get("top_violin") is not None else 6)
    if top_violin > 0:
        try:
            fig = _fig_target_violin(prep, outdir, top_violin)
            if fig:
                figures.append(fig)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"目标基因小提琴图绘制失败：{exc}")
    try:
        fig = _fig_corr_heat(corr, outdir)
        if fig:
            figures.append(fig)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"相关热图绘制失败：{exc}")

    cell_counts = {str(k): int(v) for k, v in zip(frame["cell_type"], frame["n_cells"])}
    bullets = [
        f"参考：{prep['ref_name']}（{prep['n_obs']} 细胞 × {prep['n_vars']} 基因，"
        f"{len(frame)} 个细胞类型；细胞类型列 {prep['ct_col']}，来源 {prep['ref_source']}）",
        f"签名来源：{de_label}；命中参考 {len(prep['matched'])}/{len(prep['sig_genes'])} 个基因",
        f"最高得分亚群：{top['cell_type']}（score {top['score']:.2f}，{int(top['n_cells'])} 细胞，"
        f"占 {top['pct_cells']:.1f}%）",
    ]
    if corr is not None and not corr.empty:
        # 给出每个 bulk 样本最相关的细胞类型（可解释映射）
        best = []
        for sample in [c for c in corr.columns if c != "cell_type"]:
            series = pd.to_numeric(corr[sample], errors="coerce")
            if series.notna().any():
                best.append((sample, str(corr["cell_type"].iloc[int(series.idxmax())]), float(series.max())))
        if best:
            shown = "；".join(f"{s}→{ct}(r={r:.2f})" for s, ct, r in best[:4])
            bullets.append(f"bulk 样本最相关亚群（前 4）：{shown}")

    notes.extend([
        f"打分方法：{prep['method']}（mean=签名基因跨细胞 z-score 的亚群均值；"
        "aucell=rank-AUC 简化 AUC=1−Σ名次/(k·N)；addmodule=签名均值−50 组表达量匹配随机对照集的背景均值）",
        "本步骤只读取参考的细胞类型注释，不改写参考数据；未做细胞类型注释推断",
        f"参考数据集：{prep['ref_path']}",
    ])
    if prep["backed"]:
        notes.append("参考以 backed 模式打开以控制内存；更多基因上的 marker 分析可能受限")
    metrics = {
        "ref_cells": prep["n_obs"],
        "ref_genes": prep["n_vars"],
        "n_celltypes": int(len(frame)),
        "signature_genes": int(len(prep["sig_genes"])),
        "signature_matched": int(len(prep["matched"])),
        "method": prep["method"],
        "top_celltype": str(top["cell_type"]),
        "top_score": round(float(top["score"]), 4),
        "celltype_cells": {str(k): int(v) for k, v in cell_counts.items()},
    }

    return result(
        title="单细胞映射",
        bullets=bullets,
        metrics=metrics,
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=data_out,
    )


# ---------------------------------------------------------------------------
# 分析实现
# ---------------------------------------------------------------------------

def _ana_signature_score(ctx, params=None, outdir=None):
    params = dict(params or {})
    prep = _prepare(ctx, params, outdir)
    frame = _score_frame(prep)
    out = prep["outdir"]
    warnings = list(prep["warnings"])
    expr = _mean_expression_table(prep, genes_limit=int(params.get("max_genes") or 50))
    tables = [table(
        "sc_score",
        "细胞类型 × 签名得分",
        [
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "n_cells", "label": "细胞数", "type": "number"},
            {"key": "pct_cells", "label": "占比(%)", "type": "number"},
            {"key": "score", "label": "归一化得分", "type": "number"},
            {"key": "raw_score", "label": "原始得分", "type": "number"},
            {"key": "rank", "label": "排名", "type": "number"},
        ],
        frame,
    ), table(
        "sc_gene_expr",
        "目标分子亚群表达",
        [
            {"key": "gene", "label": "分子", "type": "string"},
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "mean_expr", "label": "平均表达", "type": "number"},
            {"key": "pct_expr", "label": "表达比例(%)", "type": "number"},
        ],
        expr,
    )]
    figures = []
    try:
        figures.append(_fig_score_violin(prep, prep["scores"], out, prep["method"]))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"签名得分小提琴图绘制失败：{exc}")
    top = frame.iloc[0]
    return result(
        title="签名基因单细胞打分",
        bullets=[
            f"签名：{prep['sig_label']}（命中 {len(prep['matched'])} 个基因）",
            f"最高得分细胞类型：{top['cell_type']}（score {top['score']:.2f}，{int(top['n_cells'])} 细胞）",
            f"打分方法：{prep['method']}",
        ],
        metrics={
            "n_celltypes": int(len(frame)),
            "signature_genes": int(len(prep["sig_genes"])),
            "signature_matched": int(len(prep["matched"])),
            "method": prep["method"],
            "top_celltype": str(top["cell_type"]),
        },
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=list(prep["notes"]) + [f"签名来源：{prep['sig_label']}"],
        data_out=[bc.data_file("table_csv", export_csv(expr, out, "sc_signature_gene_expr"), "目标分子亚群表达表")],
    )


def _ana_de_map(ctx, params=None, outdir=None):
    """bulk 上调/下调签名分别映射到哪些细胞亚群，并给出可解释的方向。"""
    params = dict(params or {})
    prep = _prepare(ctx, params, outdir, require_sig=False)
    out = prep["outdir"]
    warnings = list(prep["warnings"])
    if prep["de_frame"] is None:
        raise EngineError("prerequisite", "尚未运行差异表达（de）步骤，无法做 DE→单细胞亚群映射")
    up_genes = prep["de_frame"][(prep["de_frame"]["padj"] <= 0.05) & (prep["de_frame"]["logFC"] > 0)]["gene"].tolist()
    down_genes = prep["de_frame"][(prep["de_frame"]["padj"] <= 0.05) & (prep["de_frame"]["logFC"] < 0)]["gene"].tolist()
    if not up_genes and not down_genes:
        raise EngineError("prerequisite", "差异表中没有 FDR<0.05 的显著基因，无法做 DE→单细胞映射")
    max_genes = int(params.get("max_genes") or 50)

    rows = []
    summary = {}
    for label, genes in (("上调", up_genes), ("下调", down_genes)):
        if not genes:
            continue
        matched, missing = _match_signature(genes[:max_genes], prep["ref_genes"])
        if len(matched) == 0:
            warnings.append(f"{label}签名基因在参考中都没找到，跳过该方向")
            continue
        # 借用 prep 的 X / celltypes 重新定位列
        x_lookup = {str(g): i for i, g in enumerate(prep["x_genes"])}
        rematched = []
        for bulk_gene, old_idx in matched:
            pos = x_lookup.get(str(prep["ref_genes"][old_idx]))
            if pos is not None:
                rematched.append((bulk_gene, int(pos)))
        sub = dict(prep)
        sub["matched"] = rematched
        sub["sig_genes"] = genes[:max_genes]
        scores = _score_cells(sub, prep["method"], warnings)
        frame = _aggregate_celltypes(scores, prep["celltypes"])
        sub["score_frame"] = frame
        sub["scores"] = scores
        top = frame.iloc[0]
        summary[label] = {"top_celltype": str(top["cell_type"]), "top_score": float(top["score"]),
                          "matched": len(rematched), "genes": len(genes)}
        for _, row in frame.iterrows():
            rows.append({
                "direction": label,
                "cell_type": str(row["cell_type"]),
                "score": round(float(row["score"]), 4),
                "raw_score": round(float(row["raw_score"]), 4),
                "n_cells": int(row["n_cells"]),
            })
    if not rows:
        raise EngineError("prerequisite", "上下调签名基因在参考中都未命中，无法映射")
    frame = pd.DataFrame(rows)
    tables = [table(
        "sc_de_map",
        "DE 方向 × 细胞亚群得分",
        [
            {"key": "direction", "label": "DE 方向", "type": "string"},
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "score", "label": "归一化得分", "type": "number"},
            {"key": "raw_score", "label": "原始得分", "type": "number"},
            {"key": "n_cells", "label": "细胞数", "type": "number"},
        ],
        frame,
    )]
    # 分组条形图
    figures = []
    try:
        directions = [d for d in ("上调", "下调") if d in summary]
        # 图内文字必须英文（契约 §2.5），中文方向名只出现在表格/摘要里
        direction_en = {"上调": "Up", "下调": "Down"}
        types = list(dict.fromkeys(frame["cell_type"].tolist()))
        fig, ax = bc.new_fig(8.4, 4.6)
        x = np.arange(len(types))
        width = 0.8 / max(1, len(directions))
        colors = palette(len(directions))
        for i, direction in enumerate(directions):
            sub = frame[frame["direction"] == direction].set_index("cell_type")
            vals = [float(sub["score"].get(t, 0.0)) for t in types]
            ax.bar(x + i * width - (len(directions) - 1) * width / 2.0, vals, width,
                   label=direction_en.get(direction, str(direction)),
                   color=colors[i], alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(types, rotation=30, ha="right")
        ax.set_ylabel("Normalized score")
        ax.set_title("Bulk DE direction mapped to single-cell subtypes")
        ax.legend()
        figures.append(save_fig(fig, out, "sc_de_map", title="DE direction mapped to subtypes"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"DE 映射图绘制失败：{exc}")
    bullets = []
    for label, info in summary.items():
        bullets.append(f"{label}签名 → 最相关亚群：{info['top_celltype']}（score {info['top_score']:.2f}，"
                       f"命中 {info['matched']} 个基因）")
    return result(
        title="DE→单细胞亚群映射",
        bullets=bullets,
        metrics={k: v for k, v in summary.items()},
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=[f"映射依据：bulk {label}签名在参考细胞上的 {prep['method']} 得分"
               for label in summary] + ["得分越高表示该亚群越富集该方向的 bulk 差异程序"],
    )


def _ana_celltype_prop(ctx, params=None, outdir=None):
    params = dict(params or {})
    prep = _prepare(ctx, params, outdir, require_sig=False)
    out = prep["outdir"]
    frame = _aggregate_celltypes(np.zeros(prep["n_obs"]), prep["celltypes"])
    # 占比即细胞数占比，得分无意义 → 重新整理
    frame = frame[["cell_type", "n_cells", "pct_cells"]].copy()
    frame["rank"] = np.arange(1, len(frame) + 1)
    tables = [table(
        "sc_celltype_prop",
        "细胞类型占比",
        [
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "n_cells", "label": "细胞数", "type": "number"},
            {"key": "pct_cells", "label": "占比(%)", "type": "number"},
            {"key": "rank", "label": "排名", "type": "number"},
        ],
        frame.sort_values("n_cells", ascending=False),
        note=f"参考 {prep['n_obs']} 个细胞；细胞类型列 {prep['ct_col']}",
    )]
    figures = [_fig_prop_bar(frame, out)]
    return result(
        title="参考细胞类型构成",
        bullets=[
            f"参考 {prep['ref_name']}：{prep['n_obs']} 细胞 × {prep['n_vars']} 基因，{len(frame)} 个细胞类型",
            f"最大亚群：{frame.sort_values('n_cells', ascending=False).iloc[0]['cell_type']}"
            f"（{int(frame.sort_values('n_cells', ascending=False).iloc[0]['n_cells'])} 细胞）",
        ],
        metrics={"n_celltypes": int(len(frame)), "n_cells": int(prep["n_obs"]),
                 "proportions": {str(r["cell_type"]): float(r["pct_cells"]) for _, r in frame.iterrows()}},
        tables=tables,
        figures=figures,
        warnings=list(prep["warnings"]),
        notes=[f"细胞类型列：{prep['ct_col']}（只读参考注释）"],
    )


def _ana_correlation(ctx, params=None, outdir=None):
    params = dict(params or {})
    prep = _prepare(ctx, params, outdir)
    out = prep["outdir"]
    warnings = list(prep["warnings"])
    _score_frame(prep)
    corr, n_common = _celltype_bulk_corr(prep, ctx, warnings)
    if corr is None:
        raise EngineError(
            "prerequisite",
            f"签名基因与 bulk 表达矩阵共同基因不足（{n_common} 个），无法计算细胞类型×样本相关；"
            "请确认参考与 bulk 基因命名一致",
        )
    columns = [{"key": "cell_type", "label": "细胞类型", "type": "string"}] + [
        {"key": str(s), "label": str(s), "type": "number"} for s in corr.columns if s != "cell_type"
    ]
    figures = []
    try:
        fig = _fig_corr_heat(corr, out)
        if fig:
            figures.append(fig)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"相关热图绘制失败：{exc}")
    # 每个样本最相关细胞类型
    best_rows = []
    for sample in [c for c in corr.columns if c != "cell_type"]:
        series = pd.to_numeric(corr[sample], errors="coerce")
        if series.notna().any():
            idx = int(series.idxmax())
            best_rows.append({"sample": str(sample), "best_cell_type": str(corr["cell_type"].iloc[idx]),
                              "r": round(float(series.max()), 4)})
    best = pd.DataFrame(best_rows)
    tables = [table("sc_celltype_bulk_corr", "细胞类型 × bulk 样本相关", columns, corr,
                    note=f"基于 {n_common} 个共同签名基因的 Pearson 相关")]
    if not best.empty:
        tables.append(table(
            "sc_best_match",
            "每个 bulk 样本最相关亚群",
            [
                {"key": "sample", "label": "bulk 样本", "type": "string"},
                {"key": "best_cell_type", "label": "最相关亚群", "type": "string"},
                {"key": "r", "label": "Pearson r", "type": "number"},
            ],
            best,
        ))
    bullets = [f"共同签名基因 {n_common} 个；细胞类型 {len(corr)} 个 × bulk 样本 {len(corr.columns) - 1} 个"]
    bullets += [f"样本 {r['sample']} 最相关亚群：{r['best_cell_type']}（r={r['r']:.2f}）" for _, r in best.head(5).iterrows()]
    return result(
        title="细胞类型 × bulk 样本相关",
        bullets=bullets,
        metrics={"n_common_genes": int(n_common), "n_celltypes": int(len(corr)),
                 "n_samples": int(len(corr.columns) - 1)},
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=["相关性为签名基因上的表达向量 Pearson 相关；只反映表达模式相似度，不代表因果"],
        data_out=[bc.data_file("table_csv", export_csv(corr, out, "sc_celltype_bulk_corr"), "相关矩阵")],
    )


def _ref_top_markers(prep, top_n: int = 25) -> dict[str, list[str]]:
    X = prep["X"]
    celltypes = prep["celltypes"]
    genes = prep["x_genes"]
    types = [str(t) for t in prep["score_frame"]["cell_type"].tolist()]
    means = {}
    for ct in types:
        mask = celltypes == ct
        means[ct] = X[mask].mean(axis=0) if mask.any() else np.zeros(X.shape[1])
    mat = np.vstack([means[ct] for ct in types])  # types × genes
    out: dict[str, list[str]] = {}
    for i, ct in enumerate(types):
        others = np.delete(mat, i, axis=0)
        base = others.mean(axis=0) if others.shape[0] else np.zeros(mat.shape[1])
        spec = mat[i] - base
        order = np.argsort(-spec)[:top_n]
        out[ct] = [str(genes[j]) for j in order if spec[j] > 0]
    return out


def _ana_marker_consistency(ctx, params=None, outdir=None):
    from scipy import stats

    params = dict(params or {})
    prep = _prepare(ctx, params, outdir, require_sig=False)
    out = prep["outdir"]
    warnings = list(prep["warnings"])
    de_frame = prep["de_frame"]
    if de_frame is None:
        raise EngineError("prerequisite", "尚未运行差异表达（de）步骤，无法做 marker 一致性分析")
    bulk_sig = set(de_frame[(de_frame["padj"] <= 0.05) & (de_frame["logFC"].abs() > 0)]["gene"].astype(str))
    if not bulk_sig:
        raise EngineError("prerequisite", "差异表中没有 FDR<0.05 的显著基因，无法做 marker 一致性分析")
    prep["score_frame"] = _aggregate_celltypes(np.zeros(prep["n_obs"]), prep["celltypes"])
    markers = _ref_top_markers(prep, top_n=int(params.get("top_n") or 25))
    universe = set(str(g) for g in prep["x_genes"])
    n_bulk = len(bulk_sig & universe)
    rows = []
    for ct, genes in markers.items():
        gl = set(genes)
        overlap = sorted(gl & bulk_sig)
        union = gl | (bulk_sig & universe)
        pval = float(stats.hypergeom.sf(len(overlap) - 1, max(len(universe), 1),
                                        max(n_bulk, 1), max(len(gl), 1))) if gl and n_bulk else 1.0
        rows.append({
            "cell_type": ct,
            "n_markers": len(gl),
            "n_overlap": len(overlap),
            "jaccard": round(len(overlap) / len(union), 4) if union else 0.0,
            "p_value": float(min(1.0, max(0.0, pval))),
            "overlap_genes": ", ".join(overlap[:10]) if overlap else "",
            "top_markers": ", ".join(genes[:10]),
        })
    frame = pd.DataFrame(rows).sort_values(["n_overlap", "jaccard"], ascending=[False, False]).reset_index(drop=True)
    if frame["n_markers"].sum() == 0:
        raise EngineError("prerequisite", "参考中无法计算各亚群 marker（X 为空或无差异表达）")
    padj = bc.bh_fdr(frame["p_value"].to_numpy())
    frame["padj"] = np.round(padj, 6)
    tables = [table(
        "sc_marker_consistency",
        "参考 marker 与 bulk DE 一致性",
        [
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "n_markers", "label": "参考 marker 数", "type": "number"},
            {"key": "n_overlap", "label": "与 bulk DE 重叠数", "type": "number"},
            {"key": "jaccard", "label": "Jaccard", "type": "number"},
            {"key": "p_value", "label": "超几何 p", "type": "number"},
            {"key": "padj", "label": "FDR", "type": "number"},
            {"key": "overlap_genes", "label": "重叠基因", "type": "string"},
            {"key": "top_markers", "label": "参考 top marker", "type": "string"},
        ],
        frame,
        note=f"参考 marker = 该亚群相对其它亚群平均表达最高的 {int(params.get('top_n') or 25)} 个基因",
    )]
    figures = []
    try:
        fig, ax = bc.new_fig(7.8, max(3.0, 0.42 * len(frame) + 1.6))
        colors = ["#e15759" if v > 0 else "#bab0ac" for v in frame["n_overlap"]]
        ax.barh(frame["cell_type"][::-1], frame["n_overlap"][::-1], color=colors[::-1])
        ax.set_xlabel("Overlap with bulk DE genes")
        ax.set_title("Single-cell markers overlapping bulk DE")
        figures.append(save_fig(fig, out, "sc_marker_consistency", title="Marker consistency"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"marker 一致性图绘制失败：{exc}")
    best = frame.iloc[0]
    return result(
        title="参考 marker 与 bulk DE 一致性",
        bullets=[
            f"bulk 显著基因 {len(bulk_sig)} 个，其中 {n_bulk} 个在参考基因集中",
            f"一致性最好亚群：{best['cell_type']}（重叠 {int(best['n_overlap'])} 个，"
            f"Jaccard {best['jaccard']:.3f}，FDR {best['padj']:.2g}）",
        ],
        metrics={"bulk_sig_genes": int(len(bulk_sig)), "best_celltype": str(best["cell_type"]),
                 "best_overlap": int(best["n_overlap"])},
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=["超几何检验的 background 为参考中已载入的基因集；FDR 为跨细胞类型 BH 校正"],
    )


def _ana_gene_violin(ctx, params=None, outdir=None):
    params = dict(params or {})
    genes = _as_gene_list(params.get("genes"))
    prep = _prepare(ctx, params, outdir, require_sig=not genes, genes_override=genes or None)
    out = prep["outdir"]
    warnings = list(prep["warnings"])
    prep["score_frame"] = _aggregate_celltypes(np.zeros(prep["n_obs"]), prep["celltypes"])
    if not prep["matched"]:
        raise EngineError(
            "prerequisite",
            "没有可用于小提琴图的基因：请用 genes 参数指定参考中存在的基因，或先跑 de 步骤",
        )
    top_n = int(params.get("top_violin") or len(prep["matched"]))
    figures = []
    fig = _fig_target_violin(prep, out, max(1, top_n))
    if fig:
        figures.append(fig)
    expr = _mean_expression_table(prep, genes_limit=len(prep["matched"]))
    tables = [table(
        "sc_gene_expr",
        "目标分子亚群表达",
        [
            {"key": "gene", "label": "分子", "type": "string"},
            {"key": "cell_type", "label": "细胞类型", "type": "string"},
            {"key": "mean_expr", "label": "平均表达", "type": "number"},
            {"key": "pct_expr", "label": "表达比例(%)", "type": "number"},
        ],
        expr,
    )]
    bullets = [f"绘制 {len(prep['matched'])} 个基因在 {len(prep['score_frame'])} 个细胞类型上的表达分布"]
    if not expr.empty:
        best = expr.sort_values("mean_expr", ascending=False).iloc[0]
        bullets.append(f"平均表达最高：{best['gene']} @ {best['cell_type']}（{best['mean_expr']:.2f}）")
    return result(
        title="目标基因亚群小提琴图",
        bullets=bullets,
        metrics={"genes": int(len(prep["matched"])), "n_celltypes": int(len(prep["score_frame"]))},
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=["表达值按参考原始尺度（必要时 log1p）；细胞类型来自参考 obs"],
    )


# ---------------------------------------------------------------------------
# 注册（§2.5）
# ---------------------------------------------------------------------------

_SCORE_PARAM = {
    "key": "score_method",
    "label": "打分方法",
    "type": "select",
    "default": "mean",
    "options": [
        {"value": "mean", "label": "mean（签名基因 z-score 均值）"},
        {"value": "aucell", "label": "AUCell（rank-AUC）"},
        {"value": "addmodule", "label": "AddModuleScore（对照校正）"},
    ],
}
_CELLTYPE_PARAM = {"key": "cell_type_col", "label": "细胞类型列（obs）", "type": "text", "default": ""}
_MAXGENES_PARAM = {"key": "max_genes", "label": "最多签名基因数", "type": "number", "default": 50, "min": 3, "max": 500, "step": 1}
_REF_PARAM = {"key": "ref_path", "label": "单细胞参考 h5ad", "type": "file", "default": ""}
_TOPVIOLIN_PARAM = {"key": "top_violin", "label": "小提琴图基因数", "type": "number", "default": 6, "min": 0, "max": 30, "step": 1}
_GENES_PARAM = {"key": "genes", "label": "基因列表", "type": "genes", "default": []}


STEPS = [
    {
        "id": "scmap",
        "title": "单细胞映射",
        "icon": "singlecell",
        "requires": ["de"],
        "params": [
            _REF_PARAM,
            {
                "key": "gene_set",
                "label": "签名来源",
                "type": "select",
                "default": "de_up",
                "options": [
                    {"value": "de_up", "label": "DE 上调基因"},
                    {"value": "de_down", "label": "DE 下调基因"},
                    {"value": "de_both", "label": "DE 上调+下调"},
                    {"value": "custom", "label": "自定义基因列表"},
                ],
            },
            {"key": "genes", "label": "自定义基因（gene_set=custom）", "type": "genes", "default": []},
            _CELLTYPE_PARAM,
            _SCORE_PARAM,
            _MAXGENES_PARAM,
            _TOPVIOLIN_PARAM,
        ],
        "run": run_scmap,
    }
]


ANALYSES = [
    {
        "id": "sc_signature_score",
        "title": "单细胞签名打分",
        "group": "单细胞映射",
        "method": "h5ad 参考 + 签名基因 mean/aucell/addmodule 打分（只读参考细胞类型注释）",
        "needs": ["de", "sc"],
        "params": [
            _REF_PARAM,
            {"key": "gene_set", "label": "签名来源", "type": "select", "default": "de_up",
             "options": [{"value": "de_up", "label": "DE 上调"}, {"value": "de_down", "label": "DE 下调"},
                         {"value": "de_both", "label": "DE 上下调"}, {"value": "custom", "label": "自定义"}]},
            _GENES_PARAM,
            _CELLTYPE_PARAM,
            _SCORE_PARAM,
            _MAXGENES_PARAM,
        ],
        "run": _ana_signature_score,
    },
    {
        "id": "sc_de_map",
        "title": "DE 方向→亚群映射",
        "group": "单细胞映射",
        "method": "bulk DE 上调/下调签名分别打分到参考亚群，给出方向性映射（不推断新细胞类型）",
        "needs": ["de", "sc"],
        "params": [_REF_PARAM, _CELLTYPE_PARAM, _SCORE_PARAM, _MAXGENES_PARAM],
        "run": _ana_de_map,
    },
    {
        "id": "sc_celltype_prop",
        "title": "参考细胞类型构成",
        "group": "单细胞映射",
        "method": "统计参考 h5ad 的 obs 细胞类型列（不修改参考）",
        "needs": ["sc"],
        "params": [_REF_PARAM, _CELLTYPE_PARAM],
        "run": _ana_celltype_prop,
    },
    {
        "id": "sc_correlation",
        "title": "细胞类型×bulk 相关",
        "group": "单细胞映射",
        "method": "签名基因上亚群平均表达向量与 bulk 样本表达向量的 Pearson 相关",
        "needs": ["de", "sc"],
        "params": [_REF_PARAM, {"key": "gene_set", "label": "签名来源", "type": "select", "default": "de_both",
                                "options": [{"value": "de_up", "label": "DE 上调"}, {"value": "de_down", "label": "DE 下调"},
                                            {"value": "de_both", "label": "DE 上下调"}, {"value": "custom", "label": "自定义"}]},
                   _GENES_PARAM, _CELLTYPE_PARAM, _MAXGENES_PARAM],
        "run": _ana_correlation,
    },
    {
        "id": "sc_marker_consistency",
        "title": "参考 marker 与 bulk DE 一致性",
        "group": "单细胞映射",
        "method": "参考各亚群 top 表达 marker 与 bulk 显著差异基因的重叠（超几何检验 + BH）",
        "needs": ["de", "sc"],
        "params": [_REF_PARAM, _CELLTYPE_PARAM, {"key": "top_n", "label": "每类 marker 数", "type": "number",
                                                 "default": 25, "min": 5, "max": 200, "step": 1}],
        "run": _ana_marker_consistency,
    },
    {
        "id": "sc_gene_violin",
        "title": "目标基因亚群小提琴图",
        "group": "单细胞映射",
        "method": "指定基因在参考各细胞类型中的表达分布（violin）",
        "needs": ["sc"],
        "params": [_REF_PARAM, _GENES_PARAM, _CELLTYPE_PARAM, _TOPVIOLIN_PARAM],
        "run": _ana_gene_violin,
    },
]
