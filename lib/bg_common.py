"""
dsh-geobulk — 共享工具层（所有步骤/分析模块都依赖本文件）

内容：
  * EngineError：带 error code 的异常（引擎转成 {ok:false,error:{code,message}}）
  * result()/table()：§2.5 结果信封的构造器（统一 JSON 序列化、500 行截断、NaN→None）
  * 绘图：new_fig()/save_fig()/palette() —— Agg 后端、白底、英文文字、dpi≈140
  * 统计：bh_fdr()、welch_de()（差异分析的最小公共实现）
  * 导出：export_csv()、jsonable()

约定：本文件只依赖标准库 + numpy/pandas/scipy/matplotlib（环境已装，禁止新增依赖）。
"""

from __future__ import annotations

import json
import math
import os
import re
import time

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# matplotlib：必须在 import pyplot 之前定好后端与缓存目录
# ---------------------------------------------------------------------------

def _setup_matplotlib() -> None:
    os.environ.setdefault("MPLBACKEND", "Agg")
    cache = os.environ.get("MPLCONFIGDIR")
    if not cache:
        base = os.environ.get("DS_BG_WORKDIR") or os.environ.get("DS_SC_WORKDIR") or "/tmp"
        cache = os.path.join(base, "mpl-cache")
    try:
        os.makedirs(cache, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = cache
    except OSError:
        os.environ["MPLCONFIGDIR"] = "/tmp"
    import matplotlib

    matplotlib.use("Agg", force=True)


_setup_matplotlib()

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

# 统一风格：白底、细网格、无上/右框线、英文文字
_RC = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.edgecolor": "#c8ccd4",
    "axes.labelcolor": "#2b2f36",
    "axes.titlesize": 11,
    "axes.titleweight": "medium",
    "axes.labelsize": 10,
    "axes.grid": True,
    "grid.color": "#e8eaee",
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": "#5a6068",
    "ytick.color": "#5a6068",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "font.size": 10,
    "figure.dpi": 110,
    "savefig.dpi": 140,
    "figure.autolayout": False,
}


def _apply_rc() -> None:
    for key, value in _RC.items():
        try:
            plt.rcParams[key] = value
        except Exception:  # noqa: BLE001 - 不同 matplotlib 版本容错
            pass
    # 尽量挑一个存在的无衬线字体，避免缺字体告警
    try:
        names = {f.name for f in font_manager.fontManager.ttflist}
        for cand in ("Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"):
            if cand in names:
                plt.rcParams["font.family"] = "sans-serif"
                plt.rcParams["font.sans-serif"] = [cand, "DejaVu Sans"]
                break
    except Exception:  # noqa: BLE001
        pass
    try:
        plt.rcParams["axes.unicode_minus"] = False
    except Exception:  # noqa: BLE001
        pass


_apply_rc()

# ---------------------------------------------------------------------------
# 异常
# ---------------------------------------------------------------------------

VALID_CODES = {
    "bad_request",
    "not_found",
    "prerequisite",
    "unsupported",
    "ambiguous",
    "internal",
    "timeout",
    "no_r",
    "no_python",
}


class EngineError(Exception):
    """引擎层可读错误：code 决定 HTTP 状态与前端提示。"""

    def __init__(self, code: str, message: str, **details):
        code = code if code in VALID_CODES else "internal"
        super().__init__(message)
        self.code = code
        self.message = str(message)
        self.details = details or {}

    def to_dict(self) -> dict:
        out = {"code": self.code, "message": self.message}
        for key, value in self.details.items():
            out[key] = jsonable(value)
        return out

    def __str__(self) -> str:  # pragma: no cover - 仅调试
        return f"[{self.code}] {self.message}"


def prerequisite(message: str, **details) -> EngineError:
    return EngineError("prerequisite", message, **details)


def bad_request(message: str, **details) -> EngineError:
    return EngineError("bad_request", message, **details)


# ---------------------------------------------------------------------------
# JSON 序列化 / 结果信封
# ---------------------------------------------------------------------------

_NAN = float("nan")


def jsonable(value):
    """把 numpy / pandas / NaN / Timestamp 递归转成纯 JSON 可序列化对象。"""
    if value is None:
        return None
    if isinstance(value, (str, bool)):
        return value
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (int,)):
        return int(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        f = float(value)
        if math.isnan(f) or math.isinf(f):
            return None
        return f
    if isinstance(value, (np.ndarray,)):
        return [jsonable(v) for v in value.tolist()]
    if isinstance(value, (pd.DataFrame,)):
        return [jsonable(r) for r in value.to_dict(orient="records")]
    if isinstance(value, (pd.Series,)):
        return [jsonable(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [jsonable(v) for v in value]
    if hasattr(value, "isoformat"):
        try:
            return value.isoformat()
        except Exception:  # noqa: BLE001
            return str(value)
    if hasattr(value, "item"):
        try:
            return jsonable(value.item())
        except Exception:  # noqa: BLE001
            return str(value)
    return str(value)


def _norm_columns(columns) -> list[dict]:
    out = []
    for item in columns or []:
        if isinstance(item, str):
            out.append({"key": item, "label": item, "type": "string"})
        elif isinstance(item, dict):
            col = dict(item)
            col.setdefault("key", str(col.get("key", "")))
            col.setdefault("label", col.get("label") or col["key"])
            col.setdefault("type", "string")
            out.append(col)
    return out


DEFAULT_TABLE_ROWS = 500


def table(tid: str, title: str, columns, rows, note: str | None = None, limit: int = DEFAULT_TABLE_ROWS) -> dict:
    """构造 §2.5 的 table 对象。

    rows 支持 DataFrame / list[dict] / list[list]（按 columns 顺序）。
    超过 limit 行时截断并标记 truncated=True（完整内容由调用方写 CSV 到 data_out）。
    """
    if isinstance(rows, pd.DataFrame):
        records = rows.to_dict(orient="records")
    elif rows is None:
        records = []
    elif isinstance(rows, np.ndarray):
        keys = [c["key"] for c in _norm_columns(columns)]
        records = [dict(zip(keys, r)) for r in rows.tolist()]
    else:
        records = []
        for row in rows:
            if isinstance(row, dict):
                records.append(row)
            elif isinstance(row, (list, tuple)):
                keys = [c["key"] for c in _norm_columns(columns)]
                records.append(dict(zip(keys, row)))
            else:
                records.append({"value": row})
    total = len(records)
    truncated = total > limit
    shown = records[:limit]
    out = {
        "id": str(tid),
        "title": str(title),
        "columns": _norm_columns(columns),
        "rows": [jsonable(r) for r in shown],
        "total": total,
    }
    if truncated:
        out["truncated"] = True
    if note:
        out["note"] = str(note)
    return out


def figure(fid: str, title: str, file: str, width: int | None = None, height: int | None = None) -> dict:
    """构造 §2.5 的 figure 对象（若给了宽高就带上）。"""
    item = {"id": str(fid), "title": str(title), "file": str(file)}
    if width and height:
        item["width"] = int(width)
        item["height"] = int(height)
    return item


def result(
    title: str | None = None,
    bullets=None,
    metrics=None,
    tables=None,
    figures=None,
    warnings=None,
    notes=None,
    sample_columns=None,
    sample_values=None,
    data_out=None,
    **extra,
) -> dict:
    """构造 §2.5 结果信封。所有模块统一走这里，保证字段齐全。"""
    summary = {"title": title or "", "bullets": [str(b) for b in (bullets or [])], "metrics": jsonable(metrics or {})}
    out = {
        "summary": summary,
        "tables": [t for t in (tables or []) if t],
        "figures": [f for f in (figures or []) if f],
        "warnings": [str(w) for w in (warnings or [])],
        "notes": [str(n) for n in (notes or [])],
    }
    if sample_columns:
        out["sample_columns"] = list(sample_columns)
    if sample_values:
        out["sample_values"] = {str(k): jsonable(v) for k, v in sample_values.items()}
    if data_out:
        out["data_out"] = [d for d in data_out if d]
    for key, value in extra.items():
        out[key] = jsonable(value)
    return out


def data_file(kind: str, path: str, label: str | None = None) -> dict:
    return {"kind": str(kind), "path": str(path), "label": label or os.path.basename(str(path))}


# ---------------------------------------------------------------------------
# 目录 / 文件
# ---------------------------------------------------------------------------

def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def slug(text: str, limit: int = 48) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "-", str(text or "")).strip("-")
    return (s[:limit] or "out").lower()


def export_csv(df: pd.DataFrame, outdir: str, name: str, index: bool = True) -> str:
    """把 DataFrame 写进 outdir，返回绝对路径。"""
    ensure_dir(outdir)
    filename = name if name.lower().endswith(".csv") else f"{name}.csv"
    path = os.path.join(outdir, filename)
    frame = df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)
    frame.to_csv(path, index=index, encoding="utf-8")
    return path


def export_json(obj, outdir: str, name: str) -> str:
    ensure_dir(outdir)
    filename = name if name.lower().endswith(".json") else f"{name}.json"
    path = os.path.join(outdir, filename)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(jsonable(obj), handle, ensure_ascii=False, indent=2)
    return path


def export_tsv(df: pd.DataFrame, outdir: str, name: str, index: bool = True) -> str:
    ensure_dir(outdir)
    filename = name if name.lower().endswith(".tsv") else f"{name}.tsv"
    path = os.path.join(outdir, filename)
    frame = df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)
    frame.to_csv(path, sep="\t", index=index, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# 绘图
# ---------------------------------------------------------------------------

PALETTE = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948",
    "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#1f77b4", "#ff7f0e",
    "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf", "#a0cbe8", "#ffbe7d", "#8cd17d", "#b6992d",
]

GRADIENT = ["#3b4cc0", "#6a9fdb", "#9ec6e8", "#f7f7f7", "#f0a35e", "#b40426"]


def palette(n: int) -> list[str]:
    n = max(1, int(n))
    if n <= len(PALETTE):
        return PALETTE[:n]
    out = []
    for i in range(n):
        out.append(PALETTE[i % len(PALETTE)])
    return out


def new_fig(width: float = 9.0, height: float = 5.2, **kwargs):
    """新建图表：返回 (fig, ax)。kwargs 透传给 plt.subplots。"""
    fig, ax = plt.subplots(figsize=(width, height), **kwargs)
    try:
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
    except Exception:  # noqa: BLE001
        pass
    return fig, ax


def save_fig(fig, outdir: str, name: str, title: str | None = None, dpi: int = 140, tight: bool = True) -> dict:
    """保存 PNG 并返回 §2.5 的 figure 对象。"""
    ensure_dir(outdir)
    filename = name if name.lower().endswith(".png") else f"{name}.png"
    path = os.path.join(outdir, filename)
    try:
        if tight:
            fig.tight_layout()
    except Exception:  # noqa: BLE001
        pass
    if title:
        # 已有 axes 标题时不再加总标题，避免图上出现两行重复标题
        has_axes_title = False
        try:
            has_axes_title = any(bool(ax.get_title()) for ax in fig.axes)
        except Exception:  # noqa: BLE001
            has_axes_title = False
        if not has_axes_title:
            try:
                fig.suptitle(title, fontsize=11, y=0.995)
            except Exception:  # noqa: BLE001
                pass
    fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    width = height = None
    try:
        width, height = fig.get_size_inches() * dpi
        width, height = int(width), int(height)
    except Exception:  # noqa: BLE001
        pass
    plt.close(fig)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise EngineError("internal", f"绘图失败：{filename} 未生成")
    return {"id": slug(os.path.splitext(filename)[0]), "title": title or os.path.splitext(filename)[0], "file": path,
            "width": width, "height": height}


def close_fig(fig) -> None:
    try:
        plt.close(fig)
    except Exception:  # noqa: BLE001
        pass


def png_size(path: str) -> tuple[int | None, int | None]:
    try:
        from PIL import Image

        with Image.open(path) as img:
            return int(img.width), int(img.height)
    except Exception:  # noqa: BLE001
        return None, None


def style_legend(ax, **kwargs) -> None:
    try:
        ax.legend(frameon=False, **kwargs)
    except Exception:  # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# 统计
# ---------------------------------------------------------------------------

def bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg FDR（对 NaN 容错）。"""
    arr = np.asarray(list(pvals) if not isinstance(pvals, np.ndarray) else pvals, dtype=float)
    n = arr.size
    out = np.full(n, np.nan)
    mask = ~np.isnan(arr)
    if not mask.any():
        return out
    idx = np.where(mask)[0]
    order = np.argsort(arr[mask])
    ranked = arr[mask][order]
    m = ranked.size
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(m)
    tmp[order] = q
    out[idx] = tmp
    return out


def _safe_p(t: float, df: float) -> float:
    """双侧 t 检验 p 值；不依赖 scipy 失败。"""
    try:
        from scipy import stats

        return float(2 * stats.t.sf(abs(t), max(df, 1e-6)))
    except Exception:  # noqa: BLE001
        # 正态近似回退
        return float(2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2)))))


def welch_de(
    X,
    labels,
    case: str,
    control: str,
    genes=None,
    *,
    lfc: float = 0.0,
    padj_cut: float = 0.05,
) -> pd.DataFrame:
    """Welch t 检验 + BH FDR（最小公共实现，供没有 limma 时回退）。

    X: (n_samples, n_genes) ndarray（已 log 尺度最佳）；labels: 每个样本的分组（与 X 行对齐）；
    case/control: 分组名。返回列：gene, logFC, AveExpr, t, P.Value, adj.P.Val, B, regulation。
    """
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise bad_request("表达矩阵维度不正确（应为 样本×基因）")
    labels = [str(v) for v in labels]
    samples = X.shape[0]
    if genes is None:
        genes = [f"GENE{i + 1}" for i in range(X.shape[1])]
    genes = [str(g) for g in genes]
    if len(labels) != samples:
        raise bad_request(f"分组数量（{len(labels)}）与样本数（{samples}）不一致")
    case_mask = np.array([lab == str(case) for lab in labels], dtype=bool)
    ctrl_mask = np.array([lab == str(control) for lab in labels], dtype=bool)
    n_case = int(case_mask.sum())
    n_ctrl = int(ctrl_mask.sum())
    if n_case < 2 or n_ctrl < 2:
        raise prerequisite(f"分组「{case}」({n_case}) 或「{control}」({n_ctrl}) 样本不足 2 个，无法做差异检验")
    A = X[case_mask]
    B = X[ctrl_mask]
    mean_a = A.mean(axis=0)
    mean_b = B.mean(axis=0)
    logfc = mean_a - mean_b
    var_a = A.var(axis=0, ddof=1)
    var_b = B.var(axis=0, ddof=1)
    se = np.sqrt(var_a / n_case + var_b / n_ctrl)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, logfc / se, 0.0)
    df = (var_a / n_case + var_b / n_ctrl) ** 2 / (
        (var_a / n_case) ** 2 / max(n_case - 1, 1) + (var_b / n_ctrl) ** 2 / max(n_ctrl - 1, 1)
    )
    df = np.where(np.isfinite(df) & (df > 0), df, 1.0)
    p = np.array([_safe_p(ti, di) for ti, di in zip(t, df)])
    p = np.where(np.isfinite(p), p, 1.0)
    adj = bh_fdr(p)
    ave = X.mean(axis=0)
    frame = pd.DataFrame(
        {
            "gene": genes,
            "logFC": logfc,
            "AveExpr": ave,
            "t": t,
            "P.Value": p,
            "adj.P.Val": adj,
            "B": np.zeros_like(t),
        }
    )
    frame["regulation"] = np.where(
        (frame["adj.P.Val"] <= padj_cut) & (frame["logFC"].abs() >= max(lfc, 0)),
        np.where(frame["logFC"] > 0, "up", "down"),
        "ns",
    )
    frame = frame.sort_values(["adj.P.Val", "P.Value"], ascending=[True, True]).reset_index(drop=True)
    return frame


def zscore_rows(mat) -> np.ndarray:
    """按行 z-score（用于热图），全零行置 0。"""
    arr = np.asarray(mat, dtype=np.float64)
    mean = np.nanmean(arr, axis=1, keepdims=True)
    std = np.nanstd(arr, axis=1, keepdims=True)
    std = np.where(std > 0, std, 1.0)
    out = (arr - mean) / std
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def min_max_scale(arr) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    lo, hi = np.nanmin(a), np.nanmax(a)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo <= 0:
        return np.zeros_like(a)
    return (a - lo) / (hi - lo)


def aucell_score(expr, gene_idx, *, n_genes: int | None = None) -> np.ndarray:
    """AUCell 简化版：每个样本按基因表达排序，取 top n_genes 的恢复曲线 AUC。"""
    arr = np.asarray(expr, dtype=np.float64)
    n_samples, n_expr = arr.shape
    if n_expr == 0:
        return np.zeros(n_samples)
    k = int(n_genes or max(1, round(0.05 * n_expr)))
    idx = np.asarray(list(gene_idx), dtype=int)
    idx = idx[(idx >= 0) & (idx < n_expr)]
    if idx.size == 0:
        return np.zeros(n_samples)
    out = np.zeros(n_samples)
    ranks = np.argsort(np.argsort(-arr, axis=1), axis=1)  # 0 = 最高表达
    hits = np.isin(ranks, idx)
    for i in range(n_samples):
        top = hits[i, :k]
        if top.size == 0:
            continue
        # 恢复曲线下面积 / 最大值
        curve = np.cumsum(top) / max(1, top.sum())
        out[i] = float(np.trapz(curve) / max(len(curve), 1))
    return np.clip(out, 0, 1)


def mean_top_score(expr, gene_idx) -> np.ndarray:
    arr = np.asarray(expr, dtype=np.float64)
    idx = np.asarray(list(gene_idx), dtype=int)
    idx = idx[(idx >= 0) & (idx < arr.shape[1])]
    if idx.size == 0:
        return np.zeros(arr.shape[0])
    return arr[:, idx].mean(axis=1)


IMMUNE_CELL_COLORS = {
    "B cells": "#4e79a7",
    "CD8 T cells": "#f28e2b",
    "CD4 T cells": "#e15759",
    "NK cells": "#76b7b2",
    "Monocytes": "#59a14f",
    "Macrophages": "#edc948",
    "Dendritic cells": "#b07aa1",
    "Neutrophils": "#ff9da7",
    "Tregs": "#9c755f",
    "Mast cells": "#bab0ac",
    "Endothelial": "#1f77b4",
    "Fibroblasts": "#8c564b",
    "CD4 memory T cells": "#9467bd",
    "CD8 naive T cells": "#e377c2",
    "pDC": "#7f7f7f",
    "iDC": "#bcbd22",
    "aDC": "#17becf",
    "Tem": "#a0cbe8",
    "Tcm": "#ffbe7d",
    "Th1 cells": "#8cd17d",
    "Th2 cells": "#b6992d",
    "Th17 cells": "#86b1d6",
    "Tfh": "#c49c94",
    "Memory B cells": "#d37295",
    "Naive B cells": "#6b9bd2",
    "Plasma cells": "#b4a0c8",
    "Eosinophils": "#f1ce63",
    "Activated DC": "#8c8c8c",
    "Gamma delta T cells": "#a1d99b",
    "Regulatory T cells": "#dbdb8d",
}


def color_for(name: str, index: int = 0) -> str:
    if name in IMMUNE_CELL_COLORS:
        return IMMUNE_CELL_COLORS[name]
    return PALETTE[index % len(PALETTE)]


# ---------------------------------------------------------------------------
# 时间
# ---------------------------------------------------------------------------

class Timer:
    def __init__(self):
        self.t0 = time.time()

    @property
    def ms(self) -> int:
        return int((time.time() - self.t0) * 1000)

    def __enter__(self):
        return self

    def __exit__(self, *exc):  # noqa: D401
        return False
