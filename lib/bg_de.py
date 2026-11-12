#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 差异表达模块（步骤 ``de`` + 差异类分析）。

实现策略（CONTRACT.md §8.2）：
  * ``method=limma``  → R ``limma::lmFit + eBayes + topTable``（官方实现；支持 ``covars`` 协变量与 ``paired`` 配对）
  * ``method=deseq2`` → R ``DESeq2``（需要原始 counts；负二项 GLM + Wald 检验）
  * ``method=welch``  → 纯 Python Welch t 检验 + BH FDR
  * R 不可用 / R 执行失败 → 自动回退 Welch（lite），并写 ``notes`` / ``warnings``

统一输出列：``gene, logFC, AveExpr, t, P.Value, adj.P.Val, B, regulation``
缓存：返回 ``de_store``，registry 会写入 ``dataset.results["de"]`` → 下游 ``ctx["de"]``。

对外接口
--------
* ``STEPS``：``de`` 步骤（registry 自动注册）
* ``ANALYSES``：``de_multi`` / ``de_venn`` / ``de_summary`` / ``wgcna_module``
* ``get_de(ctx)``：给其它模块读取最近的差异结果（DataFrame）

本文件只依赖标准库 + numpy/pandas/scipy + bg_common/rbridge（不新增第三方依赖）。
"""

from __future__ import annotations

import math
import os
import re

import numpy as np
import pandas as pd

from bg_common import (
    EngineError,
    bad_request,
    bh_fdr,
    data_file,
    ensure_dir,
    export_csv,
    jsonable,
    new_fig,
    palette,
    prerequisite,
    result,
    save_fig,
    table,
    welch_de,
    zscore_rows,
)
import rbridge

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

DE_COLUMNS = ["gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "regulation"]

DE_COLUMNS_META = [
    {"key": "gene", "label": "基因", "type": "string"},
    {"key": "logFC", "label": "log2FC", "type": "number"},
    {"key": "AveExpr", "label": "平均表达", "type": "number"},
    {"key": "t", "label": "t 统计量", "type": "number"},
    {"key": "P.Value", "label": "P 值", "type": "number"},
    {"key": "adj.P.Val", "label": "FDR", "type": "number"},
    {"key": "B", "label": "B（log-odds）", "type": "number"},
    {"key": "regulation", "label": "调控", "type": "string"},
]

METHODS = ["limma", "deseq2", "welch"]
METHOD_LABELS = {"limma": "limma（R 官方）", "deseq2": "DESeq2（R 官方）", "welch": "Welch t 检验（Python）"}

# 上/下调关键字（用于自动推断 case/control）
_CASE_WORDS = (
    "case", "tumor", "tumour", "cancer", "disease", "diseased", "patient", "treated", "treatment",
    "exp", "experiment", "mutant", "mutation", "ko", "knockout", "overexpression", "oe", "high",
    "resistant", "resist", "responder", "trt", "drug", "stim", "stimulated", "activated", "lesion",
)
_CTRL_WORDS = (
    "control", "ctrl", "normal", "healthy", "adjacent", "non-tumor", "nontumor", "wt", "wildtype",
    "wild-type", "sham", "vehicle", "untreated", "mock", "low", "sensitive", "non-responder",
    "unstimulated", "naive", "nc", "vector", "scramble",
)

ORDER_OPTIONS = [
    {"value": "padj", "label": "按 FDR"},
    {"value": "pvalue", "label": "按 P 值"},
    {"value": "logfc", "label": "按 |logFC|"},
    {"value": "aveexpr", "label": "按平均表达"},
    {"value": "none", "label": "不排序"},
]


# ---------------------------------------------------------------------------
# 参数归一化（UI 可能把数字/布尔传成字符串）
# ---------------------------------------------------------------------------

def _p_str(value, default: str = "") -> str:
    if value is None:
        return default
    text = str(value).strip()
    return text if text else default


def _p_float(value, default: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    if not math.isfinite(number):
        return float(default)
    return number


def _p_int(value, default: int) -> int:
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return int(default)


def _p_bool(value, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "y", "on", "是"):
        return True
    if text in ("0", "false", "no", "n", "off", "否", ""):
        return False
    return bool(default)


def _p_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set, np.ndarray, pd.Index)):
        return [str(v).strip() for v in value if str(v).strip()]
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in re.split(r"[,;|\s]+", text) if item.strip()]


def _norm_params(params: dict, ctx: dict | None = None) -> dict:
    params = dict(params or {})
    group_col = _p_str(params.get("group_col"))
    if not group_col and ctx:
        group_col = _p_str(ctx.get("group_key"))
    method = _p_str(params.get("method"), "limma").lower()
    if method in ("deseq", "deseq2_r", "negative_binomial"):
        method = "deseq2"
    if method in ("t", "ttest", "t-test", "student", "python"):
        method = "welch"
    if method not in METHODS:
        method = "limma"
    order = _p_str(params.get("order"), "padj").lower()
    if order not in {item["value"] for item in ORDER_OPTIONS}:
        order = "padj"
    return {
        "group_col": group_col,
        "case": _p_str(params.get("case")),
        "control": _p_str(params.get("control")),
        "method": method,
        "covars": _p_list(params.get("covars") or params.get("covariates")),
        "lfc": max(0.0, _p_float(params.get("lfc"), 0.0)),
        "padj": min(1.0, max(1e-8, _p_float(params.get("padj"), 0.05))),
        "paired": _p_bool(params.get("paired"), False),
        "top_n": max(1, min(100, _p_int(params.get("top_n"), 20))),
        "order": order,
        "timeout": max(60, _p_int(params.get("timeout"), 900)),
    }


# ---------------------------------------------------------------------------
# ctx 读取 / 差异结果读取
# ---------------------------------------------------------------------------

def _matrix(ctx: dict) -> tuple[np.ndarray, list[str], list[str]]:
    X = np.asarray(ctx.get("X"), dtype=np.float64)
    genes = [str(g) for g in (ctx.get("genes") or [])]
    samples = [str(s) for s in (ctx.get("samples") or [])]
    if X.ndim != 2 or X.size == 0:
        raise prerequisite("当前表达矩阵为空，无法做差异分析（请先导入数据并完成标准化）")
    if X.shape[0] != len(samples) or X.shape[1] != len(genes):
        raise EngineError(
            "internal",
            f"表达矩阵维度（{X.shape[0]}×{X.shape[1]}）与样本/基因数（{len(samples)}×{len(genes)}）不一致",
        )
    if len(samples) < 4:
        raise prerequisite(f"只有 {len(samples)} 个样本，至少需要 4 个（每组 ≥2）才能做差异分析")
    return X, genes, samples


def _info(ctx: dict, samples: list[str]) -> pd.DataFrame:
    info = ctx.get("sample_info")
    if not isinstance(info, pd.DataFrame) or info.empty:
        raise prerequisite("缺少样本信息表（sample_info），无法确定分组；请先在「样本」里导入分组表")
    info = info.copy()
    info.index = [str(i) for i in info.index]
    missing = [s for s in samples if s not in info.index]
    if missing:
        raise prerequisite(
            "样本信息表缺少 %d 个样本（如 %s），请检查样本名是否一致" % (len(missing), "、".join(missing[:3]))
        )
    return info.loc[samples]


def _de_frame_from(obj) -> pd.DataFrame | None:
    """把 de_store / DataFrame / records 统一成 DataFrame（找不到返回 None）。"""
    if obj is None:
        return None
    if isinstance(obj, pd.DataFrame):
        frame = obj.copy()
    elif isinstance(obj, dict):
        candidate = obj.get("table")
        if isinstance(candidate, pd.DataFrame):
            frame = candidate.copy()
        elif isinstance(candidate, list) and candidate and isinstance(candidate[0], dict):
            frame = pd.DataFrame(candidate)
        else:
            return None
    elif isinstance(obj, list) and obj and isinstance(obj[0], dict):
        frame = pd.DataFrame(obj)
    else:
        return None
    if frame.empty:
        return None
    if "gene" not in frame.columns:
        index = frame.index
        if getattr(index, "name", None) == "gene" or (
            len(index) and str(index[0]) not in ("0", "1", "2") and "logFC" in frame.columns
        ):
            frame = frame.reset_index().rename(columns={frame.index.name or "index": "gene"})
    if "gene" not in frame.columns:
        return None
    frame["gene"] = frame["gene"].astype(str)
    return frame


def get_de(ctx: dict) -> pd.DataFrame:
    """读取最近的差异表达结果（给其它模块复用）。

    ``ctx["de"]`` 来自 ``dataset.results["de"]``（de_store），兼容：
      * dict（含 ``table`` DataFrame 或 records）
      * 直接是 DataFrame
    返回的 DataFrame 至少有 ``gene`` 与 ``logFC`` 列；缺结果抛 ``prerequisite``。
    """
    if not isinstance(ctx, dict):
        raise prerequisite("ctx 无效，无法读取差异表达结果")
    candidates = []
    for key in ("de", "de_result"):
        if ctx.get(key) is not None:
            candidates.append(ctx.get(key))
    results = ctx.get("results") if isinstance(ctx.get("results"), dict) else {}
    if results.get("de") is not None:
        candidates.append(results.get("de"))
    for obj in candidates:
        frame = _de_frame_from(obj)
        if frame is not None and not frame.empty:
            return frame
    raise prerequisite("还没有差异表达结果，请先运行「差异表达」(de) 步骤")


def de_gene_lists(
    frame: pd.DataFrame,
    *,
    direction: str = "both",
    padj: float = 0.05,
    lfc: float = 0.0,
) -> tuple[list[str], dict]:
    """按方向取差异基因列表（优先用 regulation 列，其次重算阈值）。"""
    frame = _de_frame_from(frame) if not isinstance(frame, pd.DataFrame) else frame.copy()
    if frame is None or frame.empty:
        return [], {"n_up": 0, "n_down": 0, "n_sig": 0}
    if "regulation" in frame.columns:
        regulation = frame["regulation"].astype(str).str.lower().values
    else:
        logfc = pd.to_numeric(frame.get("logFC"), errors="coerce").fillna(0).values
        qvalue = pd.to_numeric(frame.get("adj.P.Val"), errors="coerce").fillna(1).values
        sig = qvalue <= padj
        regulation = np.where(sig & (logfc >= lfc), "up", np.where(sig & (logfc <= -lfc), "down", "ns"))
    up = frame.loc[regulation == "up", "gene"].astype(str).tolist()
    down = frame.loc[regulation == "down", "gene"].astype(str).tolist()
    direction = str(direction or "both").lower()
    if direction == "up":
        genes = up
    elif direction == "down":
        genes = down
    elif direction == "all":
        genes = frame["gene"].astype(str).tolist()
    else:
        genes = list(dict.fromkeys(up + down))
    return genes, {"n_up": len(up), "n_down": len(down), "n_sig": len(up) + len(down)}


# ---------------------------------------------------------------------------
# 分组解析
# ---------------------------------------------------------------------------

def _labels(ctx: dict, info: pd.DataFrame, samples: list[str], group_col: str) -> list[str]:
    if group_col not in info.columns:
        columns = "、".join(str(c) for c in info.columns[:12])
        raise prerequisite(f"样本信息里没有分组列「{group_col}」（可用列：{columns}）")
    series = info[group_col]
    out: list[str] = []
    for sample in samples:
        value = series.loc[sample]
        if isinstance(value, pd.Series):
            value = value.iloc[0]
        if value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value):
            out.append("NA")
        else:
            out.append(str(value).strip())
    return out


def _resolve_groups(ctx: dict, params: dict) -> tuple[str, list[str], list[str], list[str]]:
    _, _, samples = _matrix(ctx)
    info = _info(ctx, samples)
    group_col = _p_str(params.get("group_col")) or _p_str(ctx.get("group_key"))
    if not group_col:
        raise prerequisite("请指定分组列（group_col），或先在样本表里设置分组")
    labels = _labels(ctx, info, samples, group_col)
    uniq = [value for value in dict.fromkeys(labels) if value not in ("NA", "nan", "None", "")]
    if len(uniq) < 2:
        raise prerequisite(f"分组列「{group_col}」只有 {len(uniq)} 个有效取值（{uniq or '空'}），无法做两组比较")
    return group_col, labels, uniq, samples


def _resolve_case_control(params: dict, uniq: list[str], ctx: dict, group_col: str) -> tuple[str, str, list[str]]:
    notes: list[str] = []
    case = _p_str(params.get("case"))
    control = _p_str(params.get("control"))
    if case and case not in uniq:
        raise prerequisite(f"实验组「{case}」不在分组列「{group_col}」里（可选：{'、'.join(uniq)}）")
    if control and control not in uniq:
        raise prerequisite(f"对照组「{control}」不在分组列「{group_col}」里（可选：{'、'.join(uniq)}）")
    if case and control and case == control:
        raise bad_request(f"case 与 control 相同（{case}），请选择不同的两组")
    if case and control:
        return case, control, notes
    if len(uniq) > 2:
        raise prerequisite(
            f"分组列「{group_col}」有 {len(uniq)} 个取值（{'、'.join(uniq)}），请显式指定 case 与 control"
        )
    first, second = uniq[0], uniq[1]
    if not case and not control:
        inferred = _infer_pair(first, second)
        if inferred:
            control, case = inferred
            notes.append(f"未指定对比组，按关键字推断：{case} vs {control}")
        else:
            order = [str(g) for g in (ctx.get("groups") or [])]
            if len(order) >= 2 and set(order[:2]) <= set(uniq):
                case, control = order[0], order[1]
            else:
                case, control = first, second
            notes.append(f"未指定对比组，默认按分组顺序取：{case} vs {control}")
    elif not case:
        case = second if control == first else first
        notes.append(f"未指定实验组，自动选择：{case}")
    else:
        control = second if case == first else first
        notes.append(f"未指定对照组，自动选择：{control}")
    return case, control, notes


def _infer_pair(first: str, second: str) -> tuple[str, str] | None:
    def score(text: str, words) -> int:
        low = text.lower()
        return sum(1 for word in words if word in low)

    s1_case, s1_ctrl = score(first, _CASE_WORDS), score(first, _CTRL_WORDS)
    s2_case, s2_ctrl = score(second, _CASE_WORDS), score(second, _CTRL_WORDS)
    if s1_case > s1_ctrl and s2_ctrl > s2_case:
        return first, second  # (case, control)
    if s2_case > s2_ctrl and s1_ctrl > s1_case:
        return second, first
    return None


def _find_pair_col(info: pd.DataFrame, samples: list[str], case_mask, ctrl_mask) -> str:
    pattern = re.compile(r"patient|subject|pair|individual|donor|person|animal|case_id|block", re.I)
    for col in info.columns:
        if not pattern.search(str(col)):
            continue
        values = [str(info.loc[sample, col]) for sample in samples]
        table_map: dict[str, list[int]] = {}
        for i, value in enumerate(values):
            table_map.setdefault(value, []).append(i)
        usable = True
        usable_pairs = 0
        for value, idx in table_map.items():
            if len(idx) > 2 or len(idx) == 1:
                usable = False
                break
            kinds = {bool(case_mask[i]) for i in idx}
            if len(kinds) != 2:
                usable = False
                break
            usable_pairs += 1
        if usable and usable_pairs >= 2:
            return str(col)
    return ""


# ---------------------------------------------------------------------------
# 矩阵尺度转换 / counts
# ---------------------------------------------------------------------------

def _to_log_scale(X: np.ndarray, ctx: dict, method: str) -> tuple[np.ndarray, list[str]]:
    """非 lognorm 的矩阵转 log2 尺度（limma/welch 需要近似正态的表达值）。"""
    if method == "deseq2":
        return X, []
    kind = str(ctx.get("X_kind") or "unknown").lower()
    notes: list[str] = []
    if kind == "lognorm":
        return X, notes
    if kind in ("counts", "tpm", "fpkm", "cpm", "rpm"):
        notes.append(f"矩阵类型为 {kind}，已做 log2(x+1) 转换后再检验")
        return np.log2(np.clip(X, 0, None) + 1.0), notes
    # unknown：按数值特征判断
    finite = X[np.isfinite(X)]
    if finite.size:
        vmax = float(np.nanmax(finite))
        vmin = float(np.nanmin(finite))
        integral = bool(np.allclose(finite[: min(finite.size, 5000)], np.round(finite[: min(finite.size, 5000)])))
        if vmin >= 0 and (vmax > 50 or integral):
            notes.append(f"未标注矩阵类型，按数值特征（max={vmax:.0f}）判断为原始尺度，已做 log2(x+1) 转换")
            return np.log2(np.clip(X, 0, None) + 1.0), notes
    return X, notes


def _counts_matrix(ctx: dict, X: np.ndarray) -> np.ndarray | None:
    counts = ctx.get("counts")
    if counts is not None:
        arr = np.asarray(counts, dtype=np.float64)
        if arr.ndim == 2 and arr.shape == X.shape:
            return arr
    if str(ctx.get("X_kind") or "").lower() == "counts":
        return X
    return None


def _log_cpm(counts: np.ndarray) -> np.ndarray:
    lib = np.nansum(counts, axis=1, keepdims=True)
    lib = np.where(lib > 0, lib, 1.0)
    return np.log2(np.clip(counts, 0, None) / lib * 1e6 + 1.0)


# ---------------------------------------------------------------------------
# R 侧调用（limma / DESeq2）
# ---------------------------------------------------------------------------

def _r_diff(
    ctx: dict,
    params: dict,
    outdir: str,
    *,
    expr: np.ndarray,
    genes: list[str],
    samples: list[str],
    info: pd.DataFrame,
    method: str,
    pair_col: str,
    case: str,
    control: str,
    group_col: str,
) -> dict:
    """调 R（bg_de.R）。返回 {"frame": DataFrame|None, "failed": bool, "error": str, ...}。"""
    payload = {
        "expr": expr,
        "genes": genes,
        "samples": samples,
        "info": info,
        "params": {
            "method": method,
            "group_col": group_col,
            "case": case,
            "control": control,
            "covars": list(params.get("covars") or []),
            "paired": bool(params.get("paired")),
            "pair_col": pair_col,
            "lfc": params.get("lfc"),
            "padj": params.get("padj"),
            "order": params.get("order"),
            "top_n": params.get("top_n"),
            "matrix_kind": "counts" if method == "deseq2" else "log",
        },
    }
    try:
        outgoing = rbridge.run_r("bg_de.R", payload, outdir, timeout=int(params.get("timeout") or 900))
    except EngineError as error:
        return {"frame": None, "failed": True, "error": error.message, "code": error.code}
    tables = outgoing.get("tables") or {}
    frame = tables.get("de")
    if frame is None or not isinstance(frame, pd.DataFrame) or frame.empty:
        return {
            "frame": None,
            "failed": True,
            "error": "R 未返回差异结果表（tables/de.csv 缺失或为空）",
            "code": "internal",
        }
    return {
        "frame": frame,
        "failed": False,
        "warnings": list(outgoing.get("warnings") or []),
        "metrics": ((outgoing.get("result") or {}).get("summary") or {}).get("metrics") or {},
        "figures": list(outgoing.get("figures") or []),
        "ms": outgoing.get("ms"),
    }


# ---------------------------------------------------------------------------
# Python 回退实现
# ---------------------------------------------------------------------------

def _welch_paired(
    X: np.ndarray,
    genes: list[str],
    samples: list[str],
    case_mask: np.ndarray,
    ctrl_mask: np.ndarray,
    pair_col: str,
    info: pd.DataFrame,
    params: dict,
) -> pd.DataFrame:
    """配对 t 检验（按 pair_col 对齐后取组内差）。"""
    pairs: list[tuple[int, int]] = []
    case_idx = [i for i in range(len(samples)) if case_mask[i]]
    ctrl_idx = [i for i in range(len(samples)) if ctrl_mask[i]]
    for c in case_idx:
        key = str(info.loc[samples[c], pair_col])
        for k in ctrl_idx:
            if str(info.loc[samples[k], pair_col]) == key:
                pairs.append((c, k))
                break
    if len(pairs) < 2:
        return pd.DataFrame()
    ci = np.array([p[0] for p in pairs])
    ki = np.array([p[1] for p in pairs])
    diff = X[ci] - X[ki]
    n = diff.shape[0]
    mean_d = diff.mean(axis=0)
    sd_d = diff.std(axis=0, ddof=1)
    se = sd_d / math.sqrt(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, mean_d / se, 0.0)
    try:
        from scipy import stats

        pvals = 2 * stats.t.sf(np.abs(tstat), max(n - 1, 1))
    except Exception:  # noqa: BLE001
        pvals = np.array([2 * (1 - 0.5 * (1 + math.erf(abs(v) / math.sqrt(2)))) for v in tstat])
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    frame = pd.DataFrame(
        {
            "gene": genes,
            "logFC": mean_d,
            "AveExpr": X.mean(axis=0),
            "t": tstat,
            "P.Value": pvals,
            "adj.P.Val": bh_fdr(pvals),
            "B": np.zeros(len(genes)),
        }
    )
    return frame


def _ols_fallback(
    X: np.ndarray,
    genes: list[str],
    labels: np.ndarray,
    case: str,
    control: str,
    case_mask: np.ndarray,
    ctrl_mask: np.ndarray,
    info: pd.DataFrame,
    samples: list[str],
    covars: list[str],
) -> pd.DataFrame:
    """协变量校正后的组间比较：每基因 OLS（group + covars），取 group 系数。"""
    used = [s for s in samples if case_mask[samples.index(s)] or ctrl_mask[samples.index(s)]]
    keep_mask = np.array([case_mask[i] or ctrl_mask[i] for i in range(len(samples))])
    y = X[keep_mask]
    sub_labels = [labels[i] for i in range(len(samples)) if keep_mask[i]]
    n = y.shape[0]
    columns = [np.ones(n)]
    group_col = np.array([1.0 if lab == case else 0.0 for lab in sub_labels])
    columns.append(group_col)
    names = ["intercept", "group"]
    for covar in covars:
        if covar not in info.columns:
            continue
        values = np.array([str(info.loc[samples[i], covar]) for i in range(len(samples)) if keep_mask[i]])
        uniques = list(dict.fromkeys(values.tolist()))  # noqa: F841 - 仅用于判断
        if len(set(values.tolist())) < 2:
            continue
        for level in sorted(set(values.tolist()))[1:]:
            columns.append((values == level).astype(float))
            names.append(f"{covar}={level}")
    Z = np.column_stack(columns)
    rank = np.linalg.matrix_rank(Z)
    if rank < Z.shape[1] or n - rank <= 0:
        # 设计矩阵退化：丢弃协变量重新拟合
        Z = np.column_stack([np.ones(n), group_col])
        names = ["intercept", "group"]
        rank = np.linalg.matrix_rank(Z)
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
    resid = y - Z @ beta
    dof = max(n - rank, 1)
    sigma2 = (resid ** 2).sum(axis=0) / dof
    xtx_inv = np.linalg.pinv(Z.T @ Z)
    coef = beta[1]
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[1, 1], 1e-300))
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, coef / se, 0.0)
    try:
        from scipy import stats

        pvals = 2 * stats.t.sf(np.abs(tstat), dof)
    except Exception:  # noqa: BLE001
        pvals = np.array([2 * (1 - 0.5 * (1 + math.erf(abs(v) / math.sqrt(2)))) for v in tstat])
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    return pd.DataFrame(
        {
            "gene": genes,
            "logFC": coef,
            "AveExpr": X.mean(axis=0),
            "t": tstat,
            "P.Value": pvals,
            "adj.P.Val": bh_fdr(pvals),
            "B": np.zeros(len(genes)),
        }
    )


def _welch_fallback(
    X: np.ndarray,
    genes: list[str],
    labels: list[str],
    samples: list[str],
    case: str,
    control: str,
    info: pd.DataFrame,
    params: dict,
    pair_col: str,
) -> tuple[pd.DataFrame, list[str]]:
    notes: list[str] = []
    case_mask = np.array([lab == case for lab in labels])
    ctrl_mask = np.array([lab == control for lab in labels])
    if params.get("paired") and pair_col:
        frame = _welch_paired(X, genes, samples, case_mask, ctrl_mask, pair_col, info, params)
        if not frame.empty:
            notes.append(f"配对 t 检验（配对列：{pair_col}）")
            return frame, notes
        notes.append("配对列无法组成配对，改为非配对检验")
    covars = [c for c in (params.get("covars") or []) if c != params.get("group_col")]
    if covars:
        frame = _ols_fallback(X, genes, labels, case, control, case_mask, ctrl_mask, info, samples, covars)
        if not frame.empty:
            notes.append(f"协变量校正 OLS（协变量：{'、'.join(covars)}）")
            return frame, notes
    frame = welch_de(X, labels, case, control, genes=genes, lfc=params.get("lfc", 0.0), padj_cut=params.get("padj", 0.05))
    notes.append("Welch t 检验 + BH FDR")
    return frame, notes


# ---------------------------------------------------------------------------
# 结果标准化 / 排序
# ---------------------------------------------------------------------------

def _standardize(frame: pd.DataFrame, lfc: float, padj_cut: float) -> pd.DataFrame:
    frame = frame.copy()
    rename = {
        "log2FoldChange": "logFC",
        "logFC": "logFC",
        "pvalue": "P.Value",
        "PValue": "P.Value",
        "p.value": "P.Value",
        "padj": "adj.P.Val",
        "FDR": "adj.P.Val",
        "qvalue": "adj.P.Val",
        "stat": "t",
        "t": "t",
        "baseMean": "AveExpr",
        "AveExpr": "AveExpr",
        "B": "B",
        "gene": "gene",
        "Gene": "gene",
        "gene_id": "gene",
    }
    frame = frame.rename(columns={k: v for k, v in rename.items() if k in frame.columns})
    if "gene" not in frame.columns:
        if frame.index.name in ("gene", "Gene", "rownames") or not pd.api.types.is_integer_dtype(frame.index):
            frame = frame.reset_index()
            first = frame.columns[0]
            if first != "gene":
                frame = frame.rename(columns={first: "gene"})
        else:
            raise EngineError("internal", "差异结果缺少 gene 列")
    for column in DE_COLUMNS:
        if column not in frame.columns:
            if column == "regulation":
                continue
            if column == "B":
                frame["B"] = 0.0
            else:
                frame[column] = np.nan
    frame["gene"] = frame["gene"].astype(str)
    for column in ("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["B"] = frame["B"].fillna(0.0)
    # 同一基因去重（保留最小 P 值）
    if frame["gene"].duplicated().any():
        frame = (
            frame.sort_values(["P.Value"], na_position="last")
            .drop_duplicates("gene", keep="first")
            .reset_index(drop=True)
        )
    up = (frame["adj.P.Val"] <= padj_cut) & (frame["logFC"] >= max(lfc, 0)) & frame["logFC"].notna()
    down = (frame["adj.P.Val"] <= padj_cut) & (frame["logFC"] <= -max(lfc, 1e-9)) & frame["logFC"].notna()
    up = up & ~(lfc <= 0 and frame["logFC"] <= 0)
    down = down & ~(lfc <= 0 and frame["logFC"] >= 0)
    frame["regulation"] = np.where(up, "up", np.where(down, "down", "ns"))
    return frame[DE_COLUMNS]


def _sort_frame(frame: pd.DataFrame, order: str) -> pd.DataFrame:
    order = str(order or "padj").lower()
    if order == "pvalue":
        frame = frame.sort_values(["P.Value", "adj.P.Val"], ascending=[True, True], na_position="last")
    elif order == "logfc":
        frame = frame.assign(_abs=frame["logFC"].abs()).sort_values(["_abs", "P.Value"], ascending=[False, True], na_position="last")
        frame = frame.drop(columns=["_abs"])
    elif order == "aveexpr":
        frame = frame.sort_values(["AveExpr"], ascending=[False], na_position="last")
    elif order == "none":
        pass
    else:
        frame = frame.sort_values(["adj.P.Val", "P.Value"], ascending=[True, True], na_position="last")
    frame = frame.loc[:, ~frame.columns.duplicated()]
    return frame.reset_index(drop=True)


# ---------------------------------------------------------------------------
# 核心计算
# ---------------------------------------------------------------------------

def _compute_de(ctx: dict, params: dict, outdir: str) -> dict:
    """差异分析核心：返回统一结构的 dict（供 run_de 与各分析复用）。"""
    p = _norm_params(params, ctx)
    X, genes, samples = _matrix(ctx)
    info_all = _info(ctx, samples)
    group_col, labels, uniq, samples = _resolve_groups(ctx, p)
    case, control, cnotes = _resolve_case_control(p, uniq, ctx, group_col)
    labels_arr = np.array(labels)
    case_mask = labels_arr == case
    ctrl_mask = labels_arr == control
    n_case = int(case_mask.sum())
    n_ctrl = int(ctrl_mask.sum())
    if n_case < 2 or n_ctrl < 2:
        raise prerequisite(
            f"分组「{case}」({n_case}) 或「{control}」({n_ctrl}) 的有效样本不足 2 个，无法做差异检验"
        )
    pair_col = ""
    if p["paired"]:
        pair_col = _find_pair_col(info_all, samples, case_mask, ctrl_mask)
    notes: list[str] = list(cnotes)
    warnings: list[str] = []
    method = p["method"]
    effective = method
    lite = False
    r_metrics: dict = {}
    covars = [c for c in p["covars"] if c in info_all.columns and c != group_col]
    dropped = [c for c in p["covars"] if c not in info_all.columns]
    if dropped:
        warnings.append(f"协变量 {('、'.join(dropped))} 不在样本信息里，已忽略")
    params_used = dict(p)
    params_used["covars"] = covars
    params_used["pair_col"] = pair_col

    if method == "deseq2":
        counts = _counts_matrix(ctx, X)
        if counts is None:
            raise EngineError(
                "prerequisite",
                "DESeq2 需要原始 counts（整数计数矩阵）；当前矩阵类型为 "
                f"{ctx.get('X_kind') or 'unknown'}。请改用 limma/welch，或导入原始 counts 后先跑「标准化」再试。",
            )
        counts = np.clip(np.nan_to_num(counts, nan=0.0), 0, None)
        if not np.allclose(counts[: min(counts.size, 5000)], np.round(counts[: min(counts.size, 5000)])):
            warnings.append("counts 含非整数，已四舍五入后再交给 DESeq2")
        prepared = _r_diff(
            ctx, params_used, outdir, expr=counts, genes=genes, samples=samples, info=info_all,
            method="deseq2", pair_col=pair_col, case=case, control=control, group_col=group_col,
        )
        if not prepared["failed"]:
            frame = prepared["frame"]
            r_metrics = prepared.get("metrics") or {}
            warnings.extend(prepared.get("warnings") or [])
            notes.append("DESeq2 负二项 GLM + Wald 检验（R 官方实现）")
            frame = frame.copy()
            frame["B"] = 0.0
        else:
            effective, lite = "welch", True
            log_cpm = _log_cpm(counts)
            frame, extra_notes = _welch_fallback(
                log_cpm, genes, labels, samples, case, control, info_all, params_used, pair_col
            )
            notes.extend(extra_notes)
            warnings.append(f"DESeq2（R）不可用：{prepared.get('error')}；已回退 Welch t 检验（lite）")
    elif method == "limma":
        Xlog, log_notes = _to_log_scale(X, ctx, "limma")
        notes.extend(log_notes)
        prepared = _r_diff(
            ctx, params_used, outdir, expr=Xlog, genes=genes, samples=samples, info=info_all,
            method="limma", pair_col=pair_col, case=case, control=control, group_col=group_col,
        )
        if not prepared["failed"]:
            frame = prepared["frame"]
            r_metrics = prepared.get("metrics") or {}
            warnings.extend(prepared.get("warnings") or [])
            notes.append("limma::lmFit + eBayes + topTable（R 官方实现）")
        else:
            effective, lite = "welch", True
            frame, extra_notes = _welch_fallback(
                Xlog, genes, labels, samples, case, control, info_all, params_used, pair_col
            )
            notes.extend(extra_notes)
            warnings.append(f"limma（R）不可用：{prepared.get('error')}；已回退 Welch t 检验（lite）")
    else:
        Xlog, log_notes = _to_log_scale(X, ctx, "welch")
        notes.extend(log_notes)
        frame, extra_notes = _welch_fallback(
            Xlog, genes, labels, samples, case, control, info_all, params_used, pair_col
        )
        notes.extend(extra_notes)
        effective = "welch"

    frame = _standardize(frame, p["lfc"], p["padj"])
    frame = _sort_frame(frame, p["order"])
    n_up = int((frame["regulation"] == "up").sum())
    n_down = int((frame["regulation"] == "down").sum())
    return {
        "frame": frame,
        "params": p,
        "params_used": params_used,
        "method_requested": method,
        "method": effective,
        "lite": lite,
        "group_col": group_col,
        "case": case,
        "control": control,
        "labels": labels,
        "groups": uniq,
        "samples": samples,
        "genes": genes,
        "n_case": n_case,
        "n_ctrl": n_ctrl,
        "n_up": n_up,
        "n_down": n_down,
        "n_sig": n_up + n_down,
        "n_genes": int(len(frame)),
        "notes": notes,
        "warnings": warnings,
        "r_metrics": r_metrics,
        "pair_col": pair_col,
        "X_log": None,
    }


# ---------------------------------------------------------------------------
# 绘图
# ---------------------------------------------------------------------------

def _group_colors(labels) -> dict:
    names = list(dict.fromkeys([str(v) for v in labels]))
    colors = palette(len(names))
    return {name: colors[i % len(colors)] for i, name in enumerate(names)}


def _fig_volcano(frame: pd.DataFrame, params: dict, outdir: str, title_suffix: str = "") -> dict:
    lfc_cut = float(params.get("lfc") or 0.0)
    if lfc_cut <= 0:
        lfc_cut = 1.0
    padj_cut = float(params.get("padj") or 0.05)
    top_n = int(params.get("top_n") or 20)
    frame = frame.copy()
    frame["neglog"] = -np.log10(np.clip(pd.to_numeric(frame["adj.P.Val"], errors="coerce").fillna(1).values, 1e-300, 1))
    sig_up = frame["regulation"] == "up"
    sig_down = frame["regulation"] == "down"
    others = ~(sig_up | sig_down)
    fig, ax = new_fig(7.9, 6.0)
    ax.scatter(frame.loc[others, "logFC"], frame.loc[others, "neglog"], s=9, color="#b9bfc8", alpha=0.55, label="ns")
    ax.scatter(frame.loc[sig_up, "logFC"], frame.loc[sig_up, "neglog"], s=15, color="#c0392b", alpha=0.82,
               label=f"up ({int(sig_up.sum())})")
    ax.scatter(frame.loc[sig_down, "logFC"], frame.loc[sig_down, "neglog"], s=15, color="#2e6da4", alpha=0.82,
               label=f"down ({int(sig_down.sum())})")
    ax.axhline(-math.log10(max(padj_cut, 1e-300)), color="#8a919c", linestyle="--", linewidth=0.8)
    for value in (lfc_cut, -lfc_cut):
        ax.axvline(value, color="#8a919c", linestyle="--", linewidth=0.8)
    hits = frame.loc[sig_up | sig_down].nsmallest(top_n, "adj.P.Val")
    if hits.empty:
        hits = frame.nsmallest(min(top_n, 10), "adj.P.Val")
    for _, row in hits.iterrows():
        if pd.isna(row["logFC"]) or pd.isna(row["neglog"]):
            continue
        ax.annotate(str(row["gene"]), (row["logFC"], row["neglog"]), fontsize=6.2, alpha=0.9,
                    xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(adjusted P)")
    ax.set_title("Volcano plot" + (f" — {title_suffix}" if title_suffix else ""))
    ax.legend(loc="upper left", fontsize=8)
    return save_fig(fig, outdir, "volcano.png", "Volcano plot")


def _fig_ma(frame: pd.DataFrame, params: dict, outdir: str, title_suffix: str = "") -> dict:
    frame = frame.copy()
    padj_cut = float(params.get("padj") or 0.05)
    sig = (pd.to_numeric(frame["adj.P.Val"], errors="coerce").fillna(1) <= padj_cut)
    fig, ax = new_fig(7.6, 5.4)
    ax.scatter(frame["AveExpr"], frame["logFC"], s=9, color="#7f8c9b", alpha=0.45, label="ns" if sig.any() else "all")
    ax.scatter(frame.loc[sig, "AveExpr"], frame.loc[sig, "logFC"], s=12, color="#c0392b", alpha=0.72,
               label=f"FDR<{padj_cut:g} ({int(sig.sum())})")
    ax.axhline(0, color="#2b2f36", linewidth=0.8)
    ax.set_xlabel("Average log2 expression")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot" + (f" — {title_suffix}" if title_suffix else ""))
    ax.legend(fontsize=8)
    return save_fig(fig, outdir, "ma_plot.png", "MA plot")


def _fig_de_heatmap(ctx: dict, frame: pd.DataFrame, outdir: str, params: dict) -> dict | None:
    genes = [str(g) for g in (ctx.get("genes") or [])]
    samples = [str(s) for s in (ctx.get("samples") or [])]
    if not genes or not samples:
        return None
    wanted: list[str] = []
    sig = frame[frame["regulation"].isin(["up", "down"])].sort_values("adj.P.Val")
    ordered = list(sig["gene"].astype(str)) if not sig.empty else list(frame.sort_values("adj.P.Val")["gene"].astype(str))
    gene_index = {g: i for i, g in enumerate(genes)}
    for gene in ordered:
        if gene in gene_index and gene not in wanted:
            wanted.append(gene)
        if len(wanted) >= 50:
            break
    if not wanted:
        return None
    X = np.asarray(ctx.get("X"), dtype=np.float64)
    mat = X[:, [gene_index[g] for g in wanted]].T
    z = zscore_rows(mat)
    labels = np.array([str(v) for v in _ctx_labels(ctx)])
    group_names = list(dict.fromkeys(labels.tolist())) if labels.size else []
    if labels.size == len(samples) and group_names:
        order = np.argsort([group_names.index(str(v)) for v in labels])
    else:
        order = np.arange(len(samples))
    z = z[:, order]
    row_labels = [str(g) for g in wanted]
    fig, ax = new_fig(9.2, max(3.8, len(row_labels) * 0.15 + 1.6))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, interpolation="nearest")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=5.6)
    step = max(1, len(samples) // 60)
    ticks = list(range(0, len(samples), step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([samples[i] for i in order[ticks]], rotation=90, fontsize=6)
    ax.grid(False)
    if labels.size == len(samples):
        colors = _group_colors(labels)
        for tick, idx in zip(ax.get_xticklabels(), order[ticks]):
            tick.set_color(colors.get(str(labels[idx]), "#2b2f36"))
    ax.set_title("Top differential genes (row z-score)")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    return save_fig(fig, outdir, "de_heatmap.png", "DE heatmap")


def _ctx_labels(ctx: dict) -> list[str]:
    info = ctx.get("sample_info")
    samples = [str(s) for s in (ctx.get("samples") or [])]
    key = ctx.get("group_key")
    if isinstance(info, pd.DataFrame) and key and key in info.columns:
        out = []
        for sample in samples:
            try:
                value = info.loc[sample, key]
            except Exception:  # noqa: BLE001
                return ["all"] * len(samples)
            if isinstance(value, pd.Series):
                value = value.iloc[0]
            out.append(str(value))
        return out
    return ["all"] * len(samples)


def _de_figures(ctx: dict, frame: pd.DataFrame, params: dict, outdir: str) -> tuple[list[dict], list[str]]:
    warnings: list[str] = []
    figures: list[dict] = []
    suffix = ""
    try:
        figures.append(_fig_volcano(frame, params, outdir, suffix))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"火山图绘制失败：{error}")
    try:
        figures.append(_fig_ma(frame, params, outdir, suffix))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"MA 图绘制失败：{error}")
    try:
        heatmap = _fig_de_heatmap(ctx, frame, outdir, params)
        if heatmap:
            figures.append(heatmap)
    except Exception as error:  # noqa: BLE001
        warnings.append(f"差异基因热图绘制失败：{error}")
    return figures, warnings


# ---------------------------------------------------------------------------
# 表格 / 输出
# ---------------------------------------------------------------------------

def _de_tables(frame: pd.DataFrame, outdir: str) -> tuple[dict, list[dict], str, str]:
    ensure_dir(outdir)
    all_path = export_csv(frame[DE_COLUMNS], outdir, "de_all", index=False)
    sig = frame[frame["regulation"].isin(["up", "down"])]
    sig_path = export_csv(sig[DE_COLUMNS], outdir, "de_sig", index=False)
    note = f"共 {len(frame)} 个基因；显著 {len(sig)} 个（完整表：de_all.csv / de_sig.csv）"
    table_obj = table("de", "差异基因（前 500 行）", DE_COLUMNS_META, frame[DE_COLUMNS], note=note)
    data_out = [
        data_file("de_csv", all_path, "全部基因差异结果"),
        data_file("de_sig_csv", sig_path, "显著差异基因"),
    ]
    return table_obj, data_out, all_path, sig_path


def _de_bullets(ctx: dict, computed: dict, frame: pd.DataFrame) -> list[str]:
    p = computed["params"]
    method_label = METHOD_LABELS.get(computed["method"], computed["method"])
    if computed["lite"] and computed["method"] == "welch" and computed["method_requested"] != "welch":
        method_label += "（回退）"
    bullets = [
        f"{method_label}：{computed['case']} vs {computed['control']}（{computed['n_case']} vs {computed['n_ctrl']} 样本）",
        "FDR<{:g} 且 |logFC|≥{:g}：共 {} 个差异基因（上调 {} / 下调 {}）".format(
            p["padj"], p["lfc"], computed["n_sig"], computed["n_up"], computed["n_down"]
        ),
    ]
    top = frame[frame["regulation"] != "ns"].head(5)
    if not top.empty:
        pieces = [f"{row['gene']}({row['logFC']:+.2f})" for _, row in top.iterrows()]
        bullets.append("最显著基因：" + "、".join(pieces))
    if computed.get("lite"):
        bullets.append("注意：R 官方实现不可用，结果来自 Python 回退实现（lite）")
    return bullets


def _de_metrics(computed: dict) -> dict:
    p = computed["params"]
    return {
        "method": computed["method"],
        "method_requested": computed["method_requested"],
        "lite": bool(computed["lite"]),
        "n_genes": computed["n_genes"],
        "n_case": computed["n_case"],
        "n_control": computed["n_ctrl"],
        "n_sig": computed["n_sig"],
        "n_up": computed["n_up"],
        "n_down": computed["n_down"],
        "lfc": p["lfc"],
        "padj": p["padj"],
        "paired": bool(p["paired"]),
        "covars": list(p.get("covars") or []),
    }


# ---------------------------------------------------------------------------
# 步骤入口
# ---------------------------------------------------------------------------

def run_de(ctx: dict, params: dict, outdir: str) -> dict:
    """步骤 ``de``：差异表达（limma / DESeq2 / Welch）。"""
    ensure_dir(outdir)
    computed = _compute_de(ctx, params, outdir)
    frame = computed["frame"]
    p = computed["params"]
    if computed.get("lite") and computed["method_requested"] == "limma":
        p["lfc"] = p["lfc"]

    figures, figure_warnings = _de_figures(ctx, frame, computed["params"], outdir)
    table_obj, data_out, all_path, sig_path = _de_tables(frame, outdir)
    warnings = list(computed["warnings"]) + figure_warnings
    notes = list(computed["notes"])
    if computed["method"] == "welch" and computed["method_requested"] != "welch":
        notes.append("limma/DESeq2 官方实现不可用 → Python Welch 回退（lite），结果仅供探索")
    if computed["method"] == "deseq2" and computed["method_requested"] == "deseq2":
        notes.append("B 列对 DESeq2 不适用，统一填 0")

    store = {
        "group_col": computed["group_col"],
        "case": computed["case"],
        "control": computed["control"],
        "method": computed["method"],
        "method_requested": computed["method_requested"],
        "lite": bool(computed["lite"]),
        "table": frame,
        "params": {
            "lfc": p["lfc"],
            "padj": p["padj"],
            "paired": bool(p["paired"]),
            "pair_col": computed.get("pair_col") or "",
            "covars": list(p.get("covars") or []),
            "order": p["order"],
        },
        "labels": list(computed["labels"]),
        "samples": list(computed["samples"]),
        "groups": list(computed["groups"]),
        "n_case": computed["n_case"],
        "n_ctrl": computed["n_ctrl"],
        "n_up": computed["n_up"],
        "n_down": computed["n_down"],
        "n_sig": computed["n_sig"],
        "csv_all": all_path,
        "csv_sig": sig_path,
    }
    payload = result(
        title="差异表达",
        bullets=_de_bullets(ctx, computed, frame),
        metrics=_de_metrics(computed),
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=data_out,
    )
    payload["de_store"] = store
    return payload


# ---------------------------------------------------------------------------
# 分析：多组两两差异 / Venn / 摘要 / WGCNA
# ---------------------------------------------------------------------------

def _parse_contrasts(text: str, groups: list[str]) -> list[tuple[str, str]]:
    """解析 "case:control, case2:control2"；为空时返回所有两两组合（后者为实验组）。"""
    pairs: list[tuple[str, str]] = []
    for chunk in re.split(r"[,;\n]+", str(text or "")):
        chunk = chunk.strip()
        if not chunk:
            continue
        for sep in (":", "/", " vs ", " VS ", "-"):
            if sep in chunk:
                left, right = chunk.split(sep, 1)
                left, right = left.strip(), right.strip()
                if left and right:
                    pairs.append((left, right))
                break
    if pairs:
        return pairs
    clean = [g for g in groups if g]
    return [(clean[i], clean[j]) for i in range(len(clean)) for j in range(i)]


def analysis_de_multi(ctx: dict, params: dict, outdir: str) -> dict:
    """多组两两差异（一次跑完所有对比，输出合并表与对比热图）。"""
    ensure_dir(outdir)
    p = _norm_params(params, ctx)
    if not p["group_col"]:
        group_col, labels, uniq, samples = _resolve_groups(ctx, p)
    else:
        group_col, labels, uniq, samples = _resolve_groups(ctx, p)
    pairs = _parse_contrasts(params.get("contrasts"), uniq)
    if not pairs:
        raise prerequisite("分组只有一个取值，无法做多组差异（至少 2 组）")
    pairs = pairs[:6]
    for case, control in pairs:
        if case not in uniq or control not in uniq:
            raise prerequisite(f"对比「{case} vs {control}」的组名不在分组列里（可选：{'、'.join(uniq)}）")
    method = _p_str(params.get("method"), "limma").lower()
    frames: list[pd.DataFrame] = []
    warnings: list[str] = []
    notes: list[str] = []
    method_used = method
    for case, control in pairs:
        sub_params = dict(params)
        sub_params.update({"group_col": group_col, "case": case, "control": control, "method": method})
        computed = _compute_de(ctx, sub_params, outdir)
        frame = computed["frame"].copy()
        frame.insert(0, "contrast", f"{case} vs {control}")
        frame.insert(1, "case", case)
        frame.insert(2, "control", control)
        frames.append(frame)
        method_used = computed["method"]
        warnings.extend(computed["warnings"])
        notes.extend(computed["notes"])
        if computed["lite"]:
            method_used = "welch"
    merged = pd.concat(frames, ignore_index=True)
    columns = ["contrast", "gene", "logFC", "P.Value", "adj.P.Val", "regulation"]
    merged = merged[columns]
    csv_path = export_csv(merged, outdir, "de_multi_all", index=False)
    counts = (
        merged.assign(sig=merged["regulation"].isin(["up", "down"]))
        .groupby(["contrast", "regulation"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for column in ("up", "down", "ns"):
        if column not in counts.columns:
            counts[column] = 0
    counts["n_sig"] = counts["up"] + counts["down"]
    counts_path = export_csv(counts, outdir, "de_multi_counts", index=False)

    figures: list[dict] = []
    try:
        fig, ax = new_fig(8.4, 4.6)
        positions = np.arange(len(counts))
        ax.bar(positions - 0.2, counts["up"], width=0.4, color="#c0392b", alpha=0.85, label="up")
        ax.bar(positions + 0.2, counts["down"], width=0.4, color="#2e6da4", alpha=0.85, label="down")
        ax.set_xticks(positions)
        ax.set_xticklabels([str(c) for c in counts["contrast"]], rotation=20, ha="right", fontsize=8)
        ax.set_ylabel("Number of DE genes")
        ax.set_title("Differential genes per contrast")
        ax.legend(fontsize=8)
        figures.append(save_fig(fig, outdir, "de_multi_counts.png", "DE per contrast"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"对比计数图绘制失败：{error}")

    try:
        top_genes = (
            merged[merged["regulation"] != "ns"]
            .groupby("gene")["adj.P.Val"]
            .min()
            .nsmallest(30)
            .index.tolist()
        )
        if not top_genes:
            top_genes = merged.groupby("gene")["adj.P.Val"].min().nsmallest(20).index.tolist()
        pivot = merged[merged["gene"].isin(top_genes)].pivot_table(
            index="gene", columns="contrast", values="logFC", aggfunc="mean"
        )
        pivot = pivot.loc[top_genes]
        fig, ax = new_fig(8.6, max(4.0, 0.18 * len(pivot) + 1.6))
        im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
        ax.set_yticks(range(len(pivot)))
        ax.set_yticklabels([str(g) for g in pivot.index], fontsize=6)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([str(c) for c in pivot.columns], rotation=20, ha="right", fontsize=8)
        ax.grid(False)
        ax.set_title("log2FC across contrasts (top genes)")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        figures.append(save_fig(fig, outdir, "de_multi_heatmap.png", "logFC heatmap"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"对比热图绘制失败：{error}")

    bullets = [
        f"共 {len(pairs)} 个对比：" + "、".join(f"{c} vs {k}" for c, k in pairs),
        f"方法 {method_used}；合并后 {len(merged)} 行差异结果",
    ]
    if len(counts):
        detail = "；".join(f"{row['contrast']}：↑{int(row['up'])}/↓{int(row['down'])}" for _, row in counts.iterrows())
        bullets.append("差异基因数：" + detail)
    table_obj = table(
        "de_multi", "多对比差异结果（前 500 行）",
        [{"key": "contrast", "label": "对比", "type": "string"}, {"key": "gene", "label": "基因", "type": "string"},
         {"key": "logFC", "label": "log2FC", "type": "number"}, {"key": "P.Value", "label": "P 值", "type": "number"},
         {"key": "adj.P.Val", "label": "FDR", "type": "number"}, {"key": "regulation", "label": "调控", "type": "string"}],
        merged, note="完整表见 de_multi_all.csv",
    )
    payload = result(
        title="多组差异对比",
        bullets=bullets,
        metrics={
            "method": method_used,
            "n_contrasts": len(pairs),
            "n_rows": int(len(merged)),
            "pairs": [f"{c} vs {k}" for c, k in pairs],
        },
        tables=[table_obj, table("de_multi_counts", "各对比差异基因数", [{"key": c, "label": c, "type": "number" if c != "contrast" else "string"} for c in counts.columns], counts)],
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=[data_file("de_csv", csv_path, "多对比差异结果"), data_file("de_csv", counts_path, "各对比差异基因数")],
    )
    payload["store"] = {
        "de_multi": {
            "contrasts": [f"{c} vs {k}" for c, k in pairs],
            "table": merged.to_dict(orient="records"),
            "csv": csv_path,
            "method": method_used,
        }
    }
    return payload


def _venn_counts(sets: list[set], names: list[str]) -> pd.DataFrame:
    """枚举 2/3 组的所有交集区域计数。"""
    n = len(sets)
    rows = []
    all_genes = set()
    for s in sets:
        all_genes |= s
    for gene in sorted(all_genes):
        membership = tuple(i for i in range(n) if gene in sets[i])
        rows.append(membership)
    from collections import Counter

    counter = Counter(rows)
    out = []
    for membership, count in sorted(counter.items(), key=lambda item: (-len(item[0]), item[0])):
        label = " & ".join(names[i] for i in membership) if membership else "none"
        out.append({"region": label or "none", "n_genes": count, "membership": ",".join(str(i + 1) for i in membership)})
    return pd.DataFrame(out)


def analysis_de_venn(ctx: dict, params: dict, outdir: str) -> dict:
    """多对比差异基因 Venn/交集分析（2–3 个对比）。"""
    ensure_dir(outdir)
    p = _norm_params(params, ctx)
    group_col, labels, uniq, samples = _resolve_groups(ctx, p)
    direction = _p_str(params.get("direction"), "both").lower()
    if direction not in ("up", "down", "both"):
        direction = "both"
    pairs = _parse_contrasts(params.get("contrasts"), uniq)[:3]
    if not pairs:
        raise prerequisite("至少需要 1 个对比（分组需 ≥2 组）")
    method = _p_str(params.get("method"), "welch").lower()
    if method not in METHODS:
        method = "welch"
    names: list[str] = []
    sets: list[set] = []
    frames: list[pd.DataFrame] = []
    warnings: list[str] = []
    notes: list[str] = []
    for case, control in pairs:
        sub_params = dict(params)
        sub_params.update({"group_col": group_col, "case": case, "control": control, "method": method,
                           "lfc": p["lfc"], "padj": p["padj"]})
        computed = _compute_de(ctx, sub_params, outdir)
        frame = computed["frame"]
        frames.append(frame.assign(contrast=f"{case} vs {control}"))
        warnings.extend(computed["warnings"])
        notes.extend(computed["notes"])
        genes, _stats = de_gene_lists(frame, direction=direction, padj=p["padj"], lfc=p["lfc"])
        names.append(f"{case} vs {control}")
        sets.append(set(genes))
    if len(sets) < 2:
        raise prerequisite("Venn 分析需要至少 2 个对比（分组取值 ≥3 组，或指定 contrasts）")
    region = _venn_counts(sets, names)
    merged = pd.concat([frame[["contrast", "gene", "logFC", "adj.P.Val", "regulation"]] for frame in frames], ignore_index=True)
    merged = merged[merged["regulation"].isin(["up", "down"])]
    csv_path = export_csv(merged, outdir, "de_venn_genes", index=False)
    region_path = export_csv(region, outdir, "de_venn_regions", index=False)

    figures: list[dict] = []
    try:
        fig, ax = new_fig(6.4, 5.6)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis("off")
        colors = palette(len(sets))
        if len(sets) == 2:
            circles = [((4.3, 5.0), 3.0), ((5.7, 5.0), 3.0)]
        else:
            circles = [((4.2, 5.6), 2.7), ((5.8, 5.6), 2.7), ((5.0, 4.0), 2.7)]
        from matplotlib.patches import Circle

        for (center, radius), name, color in zip(circles, names, colors):
            ax.add_patch(Circle(center, radius, color=color, alpha=0.22, linewidth=1.2, edgecolor=color))
        for label, (x, y) in _venn_label_points(len(sets), region, names):
            ax.text(x, y, label, ha="center", va="center", fontsize=8.5)
        ax.set_title("DE gene overlap between contrasts")
        figures.append(save_fig(fig, outdir, "de_venn.png", "DE Venn"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"Venn 图绘制失败：{error}")

    bullets = [f"{len(sets)} 个对比（{'、'.join(names)}）的差异基因交集（方向={direction}）"]
    common = set.intersection(*sets) if sets else set()
    bullets.append(f"共有差异基因 {len(common)} 个；各对比差异基因数：" +
                   "、".join(f"{n}={len(s)}" for n, s in zip(names, sets)))
    table_obj = table(
        "de_venn", "差异基因交集区域",
        [{"key": "region", "label": "交集区域", "type": "string"},
         {"key": "n_genes", "label": "基因数", "type": "number"},
         {"key": "membership", "label": "成员编号", "type": "string"}],
        region, note="完整基因列表见 de_venn_genes.csv",
    )
    payload = result(
        title="差异基因 Venn",
        bullets=bullets,
        metrics={"n_sets": len(sets), "n_common": len(common), "direction": direction, "method": method,
                 "set_sizes": [len(s) for s in sets]},
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=[data_file("de_csv", csv_path, "差异基因（按对比）"), data_file("de_csv", region_path, "交集区域计数")],
    )
    payload["store"] = {
        "de_venn": {
            "contrasts": names,
            "direction": direction,
            "common_genes": sorted(common),
            "seed": 0,
            "table": region.to_dict(orient="records"),
        }
    }
    return payload


def _venn_label_points(n: int, region: pd.DataFrame, names: list[str]):
    """给 2/3 组 Venn 图分配文字坐标（与交集计数对应）。"""
    lookup = {row["membership"]: int(row["n_genes"]) for _, row in region.iterrows()}
    points = []
    if n == 2:
        mapping = {"1": (3.4, 5.0), "2": (6.6, 5.0), "1,2": (5.0, 5.0)}
    else:
        mapping = {"1": (3.2, 6.0), "2": (6.8, 6.0), "3": (5.0, 3.0),
                   "1,2": (5.0, 6.6), "1,3": (3.9, 4.3), "2,3": (6.1, 4.3), "1,2,3": (5.0, 5.2)}
    for membership, position in mapping.items():
        count = lookup.get(membership)
        if count is None:
            continue
        points.append((f"{count}", position))
    for i, name in enumerate(names):
        if n == 2:
            points.append((name, (2.4, 8.6) if i == 0 else (7.6, 8.6)))
        else:
            points.append((name, [(2.6, 9.0), (7.4, 9.0), (5.0, 1.6)][i]))
    return points


def analysis_de_summary(ctx: dict, params: dict, outdir: str) -> dict:
    """当前差异结果的阈值概览 + top 基因。"""
    ensure_dir(outdir)
    frame = get_de(ctx)
    store = ctx.get("de") if isinstance(ctx.get("de"), dict) else {}
    padj_default = float((store.get("params") or {}).get("padj") or params.get("padj") or 0.05)
    lfc_default = float((store.get("params") or {}).get("lfc") or params.get("lfc") or 0.0)
    thresholds = [0.0, 0.5, 1.0, 1.5, 2.0]
    rows = []
    for lfc_cut in thresholds:
        genes, stats = de_gene_lists(frame, direction="both", padj=padj_default, lfc=lfc_cut)
        rows.append({"lfc_cut": lfc_cut, "n_sig": stats["n_sig"], "n_up": stats["n_up"], "n_down": stats["n_down"]})
    grid = pd.DataFrame(rows)
    top_n = max(1, min(100, _p_int(params.get("top_n"), 25)))
    top_cols = ["gene", "logFC", "AveExpr", "P.Value", "adj.P.Val", "regulation"]
    top = frame[frame["regulation"].isin(["up", "down"])].copy()
    if top.empty:
        top = frame.copy()
    top = top.nsmallest(top_n, "adj.P.Val")[top_cols]
    grid_path = export_csv(grid, outdir, "de_summary_thresholds", index=False)
    top_path = export_csv(top, outdir, "de_summary_top", index=False)
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        fig, ax = new_fig(7.4, 4.4)
        positions = np.arange(len(grid))
        ax.bar(positions, grid["n_up"], width=0.55, color="#c0392b", alpha=0.85, label="up")
        ax.bar(positions, grid["n_down"], width=0.55, color="#2e6da4", alpha=0.85, label="down",
               bottom=grid["n_up"])
        ax.set_xticks(positions)
        ax.set_xticklabels([f"|logFC|≥{v:g}" for v in grid["lfc_cut"]], fontsize=8)
        ax.set_ylabel("Number of DE genes")
        ax.set_title(f"DE genes vs logFC threshold (FDR<{padj_default:g})")
        ax.legend(fontsize=8)
        figures.append(save_fig(fig, outdir, "de_summary_counts.png", "DE threshold summary"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"阈值统计图绘制失败：{error}")
    try:
        fig, ax = new_fig(7.2, max(3.6, 0.22 * len(top) + 1.2))
        colors = ["#c0392b" if r == "up" else "#2e6da4" for r in top["regulation"]]
        ax.barh(np.arange(len(top))[::-1], top["logFC"].values, color=colors, alpha=0.85)
        ax.set_yticks(np.arange(len(top))[::-1])
        ax.set_yticklabels([str(g) for g in top["gene"]], fontsize=6.4)
        ax.axvline(0, color="#2b2f36", linewidth=0.8)
        ax.set_xlabel("log2 fold change")
        ax.set_title("Top differential genes")
        figures.append(save_fig(fig, outdir, "de_summary_top.png", "Top DE genes"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"top 基因图绘制失败：{error}")
    bullets = [
        f"当前差异结果：{len(frame)} 个基因，显著 {int((frame['regulation'] != 'ns').sum())} 个"
        f"（上调 {int((frame['regulation'] == 'up').sum())} / 下调 {int((frame['regulation'] == 'down').sum())}）",
        "FDR<{:g} 时按 |logFC| 阈值统计：".format(padj_default)
        + "，".join(f"≥{row['lfc_cut']:g}: {row['n_sig']}" for _, row in grid.iterrows()),
    ]
    table_obj = table(
        "de_summary", "阈值统计",
        [{"key": "lfc_cut", "label": "|logFC| 阈值", "type": "number"}, {"key": "n_sig", "label": "差异基因数", "type": "number"},
         {"key": "n_up", "label": "上调", "type": "number"}, {"key": "n_down", "label": "下调", "type": "number"}],
        grid,
    )
    payload = result(
        title="差异结果摘要",
        bullets=bullets,
        metrics={"n_genes": int(len(frame)), "padj": padj_default, "lfc": lfc_default},
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        data_out=[data_file("de_csv", grid_path, "阈值统计"), data_file("de_csv", top_path, "top 差异基因")],
    )
    payload["store"] = {
        "de_summary": {
            "thresholds": grid.to_dict(orient="records"),
            "top_genes": top["gene"].astype(str).tolist(),
            "padj": padj_default,
        }
    }
    return payload


# ---- WGCNA ---------------------------------------------------------------

_WGCNA_COLORS = ["turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple"]


def _wgcna_python(
    X: np.ndarray,
    genes: list[str],
    samples: list[str],
    info: pd.DataFrame,
    params: dict,
) -> dict:
    """相关矩阵 + 层次聚类 + 模块特征基因（WGCNA 的确定性回退）。"""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    max_genes = max(200, min(4000, _p_int(params.get("top_n"), 2000)))
    variances = np.nanvar(X, axis=0)
    order = np.argsort(-np.nan_to_num(variances, nan=0.0))[:max_genes]
    order = [int(i) for i in order if np.isfinite(variances[i]) and variances[i] > 0]
    if len(order) < 20:
        raise prerequisite("可用于 WGCNA 的高变基因不足 20 个")
    sub = X[:, order]
    sub_genes = [genes[i] for i in order]
    mean = sub.mean(axis=0, keepdims=True)
    std = sub.std(axis=0, ddof=1, keepdims=True)
    std = np.where(std > 0, std, 1.0)
    z = (sub - mean) / std
    corr = np.corrcoef(z, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0)
    dist = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    condensed = squareform(dist, checks=False)
    tree = linkage(condensed, method="average")
    min_size = max(10, _p_int(params.get("min_module_size"), 30))
    modules = None
    for cut in (0.5, 0.6, 0.45, 0.7, 0.8, 0.9, 1.1):
        labels = fcluster(tree, t=cut, criterion="distance")
        sizes = pd.Series(labels).value_counts()
        kept = [label for label, size in sizes.items() if size >= min_size]
        if kept:
            modules = labels
            break
    if modules is None:
        modules = np.ones(len(sub_genes), dtype=int)
    # 小模块并入 grey（0）
    sizes = pd.Series(modules).value_counts()
    kept = {label for label, size in sizes.items() if size >= min_size}
    labels_out = np.array([label if label in kept else 0 for label in modules])
    module_ids = sorted(set(labels_out.tolist()) - {0})
    module_names = {}
    for idx, label in enumerate(module_ids):
        module_names[label] = _WGCNA_COLORS[idx % len(_WGCNA_COLORS)]
    eigengenes = []
    for label in module_ids:
        rows = labels_out == label
        block = z[:, rows]
        u, s, vt = np.linalg.svd(block, full_matrices=False)
        pc1 = u[:, 0] * s[0]
        if np.corrcoef(pc1, block.mean(axis=1))[0, 1] < 0:
            pc1 = -pc1
        eigengenes.append(pc1)
    eigengene_matrix = np.column_stack(eigengenes) if eigengenes else np.zeros((len(samples), 0))
    kme = np.zeros(len(sub_genes))
    for j, label in enumerate(module_ids):
        rows = labels_out == label
        block = z[:, rows]
        eg = eigengene_matrix[:, j]
        denom = block.std(axis=0, ddof=1) * eg.std(ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            values = np.where(denom > 0, (block * (eg - eg.mean())[:, None]).mean(axis=0) / denom, 0.0)
        kme[rows] = np.nan_to_num(values, nan=0.0)
    module_table = pd.DataFrame(
        {
            "gene": sub_genes,
            "module": [("grey" if label == 0 else f"M{label}") for label in labels_out],
            "color": ["grey" if label == 0 else module_names.get(label, "grey") for label in labels_out],
            "kME": kme,
        }
    ).sort_values(["module", "kME"], ascending=[True, False])
    eigengene_table = pd.DataFrame(
        np.column_stack([np.ones(len(samples)), eigengene_matrix]),
        columns=["sample"] + [f"M{label}" for label in module_ids],
    )
    eigengene_table["sample"] = samples
    module_summary = []
    for label in module_ids:
        rows = labels_out == label
        members = module_table[module_table["module"] == f"M{label}"].sort_values("kME", ascending=False)
        module_summary.append(
            {
                "module": f"M{label}",
                "color": module_names.get(label, "grey"),
                "n_genes": int(rows.sum()),
                "hub_genes": "、".join(members["gene"].astype(str).head(8)),
                "mean_kME": float(members["kME"].mean()) if len(members) else 0.0,
            }
        )
    trait_table = None
    group_key = ctx_group_key(ctx)
    if group_key and group_key in info.columns and eigengene_matrix.size:
        values = np.array([str(info.loc[sample, group_key]) for sample in samples])
        levels = list(dict.fromkeys(values.tolist()))
        if len(levels) >= 2:
            rows = []
            for j, label in enumerate(module_ids):
                eg = eigengene_matrix[:, j]
                for level in levels[1:]:
                    trait = (values == level).astype(float)
                    if trait.std() == 0:
                        continue
                    r = float(np.corrcoef(eg, trait)[0, 1])
                    try:
                        from scipy import stats

                        p = float(2 * stats.t.sf(abs(r) * math.sqrt(max(len(eg) - 2, 1)) / math.sqrt(max(1 - r * r, 1e-12)), max(len(eg) - 2, 1)))
                    except Exception:  # noqa: BLE001
                        p = float("nan")
                    rows.append({"module": f"M{label}", "trait": f"{group_key}={level}", "cor": r, "pvalue": p})
            trait_table = pd.DataFrame(rows)
    return {
        "modules": module_table,
        "eigengenes": eigengene_table,
        "module_summary": pd.DataFrame(module_summary),
        "traits": trait_table,
        "n_genes": len(sub_genes),
        "method": "correlation",
    }


def ctx_group_key(ctx: dict) -> str:
    return str(ctx.get("group_key") or "")


def _wgcna_figures(
    modules: pd.DataFrame,
    eigengenes: pd.DataFrame,
    module_summary: pd.DataFrame,
    traits: pd.DataFrame | None,
    samples: list[str],
    group_key: str,
    info: pd.DataFrame | None,
    outdir: str,
    warnings: list[str],
) -> list[dict]:
    figures: list[dict] = []
    try:
        fig, ax = new_fig(8.0, 4.2)
        summary = module_summary.sort_values("n_genes", ascending=False)
        ax.bar(range(len(summary)), summary["n_genes"],
               color=[_color_hex(c) for c in summary["color"]], alpha=0.9)
        ax.set_xticks(range(len(summary)))
        ax.set_xticklabels([str(m) for m in summary["module"]], rotation=0, fontsize=8)
        ax.set_ylabel("Number of genes")
        ax.set_title("WGCNA module sizes")
        figures.append(save_fig(fig, outdir, "wgcna_modules.png", "WGCNA modules"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"WGCNA 模块大小图绘制失败：{error}")
    try:
        me_columns = [c for c in eigengenes.columns if c != "sample"]
        if me_columns:
            mat = eigengenes[me_columns].to_numpy(dtype=float)
            z = zscore_rows(mat.T)
            order = None
            if info is not None and group_key and group_key in info.columns:
                labels = [str(info.loc[s, group_key]) for s in samples]
                levels = list(dict.fromkeys(labels))
                order = np.argsort([levels.index(v) for v in labels])
                z = z[:, order]
            fig, ax = new_fig(max(6.4, 0.22 * len(samples) + 2.2), max(3.0, 0.5 * len(me_columns) + 1.6))
            im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
            ax.set_yticks(range(len(me_columns)))
            ax.set_yticklabels(me_columns, fontsize=8)
            ticks = list(range(0, len(samples), max(1, len(samples) // 40)))
            ax.set_xticks(ticks)
            ax.set_xticklabels([samples[i] for i in (order[ticks] if order is not None else ticks)], rotation=90, fontsize=6)
            ax.grid(False)
            ax.set_title("Module eigengenes (z-score)")
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
            figures.append(save_fig(fig, outdir, "wgcna_eigengenes.png", "Module eigengenes"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"模块特征基因热图绘制失败：{error}")
    try:
        if traits is not None and not traits.empty:
            pivot = traits.pivot_table(index="module", columns="trait", values="cor")
            fig, ax = new_fig(max(5.6, 1.1 * len(pivot.columns) + 2.0), max(3.0, 0.5 * len(pivot) + 1.6))
            im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
            ax.set_yticks(range(len(pivot)))
            ax.set_yticklabels([str(v) for v in pivot.index], fontsize=8)
            ax.set_xticks(range(len(pivot.columns)))
            ax.set_xticklabels([str(c) for c in pivot.columns], rotation=25, ha="right", fontsize=7.5)
            for i in range(len(pivot)):
                for j in range(len(pivot.columns)):
                    value = pivot.values[i, j]
                    if np.isfinite(value):
                        ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6.5)
            ax.grid(False)
            ax.set_title("Module–trait correlation")
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
            figures.append(save_fig(fig, outdir, "wgcna_traits.png", "Module-trait correlation"))
    except Exception as error:  # noqa: BLE001
        warnings.append(f"模块-性状相关图绘制失败：{error}")
    return figures


def _color_hex(name: str) -> str:
    mapping = {
        "turquoise": "#40b7ae", "blue": "#3b7dd8", "brown": "#a0522d", "yellow": "#e5c845",
        "green": "#4f9d4f", "red": "#d0453e", "black": "#3a3f46", "pink": "#e08fb0",
        "magenta": "#c04fbf", "purple": "#8a5bb8", "grey": "#b9bfc8",
    }
    return mapping.get(str(name).lower(), "#7f8c9b")


def _r_wgcna(ctx: dict, params: dict, outdir: str, X: np.ndarray, genes: list[str], samples: list[str],
             info: pd.DataFrame) -> dict | None:
    top_n = max(200, min(6000, _p_int(params.get("top_n"), 2000)))
    variances = np.nanvar(X, axis=0)
    order = np.argsort(-np.nan_to_num(variances, nan=0.0))[:top_n]
    order = [int(i) for i in order if np.isfinite(variances[i]) and variances[i] > 0]
    if len(order) < 20:
        return None
    payload = {
        "expr": X[:, order],
        "genes": [genes[i] for i in order],
        "samples": samples,
        "info": info,
        "params": {
            "mode": "wgcna",
            "min_module_size": _p_int(params.get("min_module_size"), 30),
            "soft_power": _p_int(params.get("soft_power"), 0),
            "group_col": ctx_group_key(ctx),
            "group_case": _p_str(params.get("case")),
            "group_control": _p_str(params.get("control")),
        },
    }
    try:
        outgoing = rbridge.run_r("bg_enrich.R", payload, outdir, timeout=900)
    except EngineError as error:
        return {"failed": True, "error": error.message}
    tables = outgoing.get("tables") or {}
    modules = tables.get("wgcna_modules")
    eigengenes = tables.get("wgcna_eigengenes")
    if modules is None or eigengenes is None or modules.empty or eigengenes.empty:
        return {"failed": True, "error": "R WGCNA 未返回模块或特征基因表"}
    traits = tables.get("wgcna_traits")
    summary = tables.get("wgcna_summary")
    if summary is None or summary.empty:
        summary = (
            modules.groupby("module")
            .size()
            .reset_index(name="n_genes")
            .assign(color=lambda frame: frame["module"].map(lambda m: str(m).lower()))
        )
    return {
        "failed": False,
        "modules": modules,
        "eigengenes": eigengenes,
        "module_summary": summary,
        "traits": traits,
        "method": "wgcna",
        "warnings": list(outgoing.get("warnings") or []),
    }


def analysis_wgcna_module(ctx: dict, params: dict, outdir: str) -> dict:
    """WGCNA 共表达模块（R 官方实现优先；回退相关矩阵层次聚类）。"""
    ensure_dir(outdir)
    X, genes, samples = _matrix(ctx)
    info = _info(ctx, samples)
    group_key = ctx_group_key(ctx)
    warnings: list[str] = []
    notes: list[str] = []
    method_pref = _p_str(params.get("method"), "auto").lower()
    result_data = None
    if method_pref != "correlation":
        result_data = _r_wgcna(ctx, params, outdir, X, genes, samples, info)
        if result_data and not result_data.get("failed"):
            notes.append("WGCNA::blockwiseModules（R 官方实现）")
            warnings.extend(result_data.get("warnings") or [])
        else:
            reason = (result_data or {}).get("error") or "R 不可用"
            warnings.append(f"WGCNA（R）不可用：{reason}；已回退相关矩阵层次聚类（lite）")
            notes.append("回退实现：基因相关矩阵 + 层次聚类 + 模块特征基因 PC1")
            result_data = None
    if result_data is None:
        if method_pref == "wgcna":
            raise EngineError("no_r", "WGCNA（R）不可用，且 method=wgcna 不允许回退")
        py = _wgcna_python(X, genes, samples, info, params)
        py["warnings"] = []
        result_data = {"failed": False, **py}
        notes.append("回退实现：基因相关矩阵 + 层次聚类 + 模块特征基因 PC1")
    modules = result_data["modules"]
    eigengenes = result_data["eigengenes"]
    module_summary = result_data["module_summary"]
    traits = result_data.get("traits")
    method_used = result_data.get("method") or ("wgcna" if method_pref != "correlation" else "correlation")
    modules_path = export_csv(modules, outdir, "wgcna_gene_modules", index=False)
    summary_path = export_csv(module_summary, outdir, "wgcna_module_summary", index=False)
    eigengenes_path = export_csv(eigengenes, outdir, "wgcna_eigengenes", index=False)
    data_out = [
        data_file("wgcna_csv", modules_path, "基因-模块归属"),
        data_file("wgcna_csv", summary_path, "模块汇总"),
        data_file("wgcna_csv", eigengenes_path, "模块特征基因"),
    ]
    if traits is not None and not traits.empty:
        traits_path = export_csv(traits, outdir, "wgcna_traits", index=False)
        data_out.append(data_file("wgcna_csv", traits_path, "模块-性状相关"))
    figures = _wgcna_figures(modules, eigengenes, module_summary, traits, samples, group_key, info, outdir, warnings)
    module_names = [str(m) for m in module_summary["module"].tolist()]
    bullets = [
        f"构建 {len(modules)} 个基因的共表达网络，识别 {len(module_names)} 个模块（{method_used}）",
        "模块大小：" + "、".join(
            f"{row['module']}={int(row['n_genes'])}" for _, row in module_summary.sort_values("n_genes", ascending=False).head(10).iterrows()
        ),
    ]
    if traits is not None and not traits.empty:
        strongest = traits.reindex(traits["cor"].abs().sort_values(ascending=False).index).head(3)
        bullets.append("最强模块-性状相关：" + "、".join(
            f"{row['module']}~{row['trait']} r={row['cor']:.2f}" for _, row in strongest.iterrows()))
    table_obj = table(
        "wgcna_modules", "模块汇总",
        [{"key": "module", "label": "模块", "type": "string"}, {"key": "color", "label": "颜色", "type": "string"},
         {"key": "n_genes", "label": "基因数", "type": "number"}, {"key": "hub_genes", "label": "枢纽基因", "type": "string"},
         {"key": "mean_kME", "label": "平均 kME", "type": "number"}],
        module_summary,
    )
    payload = result(
        title="WGCNA 共表达模块",
        bullets=bullets,
        metrics={
            "method": method_used,
            "n_genes": int(result_data.get("n_genes") or len(modules)),
            "n_modules": int(len(module_summary)),
            "lite": bool(method_used != "wgcna"),
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=data_out,
    )
    payload["store"] = {
        "wgcna": {
            "method": method_used,
            "modules": modules.head(4000).to_dict(orient="records"),
            "module_summary": module_summary.to_dict(orient="records"),
            "eigengene_samples": eigengenes["sample"].astype(str).tolist(),
            "eigengene_columns": [c for c in eigengenes.columns if c != "sample"],
            "eigengene_values": eigengenes.drop(columns=["sample"]).to_numpy(dtype=float).tolist(),
            "traits": (traits.to_dict(orient="records") if traits is not None and not traits.empty else []),
        }
    }
    return payload


# ---------------------------------------------------------------------------
# 注册表
# ---------------------------------------------------------------------------

_DE_STEP_PARAMS = [
    {"key": "group_col", "label": "分组列", "type": "column", "default": None},
    {"key": "case", "label": "实验组", "type": "text", "default": ""},
    {"key": "control", "label": "对照组", "type": "text", "default": ""},
    {"key": "method", "label": "方法", "type": "select", "default": "limma",
     "options": [{"value": "limma", "label": "limma（R 官方，推荐）"},
                 {"value": "deseq2", "label": "DESeq2（R 官方，需 counts）"},
                 {"value": "welch", "label": "Welch t 检验（Python）"}]},
    {"key": "covars", "label": "协变量", "type": "columns", "default": []},
    {"key": "lfc", "label": "logFC 阈值", "type": "number", "default": 0, "min": 0, "max": 5, "step": 0.1},
    {"key": "padj", "label": "FDR 阈值", "type": "number", "default": 0.05, "min": 0.000001, "max": 1, "step": 0.01},
    {"key": "paired", "label": "配对设计", "type": "boolean", "default": False},
    {"key": "top_n", "label": "图标注重基因数", "type": "number", "default": 20, "min": 1, "max": 100, "step": 1},
    {"key": "order", "label": "结果排序", "type": "select", "default": "padj", "options": ORDER_OPTIONS},
]

_DE_ANALYSIS_PARAMS = [
    {"key": "group_col", "label": "分组列", "type": "column", "default": None},
    {"key": "contrasts", "label": "对比（case:control，逗号分隔）", "type": "text", "default": ""},
    {"key": "method", "label": "方法", "type": "select", "default": "limma",
     "options": [{"value": "limma", "label": "limma（R）"},
                 {"value": "deseq2", "label": "DESeq2（R，需 counts）"},
                 {"value": "welch", "label": "Welch（Python）"}]},
    {"key": "lfc", "label": "logFC 阈值", "type": "number", "default": 0, "min": 0, "max": 5, "step": 0.1},
    {"key": "padj", "label": "FDR 阈值", "type": "number", "default": 0.05, "min": 0.000001, "max": 1, "step": 0.01},
]

STEPS = [
    {
        "id": "de",
        "title": "差异表达",
        "icon": "de",
        "requires": ["normalize"],
        "params": _DE_STEP_PARAMS,
        "run": run_de,
    }
]

ANALYSES = [
    {
        "id": "de_multi",
        "title": "多组差异对比",
        "group": "差异分析",
        "method": "逐对比 limma/DESeq2/Welch + 合并表（R 官方优先，Python 回退）",
        "needs": ["normalize", "groups"],
        "params": _DE_ANALYSIS_PARAMS,
        "run": analysis_de_multi,
    },
    {
        "id": "de_venn",
        "title": "差异基因 Venn",
        "group": "差异分析",
        "method": "多对比差异基因交集 + matplotlib Venn（≤3 组）",
        "needs": ["de"],
        "params": _DE_ANALYSIS_PARAMS + [
            {"key": "direction", "label": "方向", "type": "select", "default": "both",
             "options": [{"value": "both", "label": "上调+下调"}, {"value": "up", "label": "仅上调"},
                         {"value": "down", "label": "仅下调"}]},
        ],
        "run": analysis_de_venn,
    },
    {
        "id": "de_summary",
        "title": "差异结果摘要",
        "group": "差异分析",
        "method": "当前 DE 表阈值统计 + top 基因条形图（确定性）",
        "needs": ["de"],
        "params": [
            {"key": "top_n", "label": "top 基因数", "type": "number", "default": 25, "min": 1, "max": 100, "step": 1},
        ],
        "run": analysis_de_summary,
    },
    {
        "id": "wgcna_module",
        "title": "WGCNA 共表达模块",
        "group": "基因集合分析",
        "method": "WGCNA::blockwiseModules（R 官方）；回退相关矩阵层次聚类 + 模块特征基因",
        "needs": ["normalize"],
        "params": [
            {"key": "top_n", "label": "参与建模的高变基因数", "type": "number", "default": 2000, "min": 200, "max": 6000, "step": 100},
            {"key": "min_module_size", "label": "最小模块基因数", "type": "number", "default": 30, "min": 10, "max": 500, "step": 5},
            {"key": "soft_power", "label": "软阈值（0=自动）", "type": "number", "default": 0, "min": 0, "max": 30, "step": 1},
            {"key": "method", "label": "实现", "type": "select", "default": "auto",
             "options": [{"value": "auto", "label": "自动（R 优先）"}, {"value": "correlation", "label": "仅 Python 回退"}]},
        ],
        "run": analysis_wgcna_module,
    },
]

__all__ = [
    "STEPS", "ANALYSES", "get_de", "de_gene_lists", "run_de",
    "analysis_de_multi", "analysis_de_venn", "analysis_de_summary", "analysis_wgcna_module",
    "DE_COLUMNS", "DE_COLUMNS_META", "METHODS", "METHOD_LABELS",
]
