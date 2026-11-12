"""
dsh-geobulk — 清洗模块（步骤 qc / filter / normalize / batch / pca + 基础分析）

实现（全部纯 Python，保证任何环境可用；方法名诚实标注）：
  * qc        : 文库大小 / 检出基因数 / 中位表达 / 头号基因占比 + 离群样本（z>3）
  * filter    : 低表达基因过滤（min_expr × min_samples_frac）+ 可选剔除离群样本
  * normalize : log2(CPM+1) / log2(x+1) / 分位数标准化 / 中位数中心化 / TMM（简化）
  * batch     : ComBat（简化经验贝叶斯位置-尺度）或 removeBatchEffect（组均值中心化）
  * pca       : SVD 主成分 + 样本相关性 + 层次聚类

注册：STEPS=[qc,filter,normalize,batch,pca]；ANALYSES 8 个基础分析。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from bg_common import (
    EngineError,
    bh_fdr,
    data_file,
    export_csv,
    new_fig,
    palette,
    result,
    save_fig,
    table,
    zscore_rows,
)

BATCH_HINTS = ("batch", "batch_id", "cohort", "plate", "run", "library_batch", "dataset", "source")


# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------

def _X(ctx) -> np.ndarray:
    arr = np.asarray(ctx["X"], dtype=np.float64)
    if arr.ndim != 2 or arr.size == 0:
        raise EngineError("prerequisite", "表达矩阵为空，请先导入数据")
    return arr


def _groups(ctx, group_col=None) -> np.ndarray:
    info = ctx.get("sample_info")
    key = group_col or ctx.get("group_key")
    if info is not None and key and key in info.columns:
        return info[key].astype(str).values
    return np.array(["all"] * len(ctx["samples"]))


def _color_map(labels) -> dict:
    uniq = list(dict.fromkeys([str(x) for x in labels]))
    return {name: palette(len(uniq))[i] for i, name in enumerate(uniq)}


def _pca_scores(X: np.ndarray, n_pcs: int = 10):
    arr = np.asarray(X, dtype=np.float64)
    arr = arr - arr.mean(axis=0, keepdims=True)
    # 按列标准化会放大低表达基因噪声，这里只中心化
    U, S, Vt = np.linalg.svd(arr, full_matrices=False)
    scores = U * S
    var = (S ** 2) / max((S ** 2).sum(), 1e-12)
    n = min(n_pcs, scores.shape[1])
    return scores[:, :n], var[:n]


def _batch_column(ctx, requested=None) -> str | None:
    info = ctx.get("sample_info")
    if info is None:
        return None
    if requested and requested in info.columns:
        return str(requested)
    lower = {str(c).lower(): str(c) for c in info.columns}
    for hint in BATCH_HINTS:
        if hint in lower:
            col = lower[hint]
            if info[col].nunique(dropna=True) > 1:
                return col
    return None


def _pca_figure(scores, var, labels, outdir: str, name: str, title: str, color_by: str = "group"):
    colors = _color_map(labels)
    fig, ax = new_fig(7.4, 5.4)
    for name_value in dict.fromkeys([str(v) for v in labels]):
        mask = np.array([str(v) == name_value for v in labels])
        ax.scatter(scores[mask, 0], scores[mask, 1], s=44, alpha=0.86, color=colors[name_value],
                   edgecolor="white", linewidth=0.6, label=str(name_value))
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)" if len(var) else "PC1")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)" if len(var) > 1 else "PC2")
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)
    return save_fig(fig, outdir, name, title)


def _dist_figure(X: np.ndarray, samples, outdir: str, name: str, title: str) -> dict:
    fig, ax = new_fig(9.0, 4.8)
    values = []
    for i in range(X.shape[0]):
        row = X[i][np.isfinite(X[i])]
        values.append(row)
    try:
        bp = ax.boxplot(values, patch_artist=True, widths=0.6, showfliers=False,
                        medianprops={"color": "#2b2f36", "linewidth": 1.0})
        for patch in bp["boxes"]:
            patch.set_facecolor("#4e79a7")
            patch.set_alpha(0.7)
            patch.set_edgecolor("#8a919c")
    except Exception:  # noqa: BLE001
        for i, row in enumerate(values):
            ax.scatter(np.full(len(row), i + 1), row, s=2, alpha=0.4)
    step = max(1, len(samples) // 20)
    ax.set_xticks(range(1, len(values) + 1, step))
    ax.set_xticklabels([str(samples[i]) for i in range(0, len(samples), step)], rotation=90, fontsize=6)
    ax.set_ylabel("Expression")
    ax.set_title(title)
    return save_fig(fig, outdir, name, title)


# ---------------------------------------------------------------------------
# 步骤：qc
# ---------------------------------------------------------------------------

def run_qc(ctx, params, outdir):
    params = dict(params or {})
    X = _X(ctx)
    samples = [str(s) for s in ctx["samples"]]
    lib = X.sum(axis=1)
    detect = (X > 0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        median_expr = np.array([np.median(row[row > 0]) if (row > 0).any() else 0.0 for row in X])
        pct_top = np.array([row.max() / row.sum() if row.sum() > 0 else 0.0 for row in X])

    def zscores(values):
        std = np.nanstd(values)
        return np.zeros_like(values) if std == 0 else (values - np.nanmean(values)) / std

    z_flags = np.vstack([np.abs(zscores(lib)), np.abs(zscores(detect)), np.abs(zscores(median_expr))])
    outlier = (np.nanmax(z_flags, axis=0) > 3).astype(int)
    labels = _groups(ctx, params.get("group_col"))

    frame = pd.DataFrame({
        "sample": samples,
        "group": labels,
        "lib_size": lib.round(1),
        "detect_genes": detect,
        "median_expr": np.round(median_expr, 3),
        "pct_top_gene": np.round(pct_top, 4),
        "qc_outlier": np.where(outlier == 1, "yes", "no"),
    })
    figures = []
    colors = _color_map(labels)
    groups = list(dict.fromkeys([str(v) for v in labels]))
    fig, axes = new_fig(11.0, 4.6, ncols=2)
    axes = np.atleast_1d(axes)
    for ax, values, title in zip(axes, [lib, detect], ["Library size", "Detected genes"]):
        data = [values[np.array([str(v) for v in labels]) == g] for g in groups]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        medianprops={"color": "#2b2f36", "linewidth": 1.2})
        for patch, g in zip(bp["boxes"], groups):
            patch.set_facecolor(colors[g])
            patch.set_alpha(0.75)
            patch.set_edgecolor("#8a919c")
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups)
        ax.set_xlabel("Group")
        ax.set_title(title)
    figures.append(save_fig(fig, outdir, "qc_box.png", "Sample QC"))
    figures.append(_dist_figure(X, samples, outdir, "qc_density.png", "Expression distribution per sample"))
    # 离群样本图
    fig, ax = new_fig(8.8, 4.2)
    ax.scatter(lib, detect, c=["#c0392b" if o else "#4e79a7" for o in outlier], s=42, alpha=0.85,
               edgecolor="white", linewidth=0.5)
    for i, sample in enumerate(samples):
        if outlier[i]:
            ax.annotate(sample, (lib[i], detect[i]), fontsize=6.5, color="#c0392b")
    ax.set_xlabel("Library size")
    ax.set_ylabel("Detected genes")
    ax.set_title("QC scatter (red = outlier)")
    figures.append(save_fig(fig, outdir, "qc_scatter.png", "QC scatter"))

    tables = [table("qc", "样本质控指标", list(frame.columns), frame)]
    path = export_csv(frame, outdir, "qc_metrics", index=False)
    n_outlier = int(outlier.sum())
    bullets = [
        f"{len(samples)} 个样本：文库大小中位数 {np.median(lib):.0f}，检出基因中位数 {int(np.median(detect))}",
        f"离群样本 {n_outlier} 个" + ("：" + "、".join(s for s, o in zip(samples, outlier) if o) if n_outlier else ""),
    ]
    payload = result(
        title="样本质控",
        bullets=bullets,
        metrics={
            "n_samples": int(len(samples)),
            "library_median": float(np.median(lib)),
            "detect_median": float(np.median(detect)),
            "n_outliers": n_outlier,
            "value_kind": ctx.get("X_kind"),
        },
        tables=tables,
        figures=figures,
        notes=[f"当前矩阵值类型：{ctx.get('X_kind')}；离群判定：文库大小/检出基因/中位表达任一 |z|>3"],
        data_out=[data_file("qc_csv", path, "质控指标表")],
        sample_columns=[
            {"key": "lib_size", "kind": "numeric"},
            {"key": "detect_genes", "kind": "numeric"},
            {"key": "median_expr", "kind": "numeric"},
            {"key": "pct_top_gene", "kind": "numeric"},
            {"key": "qc_outlier", "kind": "categorical"},
        ],
        sample_values={
            "lib_size": lib.tolist(),
            "detect_genes": detect.tolist(),
            "median_expr": np.round(median_expr, 4).tolist(),
            "pct_top_gene": np.round(pct_top, 5).tolist(),
            "qc_outlier": np.where(outlier == 1, "yes", "no").tolist(),
        },
    )
    return payload


# ---------------------------------------------------------------------------
# 步骤：filter
# ---------------------------------------------------------------------------

def run_filter(ctx, params, outdir):
    params = dict(params or {})
    X = _X(ctx)
    genes = [str(g) for g in ctx["genes"]]
    samples = [str(s) for s in ctx["samples"]]
    info = ctx.get("sample_info")
    info = info.copy() if isinstance(info, pd.DataFrame) else pd.DataFrame(index=pd.Index(samples, name="sample"))
    counts = ctx.get("counts")
    counts = None if counts is None else np.asarray(counts, dtype=np.float32)

    min_expr = float(params.get("min_expr") or 0)
    frac = float(params.get("min_samples_frac") or 0.0)
    if frac > 1:
        frac = frac / 100.0
    n_required = max(1, int(np.ceil(frac * len(samples))))
    keep_genes = ((X >= min_expr).sum(axis=0) >= n_required) if len(samples) else np.ones(len(genes), bool)
    if not keep_genes.any():
        raise EngineError("prerequisite", f"过滤后没有基因剩余（min_expr={min_expr}，样本比例={frac}），请放宽阈值")
    removed_genes = int((~keep_genes).sum())

    remove_outliers = bool(params.get("remove_outliers"))
    keep_samples = np.ones(len(samples), bool)
    removed_samples: list[str] = []
    if remove_outliers:
        if "qc_outlier" in info.columns:
            keep_samples = (info["qc_outlier"].astype(str).str.lower() != "yes").values
        else:
            lib = X.sum(axis=1)
            std = np.nanstd(lib)
            z = np.zeros_like(lib) if std == 0 else (lib - np.nanmean(lib)) / std
            keep_samples = np.abs(z) <= 3
        removed_samples = [s for s, keep in zip(samples, keep_samples) if not keep]
        if keep_samples.sum() < 4:
            keep_samples = np.ones(len(samples), bool)
            removed_samples = []
    X_new = np.asarray(X[np.ix_(keep_samples, keep_genes)], dtype=np.float32)
    genes_new = [g for g, keep in zip(genes, keep_genes) if keep]
    samples_new = [s for s, keep in zip(samples, keep_samples) if keep]
    info_new = info.loc[[s for s, keep in zip(samples, keep_samples) if keep]]

    figures = []
    fig, ax = new_fig(8.8, 4.4)
    ax.hist(np.log1p(np.abs(X).mean(axis=0)), bins=60, color="#b9bfc8", alpha=0.8, label="before")
    ax.hist(np.log1p(np.abs(X_new).mean(axis=0)), bins=60, color="#4e79a7", alpha=0.7, label="after")
    ax.set_xlabel("log1p(mean expression)")
    ax.set_ylabel("Genes")
    ax.set_title("Gene filtering")
    ax.legend(fontsize=8)
    figures.append(save_fig(fig, outdir, "filter_genes.png", "Gene filtering"))

    tables = [table("filter_summary", "过滤结果", ["item", "before", "after", "removed"], [
        {"item": "基因数", "before": len(genes), "after": len(genes_new), "removed": removed_genes},
        {"item": "样本数", "before": len(samples), "after": len(samples_new), "removed": len(removed_samples)},
    ])]
    payload = result(
        title="过滤",
        bullets=[
            f"基因：{len(genes)} → {len(genes_new)}（去掉 {removed_genes} 个低表达基因，阈值 {min_expr} × {n_required} 样本）",
            f"样本：{len(samples)} → {len(samples_new)}" + (f"（剔除 {'、'.join(removed_samples)}）" if removed_samples else ""),
        ],
        metrics={"n_genes_before": len(genes), "n_genes_after": len(genes_new),
                 "n_samples_before": len(samples), "n_samples_after": len(samples_new),
                 "removed_genes": removed_genes, "removed_samples": len(removed_samples)},
        tables=tables,
        figures=figures,
        notes=[f"保留条件：表达量 ≥{min_expr} 的样本数 ≥{n_required}"],
    )
    replace = {"X": X_new, "genes": genes_new, "samples": samples_new, "sample_info": info_new}
    if counts is not None:
        replace["counts"] = np.asarray(counts[np.ix_(keep_samples, keep_genes)], dtype=np.float32)
    payload["replace"] = replace
    return payload


# ---------------------------------------------------------------------------
# 步骤：normalize
# ---------------------------------------------------------------------------

def _quantile_normalize(X: np.ndarray) -> np.ndarray:
    """分位数标准化（preprocessCore::normalize.quantiles 的等价实现）。"""
    arr = np.array(X, dtype=np.float64)
    if arr.shape[0] < 2:
        return arr
    sorted_idx = np.argsort(arr, axis=1)
    sorted_arr = np.take_along_axis(arr, sorted_idx, axis=1)
    reference = sorted_arr.mean(axis=0)
    out = np.empty_like(arr)
    for i in range(arr.shape[0]):
        out[i, sorted_idx[i]] = reference
    return out


def _tmm_factors(X: np.ndarray) -> np.ndarray:
    """TMM 简化：以表达量最接近上四分位的样本为参照，做 30% 修剪的 M 值加权平均。"""
    arr = np.asarray(X, dtype=np.float64)
    lib = arr.sum(axis=1)
    if lib.min() <= 0:
        return np.ones(arr.shape[0])
    ref_idx = int(np.argsort(lib)[int(0.75 * (len(lib) - 1))])
    ref = arr[ref_idx]
    factors = np.ones(arr.shape[0])
    for i in range(arr.shape[0]):
        obs = arr[i]
        mask = (obs > 0) & (ref > 0)
        if mask.sum() < 10:
            continue
        log_ratio = np.log2(obs[mask] / lib[i]) - np.log2(ref[mask] / lib[ref_idx])
        abs_expr = 0.5 * (np.log2(obs[mask] / lib[i]) + np.log2(ref[mask] / lib[ref_idx]))
        lo, hi = np.quantile(abs_expr, [0.3, 0.7])
        keep = (abs_expr >= lo) & (abs_expr <= hi)
        if keep.sum() < 5:
            continue
        factors[i] = float(np.power(2.0, np.mean(log_ratio[keep])))
    factors = factors / np.exp(np.mean(np.log(factors)))
    return np.clip(factors, 0.1, 10.0)


def run_normalize(ctx, params, outdir):
    params = dict(params or {})
    X = _X(ctx)
    X_orig = X.copy()
    method = str(params.get("method") or "auto").lower()
    kind = str(ctx.get("X_kind") or "unknown")
    warnings: list[str] = []
    notes: list[str] = []
    log_transform = str(params.get("log_transform") or "auto").lower()

    do_log = False
    if log_transform == "yes":
        do_log = True
    elif log_transform == "auto":
        finite = X[np.isfinite(X)]
        integer_like = finite.size > 0 and bool(np.allclose(finite, np.round(finite)))
        do_log = bool(integer_like or (finite.size and finite.max() > 50))
    if kind == "lognorm" and log_transform != "yes":
        do_log = False
        notes.append("矩阵看起来已经是 log 尺度（X_kind=lognorm），本次不再重复取对数")

    if method == "auto":
        method = "cpm" if (do_log and kind == "counts") else ("quantile" if do_log else "median")

    Xn = X.astype(np.float64)
    if method == "cpm":
        lib = Xn.sum(axis=1, keepdims=True)
        lib[lib <= 0] = 1.0
        Xn = np.log2(Xn / lib * 1e6 + 1.0)
        notes.append("log2(CPM+1) 标准化")
    elif method == "tmm":
        factors = _tmm_factors(Xn)
        eff = Xn.sum(axis=1) * factors
        eff[eff <= 0] = 1.0
        Xn = np.log2(Xn / eff[:, None] * 1e6 + 1.0)
        notes.append("TMM 归一化因子（简化实现）+ log2(CPM+1)")
    elif method == "quantile":
        if do_log:
            Xn = np.log2(np.clip(Xn, 0, None) + 1.0)
        Xn = _quantile_normalize(Xn)
        notes.append("分位数标准化" + ("（先 log2 转换）" if do_log else ""))
    elif method == "median":
        if do_log:
            Xn = np.log2(np.clip(Xn, 0, None) + 1.0)
        row_medians = np.array([np.median(row[row > 0]) if (row > 0).any() else 0.0 for row in Xn])
        target = float(np.median(row_medians))
        Xn = Xn - (row_medians - target)[:, None]
        notes.append("中位数中心化" + ("（先 log2 转换）" if do_log else ""))
    elif method in ("none", "log"):
        if do_log:
            Xn = np.log2(np.clip(Xn, 0, None) + 1.0)
        notes.append("仅 log2 转换" if do_log else "不做标准化")
    else:
        raise EngineError("bad_request", f"未知标准化方法：{method}")

    X_out = np.asarray(Xn, dtype=np.float32)
    figures = [
        _dist_figure(X_orig, list(ctx["samples"]), outdir, "normalize_before.png", "Before normalization"),
        _dist_figure(X_out, list(ctx["samples"]), outdir, "normalize_after.png", "After normalization"),
    ]
    bullets = [f"标准化方法：{method}" + ("（含 log2 转换）" if do_log else "")]
    before_med = float(np.median([np.median(row[row > 0]) if (row > 0).any() else 0.0 for row in X_orig]))
    after_med = float(np.median([np.median(row[row > 0]) if (row > 0).any() else 0.0 for row in X_out]))
    bullets.append(f"样本中位表达：{before_med:.3g} → {after_med:.3g}")
    payload = result(
        title="标准化",
        bullets=bullets,
        metrics={"method": method, "log_transformed": do_log,
                 "median_before": before_med, "median_after": after_med},
        figures=figures,
        warnings=warnings,
        notes=notes + ["标准化后矩阵用于差异分析（logFC 为 log2 尺度）"],
    )
    payload["replace"] = {"X": X_out, "value_kind": "lognorm"}
    return payload


# ---------------------------------------------------------------------------
# 步骤：batch（ComBat / removeBatchEffect）
# ---------------------------------------------------------------------------

def _design(groups) -> np.ndarray:
    uniq = list(dict.fromkeys([str(g) for g in groups]))
    if len(uniq) < 2:
        return np.ones((len(groups), 1))
    columns = [np.ones(len(groups))]
    reference = uniq[0]
    for g in uniq[1:]:
        columns.append(np.array([1.0 if str(v) == g else 0.0 for v in groups]))
    return np.column_stack(columns)


def _combat(X: np.ndarray, batch: np.ndarray, groups) -> tuple[np.ndarray, dict]:
    """ComBat（简化经验贝叶斯位置-尺度校正）。X: 样本×基因（log 尺度）。"""
    Y = np.asarray(X, dtype=np.float64).copy()
    batch = np.asarray(batch).astype(str)
    batches = list(dict.fromkeys(batch.tolist()))
    if len(batches) < 2:
        raise EngineError("prerequisite", "批次列只有 1 个取值，无法去批次")
    G = _design(groups)
    beta, *_ = np.linalg.lstsq(G, Y, rcond=None)
    resid = Y - G @ beta
    dof = max(Y.shape[0] - G.shape[1], 1)
    sigma2 = np.maximum((resid ** 2).sum(axis=0) / dof, 1e-8)
    z = resid / np.sqrt(sigma2)
    gamma = np.zeros((len(batches), Y.shape[1]))
    delta2 = np.ones((len(batches), Y.shape[1]))
    counts = []
    for i, b in enumerate(batches):
        mask = batch == b
        counts.append(int(mask.sum()))
        if mask.sum() > 1:
            gamma[i] = z[mask].mean(axis=0)
            delta2[i] = np.maximum(z[mask].var(axis=0, ddof=1), 1e-6)
    # 经验贝叶斯：gamma 向 0 收缩，delta2 向合并方差收缩
    tau2 = np.maximum(gamma.var(axis=0, ddof=0), 1e-6)
    pooled = (delta2 * np.array(counts)[:, None]).sum(axis=0) / max(sum(counts), 1)
    n_prior = 5.0
    gamma_star = np.zeros_like(gamma)
    delta2_star = np.zeros_like(delta2)
    for i in range(len(batches)):
        weight = tau2 / (tau2 + delta2[i] / max(counts[i], 1))
        gamma_star[i] = gamma[i] * weight
        delta2_star[i] = (counts[i] * delta2[i] + n_prior * pooled) / (counts[i] + n_prior)
    corrected_z = z.copy()
    for i, b in enumerate(batches):
        mask = batch == b
        corrected_z[mask] = (z[mask] - gamma_star[i]) / np.sqrt(np.maximum(delta2_star[i], 1e-6))
    corrected = corrected_z * np.sqrt(sigma2) + G @ beta
    metrics = {
        "n_batches": len(batches),
        "batch_sizes": {b: int(c) for b, c in zip(batches, counts)},
        "mean_abs_shift": float(np.mean(np.abs(gamma_star))),
        "mean_scale": float(np.mean(np.sqrt(delta2_star))),
    }
    return corrected, metrics


def _remove_batch_effect(X: np.ndarray, batch: np.ndarray, groups) -> tuple[np.ndarray, dict]:
    """limma::removeBatchEffect 的等价简化：去掉批次均值（保留分组效应）。"""
    Y = np.asarray(X, dtype=np.float64).copy()
    batch = np.asarray(batch).astype(str)
    G = _design(groups)
    beta, *_ = np.linalg.lstsq(G, Y, rcond=None)
    resid = Y - G @ beta
    out = resid.copy()
    for b in dict.fromkeys(batch.tolist()):
        mask = batch == b
        out[mask] -= resid[mask].mean(axis=0)
    out = out + G @ beta
    metrics = {"n_batches": len(set(batch.tolist())),
               "mean_abs_shift": float(np.mean(np.abs(resid.mean(axis=0))))}
    return out, metrics


def _batch_r2(scores: np.ndarray, labels) -> float:
    """批次对前几个 PC 的解释度（1 - 组内方差/总方差）。"""
    try:
        labels = np.asarray([str(v) for v in labels])
        total = np.var(scores, axis=0).sum()
        if total <= 0:
            return 0.0
        within = 0.0
        for value in dict.fromkeys(labels.tolist()):
            subset = scores[labels == value]
            if len(subset) > 1:
                within += np.var(subset, axis=0).sum() * len(subset)
        return float(max(0.0, 1.0 - within / (total * len(labels))))
    except Exception:  # noqa: BLE001
        return 0.0


def run_batch(ctx, params, outdir):
    params = dict(params or {})
    X = _X(ctx)
    batch_col = _batch_column(ctx, params.get("batch_key"))
    if not batch_col:
        raise EngineError(
            "prerequisite",
            "样本信息里找不到批次列（batch/batch_id/cohort/plate/run 等）；请在「样本」里指定或导入带批次列的表",
        )
    info = ctx["sample_info"]
    batch = info[batch_col].astype(str).values
    if len(set(batch.tolist())) < 2:
        raise EngineError("prerequisite", f"批次列 {batch_col} 只有 1 个取值，无需去批次")
    labels = _groups(ctx, params.get("group_col"))
    method = str(params.get("method") or "combat").lower()
    if method not in ("combat", "removebatcheffect", "none"):
        raise EngineError("bad_request", f"未知去批次方法：{method}")
    if method == "none":
        corrected = np.asarray(X, dtype=np.float32)
        metrics = {"method": "none"}
        notes = ["未做校正（method=none）"]
    elif method == "removebatcheffect":
        corrected_raw, metrics = _remove_batch_effect(X, batch, labels)
        corrected = np.asarray(corrected_raw, dtype=np.float32)
        metrics["method"] = "removeBatchEffect（组均值中心化）"
        notes = ["removeBatchEffect 等价简化：去掉批次均值，保留分组效应"]
    else:
        corrected_raw, metrics = _combat(X, batch, labels)
        corrected = np.asarray(corrected_raw, dtype=np.float32)
        metrics["method"] = "ComBat（简化经验贝叶斯位置-尺度）"
        notes = ["ComBat 简化实现：位置/尺度经验贝叶斯收缩（对标 sva::ComBat，未做参考批次选项）"]

    scores_before, var_before = _pca_scores(X, 3)
    scores_after, var_after = _pca_scores(corrected, 3)
    r2_before = _batch_r2(scores_before, batch)
    r2_after = _batch_r2(scores_after, batch)
    metrics.update({"batch_r2_before": round(r2_before, 4), "batch_r2_after": round(r2_after, 4)})

    figures = []
    fig, axes = new_fig(9.6, 4.8, ncols=2)
    axes = np.atleast_1d(axes)
    colors = _color_map(batch)
    for ax, scores, title in ((axes[0], scores_before, "Before correction"), (axes[1], scores_after, "After correction")):
        for value in dict.fromkeys(batch.tolist()):
            mask = batch == value
            ax.scatter(scores[mask, 0], scores[mask, 1], s=36, alpha=0.85, color=colors[value],
                       edgecolor="white", linewidth=0.5, label=str(value))
        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
    axes[0].legend(fontsize=7)
    figures.append(save_fig(fig, outdir, "batch_pca.png", "Batch correction"))
    # 分组着色（校正后）
    fig, ax = new_fig(7.0, 5.2)
    group_colors = _color_map(labels)
    for value in dict.fromkeys([str(v) for v in labels]):
        mask = np.array([str(v) == value for v in labels])
        ax.scatter(scores_after[mask, 0], scores_after[mask, 1], s=38, alpha=0.85,
                   color=group_colors[value], edgecolor="white", linewidth=0.5, label=str(value))
    ax.set_title("After correction (color = group)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(fontsize=7)
    figures.append(save_fig(fig, outdir, "batch_group.png", "After correction by group"))
    bullets = [
        f"{method}：{metrics.get('n_batches')} 个批次（{batch_col}），"
        f"批次解释度 R² {r2_before:.3f} → {r2_after:.3f}",
        "分组与批次混杂时无法区分生物学差异" if r2_before > 0.95 else f"批次与分组可区分（R²={r2_before:.3f}）",
    ]
    payload = result(
        title="去批次",
        bullets=bullets,
        metrics=metrics,
        figures=figures,
        notes=notes,
        store={"batch": {"corrected": corrected, "batch_key": batch_col,
                         "r2_before": r2_before, "r2_after": r2_after}},
    )
    payload["replace"] = {"X": corrected}
    return payload


# ---------------------------------------------------------------------------
# 步骤：pca
# ---------------------------------------------------------------------------

def run_pca(ctx, params, outdir):
    params = dict(params or {})
    X = _X(ctx)
    n_pcs = int(params.get("n_pcs") or 20)
    scores, var = _pca_scores(X, n_pcs)
    info = ctx.get("sample_info")
    color_by = str(params.get("color_by") or "group")
    labels = _groups(ctx, None if color_by == "group" else color_by)
    if color_by != "group" and info is not None and color_by in info.columns:
        labels = info[color_by].astype(str).values

    figures = []
    # 碎石图
    n_show = min(len(var), 20)
    fig, ax = new_fig(7.6, 4.4)
    ax.bar(np.arange(1, n_show + 1), var[:n_show] * 100, color="#4e79a7", alpha=0.88)
    ax.plot(np.arange(1, n_show + 1), np.cumsum(var[:n_show]) * 100, color="#c0392b", marker="o", markersize=3,
            linewidth=1.2, label="cumulative")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title("Scree plot")
    ax.legend(fontsize=8)
    figures.append(save_fig(fig, outdir, "pca_scree.png", "Scree plot"))
    # 样本散点
    figures.append(_pca_figure(scores[:, :2], var[:2], labels, outdir, "pca_samples.png",
                               f"PCA of samples (color = {color_by})"))
    # 相关性热图
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(np.asarray(X, dtype=float))
    corr = np.nan_to_num(corr, nan=0.0)
    fig, ax = new_fig(8.4, 7.0)
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(ctx["samples"])))
    ax.set_yticks(range(len(ctx["samples"])))
    ax.set_xticklabels(ctx["samples"], rotation=90, fontsize=6)
    ax.set_yticklabels(ctx["samples"], fontsize=6)
    ax.grid(False)
    ax.set_title("Sample-sample correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    figures.append(save_fig(fig, outdir, "sample_corr.png", "Sample correlation"))
    # 层次聚类
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist

        link = linkage(pdist(np.asarray(X, dtype=float), metric="correlation"), method="average")
        fig, ax = new_fig(9.6, 5.0)
        dendrogram(link, labels=list(ctx["samples"]), ax=ax, leaf_font_size=7, color_threshold=0)
        ax.set_title("Sample clustering (correlation distance)")
        ax.set_ylabel("Distance")
        ax.grid(False)
        figures.append(save_fig(fig, outdir, "sample_dendro.png", "Sample dendrogram"))
    except Exception:  # noqa: BLE001
        pass

    explained = float(var[0]) if len(var) else 0.0
    bullets = [f"PCA：{scores.shape[1]} 个主成分，PC1 解释 {explained * 100:.1f}%，PC1+PC2 共 {float(var[:2].sum()) * 100:.1f}%"]
    payload = result(
        title="降维聚类",
        bullets=bullets,
        metrics={"n_pcs": int(scores.shape[1]), "pc1_variance": round(explained, 4),
                 "pc2_variance": round(float(var[1]), 4) if len(var) > 1 else None},
        figures=figures,
        notes=["PCA 用于样本层面质控/批次检查；样本聚类使用 1 - Pearson 相关距离 + average linkage"],
        sample_columns=[{"key": "pc1", "kind": "numeric"}, {"key": "pc2", "kind": "numeric"}],
        sample_values={"pc1": scores[:, 0].tolist(), "pc2": scores[:, 1].tolist()},
    )
    return payload


# ---------------------------------------------------------------------------
# 基础分析
# ---------------------------------------------------------------------------

def run_gene_box(ctx, params, outdir):
    params = dict(params or {})
    genes_in = params.get("genes") or params.get("gene")
    if isinstance(genes_in, str):
        genes_in = [g.strip() for g in genes_in.replace(",", " ").split() if g.strip()]
    genes = [str(g) for g in (genes_in or [])]
    index = {g: i for i, g in enumerate(ctx["genes"])}
    genes = [g for g in genes if g in index]
    if not genes:
        raise EngineError("prerequisite", "请提供存在于矩阵中的基因（genes）")
    X = _X(ctx)
    labels = _groups(ctx, params.get("group_col"))
    group_names = list(dict.fromkeys([str(v) for v in labels]))
    colors = _color_map(labels)
    figures = []
    for gene in genes[:6]:
        values = X[:, index[gene]]
        data = [values[np.array([str(v) for v in labels]) == g] for g in group_names]
        fig, ax = new_fig(5.6, 4.4)
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                        medianprops={"color": "#2b2f36", "linewidth": 1.1})
        for patch, g in zip(bp["boxes"], group_names):
            patch.set_facecolor(colors[g])
            patch.set_alpha(0.78)
            patch.set_edgecolor("#8a919c")
        ax.set_xticks(range(1, len(group_names) + 1))
        ax.set_xticklabels(group_names)
        ax.set_ylabel("Expression")
        ax.set_title(f"{gene} by group")
        figures.append(save_fig(fig, outdir, f"gene_box_{gene}.png", f"{gene} expression"))
    return result(title="基因表达箱线图", bullets=[f"{len(genes)} 个基因按 {len(group_names)} 组画箱线图"],
                  metrics={"n_genes": len(genes)}, figures=figures,
                  tables=[table("gene_means", "各组均值", [{"key": "gene", "label": "基因", "type": "string"}] +
                                [{"key": g, "label": g, "type": "number"} for g in group_names],
                                [{"gene": gene, **{g: float(np.mean(X[np.array([str(v) for v in labels]) == g, index[gene]]))
                                                   for g in group_names}} for gene in genes])])


def run_gene_heatmap(ctx, params, outdir):
    params = dict(params or {})
    genes_in = params.get("genes") or params.get("gene")
    if isinstance(genes_in, str):
        genes_in = [g.strip() for g in genes_in.replace(",", " ").split() if g.strip()]
    index = {g: i for i, g in enumerate(ctx["genes"])}
    genes = [str(g) for g in (genes_in or []) if str(g) in index]
    if not genes:
        raise EngineError("prerequisite", "请提供存在于矩阵中的基因（genes）")
    X = _X(ctx)
    mat = X[:, [index[g] for g in genes]].T
    z = zscore_rows(mat)
    labels = _groups(ctx, params.get("group_col"))
    group_names = list(dict.fromkeys([str(v) for v in labels]))
    order = np.argsort([group_names.index(str(v)) for v in labels])
    fig, ax = new_fig(8.8, max(3.8, 0.24 * len(genes) + 1.6))
    im = ax.imshow(z[:, order], aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xticks(range(len(ctx["samples"])))
    ax.set_xticklabels([ctx["samples"][i] for i in order], rotation=90, fontsize=6)
    ax.grid(False)
    ax.set_title("Gene expression heatmap (row z-score)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    figure = save_fig(fig, outdir, "gene_heatmap.png", "Gene heatmap")
    return result(title="基因热图", bullets=[f"{len(genes)} 个基因 × {len(ctx['samples'])} 样本"],
                  metrics={"n_genes": len(genes)}, figures=[figure])


def run_sample_correlation(ctx, params, outdir):
    X = _X(ctx)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.nan_to_num(np.corrcoef(X))
    fig, ax = new_fig(8.4, 7.0)
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(ctx["samples"])))
    ax.set_yticks(range(len(ctx["samples"])))
    ax.set_xticklabels(ctx["samples"], rotation=90, fontsize=6)
    ax.set_yticklabels(ctx["samples"], fontsize=6)
    ax.grid(False)
    ax.set_title("Sample-sample correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    figure = save_fig(fig, outdir, "sample_corr.png", "Sample correlation")
    off = corr[~np.eye(corr.shape[0], dtype=bool)]
    return result(title="样本相关性",
                  bullets=[f"样本间 Pearson 相关：中位数 {np.median(off):.3f}，最低 {off.min():.3f}"],
                  metrics={"median": float(np.median(off)), "min": float(off.min()), "max": float(off.max())},
                  figures=[figure])


def run_sample_dendro(ctx, params, outdir):
    X = _X(ctx)
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist

    link = linkage(pdist(X, metric="correlation"), method="average")
    fig, ax = new_fig(9.6, 5.0)
    dendrogram(link, labels=list(ctx["samples"]), ax=ax, leaf_font_size=7, color_threshold=0)
    ax.set_title("Sample clustering (correlation distance)")
    ax.set_ylabel("Distance")
    ax.grid(False)
    figure = save_fig(fig, outdir, "sample_dendro.png", "Sample dendrogram")
    return result(title="样本层次聚类", bullets=["按 1 - Pearson 相关距离做 average linkage 聚类"],
                  metrics={"n_samples": len(ctx["samples"])}, figures=[figure])


def run_expr_distribution(ctx, params, outdir):
    X = _X(ctx)
    figure = _dist_figure(X, list(ctx["samples"]), outdir, "expr_dist.png", "Expression distribution")
    medians = [float(np.median(row[row > 0])) if (row > 0).any() else 0.0 for row in X]
    return result(title="表达分布", bullets=[f"样本中位表达范围 {min(medians):.3g} – {max(medians):.3g}"],
                  metrics={"median_min": min(medians), "median_max": max(medians)},
                  figures=[figure])


def run_outlier_detect(ctx, params, outdir):
    X = _X(ctx)
    samples = list(ctx["samples"])
    metrics = {
        "lib_size": X.sum(axis=1),
        "detect_genes": (X > 0).sum(axis=1),
        "median_expr": np.array([np.median(row[row > 0]) if (row > 0).any() else 0.0 for row in X]),
    }
    flags = np.zeros(len(samples), dtype=bool)
    details = {}
    for name, values in metrics.items():
        std = np.nanstd(values)
        z = np.zeros_like(values) if std == 0 else (values - np.nanmean(values)) / std
        hot = np.abs(z) > 3
        flags |= hot
        details[name] = np.where(hot, "outlier", "ok").tolist()
    frame = pd.DataFrame({"sample": samples, "lib_size": metrics["lib_size"].round(1),
                          "detect_genes": metrics["detect_genes"], "median_expr": metrics["median_expr"].round(3),
                          "outlier": np.where(flags, "yes", "no")})
    fig, ax = new_fig(8.8, 4.2)
    ax.scatter(metrics["lib_size"], metrics["detect_genes"],
               c=["#c0392b" if f else "#4e79a7" for f in flags], s=42, alpha=0.85,
               edgecolor="white", linewidth=0.5)
    for i, sample in enumerate(samples):
        if flags[i]:
            ax.annotate(sample, (metrics["lib_size"][i], metrics["detect_genes"][i]), fontsize=6.5, color="#c0392b")
    ax.set_xlabel("Library size")
    ax.set_ylabel("Detected genes")
    ax.set_title("Outlier detection")
    figure = save_fig(fig, outdir, "outliers.png", "Outlier detection")
    return result(title="离群样本检测",
                  bullets=[f"{len(samples)} 个样本中检出 {int(flags.sum())} 个离群样本" +
                           ("：" + "、".join(s for s, f in zip(samples, flags) if f) if flags.any() else "")],
                  metrics={"n_outliers": int(flags.sum())},
                  tables=[table("outliers", "离群检测明细", list(frame.columns), frame)],
                  figures=[figure],
                  sample_columns=[{"key": "outlier", "kind": "categorical"}],
                  sample_values={"outlier": np.where(flags, "yes", "no").tolist()})


def run_batch_diagnose(ctx, params, outdir):
    params = dict(params or {})
    X = _X(ctx)
    batch_col = _batch_column(ctx, params.get("batch_key"))
    if not batch_col:
        raise EngineError("prerequisite", "样本信息里没有批次列")
    info = ctx["sample_info"]
    batch = info[batch_col].astype(str).values
    scores, var = _pca_scores(X, 5)
    r2 = _batch_r2(scores, batch)
    if ctx.get("group_key"):
        crosstab = pd.crosstab(pd.Series(batch, name="batch"), pd.Series(_groups(ctx), name="group"))
        confounded = bool((crosstab.values > 0).sum(axis=1).max() == 1 and crosstab.shape[0] > 1)
    else:
        crosstab = pd.DataFrame()
        confounded = False
    fig, axes = new_fig(9.6, 4.8, ncols=2)
    axes = np.atleast_1d(axes)
    colors = _color_map(batch)
    for value in dict.fromkeys(batch.tolist()):
        mask = batch == value
        axes[0].scatter(scores[mask, 0], scores[mask, 1], s=36, alpha=0.85, color=colors[value],
                        edgecolor="white", linewidth=0.5, label=str(value))
    axes[0].set_title(f"PCA by batch (R²={r2:.3f})")
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].legend(fontsize=7)
    if not crosstab.empty:
        im = axes[1].imshow(crosstab.values, cmap="Blues", aspect="auto")
        axes[1].set_xticks(range(crosstab.shape[1]))
        axes[1].set_xticklabels(crosstab.columns, rotation=45, fontsize=7)
        axes[1].set_yticks(range(crosstab.shape[0]))
        axes[1].set_yticklabels(crosstab.index, fontsize=7)
        axes[1].set_title("Batch × group")
        axes[1].grid(False)
        fig.colorbar(im, ax=axes[1], fraction=0.03, pad=0.02)
    figure = save_fig(fig, outdir, "batch_diagnose.png", "Batch diagnosis")
    return result(title="批次效应诊断",
                  bullets=[f"批次 {batch_col} 解释前 5 个主成分的 R²={r2:.3f}" +
                           ("；批次与分组完全混杂！" if confounded else "")],
                  metrics={"batch_r2": round(r2, 4), "n_batches": int(len(set(batch.tolist()))),
                           "confounded": confounded},
                  tables=[table("batch_group", "批次 × 分组", [{"key": "batch", "label": "批次", "type": "string"}] +
                                [{"key": str(c), "label": str(c), "type": "number"} for c in crosstab.columns],
                                crosstab.reset_index().rename(columns={"index": "batch"}))] if not crosstab.empty else [],
                  figures=[figure])


STEPS = [
    {
        "id": "qc",
        "title": "样本质控",
        "requires": [],
        "params": [
            {"key": "group_col", "label": "分组列（用于分组展示）", "type": "column", "default": None},
        ],
        "run": run_qc,
    },
    {
        "id": "filter",
        "title": "过滤",
        "requires": ["qc"],
        "params": [
            {"key": "min_expr", "label": "最小表达量", "type": "number", "default": 1, "min": 0, "max": 100, "step": 0.5},
            {"key": "min_samples_frac", "label": "最低样本比例", "type": "number", "default": 0.2, "min": 0, "max": 1, "step": 0.05},
            {"key": "remove_outliers", "label": "剔除质控离群样本", "type": "boolean", "default": False},
        ],
        "run": run_filter,
    },
    {
        "id": "normalize",
        "title": "标准化",
        "requires": ["filter"],
        "params": [
            {"key": "method", "label": "方法", "type": "select", "default": "auto",
             "options": [{"value": "auto", "label": "自动"}, {"value": "cpm", "label": "log2(CPM+1)"},
                         {"value": "quantile", "label": "分位数标准化"}, {"value": "tmm", "label": "TMM"},
                         {"value": "median", "label": "中位数中心化"}, {"value": "none", "label": "不标准化"}]},
            {"key": "log_transform", "label": "log2 转换", "type": "select", "default": "auto",
             "options": [{"value": "auto", "label": "自动"}, {"value": "yes", "label": "强制转换"}, {"value": "no", "label": "不转换"}]},
        ],
        "run": run_normalize,
    },
    {
        "id": "batch",
        "title": "去批次",
        "requires": ["normalize"],
        "params": [
            {"key": "batch_key", "label": "批次列", "type": "column", "default": None},
            {"key": "group_col", "label": "分组列（保留其效应）", "type": "column", "default": None},
            {"key": "method", "label": "方法", "type": "select", "default": "combat",
             "options": [{"value": "combat", "label": "ComBat"}, {"value": "removeBatchEffect", "label": "removeBatchEffect"},
                         {"value": "none", "label": "不校正"}]},
        ],
        "run": run_batch,
    },
    {
        "id": "pca",
        "title": "降维聚类",
        "requires": ["normalize"],
        "params": [
            {"key": "n_pcs", "label": "主成分数", "type": "number", "default": 20, "min": 2, "max": 100, "step": 1},
            {"key": "color_by", "label": "着色列", "type": "column", "default": "group"},
        ],
        "run": run_pca,
    },
]

ANALYSES = [
    {"id": "sample_correlation", "title": "样本相关性热图", "group": "基本分析",
     "method": "Pearson 相关 + 热图", "needs": [], "params": [], "run": run_sample_correlation},
    {"id": "sample_dendro", "title": "样本层次聚类", "group": "基本分析",
     "method": "1 - Pearson 距离 + average linkage", "needs": [], "params": [], "run": run_sample_dendro},
    {"id": "gene_boxplot", "title": "基因表达箱线图", "group": "基本分析",
     "method": "按分组画箱线图", "needs": [], "params": [
         {"key": "genes", "label": "基因列表", "type": "genes", "default": []},
         {"key": "group_col", "label": "分组列", "type": "column", "default": None},
     ], "run": run_gene_box},
    {"id": "gene_heatmap", "title": "基因表达热图", "group": "基本分析",
     "method": "行 z-score 热图", "needs": [], "params": [
         {"key": "genes", "label": "基因列表", "type": "genes", "default": []},
         {"key": "group_col", "label": "分组列", "type": "column", "default": None},
     ], "run": run_gene_heatmap},
    {"id": "expr_distribution", "title": "表达分布", "group": "基本分析",
     "method": "逐样本箱线图", "needs": [], "params": [], "run": run_expr_distribution},
    {"id": "outlier_detect", "title": "离群样本检测", "group": "基本分析",
     "method": "文库大小/检出基因/中位表达 z>3", "needs": [], "params": [], "run": run_outlier_detect},
    {"id": "batch_diagnose", "title": "批次效应诊断", "group": "基本分析",
     "method": "PCA R² + 批次×分组列联表", "needs": [], "params": [], "run": run_batch_diagnose},
    {"id": "pca_samples", "title": "PCA 样本散点", "group": "基本分析",
     "method": "SVD 主成分散点", "needs": [], "params": [
         {"key": "color_by", "label": "着色列", "type": "column", "default": "group"},
     ], "run": run_pca},
]
