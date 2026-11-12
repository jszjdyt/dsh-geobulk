"""
dsh-geobulk — 通用出图端点（/geobulk/plot）

UI 的"按需重绘"入口：不依赖各模块内部实现，直接基于 ctx（当前矩阵 + 缓存结果）画常见图。
步骤/分析自己返回的图仍由各模块生成；本文件只负责 plot 命令的 20 种 kind。

约定：图存 <outdir>/<kind>-<hash>.png，白底/英文文字/dpi≈140；缺前置数据抛 prerequisite。
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import pandas as pd

from bg_common import (
    EngineError,
    GRADIENT,
    bh_fdr,
    color_for,
    ensure_dir,
    new_fig,
    palette,
    prerequisite,
    save_fig,
    zscore_rows,
)

PLOT_KINDS = [
    "qc_box", "qc_density", "pca_scatter", "sample_corr", "sample_dendro",
    "volcano", "ma_plot", "de_heatmap", "gene_box",
    "enrich_bar", "enrich_dot", "immune_stack", "immune_corr",
    "km_curve", "forest", "risk_plot", "roc_time", "batch_pca", "circos", "splicing_heat",
]


def _hash(*parts) -> str:
    raw = json.dumps([str(p) for p in parts], ensure_ascii=False, sort_keys=True)
    return hashlib.md5(raw.encode("utf-8")).hexdigest()[:10]


def _de_frame(ctx):
    de = ctx.get("de")
    if not de:
        raise prerequisite("还没有差异表达结果，请先运行「差异表达」")
    frame = de.get("table")
    if isinstance(frame, pd.DataFrame) and not frame.empty:
        return frame
    raise prerequisite("差异结果为空")


def _store(ctx, key):
    value = (ctx.get("results") or {}).get(key)
    if not value:
        return None
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, dict):
        frame = value.get("table")
        if isinstance(frame, pd.DataFrame):
            return frame
        if isinstance(frame, list) and frame:
            return pd.DataFrame(frame)
        for key2 in ("table", "de_table", "cells", "scores"):
            frame = value.get(key2)
            if isinstance(frame, pd.DataFrame):
                return frame
            if isinstance(frame, list) and frame and isinstance(frame[0], dict):
                return pd.DataFrame(frame)
    return None


def _groups(ctx):
    labels = ctx.get("sample_info")
    key = ctx.get("group_key")
    if labels is not None and key and key in labels.columns:
        return labels[key].astype(str).values
    return np.array(["all"] * len(ctx["samples"]))


def _color_map(labels) -> dict:
    uniq = list(dict.fromkeys([str(x) for x in labels]))
    return {name: palette(len(uniq))[i] for i, name in enumerate(uniq)}


# ---------------------------------------------------------------------------
# 各类图
# ---------------------------------------------------------------------------

def plot_qc_box(ctx, params, outdir):
    X = np.asarray(ctx["X"], dtype=float)
    if X.size == 0:
        raise prerequisite("表达矩阵为空")
    detect = (X > 0).sum(axis=1)
    lib = X.sum(axis=1)
    labels = _groups(ctx)
    colors = _color_map(labels)
    fig, axes = new_fig(11.0, 4.6, ncols=2)
    axes = np.atleast_1d(axes)
    for ax, values, title in zip(axes, [lib, detect], ["Library size", "Detected genes"]):
        order = np.argsort([str(v) for v in labels])
        data = [values[order][np.array([str(labels[i]) for i in order]) == g] for g in dict.fromkeys([str(v) for v in labels])]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        medianprops={"color": "#2b2f36", "linewidth": 1.2})
        for patch, name in zip(bp["boxes"], dict.fromkeys([str(v) for v in labels])):
            patch.set_facecolor(colors.get(name, "#4e79a7"))
            patch.set_alpha(0.75)
            patch.set_edgecolor("#8a919c")
        ax.set_xticks(range(1, len(data) + 1))
        ax.set_xticklabels(list(dict.fromkeys([str(v) for v in labels])), rotation=0)
        ax.set_title(title)
        ax.set_xlabel("Group")
    return save_fig(fig, outdir, f"qc_box-{_hash('qc_box', list(labels)[:5])}.png", "Sample QC")


def plot_qc_density(ctx, params, outdir):
    X = np.asarray(ctx["X"], dtype=float)
    fig, ax = new_fig(8.4, 4.8)
    for i in range(min(X.shape[0], 60)):
        values = X[i][np.isfinite(X[i])]
        if values.size == 0:
            continue
        try:
            import scipy.stats as st

            kde = st.gaussian_kde(values)
            grid = np.linspace(values.min(), values.max(), 200)
            ax.plot(grid, kde(grid), color="#4e79a7", alpha=0.25, linewidth=0.9)
        except Exception:  # noqa: BLE001
            break
    ax.set_title("Expression distribution per sample")
    ax.set_xlabel("Expression")
    ax.set_ylabel("Density")
    return save_fig(fig, outdir, f"qc_density-{_hash('qc_density', X.shape)}.png", "Sample density")


def _pca_components(X, n_pcs=10):
    arr = np.asarray(X, dtype=float)
    arr = arr - arr.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(arr, full_matrices=False)
    scores = U * S
    var = (S ** 2) / max((S ** 2).sum(), 1e-12)
    return scores[:, :n_pcs], var[:n_pcs]


def plot_pca_scatter(ctx, params, outdir):
    X = np.asarray(ctx["X"], dtype=float)
    if X.shape[0] < 2:
        raise prerequisite("样本数不足 2，无法做 PCA")
    scores, var = _pca_components(X, 5)
    color_by = str(params.get("color_by") or "group")
    labels = _groups(ctx)
    info = ctx.get("sample_info")
    if color_by != "group" and info is not None and color_by in info.columns:
        labels = info[color_by].astype(str).values
    colors = _color_map(labels)
    fig, ax = new_fig(7.4, 5.4)
    for name in dict.fromkeys([str(x) for x in labels]):
        mask = np.array([str(v) == name for v in labels])
        ax.scatter(scores[mask, 0], scores[mask, 1], s=46, alpha=0.85,
                   color=colors[name], edgecolor="white", linewidth=0.6, label=name)
    for i, sample in enumerate(ctx["samples"]):
        if params.get("labels") and i < len(ctx["samples"]):
            ax.annotate(str(sample), (scores[i, 0], scores[i, 1]), fontsize=6.5, alpha=0.7)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)" if len(var) > 1 else "PC2")
    ax.set_title("PCA of samples")
    ax.legend(loc="best", fontsize=8)
    return save_fig(fig, outdir, f"pca-{_hash('pca', X.shape, color_by)}.png", "PCA")


def plot_sample_corr(ctx, params, outdir):
    X = np.asarray(ctx["X"], dtype=float)
    if X.shape[0] < 2:
        raise prerequisite("样本数不足 2")
    corr = np.corrcoef(X)
    fig, ax = new_fig(8.6, 7.2)
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(ctx["samples"])))
    ax.set_yticks(range(len(ctx["samples"])))
    ax.set_xticklabels(ctx["samples"], rotation=90, fontsize=6)
    ax.set_yticklabels(ctx["samples"], fontsize=6)
    ax.grid(False)
    ax.set_title("Sample-sample correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    return save_fig(fig, outdir, f"sample_corr-{_hash('corr', X.shape)}.png", "Sample correlation")


def plot_sample_dendro(ctx, params, outdir):
    X = np.asarray(ctx["X"], dtype=float)
    if X.shape[0] < 3:
        raise prerequisite("样本数不足 3，无法聚类")
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist

        dist = pdist(X, metric="correlation")
        link = linkage(dist, method="average")
    except Exception as error:  # noqa: BLE001
        raise EngineError("internal", f"层次聚类失败：{error}")
    fig, ax = new_fig(9.6, 5.0)
    dendrogram(link, labels=list(ctx["samples"]), ax=ax, leaf_font_size=7, color_threshold=0)
    ax.set_title("Sample clustering (correlation distance)")
    ax.set_ylabel("Distance")
    ax.grid(False)
    return save_fig(fig, outdir, f"dendro-{_hash('dendro', X.shape)}.png", "Sample dendrogram")


def plot_volcano(ctx, params, outdir):
    frame = _de_frame(ctx)
    lfc_cut = float(params.get("lfc", 1.0) or 1.0)
    padj_cut = float(params.get("padj", 0.05) or 0.05)
    top_n = int(params.get("top_n", 12) or 12)
    frame = frame.copy()
    frame["neglog"] = -np.log10(np.clip(frame["adj.P.Val"].astype(float).values, 1e-300, 1))
    sig = (frame["adj.P.Val"] <= padj_cut) & (frame["logFC"].abs() >= lfc_cut)
    up = sig & (frame["logFC"] > 0)
    down = sig & (frame["logFC"] < 0)
    fig, ax = new_fig(7.8, 6.0)
    ax.scatter(frame.loc[~sig, "logFC"], frame.loc[~sig, "neglog"], s=10, color="#b9bfc8", alpha=0.55, label="ns")
    ax.scatter(frame.loc[up, "logFC"], frame.loc[up, "neglog"], s=14, color="#c0392b", alpha=0.8, label=f"up ({int(up.sum())})")
    ax.scatter(frame.loc[down, "logFC"], frame.loc[down, "neglog"], s=14, color="#2e6da4", alpha=0.8, label=f"down ({int(down.sum())})")
    ax.axhline(-np.log10(padj_cut), color="#8a919c", linestyle="--", linewidth=0.8)
    ax.axvline(lfc_cut, color="#8a919c", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_cut, color="#8a919c", linestyle="--", linewidth=0.8)
    hits = frame.loc[sig].nsmallest(top_n, "adj.P.Val")
    for _, row in hits.iterrows():
        ax.annotate(str(row["gene"]), (row["logFC"], row["neglog"]), fontsize=6.5, alpha=0.85)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(adjusted P)")
    ax.set_title("Volcano plot")
    ax.legend(loc="upper right", fontsize=8)
    return save_fig(fig, outdir, f"volcano-{_hash('volcano', lfc_cut, padj_cut)}.png", "Volcano")


def plot_ma(ctx, params, outdir):
    frame = _de_frame(ctx)
    fig, ax = new_fig(7.4, 5.4)
    ax.scatter(frame["AveExpr"], frame["logFC"], s=9, color="#7f8c9b", alpha=0.45)
    sig = frame["adj.P.Val"] <= float(params.get("padj", 0.05) or 0.05)
    ax.scatter(frame.loc[sig, "AveExpr"], frame.loc[sig, "logFC"], s=12, color="#c0392b", alpha=0.7, label="significant")
    ax.axhline(0, color="#2b2f36", linewidth=0.8)
    ax.set_xlabel("Average expression")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot")
    ax.legend(fontsize=8)
    return save_fig(fig, outdir, f"ma-{_hash('ma', frame.shape)}.png", "MA plot")


def plot_de_heatmap(ctx, params, outdir):
    frame = _de_frame(ctx)
    top_n = int(params.get("top_n", 50) or 50)
    genes = params.get("genes")
    if genes:
        wanted = [g for g in _as_list(genes) if g in ctx["genes"]]
    else:
        wanted = frame.sort_values("adj.P.Val").head(top_n)["gene"].astype(str).tolist()
    wanted = [g for g in wanted if g in ctx["genes"]]
    if not wanted:
        raise prerequisite("没有可画的差异基因")
    idx = [ctx["genes"].index(g) for g in wanted]
    mat = np.asarray(ctx["X"], dtype=float)[:, idx].T
    z = zscore_rows(mat)
    labels = _groups(ctx)
    group_names = list(dict.fromkeys([str(l) for l in labels]))
    order = np.argsort([group_names.index(str(l)) for l in labels])
    z = z[:, order]
    fig, ax = new_fig(8.8, max(4.2, len(wanted) * 0.16))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(wanted)))
    ax.set_yticklabels(wanted, fontsize=6)
    ax.set_xticks(range(len(ctx["samples"])))
    xticklabels = [ctx["samples"][i] for i in order]
    ax.set_xticklabels(xticklabels, rotation=90, fontsize=6)
    ax.grid(False)
    ax.set_title("Top differential genes (row z-score)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    return save_fig(fig, outdir, f"de_heat-{_hash('de_heat', wanted[:5], top_n)}.png", "DE heatmap")


def plot_gene_box(ctx, params, outdir):
    genes = _as_list(params.get("genes") or params.get("gene"))
    genes = [g for g in genes if g in ctx["genes"]]
    if not genes:
        raise prerequisite("请指定至少一个存在于矩阵中的基因（genes）")
    genes = genes[:8]
    X = np.asarray(ctx["X"], dtype=float)
    labels = _groups(ctx)
    group_names = list(dict.fromkeys([str(l) for l in labels]))
    colors = _color_map(labels)
    fig, axes = new_fig(max(8.0, 2.2 * len(genes)), 4.2, ncols=len(genes))
    axes = np.atleast_1d(axes)
    for ax, gene in zip(axes, genes):
        values = X[:, ctx["genes"].index(gene)]
        data = [values[np.array([str(l) for l in labels]) == g] for g in group_names]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                        medianprops={"color": "#2b2f36", "linewidth": 1.1})
        for patch, name in zip(bp["boxes"], group_names):
            patch.set_facecolor(colors.get(name, "#4e79a7"))
            patch.set_alpha(0.78)
            patch.set_edgecolor("#8a919c")
        ax.set_xticks(range(1, len(group_names) + 1))
        ax.set_xticklabels(group_names, fontsize=7.5)
        ax.set_title(gene, fontsize=9)
    fig.suptitle("Gene expression by group", fontsize=10)
    return save_fig(fig, outdir, f"genebox-{_hash('genebox', genes, group_names)}.png", "Gene boxplot")


def plot_enrich_bar(ctx, params, outdir):
    frame = _store(ctx, "enrich")
    if frame is None or frame.empty:
        raise prerequisite("还没有功能富集结果，请先运行「功能富集」")
    top_n = int(params.get("top_n", 15) or 15)
    label_col = "Description" if "Description" in frame.columns else frame.columns[min(1, len(frame.columns) - 1)]
    frame = frame.head(top_n).iloc[::-1]
    fig, ax = new_fig(9.2, max(3.6, 0.36 * len(frame) + 1.4))
    values = -np.log10(np.clip(frame["p.adjust"].astype(float).values if "p.adjust" in frame.columns else frame["pvalue"].astype(float).values, 1e-300, 1))
    ax.barh(np.arange(len(frame)), values, color="#4e79a7", alpha=0.85)
    ax.set_yticks(np.arange(len(frame)))
    ax.set_yticklabels([str(v)[:60] for v in frame[label_col]], fontsize=7.5)
    ax.set_xlabel("-log10(adjusted P)")
    ax.set_title("Enrichment (top terms)")
    return save_fig(fig, outdir, f"enrich_bar-{_hash('enrich_bar', top_n)}.png", "Enrichment bar")


def plot_enrich_dot(ctx, params, outdir):
    frame = _store(ctx, "enrich")
    if frame is None or frame.empty:
        raise prerequisite("还没有功能富集结果，请先运行「功能富集」")
    top_n = int(params.get("top_n", 15) or 15)
    frame = frame.head(top_n).iloc[::-1]
    label_col = "Description" if "Description" in frame.columns else frame.columns[min(1, len(frame.columns) - 1)]
    count = frame["Count"].astype(float).values if "Count" in frame.columns else np.ones(len(frame))
    padj = frame["p.adjust"].astype(float).values if "p.adjust" in frame.columns else frame["pvalue"].astype(float).values
    fig, ax = new_fig(9.4, max(3.6, 0.34 * len(frame) + 1.4))
    sc = ax.scatter(-np.log10(np.clip(padj, 1e-300, 1)), np.arange(len(frame)),
                    s=24 + 40 * (count / max(count.max(), 1)), c=-np.log10(np.clip(padj, 1e-300, 1)),
                    cmap="viridis", edgecolor="white", linewidth=0.5)
    ax.set_yticks(np.arange(len(frame)))
    ax.set_yticklabels([str(v)[:60] for v in frame[label_col]], fontsize=7.5)
    ax.set_xlabel("-log10(adjusted P)")
    ax.set_title("Enrichment (dot)")
    fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02, label="-log10(adj.P)")
    return save_fig(fig, outdir, f"enrich_dot-{_hash('enrich_dot', top_n)}.png", "Enrichment dot")


def plot_immune_stack(ctx, params, outdir):
    frame = _store(ctx, "immune")
    if frame is None or frame.empty:
        raise prerequisite("还没有免疫浸润结果，请先运行「免疫浸润」")
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.empty:
        raise prerequisite("免疫结果里没有数值列")
    data = numeric.iloc[:200]
    fig, ax = new_fig(10.6, 5.2)
    bottom = np.zeros(len(data))
    for i, col in enumerate(data.columns[:40]):
        values = data[col].fillna(0).astype(float).values
        if values.max() <= 1.2 and values.min() >= -0.1:
            values = values / max(values.sum(), 1e-9) * 100
        ax.bar(np.arange(len(data)), values, bottom=bottom, color=color_for(str(col), i), width=1.0, label=str(col))
        bottom += values
    ax.set_xlabel("Samples")
    ax.set_ylabel("Estimated proportion / score")
    ax.set_title("Immune infiltration")
    ax.grid(False)
    if data.columns.size <= 24:
        ax.legend(fontsize=6, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return save_fig(fig, outdir, f"immune_stack-{_hash('immune_stack', data.shape)}.png", "Immune stack")


def plot_immune_corr(ctx, params, outdir):
    frame = _store(ctx, "immune")
    if frame is None or frame.empty:
        raise prerequisite("还没有免疫浸润结果")
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.shape[1] < 2:
        raise prerequisite("免疫结果列数不足 2")
    corr = numeric.corr().values
    fig, ax = new_fig(7.6, 6.6)
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(numeric.columns)))
    ax.set_yticks(range(len(numeric.columns)))
    ax.set_xticklabels([str(c)[:16] for c in numeric.columns], rotation=90, fontsize=6.5)
    ax.set_yticklabels([str(c)[:20] for c in numeric.columns], fontsize=6.5)
    ax.grid(False)
    ax.set_title("Immune cell correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    return save_fig(fig, outdir, f"immune_corr-{_hash('immune_corr', numeric.shape)}.png", "Immune correlation")


def plot_km(ctx, params, outdir):
    stored = (ctx.get("results") or {}).get("survival")
    info = ctx.get("sample_info")
    if info is None:
        raise prerequisite("缺少样本信息")
    time_col = params.get("time_col") or (stored or {}).get("time_col")
    event_col = params.get("event_col") or (stored or {}).get("event_col")
    if (not time_col or not event_col) and info is not None:
        # 自动识别：时间列/事件列常见命名
        time_hints = ("os_time", "time", "survival_time", "dfs_time", "rfs_time", "pfs_time", "month", "months", "days")
        event_hints = ("os_event", "event", "status", "dead", "death", "censor", "vital")
        lower = {str(c).lower(): str(c) for c in info.columns}
        if not time_col:
            for hint in time_hints:
                match = next((orig for low, orig in lower.items() if hint == low or low.endswith(hint)), None)
                if match:
                    time_col = match
                    break
        if not event_col:
            for hint in event_hints:
                match = next((orig for low, orig in lower.items() if hint == low or low.endswith(hint)), None)
                if match:
                    event_col = match
                    break
    if not time_col or not event_col:
        raise prerequisite("请指定生存时间列与事件列（time_col/event_col）")
    gene = params.get("gene")
    group = params.get("group_col")
    if gene and gene in ctx["genes"]:
        values = np.asarray(ctx["X"], dtype=float)[:, ctx["genes"].index(gene)]
        cutoff = float(np.median(values))
        labels = np.where(values > cutoff, "high", "low")
        title = f"{gene} (median split)"
    elif group:
        labels = info[group].astype(str).values
        title = f"By {group}"
    else:
        labels = _groups(ctx)
        title = "By group"
    try:
        from lifelines import KaplanMeierFitter
        from lifelines.statistics import multivariate_logrank_test
    except Exception as error:  # noqa: BLE001
        raise EngineError("unsupported", f"lifelines 不可用：{error}")
    times = pd.to_numeric(info[time_col], errors="coerce").values
    events = pd.to_numeric(info[event_col], errors="coerce").fillna(0).values
    mask = np.isfinite(times) & (times > 0)
    fig, ax = new_fig(7.4, 5.4)
    colors = _color_map(labels)
    for name in dict.fromkeys([str(l) for l in labels]):
        sel = mask & (np.array([str(v) for v in labels]) == name)
        if sel.sum() < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(times[sel], events[sel], label=f"{name} (n={int(sel.sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=colors.get(name))
    p_value = None
    try:
        res = multivariate_logrank_test(times[mask], labels[mask], events[mask])
        p_value = float(res.p_value)
    except Exception:  # noqa: BLE001
        pass
    if p_value is not None:
        ax.set_title(f"{title}  log-rank P = {p_value:.3g}")
    else:
        ax.set_title(title)
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    return save_fig(fig, outdir, f"km-{_hash('km', str(gene or group), labels[:5])}.png", "Kaplan-Meier")


def plot_forest(ctx, params, outdir):
    frame = _store(ctx, "survival")
    if frame is None or frame.empty:
        raise prerequisite("还没有生存分析结果")
    if not {"HR", "gene"}.issubset(set(map(str, frame.columns))):
        raise prerequisite("生存结果缺少 HR/gene 列")
    frame = frame.head(25).iloc[::-1]
    hr = pd.to_numeric(frame["HR"], errors="coerce").values
    lo = pd.to_numeric(frame["lower"], errors="coerce").values if "lower" in frame.columns else hr * 0.8
    hi = pd.to_numeric(frame["upper"], errors="coerce").values if "upper" in frame.columns else hr * 1.2
    fig, ax = new_fig(7.6, max(3.4, 0.32 * len(frame) + 1.2))
    y = np.arange(len(frame))
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="o", color="#c0392b", ecolor="#8a919c", capsize=3, markersize=4)
    ax.axvline(1.0, color="#2b2f36", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([str(v) for v in frame["gene"]], fontsize=7.5)
    ax.set_xlabel("Hazard ratio (95% CI)")
    ax.set_title("Forest plot")
    return save_fig(fig, outdir, f"forest-{_hash('forest', frame.shape)}.png", "Forest")


def plot_risk(ctx, params, outdir):
    stored = (ctx.get("results") or {}).get("survival")
    frame = None
    if isinstance(stored, dict):
        candidate = stored.get("risk_table")
        if isinstance(candidate, pd.DataFrame):
            frame = candidate
        elif isinstance(candidate, list) and candidate:
            frame = pd.DataFrame(candidate)
    if frame is None or frame.empty:
        raise prerequisite("还没有风险评分结果，请先运行「高级预后」")
    score_col = "risk_score" if "risk_score" in frame.columns else frame.select_dtypes(include=[np.number]).columns[0]
    frame = frame.sort_values(score_col)
    fig, axes = new_fig(9.6, 4.8, ncols=2)
    axes = np.atleast_1d(axes)
    colors = np.where(frame.get("risk_group", pd.Series(["low"] * len(frame))).astype(str) == "high", "#c0392b", "#2e6da4")
    axes[0].scatter(np.arange(len(frame)), frame[score_col], s=18, c=colors)
    axes[0].axhline(float(np.median(frame[score_col])), color="#8a919c", linestyle="--", linewidth=0.8)
    axes[0].set_title("Risk score")
    axes[0].set_xlabel("Samples (sorted)")
    for col, title, color in [("time", "Survival time", "#2b2f36"), ("event", "Event", "#4e79a7")]:
        pass
    # 右图：时间/事件散点
    info = ctx.get("sample_info")
    if info is not None:
        time_col = params.get("time_col") or (stored or {}).get("time_col")
        event_col = params.get("event_col") or (stored or {}).get("event_col")
        if time_col and time_col in info.columns:
            ordered = [ctx["samples"].index(s) for s in frame.index if s in ctx["samples"]]
            if ordered:
                times = pd.to_numeric(info[time_col].iloc[ordered], errors="coerce").values
                events = pd.to_numeric(info[event_col].iloc[ordered], errors="coerce").fillna(0).values if event_col else np.zeros(len(ordered))
                axes[1].scatter(np.arange(len(times)), times, c=np.where(events > 0, "#c0392b", "#7f8c9b"), s=18)
                axes[1].set_title("Survival status")
                axes[1].set_xlabel("Samples (sorted by risk)")
                axes[1].set_ylabel("Time")
    return save_fig(fig, outdir, f"risk-{_hash('risk', len(frame))}.png", "Risk score")


def plot_roc_time(ctx, params, outdir):
    stored = (ctx.get("results") or {}).get("survival")
    roc = None
    if isinstance(stored, dict):
        roc = stored.get("roc")
    if roc is None:
        raise prerequisite("还没有时间依赖 ROC 结果，请先运行「高级预后」或「临床生存」")
    fig, ax = new_fig(6.6, 5.6)
    items = roc if isinstance(roc, list) else list(roc.items()) if isinstance(roc, dict) else []
    for i, item in enumerate(items):
        if isinstance(item, dict):
            fpr = item.get("fpr") or item.get("x")
            tpr = item.get("tpr") or item.get("y")
            label = item.get("label") or item.get("time") or f"time{i}"
            auc = item.get("auc")
        else:
            label, payload = item
            fpr = payload.get("fpr") if isinstance(payload, dict) else None
            tpr = payload.get("tpr") if isinstance(payload, dict) else None
            auc = payload.get("auc") if isinstance(payload, dict) else None
        if not fpr or not tpr:
            continue
        ax.plot(fpr, tpr, color=palette(len(items))[i], linewidth=1.5,
                label=f"{label}" + (f" AUC={float(auc):.3f}" if auc is not None else ""))
    ax.plot([0, 1], [0, 1], color="#b9bfc8", linestyle="--", linewidth=0.8)
    ax.set_xlabel("1 - Specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_title("Time-dependent ROC")
    ax.legend(fontsize=8)
    return save_fig(fig, outdir, f"roc-{_hash('roc', len(items))}.png", "Time ROC")


def plot_batch_pca(ctx, params, outdir):
    info = ctx.get("sample_info")
    batch_col = params.get("batch_key") or (info.columns[0] if info is not None and len(info.columns) else None)
    if info is None or not batch_col or batch_col not in info.columns:
        raise prerequisite("样本信息里没有批次列（batch_key）")
    X = np.asarray(ctx["X"], dtype=float)
    scores, var = _pca_components(X, 5)
    labels = info[batch_col].astype(str).values
    colors = _color_map(labels)
    fig, axes = new_fig(9.2, 4.8, ncols=2)
    axes = np.atleast_1d(axes)
    modes = [("Raw expression", X)]
    corrected = None
    if isinstance((ctx.get("results") or {}).get("batch"), dict):
        maybe = (ctx["results"]["batch"]).get("corrected")
        if isinstance(maybe, np.ndarray):
            corrected = maybe
    if corrected is not None:
        modes.append(("After correction", corrected))
    for ax, (title, mat) in zip(axes, modes):
        sc, _ = _pca_components(mat, 3)
        for name in dict.fromkeys([str(l) for l in labels]):
            mask = np.array([str(v) == name for v in labels])
            ax.scatter(sc[mask, 0], sc[mask, 1], s=30, alpha=0.8, color=colors[name], label=name, edgecolor="white", linewidth=0.5)
        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
    axes[0].legend(fontsize=7)
    return save_fig(fig, outdir, f"batchpca-{_hash('batch', batch_col)}.png", "Batch PCA")


def plot_circos(ctx, params, outdir):
    module = (ctx.get("modules") or {}).get("circ")
    if not module:
        raise prerequisite("还没有挂载 circRNA 数据")
    frame = module.get("table")
    if not isinstance(frame, pd.DataFrame):
        frame = module.get("info") if isinstance(module.get("info"), pd.DataFrame) else None
    ids = []
    if isinstance(frame, pd.DataFrame) and len(frame.columns):
        for col in frame.columns:
            if str(col).lower() in ("circ_id", "circrna", "id", "circid"):
                ids = [str(v) for v in frame[col].tolist()]
                break
        if not ids:
            ids = [str(v) for v in frame.iloc[:, 0].tolist()]
    elif module.get("genes"):
        ids = [str(v) for v in module["genes"]]
    if not ids:
        raise prerequisite("circRNA 模块里没有可解析的 ID")
    import re

    coords = []
    for cid in ids:
        match = re.search(r"chr([0-9XYM]+)[:\-_](\d+)[-_:](\d+)", cid)
        if match:
            coords.append((match.group(1), int(match.group(2)), int(match.group(3)), cid))
    if not coords:
        raise prerequisite("circRNA ID 无法解析出染色体坐标（示例：chr1:1000-2000:+）")
    chroms = list(dict.fromkeys([c[0] for c in coords]))
    fig, ax = new_fig(7.2, 7.2, subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    n_chr = len(chroms)
    for i, (chrom, start, end, cid) in enumerate(coords):
        idx = chroms.index(chrom)
        theta = (idx + 0.5) / n_chr * 2 * np.pi
        length = np.clip((end - start) / 1e6, 0.05, 3.0)
        ax.bar(theta, length + 0.4, bottom=0.6, width=0.9, color=palette(n_chr)[idx], alpha=0.8)
    ax.set_yticklabels([])
    ax.set_xticks(np.arange(n_chr) / n_chr * 2 * np.pi + np.pi / n_chr)
    ax.set_xticklabels([f"chr{c}" for c in chroms], fontsize=8)
    ax.set_title("circRNA genomic distribution")
    return save_fig(fig, outdir, f"circos-{_hash('circos', len(coords))}.png", "circRNA landscape")


def plot_splicing_heat(ctx, params, outdir):
    module = (ctx.get("modules") or {}).get("splicing")
    if not module:
        raise prerequisite("还没有挂载可变剪接数据")
    frame = module.get("table")
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise prerequisite("可变剪接数据为空")
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.empty:
        raise prerequisite("可变剪接表里没有数值列（PSI/IncLevel）")
    mat = zscore_rows(numeric.head(60).values.T)
    fig, ax = new_fig(9.2, 5.4)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_title("Splicing events (row z-score)")
    ax.set_xlabel("Events")
    ax.set_ylabel("Samples")
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    return save_fig(fig, outdir, f"splicing-{_hash('splicing', mat.shape)}.png", "Splicing heatmap")


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [v.strip() for v in value.replace(",", " ").split() if v.strip()]
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value if str(v).strip()]
    return [str(value)]


HANDLERS = {
    "qc_box": plot_qc_box,
    "qc_density": plot_qc_density,
    "pca_scatter": plot_pca_scatter,
    "sample_corr": plot_sample_corr,
    "sample_dendro": plot_sample_dendro,
    "volcano": plot_volcano,
    "ma_plot": plot_ma,
    "de_heatmap": plot_de_heatmap,
    "gene_box": plot_gene_box,
    "enrich_bar": plot_enrich_bar,
    "enrich_dot": plot_enrich_dot,
    "immune_stack": plot_immune_stack,
    "immune_corr": plot_immune_corr,
    "km_curve": plot_km,
    "forest": plot_forest,
    "risk_plot": plot_risk,
    "roc_time": plot_roc_time,
    "batch_pca": plot_batch_pca,
    "circos": plot_circos,
    "splicing_heat": plot_splicing_heat,
}


def render(ctx, kind: str, params: dict, outdir: str) -> dict:
    kind = str(kind or "")
    handler = HANDLERS.get(kind)
    if handler is None:
        raise EngineError("bad_request", f"不支持的图类型：{kind}（可选：{', '.join(PLOT_KINDS)}）")
    ensure_dir(outdir)
    figure = handler(ctx, params or {}, outdir)
    figure["kind"] = kind
    try:
        figure["bytes"] = os.path.getsize(figure["file"])
    except OSError:
        figure["bytes"] = None
    return figure
