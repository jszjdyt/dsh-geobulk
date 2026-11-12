"""
dsh-geobulk — 功能富集模块（步骤 enrich + 通路类分析）

实现策略（诚实标注）：
  * ORA：超几何检验（scipy），基因集来自 bg_genesets 内置 GMT / 在线 Enrichr 缓存；
    输出列名对齐 clusterProfiler（ID/Description/GeneRatio/BgRatio/pvalue/p.adjust/Count/geneID）。
  * GSEA：加权 Kolmogorov-Smirnov running-sum + 基因标签置换 p 值（≥1000 次，固定 seed）。
  * ssGSEA/GSVA-lite：按样本做秩次加权的富集打分。
  * 不依赖 R，保证任何环境下可用；方法与集群来源在 notes 里写明。

注册：
  STEPS   = [enrich]
  ANALYSES= ora_go/ora_kegg/ora_msigdb/gsea_go/gsea_kegg/gsea_hallmark/
            gsva_pathway/ssgsea_pathway/pathway_heatmap/enrich_summary
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
)
from bg_genesets import available_dbs, get_gene_sets

MIN_SET_SIZE = 5
MAX_SET_SIZE = 500
DEFAULT_DBS = ["GO_BP", "GO_CC", "GO_MF", "KEGG"]
DB_LABELS = {
    "GO_BP": "GO Biological Process",
    "GO_CC": "GO Cellular Component",
    "GO_MF": "GO Molecular Function",
    "KEGG": "KEGG Pathway",
    "HALLMARK": "MSigDB Hallmark",
    "REACTOME": "Reactome",
    "IMMUNE_28": "Immune cell signatures",
}


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def _de_frame(ctx) -> pd.DataFrame:
    de = ctx.get("de")
    if isinstance(de, pd.DataFrame):
        frame = de
    elif isinstance(de, dict):
        frame = de.get("table")
        if not isinstance(frame, pd.DataFrame):
            rows = de.get("rows")
            frame = pd.DataFrame(rows) if rows else None
    else:
        frame = None
    if frame is None or frame.empty:
        raise EngineError("prerequisite", "还没有差异表达结果，请先运行「差异表达」")
    frame = frame.copy()
    for col in ("logFC", "P.Value", "adj.P.Val"):
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame.dropna(subset=["gene"])
    return frame


def _case_control(ctx, params) -> tuple[str, str]:
    de = ctx.get("de") if isinstance(ctx.get("de"), dict) else {}
    case = params.get("case") or de.get("case")
    control = params.get("control") or de.get("control")
    if not case or not control:
        groups = ctx.get("groups") or []
        if len(groups) >= 2:
            case = case or groups[0]
            control = control or groups[1]
        else:
            raise EngineError("prerequisite", "请提供 case / control 分组名")
    return str(case), str(control)


def _set_genes(value) -> list[str]:
    if isinstance(value, dict):
        return [str(g) for g in (value.get("genes") or [])]
    if isinstance(value, (list, tuple, set)):
        return [str(g) for g in value]
    return []


def _set_name(key: str, value) -> str:
    if isinstance(value, dict) and value.get("name"):
        return str(value.get("name"))
    return str(key)


def _prepare_sets(db: str, workdir=None, min_size: int = MIN_SET_SIZE, max_size: int = MAX_SET_SIZE) -> dict:
    raw = get_gene_sets(db, workdir)
    if not raw:
        raise EngineError("unsupported", f"基因集库 {db} 不可用（内置与在线均未取到）")
    sets: dict[str, list[str]] = {}
    names: dict[str, str] = {}
    for key, value in raw.items():
        genes = _set_genes(value)
        if not genes:
            continue
        if len(genes) < min_size or len(genes) > max_size:
            continue
        sets[str(key)] = genes
        names[str(key)] = _set_name(key, value)
    if not sets:
        raise EngineError("unsupported", f"基因集库 {db} 里没有满足大小限制的条目")
    return {"sets": sets, "names": names, "source": str(db)}


def _universe(ctx, de: pd.DataFrame) -> set[str]:
    universe = {str(g).upper() for g in ctx.get("genes") or []}
    if not universe:
        universe = {str(g).upper() for g in de["gene"].astype(str)}
    return universe


def _gene_list(de: pd.DataFrame, direction: str, padj_cut: float, lfc_cut: float) -> list[str]:
    frame = de
    if padj_cut and "adj.P.Val" in frame.columns:
        frame = frame[frame["adj.P.Val"] <= float(padj_cut)]
    if lfc_cut and "logFC" in frame.columns:
        frame = frame[frame["logFC"].abs() >= float(lfc_cut)]
    if direction == "up":
        frame = frame[frame["logFC"] > 0]
    elif direction == "down":
        frame = frame[frame["logFC"] < 0]
    return [str(g).upper() for g in frame["gene"].astype(str)]


# ---------------------------------------------------------------------------
# 统计核心
# ---------------------------------------------------------------------------

def _hypergeom_pvalue(k: int, K: int, n: int, N: int) -> float:
    """P(X >= k)：N=背景基因数，K=通路基因数，n=查询基因数，k=命中数。"""
    if k <= 0:
        return 1.0
    try:
        from scipy.stats import hypergeom

        return float(hypergeom.sf(k - 1, N, K, n))
    except Exception:  # noqa: BLE001
        # 组合数回退
        from math import comb

        total = comb(N, n)
        if total <= 0:
            return 1.0
        acc = 0
        for i in range(k, min(K, n) + 1):
            acc += comb(K, i) * comb(N - K, n - i)
        return float(min(1.0, acc / total))


def ora(genes, universe, gene_sets: dict, names: dict, db_label: str) -> pd.DataFrame:
    query = {str(g).upper() for g in genes}
    universe = {str(g).upper() for g in universe}
    query &= universe
    if len(query) < 3:
        return pd.DataFrame()
    rows = []
    N = len(universe)
    for key, members in gene_sets.items():
        member_set = {str(g).upper() for g in members}
        in_universe = member_set & universe
        K = len(in_universe)
        if K < MIN_SET_SIZE:
            continue
        hits = query & in_universe
        k = len(hits)
        if k < 2:
            continue
        pvalue = _hypergeom_pvalue(k, K, len(query), N)
        rows.append({
            "Database": db_label,
            "ID": key,
            "Description": names.get(key, key),
            "GeneRatio": f"{k}/{len(query)}",
            "BgRatio": f"{K}/{N}",
            "Count": k,
            "pvalue": pvalue,
            "geneID": "/".join(sorted(hits)),
            "FoldEnrichment": (k / max(len(query), 1)) / (K / max(N, 1)),
        })
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    frame["p.adjust"] = bh_fdr(frame["pvalue"].values)
    frame = frame.sort_values(["p.adjust", "pvalue"]).reset_index(drop=True)
    return frame


def _es_score(ranked_genes: list[str], gene_set: set[str], weights: np.ndarray) -> float:
    n = len(ranked_genes)
    hits = np.array([1.0 if gene in gene_set else 0.0 for gene in ranked_genes])
    n_hit = hits.sum()
    if n_hit == 0 or n_hit >= n:
        return 0.0
    hit_weights = np.abs(weights) * hits
    denom = hit_weights.sum() or 1.0
    running = np.cumsum(hit_weights / denom - (1 - hits) / max(n - n_hit, 1))
    idx = int(np.argmax(np.abs(running)))
    return float(running[idx])


def gsea(de: pd.DataFrame, gene_sets: dict, names: dict, db_label: str, nperm: int = 1000, seed: int = 0,
         min_size: int = MIN_SET_SIZE, max_size: int = MAX_SET_SIZE) -> pd.DataFrame:
    if "logFC" not in de.columns:
        raise EngineError("prerequisite", "GSEA 需要差异结果的 logFC 列")
    frame = de.dropna(subset=["logFC"]).copy()
    frame["gene"] = frame["gene"].astype(str).str.upper()
    frame = frame.drop_duplicates(subset=["gene"])
    frame = frame.sort_values("logFC", ascending=False)
    ranked = frame["gene"].tolist()
    weights = pd.to_numeric(frame["logFC"], errors="coerce").fillna(0.0).values
    n = len(ranked)
    if n < 50:
        raise EngineError("prerequisite", f"排序基因只有 {n} 个，太少，无法做 GSEA")
    rng = np.random.default_rng(seed)
    rows = []
    for key, members in gene_sets.items():
        member_set = {str(g).upper() for g in members}
        size = len(member_set & set(ranked))
        if size < min_size or size > max_size:
            continue
        es = _es_score(ranked, member_set, weights)
        # 置换：随机抽取同样大小的基因集
        null = np.empty(nperm, dtype=np.float64)
        pool = np.array(ranked)
        for i in range(nperm):
            sample = rng.choice(pool, size=size, replace=False)
            null[i] = _es_score(ranked, set(sample.tolist()), weights)
        if es >= 0:
            pvalue = float((np.sum(null >= es) + 1) / (nperm + 1))
        else:
            pvalue = float((np.sum(null <= es) + 1) / (nperm + 1))
        rows.append({
            "Database": db_label,
            "ID": key,
            "Description": names.get(key, key),
            "SetSize": int(size),
            "ES": es,
            "NES": float(es / (np.std(null) + 1e-9)) if np.std(null) > 0 else 0.0,
            "pvalue": pvalue,
            "geneID": "/".join(sorted(member_set & set(ranked))[:200]),
        })
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["p.adjust"] = bh_fdr(out["pvalue"].values)
    return out.sort_values(["p.adjust", "pvalue"]).reset_index(drop=True)


def ssgsea_scores(ctx, gene_sets: dict, names: dict, max_sets: int = 60) -> pd.DataFrame:
    """按样本的 ssGSEA 简化打分（秩次加权）。"""
    X = np.asarray(ctx["X"], dtype=np.float64)
    genes = [str(g).upper() for g in ctx["genes"]]
    index = {g: i for i, g in enumerate(genes)}
    # 每个样本计算秩
    ranks = np.argsort(np.argsort(-X, axis=1), axis=1)  # 0=最高
    n_genes = X.shape[1]
    score_rows = {}
    used = 0
    for key, members in gene_sets.items():
        idx = [index[g] for g in (str(x).upper() for x in members) if g in index]
        if len(idx) < MIN_SET_SIZE:
            continue
        hits = np.zeros((X.shape[0], n_genes), dtype=bool)
        hits[:, idx] = True
        # 加权 rank 得分：命中基因秩的均值越小 → 越富集
        ranks_f = ranks + 1.0
        weights = np.where(hits, 1.0 / ranks_f, 0.0)
        score = weights.sum(axis=1)
        # 归一：减去随机期望（用全部基因的 1/rank 均值）
        baseline = (1.0 / np.arange(1, n_genes + 1)).mean()
        z = (score - len(idx) * baseline) / (np.sqrt(len(idx)) * baseline + 1e-12)
        score_rows[names.get(key, key)] = z
        used += 1
        if used >= max_sets:
            break
    if not score_rows:
        raise EngineError("prerequisite", "没有基因集能与表达矩阵匹配（基因名可能不是 symbol）")
    return pd.DataFrame(score_rows, index=list(ctx["samples"])).T


# ---------------------------------------------------------------------------
# 图
# ---------------------------------------------------------------------------

def _enrich_figures(frame: pd.DataFrame, outdir: str, prefix: str, title: str, top_n: int = 15) -> list[dict]:
    figures = []
    if frame is None or frame.empty:
        return figures
    top = frame.head(top_n).iloc[::-1]
    labels = [str(v)[:52] for v in top["Description"]]
    values = -np.log10(np.clip(top["p.adjust"].astype(float).values, 1e-300, 1))
    fig, ax = new_fig(9.2, max(3.6, 0.34 * len(top) + 1.4))
    ax.barh(np.arange(len(top)), values, color="#4e79a7", alpha=0.88)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("-log10(adjusted P)")
    ax.set_title(f"{title} (top {len(top)})")
    figures.append(save_fig(fig, outdir, f"{prefix}_bar.png", f"{title} bar"))
    counts = top["Count"].astype(float).values if "Count" in top.columns else np.ones(len(top))
    fig, ax = new_fig(9.4, max(3.6, 0.34 * len(top) + 1.4))
    sc = ax.scatter(values, np.arange(len(top)),
                    s=26 + 34 * (counts / max(counts.max(), 1)),
                    c=values, cmap="viridis", edgecolor="white", linewidth=0.5)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("-log10(adjusted P)")
    ax.set_title(f"{title} (dot)")
    fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02, label="-log10(adj.P)")
    figures.append(save_fig(fig, outdir, f"{prefix}_dot.png", f"{title} dot"))
    return figures


# ---------------------------------------------------------------------------
# 步骤
# ---------------------------------------------------------------------------

def _run_enrichment(ctx, params: dict, outdir: str):
    de = _de_frame(ctx)
    dbs = params.get("db") or DEFAULT_DBS
    if isinstance(dbs, str):
        dbs = [d for d in dbs.replace(",", " ").split() if d]
    dbs = [d for d in dbs if d in available_dbs()] or ["GO_BP"]
    method = str(params.get("method") or "ora").lower()
    pvalue_cut = float(params.get("pvalue") or 0.05)
    qvalue_cut = float(params.get("qvalue") or 0.2)
    direction = str(params.get("direction") or "both")
    lfc_cut = float(params.get("lfc") or 0)
    padj_cut = float(params.get("padj") or 0.05)
    top_n = int(params.get("top_n") or 20)
    workdir = ctx.get("workdir")

    universe = _universe(ctx, de)
    tables = []
    frames = []
    warnings = []
    notes = []
    for db in dbs:
        try:
            prepared = _prepare_sets(db, workdir)
        except EngineError as error:
            warnings.append(f"{db} 跳过：{error.message}")
            continue
        label = DB_LABELS.get(db, db)
        if method == "gsea":
            frame = gsea(de, prepared["sets"], prepared["names"], label)
            if frame is not None and not frame.empty:
                frames.append(frame)
            continue
        directions = ["up", "down"] if direction in ("both", "all") else [direction]
        parts = []
        for item in directions:
            genes = _gene_list(de, item, padj_cut, lfc_cut) if item in ("up", "down") else _gene_list(de, "all", padj_cut, lfc_cut)
            if len(genes) < 3:
                continue
            frame = ora(genes, universe, prepared["sets"], prepared["names"], f"{label} ({'up' if item == 'up' else 'down'})")
            if frame is not None and not frame.empty:
                parts.append(frame)
        if parts:
            frames.append(pd.concat(parts, ignore_index=True))
        else:
            warnings.append(f"{db}：没有条目通过筛选（可放宽 padj/lfc 或换库）")

    if frames:
        combined = pd.concat(frames, ignore_index=True).sort_values(["p.adjust", "pvalue"]).reset_index(drop=True)
    else:
        combined = pd.DataFrame(columns=["Database", "ID", "Description", "GeneRatio", "BgRatio", "Count",
                                         "pvalue", "p.adjust", "geneID", "FoldEnrichment"])
        warnings.append("未发现显著富集条目（可放宽阈值或补充基因集）")

    sig = combined[combined["p.adjust"] <= pvalue_cut] if not combined.empty else combined
    if not combined.empty:
        tables.append(table("enrich", "富集结果", list(combined.columns)[:10], combined.head(500),
                            note=f"共 {len(combined)} 条，显著（p.adjust≤{pvalue_cut}）{len(sig)} 条"))
        path = export_csv(combined, outdir, "enrich_all")
        data_out = [data_file("enrich_csv", path, "富集全表")]
    else:
        data_out = []

    figures = _enrich_figures(sig if not sig.empty else combined, outdir, "enrich",
                              DB_LABELS.get(dbs[0], dbs[0]) if len(dbs) == 1 else "All databases", top_n=min(top_n, 20))
    bullets = [
        f"{'GSEA' if method == 'gsea' else 'ORA'}：{'+'.join(dbs)}，共 {len(combined)} 条，"
        f"显著 {len(sig)} 条（p.adjust≤{pvalue_cut}）",
    ]
    if not sig.empty:
        top = sig.iloc[0]
        bullets.append(f"最显著：{top['Description']}（p.adjust={float(top['p.adjust']):.2e}，Count={int(top['Count']) if 'Count' in top else '—'}）")
    metrics = {
        "method": method,
        "databases": dbs,
        "n_terms": int(len(combined)),
        "n_significant": int(len(sig)),
        "n_genes": int(len(_gene_list(de, direction, padj_cut, lfc_cut))),
    }
    payload = result(
        title="功能富集",
        bullets=bullets,
        metrics=metrics,
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=notes + [f"富集实现：{'加权 KS + 置换（GSEA）' if method == 'gsea' else '超几何检验（ORA）'}；基因集：内置 GMT / Enrichr 缓存"],
        data_out=data_out,
    )
    payload["store"] = {
        "enrich": {
            "method": method,
            "databases": dbs,
            "table": combined.head(2000).to_dict("records"),
            "n_terms": int(len(combined)),
            "n_significant": int(len(sig)),
        }
    }
    return payload


def run_enrich(ctx, params, outdir):
    return _run_enrichment(ctx, {**(params or {}), "method": (params or {}).get("method") or "ora"}, outdir)


def run_ora(ctx, params, outdir):
    return _run_enrichment(ctx, {**(params or {}), "method": "ora"}, outdir)


def run_gsea_step(ctx, params, outdir):
    return _run_enrichment(ctx, {**(params or {}), "method": "gsea"}, outdir)


# ---------------------------------------------------------------------------
# 分析
# ---------------------------------------------------------------------------

def _single_db_analysis(db: str):
    def run(ctx, params, outdir):
        return _run_enrichment(ctx, {**(params or {}), "db": [db], "method": str((params or {}).get("method") or "ora")}, outdir)
    return run


def run_pathway_scores(ctx, params, outdir):
    """ssGSEA/GSVA-lite：样本 × 通路活性表 + 热图。"""
    db = str((params or {}).get("db") or "HALLMARK")
    prepared = _prepare_sets(db, ctx.get("workdir"), max_size=500)
    scores = ssgsea_scores(ctx, prepared["sets"], prepared["names"], max_sets=int((params or {}).get("max_pathways") or 50))
    outdir_local = outdir
    frame = scores.copy()
    frame.insert(0, "pathway", frame.index)
    path = export_csv(frame, outdir_local, f"pathway_scores_{db.lower()}")
    # 热图
    mat = scores.values
    mat = (mat - np.nanmean(mat, axis=1, keepdims=True)) / (np.nanstd(mat, axis=1, keepdims=True) + 1e-9)
    fig, ax = new_fig(9.6, 6.2)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(np.arange(len(scores.index)))
    ax.set_yticklabels([str(v)[:34] for v in scores.index], fontsize=7)
    ax.set_xticks(np.arange(len(scores.columns)))
    ax.set_xticklabels(scores.columns, rotation=90, fontsize=6.5)
    ax.grid(False)
    ax.set_title(f"Pathway activity ({db})")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    figure = save_fig(fig, outdir_local, f"pathway_heat_{db.lower()}.png", "Pathway activity")
    top_paths = scores.mean(axis=1).sort_values(ascending=False).head(5)
    return result(
        title="通路活性打分",
        bullets=[f"{db}：{len(scores.index)} 条通路 × {len(scores.columns)} 样本",
                 "活性最高：" + "、".join(f"{k}({v:.2f})" for k, v in top_paths.items())],
        metrics={"db": db, "n_pathways": int(len(scores.index)), "n_samples": int(len(scores.columns))},
        tables=[table("pathway_scores", "样本×通路活性", [{"key": "pathway", "label": "Pathway", "type": "string"}] +
                      [{"key": str(s), "label": str(s), "type": "number"} for s in scores.columns],
                      frame.head(300))],
        figures=[figure],
        warnings=[],
        notes=["ssGSEA 简化实现：按样本内秩次加权，跨样本可比"],
        data_out=[data_file("pathway_csv", path, "通路活性表")],
        store={"enrich": {"pathway_scores_db": db}},
    )


def run_pathway_heatmap(ctx, params, outdir):
    return run_pathway_scores(ctx, params, outdir)


def run_enrich_summary(ctx, params, outdir):
    payload = _run_enrichment(ctx, {**(params or {}), "method": "ora"}, outdir)
    stored = (ctx.get("results") or {}).get("enrich") or {}
    frame = pd.DataFrame(stored.get("table") or [])
    if frame.empty:
        raise EngineError("prerequisite", "还没有富集结果，请先运行「功能富集」")
    counts = frame.groupby("Database").size().reset_index(name="n_terms") if "Database" in frame.columns else pd.DataFrame()
    bullets = [f"当前富集结果：{len(frame)} 条，来自 {frame['Database'].nunique() if 'Database' in frame.columns else 1} 个库"]
    return result(
        title="富集结果汇总",
        bullets=bullets,
        metrics={"n_terms": int(len(frame))},
        tables=[table("enrich_summary", "各库条目数", [{"key": "Database", "label": "数据库", "type": "string"},
                                                     {"key": "n_terms", "label": "条目数", "type": "number"}], counts)],
        figures=payload.get("figures", [])[:1],
    )


STEPS = [
    {
        "id": "enrich",
        "title": "功能富集",
        "requires": ["de"],
        "params": [
            {"key": "db", "label": "基因集库", "type": "select_multiple", "default": DEFAULT_DBS,
             "options": [{"value": k, "label": v} for k, v in DB_LABELS.items() if k in available_dbs()]},
            {"key": "method", "label": "方法", "type": "select", "default": "ora",
             "options": [{"value": "ora", "label": "ORA 超几何"}, {"value": "gsea", "label": "GSEA 预排序"}]},
            {"key": "direction", "label": "方向", "type": "select", "default": "both",
             "options": [{"value": "both", "label": "上调+下调"}, {"value": "up", "label": "仅上调"}, {"value": "down", "label": "仅下调"}]},
            {"key": "padj", "label": "差异基因 FDR", "type": "number", "default": 0.05, "min": 0, "max": 1, "step": 0.01},
            {"key": "lfc", "label": "差异倍数阈值", "type": "number", "default": 0, "min": 0, "max": 10, "step": 0.1},
            {"key": "pvalue", "label": "富集 p 阈值", "type": "number", "default": 0.05, "min": 0, "max": 1, "step": 0.01},
            {"key": "top_n", "label": "展示条目数", "type": "number", "default": 20, "min": 5, "max": 100, "step": 5},
        ],
        "run": run_enrich,
    }
]

ANALYSES = [
    {"id": "ora_go", "title": "GO 富集（ORA）", "group": "通路分析",
     "method": "超几何检验 + 内置 GO BP/CC/MF 基因集", "needs": ["de"],
     "params": [{"key": "db", "label": "GO 分支", "type": "select", "default": "GO_BP",
                 "options": [{"value": "GO_BP", "label": "BP"}, {"value": "GO_CC", "label": "CC"}, {"value": "GO_MF", "label": "MF"}]}],
     "run": _single_db_analysis("GO_BP")},
    {"id": "ora_kegg", "title": "KEGG 通路富集", "group": "通路分析",
     "method": "超几何检验 + 内置 KEGG 基因集", "needs": ["de"], "params": [], "run": _single_db_analysis("KEGG")},
    {"id": "ora_msigdb", "title": "MSigDB 富集（Hallmark/Reactome）", "group": "通路分析",
     "method": "超几何检验 + 内置 Hallmark/Reactome 基因集", "needs": ["de"],
     "params": [{"key": "db", "label": "库", "type": "select", "default": "HALLMARK",
                 "options": [{"value": "HALLMARK", "label": "Hallmark"}, {"value": "REACTOME", "label": "Reactome"}]}],
     "run": _single_db_analysis("HALLMARK")},
    {"id": "gsea_go", "title": "GSEA（GO）", "group": "通路分析",
     "method": "加权 KS running-sum + 1000 次置换", "needs": ["de"], "params": [],
     "run": lambda ctx, params, outdir: _run_enrichment(ctx, {**params, "db": ["GO_BP"], "method": "gsea"}, outdir)},
    {"id": "gsea_kegg", "title": "GSEA（KEGG）", "group": "通路分析",
     "method": "加权 KS running-sum + 置换", "needs": ["de"], "params": [],
     "run": lambda ctx, params, outdir: _run_enrichment(ctx, {**params, "db": ["KEGG"], "method": "gsea"}, outdir)},
    {"id": "gsea_hallmark", "title": "GSEA（Hallmark）", "group": "通路分析",
     "method": "加权 KS running-sum + 置换", "needs": ["de"], "params": [],
     "run": lambda ctx, params, outdir: _run_enrichment(ctx, {**params, "db": ["HALLMARK"], "method": "gsea"}, outdir)},
    {"id": "gsva_pathway", "title": "GSVA-lite 通路活性", "group": "通路分析",
     "method": "ssGSEA 简化（秩次加权）", "needs": [],
     "params": [{"key": "db", "label": "基因集库", "type": "select", "default": "HALLMARK",
                 "options": [{"value": "HALLMARK", "label": "Hallmark"}, {"value": "KEGG", "label": "KEGG"},
                             {"value": "REACTOME", "label": "Reactome"}, {"value": "GO_BP", "label": "GO BP"}]}],
     "run": run_pathway_scores},
    {"id": "ssgsea_pathway", "title": "ssGSEA 通路打分", "group": "通路分析",
     "method": "ssGSEA 简化（秩次加权）", "needs": [], "params": [], "run": run_pathway_scores},
    {"id": "pathway_heatmap", "title": "通路活性热图", "group": "通路分析",
     "method": "样本×通路活性 z-score 热图", "needs": [], "params": [], "run": run_pathway_heatmap},
    {"id": "enrich_summary", "title": "富集结果汇总", "group": "通路分析",
     "method": "按基因集库统计条目数", "needs": ["enrich"], "params": [], "run": run_enrich_summary},
]
