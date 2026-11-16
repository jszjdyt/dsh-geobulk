"""
dsh-geobulk — 高级预后模块（步骤 prognosis + 预后类分析）

实现：Cox（lifelines）/ LASSO-Cox（lifelines L1 惩罚；R glmnet 不可用时的替代）/
      RSF-lite（bootstrap Cox 集成）。输出风险评分、风险分组、KM、时间 ROC、C-index、
      校准曲线、决策曲线（DCA）与列线图数据。

注册：
  STEPS   = [prognosis]
  ANALYSES= risk_score_cox / risk_score_lasso / risk_score_rsf / risk_heatmap /
            cindex / calibration / dca / nomogram / time_roc_multi
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from bg_common import (
    EngineError,
    data_file,
    export_csv,
    new_fig,
    palette,
    result,
    save_fig,
    table,
    zscore_rows,
)
from bg_survival import clinical, time_roc, _roc_figure


# ---------------------------------------------------------------------------
# 基因选择与模型
# ---------------------------------------------------------------------------

def _candidate_genes(ctx, params, limit: int = 10) -> list[str]:
    genes = ctx.get("genes") or []
    requested = params.get("genes") or []
    if isinstance(requested, str):
        requested = [g.strip() for g in requested.replace(",", " ").split() if g.strip()]
    if requested:
        return [g for g in requested if g in set(genes)][:limit]
    de = ctx.get("de")
    frame = None
    if isinstance(de, dict):
        frame = de.get("table")
    elif isinstance(de, pd.DataFrame):
        frame = de
    if isinstance(frame, pd.DataFrame) and not frame.empty and "gene" in frame.columns:
        frame = frame.sort_values("adj.P.Val") if "adj.P.Val" in frame.columns else frame
        picked = [str(g) for g in frame["gene"].head(limit * 3)]
        return [g for g in picked if g in set(genes)][:limit]
    # 退化：取方差最大的基因
    X = np.asarray(ctx["X"], dtype=float)
    variance = np.nanvar(X, axis=0)
    order = np.argsort(-np.nan_to_num(variance))[:limit]
    return [genes[i] for i in order]


def _gene_matrix(ctx, genes: list[str]) -> np.ndarray:
    index = {g: i for i, g in enumerate(ctx.get("genes") or [])}
    X = np.asarray(ctx["X"], dtype=float)
    cols = [index[g] for g in genes if g in index]
    return X[:, cols]


def _fit_cox(frame: pd.DataFrame, time_col: str, event_col: str, columns: list[str], penalizer: float = 0.01,
             l1_ratio: float | None = None):
    from lifelines import CoxPHFitter

    data = frame[[time_col, event_col] + columns].copy()
    data = data.dropna()
    if len(data) < 10:
        raise EngineError("prerequisite", f"可用于建模的样本只有 {len(data)} 个（至少 10 个）")
    events = float(pd.to_numeric(data[event_col], errors="coerce").fillna(0).sum())
    if events < 5:
        raise EngineError("prerequisite", f"事件数只有 {int(events)} 个，无法建立预后模型（至少 5 个）")
    kwargs = {"penalizer": penalizer}
    if l1_ratio is not None:
        kwargs["l1_ratio"] = l1_ratio
    cph = CoxPHFitter(**kwargs)
    cph.fit(data, duration_col=time_col, event_col=event_col)
    return cph, data


def _risk_frame(data: pd.DataFrame, time_col: str, event_col: str, score: np.ndarray,
                cutoff: str = "median") -> pd.DataFrame:
    frame = pd.DataFrame({
        "sample": data.index.astype(str),
        "risk_score": score,
        "time": pd.to_numeric(data[time_col], errors="coerce").values,
        "event": pd.to_numeric(data[event_col], errors="coerce").fillna(0).values,
    })
    if cutoff == "quartile":
        threshold = float(np.nanquantile(score, 0.75))
    elif cutoff == "optimal":
        threshold = float(np.nanmedian(score))
    else:
        threshold = float(np.nanmedian(score))
    frame["risk_group"] = np.where(frame["risk_score"] > threshold, "high", "low")
    frame["threshold"] = threshold
    return frame


def _km_by_risk(frame: pd.DataFrame, outdir: str, title: str, name: str) -> tuple[dict | None, float | None]:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    high = frame["risk_group"] == "high"
    if high.sum() < 3 or (~high).sum() < 3:
        return None, None
    fig, ax = new_fig(7.4, 5.4)
    for mask, label, color in ((~high, "low", "#2e6da4"), (high, "high", "#c0392b")):
        subset = frame[mask]
        KaplanMeierFitter().fit(subset["time"], subset["event"], label=f"{label} (n={len(subset)})").plot_survival_function(
            ax=ax, ci_show=False, color=color)
    p_value = float(logrank_test(frame.loc[high, "time"], frame.loc[~high, "time"],
                                 frame.loc[high, "event"], frame.loc[~high, "event"]).p_value)
    ax.set_title(f"{title}  log-rank P = {p_value:.3g}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    return save_fig(fig, outdir, name, title), p_value


def _risk_plot(frame: pd.DataFrame, outdir: str, title: str, name: str) -> dict:
    frame = frame.sort_values("risk_score").reset_index(drop=True)
    fig, axes = new_fig(9.8, 4.8, ncols=2)
    axes = np.atleast_1d(axes)
    colors = np.where(frame["risk_group"] == "high", "#c0392b", "#2e6da4")
    axes[0].scatter(np.arange(len(frame)), frame["risk_score"], c=colors, s=22)
    axes[0].axhline(float(frame["threshold"].iloc[0]), color="#8a919c", linestyle="--", linewidth=0.9)
    axes[0].set_title("Risk score (sorted)")
    axes[0].set_xlabel("Samples")
    axes[0].set_ylabel("Risk score")
    axes[1].scatter(np.arange(len(frame)), frame["time"], c=np.where(frame["event"] > 0, "#c0392b", "#7f8c9b"), s=22)
    axes[1].set_title("Survival status")
    axes[1].set_xlabel("Samples (sorted by risk)")
    axes[1].set_ylabel("Time")
    return save_fig(fig, outdir, name, title)


def _cindex(time, event, risk) -> float:
    from lifelines.utils import concordance_index

    return float(concordance_index(time, -np.asarray(risk, dtype=float), event))


# ---------------------------------------------------------------------------
# 校准 / DCA
# ---------------------------------------------------------------------------

def _baseline_survival(cph, data: pd.DataFrame, time_col: str) -> tuple[np.ndarray, np.ndarray]:
    """基线生存 S0(t)：**直接取已拟合 Cox 的 baseline_survival_**（GEO-05）。

    旧实现把「所有患者预测生存率的均值」当 S0(t)，再按 exp(lp) 幂变换一次，
    等于把协变量效应施加了两遍。这里只从拟合模型自身属性导出，取不到就明确报错，
    绝不合成一条 `exp(-0.05*t)` 之类的假曲线。
    """
    frame = getattr(cph, "baseline_survival_", None)
    if frame is None or not isinstance(frame, pd.DataFrame) or frame.empty:
        raise EngineError("internal", "拟合好的 Cox 模型没有 baseline_survival_，无法给出基线生存曲线")
    timeline = np.asarray(frame.index, dtype=float)
    baseline = np.asarray(frame.iloc[:, 0], dtype=float)
    if not np.all(np.isfinite(baseline)):
        raise EngineError("internal", "拟合好的 Cox 模型 baseline_survival_ 含非有限值")
    return timeline, baseline


def _predict_survival(cph, data: pd.DataFrame, horizons: list[float]) -> np.ndarray:
    """按已拟合模型取个体生存概率（行=horizon，列=样本）。缺必要基因时明确失败。"""
    from lifelines import CoxPHFitter  # noqa: F401  # 仅确认依赖可导入

    times = [float(h) for h in horizons if float(h) > 0]
    if not times:
        raise EngineError("prerequisite", "没有可用的预测时间点")
    try:
        surv = cph.predict_survival_function(data, times=times)
    except Exception as error:  # noqa: BLE001
        raise EngineError(
            "prerequisite",
            f"无法用已拟合的 Cox 模型做个体预测（可能缺少模型必需的基因/协变量列）：{error}",
        )
    values = np.asarray(surv.values, dtype=float)
    if values.shape[0] == len(times) and values.shape[1] == len(data):
        return values
    if values.shape[0] == len(data) and values.shape[1] == len(times):
        return values.T
    raise EngineError("internal", f"predict_survival_function 返回了意外形状：{values.shape}")


def calibration_curve(cph, data: pd.DataFrame, time_col: str, event_col: str, horizons: list[float]) -> pd.DataFrame:
    """校准曲线：预测生存概率（同一已拟合模型）vs 分箱内 KM 观测生存率。"""
    horizons = sorted({float(h) for h in horizons if float(h) > 0})
    if not horizons:
        return pd.DataFrame()
    surv = _predict_survival(cph, data, horizons)
    rows = []
    from lifelines import KaplanMeierFitter

    for row_index, horizon in enumerate(horizons):
        predicted = np.clip(surv[row_index], 0.0, 1.0)
        bins = pd.qcut(pd.Series(predicted).rank(method="first"), q=4, labels=False, duplicates="drop")
        for b in sorted(pd.Series(bins).dropna().unique()):
            mask = np.asarray(bins == b)
            if mask.sum() < 3:
                continue
            try:
                kmf = KaplanMeierFitter().fit(
                    pd.to_numeric(data[time_col], errors="coerce").values[mask],
                    pd.to_numeric(data[event_col], errors="coerce").values[mask],
                )
                observed = float(kmf.predict(horizon))
            except Exception:  # noqa: BLE001
                observed = float(np.mean(pd.to_numeric(data[time_col], errors="coerce").values[mask] > horizon))
            rows.append({"horizon": horizon, "bin": int(b), "n": int(mask.sum()),
                         "predicted": float(np.mean(predicted[mask])), "observed": observed})
    return pd.DataFrame(rows)


def decision_curve(cph, data: pd.DataFrame, time_col: str, event_col: str, horizon: float, steps: int = 25) -> pd.DataFrame:
    """生存结局的决策曲线（KM 版净获益，GEO-05）。

    旧实现把「时点前未观察到事件」的人一律记为阴性再算 TP/FP，等于把时点前删失者
    当成随访充分且无事件，净获益会系统性偏移。这里改为在「按阈值处理 / 不处理」两组内
    分别用 Kaplan-Meier 估计 t 时点事件风险，再算净获益：

        NB(pt) = [1 - S_treat(t)] - [1 - S_notreat(t)] * pt/(1-pt)

    与 dcurves/rmda 的 survival DCA 采用同一 KM 估计思路；不使用 IPCW 加权。
    """
    from lifelines import KaplanMeierFitter

    horizon = float(horizon)
    if horizon <= 0:
        return pd.DataFrame()
    surv = _predict_survival(cph, data, [horizon])
    risk = 1.0 - np.clip(surv[0], 0.0, 1.0)
    times = pd.to_numeric(data[time_col], errors="coerce").values.astype(float)
    events = pd.to_numeric(data[event_col], errors="coerce").fillna(0).values.astype(float)
    n = len(times)

    def km_event(mask: np.ndarray) -> float | None:
        if int(mask.sum()) < 3 or float(events[mask].sum()) < 1:
            return None
        try:
            return float(1.0 - KaplanMeierFitter().fit(times[mask], events[mask]).predict(horizon))
        except Exception:  # noqa: BLE001
            return None

    overall = km_event(np.ones(n, dtype=bool))
    rows = []
    for threshold in np.linspace(0.02, 0.8, steps):
        positive = risk >= threshold
        if int(positive.sum()) < 3 or int((~positive).sum()) < 3:
            continue
        event_treat = km_event(positive)
        event_untreat = km_event(~positive)
        if event_treat is None or event_untreat is None:
            continue
        odds = threshold / (1.0 - threshold)
        net_benefit = event_treat - event_untreat * odds
        treat_all = None if overall is None else overall - (1.0 - overall) * odds
        rows.append({
            "threshold": float(threshold),
            "net_benefit": float(net_benefit),
            "treat_all": None if treat_all is None else float(treat_all),
            "event_risk_treat": float(event_treat),
            "event_risk_untreat": float(event_untreat),
            "n_positive": int(positive.sum()),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 步骤
# ---------------------------------------------------------------------------

def _prognosis_payload(ctx, params, outdir, model: str | None = None) -> dict:
    params = dict(params or {})
    model = (model or params.get("model") or "cox").lower()
    data = clinical(ctx, params)
    time_col, event_col = data["time_col"], data["event_col"]
    frame = data["frame"].copy()
    genes = _candidate_genes(ctx, params, limit=int(params.get("max_genes") or 10))
    if not genes:
        raise EngineError("prerequisite", "没有可用基因（请先做差异表达，或用 params.genes 指定）")
    matrix = _gene_matrix(ctx, genes)
    gene_frame = pd.DataFrame(matrix, index=[str(item) for item in ctx.get("samples") or []], columns=genes)
    gene_frame = gene_frame.reindex([str(i) for i in frame.index])
    for gene in genes:
        frame[gene] = gene_frame[gene].values
    warnings: list[str] = []
    notes: list[str] = []
    tables: list[dict] = []
    figures: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []

    if model == "lasso":
        cph, used = _fit_cox(frame, time_col, event_col, genes, penalizer=0.1, l1_ratio=1.0)
        selected = [g for g in genes if abs(float(cph.params_.get(g, 0.0))) > 1e-6]
        if not selected:
            selected = genes[: max(1, len(genes) // 2)]
        coefs = cph.params_.reindex(selected).fillna(0.0)
        risk = np.asarray(frame[selected].values @ coefs.values, dtype=float)
        method = "LASSO-Cox（lifelines L1，penalizer=0.1）"
        model_kind = "lasso_cox"
        notes.append("lifelines 的 L1 惩罚替代 glmnet；系数为惩罚后估计")
    elif model in ("rsf", "random_forest", "ensemble"):
        # 随机特征子集 + 患者 bootstrap 的 Cox 集成（**不是**随机生存森林，GEO-06）
        rng = np.random.default_rng(int(params.get("seed") or 0))
        n_boot = max(10, int(params.get("n_boot") or 200))
        n_genes = len(genes)
        subset_size = max(1, min(n_genes, int(np.sqrt(n_genes) + 2)))
        preds: list[np.ndarray] = []
        coef_rows_all: list[pd.Series] = []
        n_rows = len(frame)
        for _ in range(n_boot):
            size = max(2, min(subset_size, n_genes))
            pick = rng.choice(n_genes, size=size, replace=False) if size < n_genes else np.arange(n_genes)
            picked = [genes[i] for i in np.atleast_1d(pick)]
            # 患者（行）bootstrap：没有这一步就只是「随机抽基因」，不构成集成学习
            rows = rng.integers(0, n_rows, n_rows)
            sample_frame = frame.iloc[rows]
            try:
                cph_b, _ = _fit_cox(sample_frame, time_col, event_col, picked, penalizer=0.05)
                preds.append(np.asarray(cph_b.predict_log_partial_hazard(frame).values, dtype=float))
                coef_rows_all.append(cph_b.params_.reindex(genes).fillna(0.0))
            except EngineError:
                continue
        if not preds:
            raise EngineError("prerequisite", "bootstrap Cox 全部失败（样本/事件太少）")
        risk = np.nanmean(np.vstack(preds), axis=0)
        selected = genes
        method = f"随机特征子集 + 患者 bootstrap 的 Cox 集成（{len(preds)} 次）"
        model_kind = "subset_bootstrap_cox_ensemble"
        notes.append(
            "这不是随机生存森林（RSF）：实现是「随机抽基因子集 + 患者有放回重采样」的 Cox 线性预测均值；"
            "如需官方 RSF 请安装 sksurv 后使用专用分析"
        )
        # 校准/DCA 需要一个有预测接口的单一模型：对集成风险再拟合一个单变量 Cox（重校准模型）
        calib_frame = frame[[time_col, event_col]].copy()
        calib_frame["_ensemble_risk"] = risk
        cph, _ = _fit_cox(calib_frame, time_col, event_col, ["_ensemble_risk"], penalizer=0.0)
        coef_frame = pd.DataFrame([
            {"gene": gene,
             "coef": float(np.mean([c.get(gene, 0.0) for c in coef_rows_all])) if coef_rows_all else 0.0,
             "selection_freq": float(np.mean([abs(c.get(gene, 0.0)) > 1e-9 for c in coef_rows_all])) if coef_rows_all else 0.0}
            for gene in genes
        ])
        notes.append("系数表为集成内各次拟合的平均系数与入选频率；校准/DCA 用「对集成风险再拟合的单变量 Cox」")
    else:
        cph, _ = _fit_cox(frame, time_col, event_col, genes, penalizer=0.01)
        selected = genes
        risk = np.asarray(cph.predict_log_partial_hazard(frame).values, dtype=float)
        method = "Cox 比例风险模型（lifelines）"
        model_kind = "cox"

    # 系数表（rsf 分支已在上面生成）
    if model in ("rsf", "random_forest", "ensemble"):
        pass
    else:
        coef_frame = pd.DataFrame([
            {"gene": gene, "coef": float(cph.params_[gene]),
             "HR": float(np.exp(cph.params_[gene])),
             "p": float(cph.summary.loc[gene, "p"]) if gene in cph.summary.index else None}
            for gene in selected if gene in cph.params_.index
        ])
    if not coef_frame.empty:
        tables.append(table("model_coef", "模型系数", list(coef_frame.columns), coef_frame))

    risk_frame = _risk_frame(frame, time_col, event_col, risk, str(params.get("cutoff") or "median"))
    cutoff = float(risk_frame["threshold"].iloc[0])
    cindex = _cindex(risk_frame["time"].values, risk_frame["event"].values, risk_frame["risk_score"].values)
    km_figure, logrank_p = _km_by_risk(risk_frame, outdir, f"Risk groups ({model})", "risk_km.png")
    if km_figure:
        figures.append(km_figure)
    figures.append(_risk_plot(risk_frame, outdir, f"Risk score ({model})", "risk_plot.png"))

    # 风险热图（基因 z-score，按风险排序）
    mat = np.column_stack([frame[g].values for g in selected])
    z = zscore_rows(mat.T)
    order = np.argsort(risk_frame["risk_score"].values)
    fig, ax = new_fig(9.6, max(3.6, 0.22 * len(selected) + 1.6))
    im = ax.imshow(z[:, order], aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(np.arange(len(selected)))
    ax.set_yticklabels(selected, fontsize=7)
    ax.set_xticks([])
    ax.set_xlabel("Samples (sorted by risk)")
    ax.set_title("Risk gene heatmap")
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    figures.append(save_fig(fig, outdir, "risk_heatmap.png", "Risk heatmap"))

    # 时间 ROC
    horizons = [float(v) for v in np.nanquantile(risk_frame["time"], [0.25, 0.5, 0.75]) if v > 0]
    rocs = time_roc(risk_frame["time"].values, risk_frame["event"].values, risk_frame["risk_score"].values, horizons[:3])
    roc_figure = _roc_figure(rocs, outdir, "Time-dependent ROC (risk score)", "risk_roc.png")
    if roc_figure:
        figures.append(roc_figure)

    # 校准 & DCA & 列线图数据（Cox 系）
    calibration = pd.DataFrame()
    dca = pd.DataFrame()
    try:
        calibration = calibration_curve(cph, frame, time_col, event_col, horizons[:2])
        if not calibration.empty:
            fig, ax = new_fig(5.6, 5.0)
            ax.plot([0, 1], [0, 1], color="#b9bfc8", linestyle="--", linewidth=0.8)
            for i, horizon in enumerate(calibration["horizon"].unique()):
                subset = calibration[calibration["horizon"] == horizon]
                ax.plot(subset["predicted"], subset["observed"], marker="o", color=palette(4)[i],
                        label=f"{horizon:g} month")
            ax.set_xlabel("Predicted survival probability")
            ax.set_ylabel("Observed survival probability")
            ax.set_title("Calibration curve")
            ax.legend(fontsize=8)
            figures.append(save_fig(fig, outdir, "calibration.png", "Calibration"))
            tables.append(table("calibration", "校准曲线点", list(calibration.columns), calibration))
    except EngineError as error:
        warnings.append(f"校准曲线跳过：{error.message}")

    try:
        if horizons:
            dca = decision_curve(cph, frame, time_col, event_col, float(horizons[-1]))
            if not dca.empty:
                fig, ax = new_fig(6.2, 4.8)
                ax.plot(dca["threshold"], dca["net_benefit"], color="#c0392b", linewidth=1.6, label="Model")
                treat_all = dca.dropna(subset=["treat_all"])
                if not treat_all.empty:
                    ax.plot(treat_all["threshold"], treat_all["treat_all"], color="#7f8c9b",
                            linestyle="--", linewidth=1.0, label="Treat all")
                ax.axhline(0, color="#2b2f36", linewidth=0.8)
                ax.set_xlabel("Risk threshold")
                ax.set_ylabel("Net benefit")
                ax.set_title(f"Decision curve (t={float(horizons[-1]):g})")
                ax.legend(fontsize=8)
                figures.append(save_fig(fig, outdir, "dca.png", "Decision curve"))
                tables.append(table("dca", "决策曲线", list(dca.columns), dca))
    except EngineError as error:
        warnings.append(f"DCA 跳过：{error.message}")

    # 列线图数据（系数 + 基线生存）：基线直接来自已拟合模型，缺失则明确失败（GEO-05）
    try:
        timeline, baseline = _baseline_survival(cph, frame, time_col)
        # 系数一律取自系数表（rsf 集成分支的 cph 只是重校准模型，不含基因项）
        coef_lookup = {row["gene"]: float(row["coef"]) for _, row in coef_frame.iterrows()} if not coef_frame.empty else {}
        nomogram_rows = [{"variable": gene, "coef": coef_lookup[gene],
                          "HR": float(np.exp(coef_lookup[gene])),
                          "points_per_unit": float(abs(coef_lookup[gene]) * 100)}
                         for gene in selected if gene in coef_lookup]
        nomogram_rows.append({"variable": "baseline", "coef": 0.0, "HR": None, "points_per_unit": None})
        nomogram = pd.DataFrame(nomogram_rows)
        tables.append(table("nomogram", "列线图数据（系数与单位分值）", list(nomogram.columns), nomogram))
        tables.append(table("baseline_survival", "基线生存 S0(t)", [
            {"key": "time", "label": "时间", "type": "number"},
            {"key": "S0", "label": "S0(t)", "type": "number"}],
            [{"time": float(t), "S0": float(s)} for t, s in zip(timeline[::10], baseline[::10])]))
    except EngineError as error:
        warnings.append(f"列线图/基线生存未生成：{error.message}")

    tables.insert(0, table("risk_score", "风险评分表", list(risk_frame.columns), risk_frame))
    path = export_csv(risk_frame, outdir, "risk_scores")
    data_out.append(data_file("risk_csv", path, "风险评分表"))

    # GEO-06：明确标注这是训练集表观性能，不是外部验证结果
    evaluation_scope = "apparent_training_set"
    bullets.append(f"{method}：{len(selected)} 个基因，C-index={cindex:.3f}（**训练集表观性能**）")
    if logrank_p is not None:
        bullets.append(f"高低风险组 log-rank P={logrank_p:.3g}（训练集表观）")
    if rocs:
        bullets.append("时间 ROC：" + "、".join(f"{r['time']:g}={r['auc']:.3f}" for r in rocs) + "（训练集表观）")
    metrics = {"model": model, "model_kind": model_kind, "n_genes": len(selected),
               "cindex": cindex, "cindex_scope": evaluation_scope,
               "logrank_p": logrank_p,
               "cutoff": cutoff, "n_samples": int(len(risk_frame)),
               "n_events": int((risk_frame["event"] > 0).sum()),
               "auc": {f"{r['time']:g}": round(r["auc"], 4) for r in rocs}}
    metrics["evaluation_scope"] = "训练集表观性能（apparent）；基因筛选/拟合/阈值/指标均在同一样本上完成，未做外部验证"

    payload = result(
        title="高级预后分析",
        bullets=bullets,
        metrics=metrics,
        tables=tables,
        figures=figures,
        warnings=warnings + list(data.get("warnings") or []),
        notes=notes + [
            "风险评分 = 线性预测值；分组阈值为中位数（可用 cutoff 调整）",
            "评估口径：全部指标来自训练集表观性能，不能当作外部验证/前瞻验证结果；临床应用前需要独立队列验证",
        ],
        data_out=data_out,
    )
    payload["store"] = {
        "survival": {
            "time_col": time_col,
            "event_col": event_col,
            "time_col_source": "prognosis",
            "risk_table": risk_frame.to_dict("records"),
            "roc": rocs,
            "cindex": cindex,
            "model": model,
            "genes": selected,
        }
    }
    return payload


def run_prognosis(ctx, params, outdir):
    return _prognosis_payload(ctx, params, outdir)


# ---------------------------------------------------------------------------
# 分析
# ---------------------------------------------------------------------------

def run_risk_cox(ctx, params, outdir):
    return _prognosis_payload(ctx, params, outdir, model="cox")


def run_risk_lasso(ctx, params, outdir):
    return _prognosis_payload(ctx, params, outdir, model="lasso")


def run_risk_rsf(ctx, params, outdir):
    return _prognosis_payload(ctx, params, outdir, model="rsf")


def run_risk_heatmap(ctx, params, outdir):
    payload = _prognosis_payload(ctx, params, outdir, model=str((params or {}).get("model") or "cox"))
    keep = [fig for fig in payload.get("figures", []) if "risk_heatmap" in str(fig.get("file", ""))]
    payload["figures"] = keep or payload.get("figures", [])[:1]
    payload["summary"]["title"] = "风险基因热图"
    return payload


def run_cindex(ctx, params, outdir):
    data = clinical(ctx, params)
    frame = data["frame"]
    time_col, event_col = data["time_col"], data["event_col"]
    genes = _candidate_genes(ctx, params, limit=int((params or {}).get("max_genes") or 8))
    matrix = _gene_matrix(ctx, genes)
    sample_index = [str(item) for item in ctx.get("samples") or []]
    frame_index = [str(i) for i in frame.index]
    rows = []
    for i, gene in enumerate(genes):
        try:
            aligned = pd.Series(matrix[:, i], index=sample_index).reindex(frame_index).values
            score = _cindex(data["time"], data["event"], aligned)
            rows.append({"model": gene, "cindex": score})
        except Exception:  # noqa: BLE001
            continue
    if not rows:
        raise EngineError("prerequisite", "无法计算任何 C-index")
    out = pd.DataFrame(rows).sort_values("cindex", ascending=False)
    fig, ax = new_fig(7.4, max(3.2, 0.3 * len(out) + 1.4))
    ax.barh(np.arange(len(out))[::-1], out["cindex"].values, color="#4e79a7")
    ax.set_yticks(np.arange(len(out))[::-1])
    ax.set_yticklabels(out["model"], fontsize=8)
    ax.set_xlim(0.4, 1.0)
    ax.axvline(0.5, color="#b9bfc8", linestyle="--", linewidth=0.8)
    ax.set_xlabel("C-index")
    ax.set_title("Concordance index (single gene)")
    figure = save_fig(fig, outdir, "cindex.png", "C-index")
    return result(title="C-index 比较", bullets=[f"{len(out)} 个基因；最高 {out.iloc[0]['model']} C-index={out.iloc[0]['cindex']:.3f}"],
                  metrics={"n_genes": int(len(out))},
                  tables=[table("cindex", "C-index", list(out.columns), out)], figures=[figure])


def run_calibration(ctx, params, outdir):
    payload = _prognosis_payload(ctx, params, outdir, model="cox")
    keep = [fig for fig in payload.get("figures", []) if "calibration" in str(fig.get("file", ""))]
    tables = [t for t in payload.get("tables", []) if t.get("id") == "calibration"]
    payload["figures"] = keep
    payload["tables"] = tables
    payload["summary"]["title"] = "校准曲线"
    return payload


def run_dca(ctx, params, outdir):
    payload = _prognosis_payload(ctx, params, outdir, model="cox")
    keep = [fig for fig in payload.get("figures", []) if "dca" in str(fig.get("file", ""))]
    tables = [t for t in payload.get("tables", []) if t.get("id") == "dca"]
    payload["figures"] = keep
    payload["tables"] = tables
    payload["summary"]["title"] = "决策曲线"
    return payload


def run_nomogram(ctx, params, outdir):
    payload = _prognosis_payload(ctx, params, outdir, model="cox")
    tables = [t for t in payload.get("tables", []) if t.get("id") in ("nomogram", "baseline_survival", "model_coef")]
    payload["tables"] = tables
    payload["figures"] = []
    payload["summary"]["title"] = "列线图数据"
    payload["summary"]["bullets"] = ["列线图数据表：每个变量的系数/HR 与单位分值，以及基线生存 S0(t)；可用 rms::nomogram 直接作图"]
    return payload


def run_time_roc_multi(ctx, params, outdir):
    payload = _prognosis_payload(ctx, params, outdir, model="cox")
    keep = [fig for fig in payload.get("figures", []) if "roc" in str(fig.get("file", ""))]
    payload["figures"] = keep
    payload["summary"]["title"] = "多时点 ROC"
    return payload


STEPS = [
    {
        "id": "prognosis",
        "title": "高级预后",
        "requires": ["survival"],
        "params": [
            {"key": "model", "label": "模型", "type": "select", "default": "cox",
             "options": [{"value": "cox", "label": "Cox"}, {"value": "lasso", "label": "LASSO-Cox"},
                         {"value": "rsf", "label": "RSF-lite"}]},
            {"key": "genes", "label": "建模基因（留空=差异基因 top10）", "type": "genes", "default": []},
            {"key": "cutoff", "label": "风险分组", "type": "select", "default": "median",
             "options": [{"value": "median", "label": "中位数"}, {"value": "quartile", "label": "上四分位"}]},
            {"key": "time_col", "label": "时间列", "type": "column", "default": None},
            {"key": "event_col", "label": "事件列", "type": "column", "default": None},
        ],
        "run": run_prognosis,
    }
]

ANALYSES = [
    {"id": "risk_score_cox", "title": "Cox 风险评分", "group": "高级预后",
     "method": "多因素 Cox（lifelines）风险评分", "needs": ["clinical"], "params": [], "run": run_risk_cox},
    {"id": "risk_score_lasso", "title": "LASSO-Cox 风险评分", "group": "高级预后",
     "method": "L1 惩罚 Cox（lifelines）", "needs": ["clinical"], "params": [], "run": run_risk_lasso},
    {"id": "risk_score_rsf", "title": "RSF-lite 风险评分", "group": "高级预后",
     "method": "bootstrap Cox 集成（sksurv 不可用时的近似）", "needs": ["clinical"], "params": [], "run": run_risk_rsf},
    {"id": "risk_heatmap", "title": "风险基因热图", "group": "高级预后",
     "method": "按风险排序的基因 z-score 热图", "needs": ["clinical"], "params": [], "run": run_risk_heatmap},
    {"id": "cindex", "title": "C-index 比较", "group": "高级预后",
     "method": "单基因 C-index（lifelines concordance_index）", "needs": ["clinical"], "params": [], "run": run_cindex},
    {"id": "calibration", "title": "校准曲线", "group": "高级预后",
     "method": "Cox 预测生存 vs KM 观测（四分箱）", "needs": ["clinical"], "params": [], "run": run_calibration},
    {"id": "dca", "title": "决策曲线（DCA）", "group": "高级预后",
     "method": "净收益 = TP/n - FP/n×pt/(1-pt)", "needs": ["clinical"], "params": [], "run": run_dca},
    {"id": "nomogram", "title": "列线图数据", "group": "高级预后",
     "method": "Cox 系数 + 基线生存 S0(t)", "needs": ["clinical"], "params": [], "run": run_nomogram},
    {"id": "time_roc_multi", "title": "多时点 ROC", "group": "高级预后",
     "method": "按时间点二分类 + 删失排除近似", "needs": ["clinical"], "params": [], "run": run_time_roc_multi},
]
