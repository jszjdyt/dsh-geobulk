"""
dsh-geobulk — 临床生存分析模块（步骤 survival + 生存类分析）

实现：lifelines（KM / log-rank / Cox），时间依赖 ROC 用「按时间点二分类 + 删失排除」近似，
      并在方法名与 notes 里诚实标注（非 IPCW 精确估计）。

注册：
  STEPS   = [survival]
  ANALYSES= km_curve / logrank_multi / cox_uni / cox_multi / cutoff_scan / time_roc /
            subgroup_forest / clinical_table
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

TIME_HINTS = ("os_time", "os.time", "time", "survival_time", "dfs_time", "rfs_time", "pfs_time",
              "month", "months", "days", "follow")
EVENT_HINTS = ("os_event", "event", "status", "dead", "death", "censor", "vital")


# ---------------------------------------------------------------------------
# 临床数据
# ---------------------------------------------------------------------------

def clinical(ctx, params: dict | None = None) -> dict:
    """取出并校验生存时间/事件列，返回 {frame, time_col, event_col, time, event}。"""
    params = params or {}
    info = ctx.get("sample_info")
    if info is None or len(info) == 0:
        raise EngineError("prerequisite", "样本信息为空，无法做生存分析")
    frame = info.copy()
    lower = {str(c).lower(): str(c) for c in frame.columns}
    time_col = params.get("time_col") or None
    event_col = params.get("event_col") or None
    if not time_col:
        for hint in TIME_HINTS:
            match = next((orig for low, orig in lower.items() if low == hint or low.endswith(hint)), None)
            if match:
                time_col = match
                break
    if not event_col:
        for hint in EVENT_HINTS:
            match = next((orig for low, orig in lower.items() if low == hint or low.endswith(hint)), None)
            if match:
                event_col = match
                break
    if not time_col or not event_col:
        raise EngineError(
            "prerequisite",
            "样本信息里找不到生存时间或事件列（可识别 os_time/event 等命名）；请先在「样本」里 join 临床表，或用 params 指定 time_col/event_col",
        )
    times = pd.to_numeric(frame[time_col], errors="coerce")
    events = pd.to_numeric(frame[event_col], errors="coerce")
    # 事件列可能是 "Dead"/"Alive" 之类的字符串
    if events.isna().all():
        text = frame[event_col].astype(str).str.lower().str.strip()
        events = text.isin(["1", "1.0", "dead", "death", "event", "yes", "true", "deceased", "progression"]).astype(float)
    keep = times.notna() & (times > 0)
    frame = frame.loc[keep]
    times = times.loc[keep].astype(float)
    events = events.loc[keep].fillna(0).astype(float)
    if len(frame) < 6:
        raise EngineError("prerequisite", f"有效生存样本只有 {len(frame)} 个，至少需要 6 个")
    n_events = int((events > 0).sum())
    if n_events < 3:
        raise EngineError("prerequisite", f"事件数只有 {n_events} 个，无法做生存分析（建议至少 5 个事件）")
    return {"frame": frame, "time_col": time_col, "event_col": event_col, "time": times.values, "event": events.values}


def align_values(ctx, values, index) -> np.ndarray:
    """把按 ctx.samples 顺序的向量对齐到（可能已过滤的）临床表索引。"""
    series = pd.Series(np.asarray(values, dtype=float), index=[str(s) for s in ctx.get("samples") or []])
    return series.reindex([str(i) for i in index]).values.astype(float)


def _split(values: np.ndarray, labels, cutoff: str, time, event) -> tuple[np.ndarray, float]:
    """按 cutoff 策略分组，返回 (0/1 标签, 阈值)。"""
    values = np.asarray(values, dtype=float)
    if cutoff == "optimal":
        best = None
        for q in np.linspace(0.15, 0.85, 15):
            threshold = float(np.nanquantile(values, q))
            group = (values > threshold).astype(int)
            if group.sum() < 3 or (1 - group).sum() < 3:
                continue
            try:
                from lifelines.statistics import logrank_test

                p = float(logrank_test(time[group == 1], time[group == 0], event[group == 1], event[group == 0]).p_value)
            except Exception:  # noqa: BLE001
                continue
            if best is None or p < best[1]:
                best = (threshold, p)
        if best is None:
            median = float(np.nanmedian(values))
            return (values > median).astype(int), median
        return (values > best[0]).astype(int), best[0]
    if cutoff == "quartile":
        q = float(np.nanquantile(values, 0.75))
        return (values > q).astype(int), q
    median = float(np.nanmedian(values))
    return (values > median).astype(int), median


def _km_plot(time, event, group, outdir: str, title: str, name: str, threshold_note: str = "") -> tuple[dict | None, float | None, pd.DataFrame]:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test

    group = np.asarray(group).astype(int)
    if group.sum() < 2 or (1 - group).sum() < 2:
        return None, None, pd.DataFrame()
    fig, ax = new_fig(7.4, 5.4)
    colors = {"low": "#2e6da4", "high": "#c0392b"}
    rows = []
    for value, label in ((0, "low"), (1, "high")):
        mask = group == value
        if mask.sum() < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(time[mask], event[mask], label=f"{label} (n={int(mask.sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=colors.get(label))
        median = kmf.median_survival_time_
        rows.append({"group": label, "n": int(mask.sum()), "events": int(event[mask].sum()),
                     "median_survival": None if not np.isfinite(median) else float(median)})
    p_value = None
    try:
        p_value = float(multivariate_logrank_test(time, group, event).p_value)
    except Exception:  # noqa: BLE001
        pass
    ax.set_title(f"{title}" + (f"  log-rank P = {p_value:.3g}" if p_value is not None else ""))
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    if threshold_note:
        ax.text(0.02, 0.03, threshold_note, transform=ax.transAxes, fontsize=7, color="#5a6068")
    figure = save_fig(fig, outdir, name, title)
    return figure, p_value, pd.DataFrame(rows)


def _cox_fit(frame: pd.DataFrame, time_col: str, event_col: str, covariates: list[str]):
    from lifelines import CoxPHFitter

    data = frame[[time_col, event_col] + covariates].copy()
    for col in covariates:
        if not pd.api.types.is_numeric_dtype(data[col]):
            codes = pd.Categorical(data[col].astype(str))
            data[col] = codes.codes.astype(float)
    data = data.dropna()
    if len(data) < max(8, len(covariates) * 3):
        raise EngineError("prerequisite", f"样本太少（{len(data)}）无法拟合 {len(covariates)} 个协变量的 Cox 模型")
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(data, duration_col=time_col, event_col=event_col)
    return cph, data


def _cox_table(cph, covariates: list[str]) -> pd.DataFrame:
    summary = cph.summary
    rows = []
    for covariate in covariates:
        if covariate not in summary.index:
            continue
        row = summary.loc[covariate]
        rows.append({
            "gene": covariate,
            "HR": float(np.exp(row["coef"])),
            "lower": float(np.exp(row["coef lower 95%"])),
            "upper": float(np.exp(row["coef upper 95%"])),
            "z": float(row["z"]),
            "p": float(row["p"]),
        })
    return pd.DataFrame(rows)


def _forest_plot(frame: pd.DataFrame, outdir: str, title: str, name: str) -> dict | None:
    if frame is None or frame.empty:
        return None
    frame = frame.head(25).iloc[::-1]
    hr = frame["HR"].astype(float).values
    lo = frame["lower"].astype(float).values
    hi = frame["upper"].astype(float).values
    fig, ax = new_fig(7.6, max(3.4, 0.32 * len(frame) + 1.2))
    y = np.arange(len(frame))
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="o", color="#c0392b", ecolor="#8a919c",
                capsize=3, markersize=4)
    ax.axvline(1.0, color="#2b2f36", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([str(v) for v in frame["gene"]], fontsize=7.5)
    ax.set_xlabel("Hazard ratio (95% CI)")
    ax.set_title(title)
    return save_fig(fig, outdir, name, title)


# ---------------------------------------------------------------------------
# 时间依赖 ROC（删失排除近似）
# ---------------------------------------------------------------------------

def time_roc(time, event, score, horizons) -> list[dict]:
    out = []
    score = np.asarray(score, dtype=float)
    for horizon in horizons:
        horizon = float(horizon)
        case = (time <= horizon) & (event > 0)
        control = time > horizon
        usable = case | control
        if case.sum() < 3 or control.sum() < 3:
            continue
        y = case[usable].astype(int)
        s = score[usable]
        # AUC = P(score_case > score_control)
        pos = s[y == 1]
        neg = s[y == 0]
        try:
            from sklearn.metrics import roc_curve, auc

            fpr, tpr, _ = roc_curve(y, s)
            auc_value = float(auc(y, s))
            points = {"fpr": [round(float(v), 4) for v in fpr], "tpr": [round(float(v), 4) for v in tpr]}
        except Exception:  # noqa: BLE001
            auc_value = float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))
            points = {"fpr": [0, 1], "tpr": [0, 1]}
        out.append({"time": horizon, "auc": auc_value, "n_case": int(case.sum()), "n_control": int(control.sum()), **points})
    return out


def _roc_figure(rocs: list[dict], outdir: str, title: str, name: str) -> dict | None:
    if not rocs:
        return None
    fig, ax = new_fig(6.6, 5.6)
    for i, item in enumerate(rocs):
        fpr = item.get("fpr") or []
        tpr = item.get("tpr") or []
        if not fpr or not tpr:
            continue
        ax.plot(fpr, tpr, color=palette(len(rocs))[i], linewidth=1.5,
                label=f"{item['time']:g} AUC={item['auc']:.3f}")
    ax.plot([0, 1], [0, 1], color="#b9bfc8", linestyle="--", linewidth=0.8)
    ax.set_xlabel("1 - Specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_title(title)
    ax.legend(fontsize=8)
    return save_fig(fig, outdir, name, title)


# ---------------------------------------------------------------------------
# 步骤
# ---------------------------------------------------------------------------

def _candidate_genes(ctx, limit: int = 8) -> list[str]:
    de = ctx.get("de")
    if isinstance(de, dict):
        frame = de.get("table")
        if isinstance(frame, pd.DataFrame) and not frame.empty and "gene" in frame.columns:
            frame = frame.sort_values("adj.P.Val") if "adj.P.Val" in frame.columns else frame
            genes = [str(g) for g in frame["gene"].head(limit * 3)]
            return [g for g in genes if g in set(ctx.get("genes") or [])][:limit]
    return []


def run_survival(ctx, params, outdir):
    params = dict(params or {})
    data = clinical(ctx, params)
    info = data["frame"]
    time = data["time"]
    event = data["event"]
    time_col = data["time_col"]
    event_col = data["event_col"]
    cutoff = str(params.get("cutoff") or "median")
    gene = str(params.get("gene") or "").strip()
    genes = ctx.get("genes") or []
    warnings: list[str] = []
    tables: list[dict] = []
    figures: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []
    metrics: dict = {"time_col": time_col, "event_col": event_col, "cutoff": cutoff, "n_samples": int(len(info))}

    # 分组值
    sample_index = [str(item) for item in ctx.get("samples") or []]
    if gene:
        if gene not in genes:
            raise EngineError("not_found", f"表达矩阵里没有基因 {gene}")
        raw = np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)]
        values = pd.Series(raw, index=sample_index).reindex([str(i) for i in info.index]).values
        group, threshold = _split(values, None, cutoff, time, event)
        title = f"{gene} ({cutoff} split, threshold={threshold:.3g})"
        figure, p_value, km_table = _km_plot(time, event, group, outdir, title, f"km_{gene}.png",
                                             f"threshold={threshold:.3g}")
        if figure:
            figures.append(figure)
        if km_table is not None and not km_table.empty:
            tables.append(table("km_stats", f"{gene} 高低组生存", list(km_table.columns), km_table))
        bullets.append(f"{gene} {cutoff} 分组：log-rank P={p_value:.3g}" if p_value is not None else f"{gene} 分组完成")
        metrics.update({"gene": gene, "threshold": float(threshold), "logrank_p": p_value})
        try:
            cox_frame = pd.DataFrame({time_col: time, event_col: event, gene: values}, index=info.index)
            cph, _ = _cox_fit(cox_frame, time_col, event_col, [gene])
            cox_table = _cox_table(cph, [gene])
            if not cox_table.empty:
                tables.append(table("cox", "单因素 Cox", list(cox_table.columns), cox_table))
                forest = _forest_plot(cox_table, outdir, f"Cox: {gene}", "forest.png")
                if forest:
                    figures.append(forest)
                hr = float(cox_table.iloc[0]["HR"])
                bullets.append(f"Cox HR={hr:.2f}（95%CI {float(cox_table.iloc[0]['lower']):.2f}–{float(cox_table.iloc[0]['upper']):.2f}，p={float(cox_table.iloc[0]['p']):.3g}）")
                metrics.update({"HR": hr, "cox_p": float(cox_table.iloc[0]["p"])})
                path = export_csv(cox_table, outdir, "cox_univariate")
                data_out.append(data_file("survival_csv", path, "单因素 Cox"))
        except EngineError as error:
            warnings.append(f"Cox 未拟合：{error.message}")
    else:
        group_col = params.get("group_col") or ctx.get("group_key")
        if not group_col or group_col not in info.columns:
            # 退化：用第一个候选基因
            candidates = _candidate_genes(ctx, 1)
            if not candidates:
                raise EngineError("prerequisite", "缺少分组列，也没有可用于分组的差异基因；请指定 gene 或 group_col")
            gene = candidates[0]
            raw = np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)]
            aligned = pd.Series(raw, index=sample_index).reindex([str(i) for i in info.index]).values
            info = info.copy()
            median = np.nanmedian(aligned)
            info["_risk_group"] = np.where(aligned > median, "high", "low")
            group_col = "_risk_group"
            warnings.append(f"没有可用分组列，已按 {gene} 中位数分组")
            title = f"{gene} median split"
        else:
            title = f"By {group_col}"
        group_text = info[group_col].astype(str).values
        # 二分类：取前两组的第二组作为 high
        uniq = [g for g in dict.fromkeys(group_text)]
        group = np.array([1 if g == uniq[-1] else 0 for g in group_text]) if len(uniq) >= 2 else np.zeros(len(group_text))
        figure, p_value, km_table = _km_plot(time, event, group, outdir, title, "km_group.png",
                                             f"groups: {uniq[0] if len(uniq) > 1 else '—'} vs {uniq[-1] if uniq else '—'}")
        if figure:
            figures.append(figure)
        if km_table is not None and not km_table.empty:
            tables.append(table("km_stats", "分组生存", list(km_table.columns), km_table))
        bullets.append(f"分组（{group_col}）log-rank P={p_value:.3g}" if p_value is not None else "分组生存分析完成")
        metrics.update({"group_col": group_col, "logrank_p": p_value})

    # Cox 多因素（协变量由参数给，否则跳过）
    covariates = params.get("covariates") or []
    if isinstance(covariates, str):
        covariates = [c.strip() for c in covariates.replace(",", " ").split() if c.strip()]
    covariates = [c for c in covariates if c in info.columns]
    if covariates:
        try:
            cph, _ = _cox_fit(info, time_col, event_col, covariates)
            cox_table = _cox_table(cph, covariates)
            if not cox_table.empty:
                tables.append(table("cox_multi", "多因素 Cox", list(cox_table.columns), cox_table))
                forest = _forest_plot(cox_table, outdir, "Multivariable Cox", "forest_multi.png")
                if forest:
                    figures.append(forest)
                bullets.append(f"多因素 Cox：{len(covariates)} 个协变量，C-index={float(cph.concordance_index_):.3f}")
                metrics.update({"cox_multi_cindex": float(cph.concordance_index_)})
        except EngineError as error:
            warnings.append(f"多因素 Cox 未拟合：{error.message}")

    # 时间依赖 ROC（若有基因/分组）
    rocs = []
    if gene:
        values = pd.Series(np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)],
                           index=sample_index).reindex([str(i) for i in info.index]).values
        horizons = sorted({float(q) for q in np.nanquantile(time, [0.25, 0.5, 0.75]) if q > 0})
        rocs = time_roc(time, event, values, horizons[:3])
        figure = _roc_figure(rocs, outdir, f"Time-dependent ROC ({gene})", "roc_time.png")
        if figure:
            figures.append(figure)
        if rocs:
            tables.append(table("time_roc", "时间依赖 ROC", [
                {"key": "time", "label": "时间", "type": "number"},
                {"key": "auc", "label": "AUC", "type": "number"},
                {"key": "n_case", "label": "事件数", "type": "number"},
                {"key": "n_control", "label": "对照数", "type": "number"},
            ], [{"time": r["time"], "auc": r["auc"], "n_case": r["n_case"], "n_control": r["n_control"]} for r in rocs]))

    if not figures:
        raise EngineError("internal", "生存分析没有生成任何图，请检查数据")

    cindex = None
    if gene:
        try:
            from lifelines.utils import concordance_index

            values = pd.Series(np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)],
                               index=sample_index).reindex([str(i) for i in info.index]).values
            cindex = float(concordance_index(time, -values, event))
            metrics["cindex"] = cindex
        except Exception:  # noqa: BLE001
            pass

    payload = result(
        title="临床生存",
        bullets=bullets,
        metrics=metrics,
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=["时间依赖 ROC 为「按时间点二分类 + 删失排除」近似（非 IPCW 精确估计）",
               f"识别到时间列 {time_col}、事件列 {event_col}"],
        data_out=data_out,
    )
    payload["store"] = {
        "survival": {
            "time_col": time_col,
            "event_col": event_col,
            "cutoff": cutoff,
            "gene": gene or None,
            "logrank_p": metrics.get("logrank_p"),
            "roc": rocs,
            "cindex": cindex,
            "risk_table": [],
        }
    }
    return payload


# ---------------------------------------------------------------------------
# 分析
# ---------------------------------------------------------------------------

def _survival_from_params(ctx, params, outdir, *, gene=None, cutoff=None):
    merged = dict(params or {})
    if gene:
        merged["gene"] = gene
    if cutoff:
        merged["cutoff"] = cutoff
    return run_survival(ctx, merged, outdir)


def run_logrank_multi(ctx, params, outdir):
    data = clinical(ctx, params)
    info = data["frame"]
    time = data["time"]
    event = data["event"]
    group_col = (params or {}).get("group_col") or ctx.get("group_key")
    if not group_col or group_col not in info.columns:
        raise EngineError("prerequisite", "需要分组列（group_col）")
    from lifelines.statistics import multivariate_logrank_test

    group = info[group_col].astype(str).values
    result_test = multivariate_logrank_test(time, group, event)
    rows = [{"group": g, "n": int((group == g).sum()), "events": int(event[group == g].sum())}
            for g in dict.fromkeys(group)]
    frame = pd.DataFrame(rows)
    fig, ax = new_fig(7.4, 5.4)
    from lifelines import KaplanMeierFitter

    for i, g in enumerate(dict.fromkeys(group)):
        mask = group == g
        if mask.sum() < 2:
            continue
        KaplanMeierFitter().fit(time[mask], event[mask], label=f"{g} (n={int(mask.sum())})").plot_survival_function(
            ax=ax, ci_show=False, color=palette(len(set(group)))[i])
    ax.set_title(f"Kaplan-Meier by {group_col}  log-rank P = {float(result_test.p_value):.3g}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    figure = save_fig(fig, outdir, "km_multi.png", "KM multi-group")
    return result(title="多组生存比较", bullets=[f"log-rank P={float(result_test.p_value):.3g}（{len(frame)} 组）"],
                  metrics={"p_value": float(result_test.p_value), "n_groups": int(len(frame))},
                  tables=[table("km_groups", "各组样本与事件", list(frame.columns), frame)],
                  figures=[figure])


def run_cox_uni(ctx, params, outdir):
    data = clinical(ctx, params)
    info = data["frame"]
    time_col, event_col = data["time_col"], data["event_col"]
    gene = str((params or {}).get("gene") or "").strip()
    genes = ctx.get("genes") or []
    if gene and gene in genes:
        candidates = [gene]
    else:
        candidates = _candidate_genes(ctx, int((params or {}).get("top_n") or 10))
    if not candidates:
        raise EngineError("prerequisite", "没有候选基因（请先跑差异表达或提供 gene）")
    X = np.asarray(ctx["X"], dtype=float)
    index = {g: i for i, g in enumerate(genes)}
    rows = []
    warnings = []
    for name in candidates:
        if name not in index:
            continue
        aligned = align_values(ctx, X[:, index[name]], info.index)
        frame = pd.DataFrame({time_col: data["time"], event_col: data["event"], name: aligned}, index=info.index)
        try:
            cph, _ = _cox_fit(frame, time_col, event_col, [name])
            cox = _cox_table(cph, [name])
            if not cox.empty:
                rows.append(cox.iloc[0].to_dict())
        except EngineError as error:
            warnings.append(f"{name}: {error.message}")
    if not rows:
        raise EngineError("prerequisite", "Cox 拟合全部失败（样本或事件太少）")
    out = pd.DataFrame(rows)
    out["p.adjust"] = bh_fdr(out["p"].values)
    out = out.sort_values("p").reset_index(drop=True)
    forest = _forest_plot(out, outdir, "Univariable Cox", "forest_uni.png")
    return result(
        title="单因素 Cox 森林图",
        bullets=[f"{len(out)} 个基因；最显著 {out.iloc[0]['gene']}（HR={out.iloc[0]['HR']:.2f}, p={out.iloc[0]['p']:.3g}）"],
        metrics={"n_genes": int(len(out))},
        tables=[table("cox_uni", "单因素 Cox", list(out.columns), out)],
        figures=[forest] if forest else [],
        warnings=warnings,
        notes=["协变量未校正的单因素 Cox"],
    )


def run_cox_multi(ctx, params, outdir):
    data = clinical(ctx, params)
    covariates = (params or {}).get("covariates") or _candidate_genes(ctx, 5)
    if isinstance(covariates, str):
        covariates = [c.strip() for c in covariates.replace(",", " ").split() if c.strip()]
    genes = ctx.get("genes") or []
    X = np.asarray(ctx["X"], dtype=float)
    frame = data["frame"].copy()
    used = [c for c in covariates if c in genes]
    for name in used:
        frame[name] = align_values(ctx, X[:, genes.index(name)], frame.index)
    if not used:
        raise EngineError("prerequisite", "没有可用协变量（基因名需在表达矩阵里）")
    cph, _ = _cox_fit(frame, data["time_col"], data["event_col"], used)
    cox = _cox_table(cph, used)
    forest = _forest_plot(cox, outdir, "Multivariable Cox", "forest_multi.png")
    return result(
        title="多因素 Cox",
        bullets=[f"{len(used)} 个变量，C-index={float(cph.concordance_index_):.3f}"],
        metrics={"cindex": float(cph.concordance_index_), "n_vars": len(used)},
        tables=[table("cox_multi", "多因素 Cox", list(cox.columns), cox)],
        figures=[forest] if forest else [],
        notes=["协变量默认取差异表达最显著基因（可用 params.covariates 指定）"],
    )


def run_cutoff_scan(ctx, params, outdir):
    data = clinical(ctx, params)
    time, event = data["time"], data["event"]
    gene = str((params or {}).get("gene") or "").strip()
    genes = ctx.get("genes") or []
    if not gene:
        candidates = _candidate_genes(ctx, 1)
        if not candidates:
            raise EngineError("prerequisite", "请提供 gene")
        gene = candidates[0]
    if gene not in genes:
        raise EngineError("not_found", f"表达矩阵里没有基因 {gene}")
    values = align_values(ctx, np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)], data["frame"].index)
    from lifelines.statistics import logrank_test

    rows = []
    for q in np.linspace(0.1, 0.9, 17):
        threshold = float(np.nanquantile(values, q))
        group = values > threshold
        if group.sum() < 3 or (~group).sum() < 3:
            continue
        p = float(logrank_test(time[group], time[~group], event[group], event[~group]).p_value)
        rows.append({"quantile": float(q), "threshold": threshold, "n_high": int(group.sum()),
                     "n_low": int((~group).sum()), "p_value": p, "neg_log10_p": -np.log10(max(p, 1e-300))})
    if not rows:
        raise EngineError("prerequisite", "无法扫描出有效阈值（样本太少）")
    frame = pd.DataFrame(rows)
    best = frame.loc[frame["p_value"].idxmin()]
    fig, ax = new_fig(7.8, 4.8)
    ax.plot(frame["threshold"], frame["neg_log10_p"], color="#4e79a7", marker="o", markersize=3)
    ax.axvline(float(best["threshold"]), color="#c0392b", linestyle="--", linewidth=1,
               label=f"best={float(best['threshold']):.3g}")
    ax.set_xlabel(f"{gene} expression")
    ax.set_ylabel("-log10(log-rank P)")
    ax.set_title("Optimal cutoff scan")
    ax.legend(fontsize=8)
    figure = save_fig(fig, outdir, f"cutoff_{gene}.png", "Cutoff scan")
    return result(
        title="最佳截断值扫描",
        bullets=[f"{gene}：最佳阈值 {float(best['threshold']):.3g}（P={float(best['p_value']):.3g}，高表达 {int(best['n_high'])} 例）"],
        metrics={"gene": gene, "best_threshold": float(best["threshold"]), "best_p": float(best["p_value"])},
        tables=[table("cutoff_scan", "阈值扫描", list(frame.columns), frame)],
        figures=[figure],
        notes=["逐分位点扫描 log-rank 统计量；最优阈值存在选择偏倚，建议外部验证"],
    )


def run_time_roc(ctx, params, outdir):
    data = clinical(ctx, params)
    gene = str((params or {}).get("gene") or "").strip()
    genes = ctx.get("genes") or []
    if not gene:
        candidates = _candidate_genes(ctx, 1)
        if not candidates:
            raise EngineError("prerequisite", "请提供 gene")
        gene = candidates[0]
    if gene not in genes:
        raise EngineError("not_found", f"表达矩阵里没有基因 {gene}")
    values = align_values(ctx, np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)], data["frame"].index)
    horizons = (params or {}).get("horizons")
    if not horizons:
        horizons = [float(v) for v in np.nanquantile(data["time"], [0.25, 0.5, 0.75])]
    elif isinstance(horizons, str):
        horizons = [float(v) for v in horizons.replace(",", " ").split()]
    rocs = time_roc(data["time"], data["event"], values, horizons)
    if not rocs:
        raise EngineError("prerequisite", "没有足够的事件/对照做时间 ROC")
    figure = _roc_figure(rocs, outdir, f"Time-dependent ROC ({gene})", f"roc_{gene}.png")
    return result(
        title="时间依赖 ROC",
        bullets=[f"{gene}：" + "、".join(f"{r['time']:g} 时点 AUC={r['auc']:.3f}" for r in rocs)],
        metrics={"gene": gene, "auc": {f"{r['time']:g}": round(r["auc"], 4) for r in rocs}},
        tables=[table("time_roc", "时间依赖 ROC", [
            {"key": "time", "label": "时间", "type": "number"},
            {"key": "auc", "label": "AUC", "type": "number"},
            {"key": "n_case", "label": "事件数", "type": "number"},
            {"key": "n_control", "label": "对照数", "type": "number"}],
            [{"time": r["time"], "auc": r["auc"], "n_case": r["n_case"], "n_control": r["n_control"]} for r in rocs])],
        figures=[figure] if figure else [],
        notes=["时间 ROC 采用按时间点二分类 + 删失排除近似"],
    )


def run_subgroup_forest(ctx, params, outdir):
    data = clinical(ctx, params)
    gene = str((params or {}).get("gene") or "").strip()
    genes = ctx.get("genes") or []
    if not gene:
        candidates = _candidate_genes(ctx, 1)
        if not candidates:
            raise EngineError("prerequisite", "请提供 gene")
        gene = candidates[0]
    if gene not in genes:
        raise EngineError("not_found", f"表达矩阵里没有基因 {gene}")
    info = data["frame"].copy()
    values = align_values(ctx, np.asarray(ctx["X"], dtype=float)[:, genes.index(gene)], info.index)
    group, _ = _split(values, None, "median", data["time"], data["event"])
    info["_gene_group"] = np.where(group == 1, "high", "low")
    subgroups = (params or {}).get("subgroups") or []
    if isinstance(subgroups, str):
        subgroups = [s.strip() for s in subgroups.replace(",", " ").split() if s.strip()]
    if not subgroups:
        subgroups = [c for c in info.columns if c.lower() in ("sex", "age", "stage", "grade", "batch")][:3]
    rows = []
    for subgroup in subgroups:
        if subgroup not in info.columns:
            continue
        for level in dict.fromkeys(info[subgroup].astype(str)):
            subset = info[info[subgroup].astype(str) == level]
            if len(subset) < 10 or subset["_gene_group"].nunique() < 2:
                continue
            frame = pd.DataFrame({data["time_col"]: pd.to_numeric(subset[data["time_col"]], errors="coerce"),
                                  data["event_col"]: pd.to_numeric(subset[data["event_col"]], errors="coerce").fillna(0),
                                  gene: values[info.index.isin(subset.index)]})
            try:
                cph, _ = _cox_fit(frame, data["time_col"], data["event_col"], [gene])
                cox = _cox_table(cph, [gene])
                if not cox.empty:
                    row = cox.iloc[0].to_dict()
                    row["subgroup"] = f"{subgroup}={level}"
                    rows.append(row)
            except EngineError:
                continue
    if not rows:
        raise EngineError("prerequisite", "没有足够的亚组可做森林图（每组需 ≥10 例）")
    out = pd.DataFrame(rows)
    fig, ax = new_fig(7.8, max(3.4, 0.34 * len(out) + 1.2))
    hr = out["HR"].astype(float).values
    lo = out["lower"].astype(float).values
    hi = out["upper"].astype(float).values
    y = np.arange(len(out))
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="o", color="#c0392b", ecolor="#8a919c", capsize=3, markersize=4)
    ax.axvline(1.0, color="#2b2f36", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(out["subgroup"].tolist(), fontsize=7.5)
    ax.set_xlabel("Hazard ratio (95% CI)")
    ax.set_title(f"Subgroup forest ({gene})")
    figure = save_fig(fig, outdir, f"forest_subgroup_{gene}.png", "Subgroup forest")
    return result(title="亚组森林图", bullets=[f"{gene}：{len(out)} 个亚组"],
                  metrics={"gene": gene, "n_subgroups": int(len(out))},
                  tables=[table("subgroup", "亚组 Cox", list(out.columns), out)],
                  figures=[figure])


def run_clinical_table(ctx, params, outdir):
    data = clinical(ctx, params)
    info = data["frame"]
    group_col = (params or {}).get("group_col") or ctx.get("group_key")
    if not group_col or group_col not in info.columns:
        raise EngineError("prerequisite", "需要分组列（group_col）")
    groups = [g for g in dict.fromkeys(info[group_col].astype(str))]
    rows = []
    for col in info.columns:
        if col == group_col:
            continue
        series = info[col]
        if pd.api.types.is_numeric_dtype(series) and series.nunique() > 8:
            row = {"variable": col, "type": "continuous"}
            for g in groups:
                values = pd.to_numeric(series[info[group_col].astype(str) == g], errors="coerce").dropna()
                row[g] = f"{values.mean():.1f} ± {values.std():.1f}" if len(values) else "—"
            from scipy import stats

            try:
                samples = [pd.to_numeric(series[info[group_col].astype(str) == g], errors="coerce").dropna() for g in groups]
                row["p"] = float(stats.kruskal(*[s for s in samples if len(s)]) .pvalue) if len(groups) > 2 else float(
                    stats.mannwhitneyu(samples[0], samples[1]).pvalue) if all(len(s) for s in samples[:2]) else None
            except Exception:  # noqa: BLE001
                row["p"] = None
            rows.append(row)
        else:
            row = {"variable": col, "type": "categorical"}
            for g in groups:
                counts = series[info[group_col].astype(str) == g].astype(str).value_counts()
                row[g] = "；".join(f"{k}:{v}" for k, v in counts.head(3).items())
            rows.append(row)
    frame = pd.DataFrame(rows)
    return result(title="临床基线表", bullets=[f"{len(frame)} 个变量 × {len(groups)} 组"],
                  metrics={"n_vars": int(len(frame)), "n_groups": len(groups)},
                  tables=[table("clinical", "基线特征", list(frame.columns), frame)],
                  figures=[])


STEPS = [
    {
        "id": "survival",
        "title": "临床生存",
        "requires": ["normalize"],
        "params": [
            {"key": "time_col", "label": "生存时间列", "type": "column", "default": None},
            {"key": "event_col", "label": "事件列", "type": "column", "default": None},
            {"key": "gene", "label": "分析基因（可选）", "type": "text", "default": ""},
            {"key": "group_col", "label": "分组列（可选）", "type": "column", "default": None},
            {"key": "cutoff", "label": "分组方式", "type": "select", "default": "median",
             "options": [{"value": "median", "label": "中位数"}, {"value": "optimal", "label": "最佳截断值"}, {"value": "quartile", "label": "上四分位"}]},
            {"key": "covariates", "label": "多因素协变量", "type": "genes", "default": []},
        ],
        "run": run_survival,
    }
]

ANALYSES = [
    {"id": "km_curve", "title": "KM 生存曲线", "group": "临床生存",
     "method": "Kaplan-Meier + log-rank（lifelines）", "needs": ["clinical"], "params": [], "run": run_survival},
    {"id": "logrank_multi", "title": "多组 log-rank", "group": "临床生存",
     "method": "多组 log-rank（lifelines）", "needs": ["clinical", "groups"], "params": [], "run": run_logrank_multi},
    {"id": "cox_uni", "title": "单因素 Cox 森林图", "group": "临床生存",
     "method": "单因素 Cox（lifelines）+ BH", "needs": ["clinical"], "params": [], "run": run_cox_uni},
    {"id": "cox_multi", "title": "多因素 Cox", "group": "临床生存",
     "method": "多因素 Cox（lifelines）", "needs": ["clinical"], "params": [], "run": run_cox_multi},
    {"id": "cutoff_scan", "title": "最佳截断值扫描", "group": "临床生存",
     "method": "逐分位点 log-rank 扫描", "needs": ["clinical"], "params": [], "run": run_cutoff_scan},
    {"id": "time_roc", "title": "时间依赖 ROC", "group": "临床生存",
     "method": "按时间点二分类 + 删失排除近似", "needs": ["clinical"], "params": [], "run": run_time_roc},
    {"id": "subgroup_forest", "title": "亚组森林图", "group": "临床生存",
     "method": "各亚组 Cox HR", "needs": ["clinical"], "params": [], "run": run_subgroup_forest},
    {"id": "clinical_table", "title": "临床基线表", "group": "临床生存",
     "method": "连续变量 mean±sd / 分类变量频数", "needs": ["clinical", "groups"], "params": [], "run": run_clinical_table},
]
