#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 生存/预后回归测试（任务单 GEO-04 / GEO-05 / GEO-06）

  GEO-04  事件列逐值清洗：混合编码/TCGA 风格/未知文本/缺失/非 0/1 数字都要有明确行为，
          未知状态不得当删失；清洗后的 time/event 必须写回 frame，让 KM/Cox/预后同源。
  GEO-05  校准/DCA/列线图概率全部来自同一个已拟合 Cox；基线取自 baseline_survival_；
          缺必需列时明确失败，不合成替代数值；DCA 用 KM 估计而不是把删失当阴性。
  GEO-06  结果必须标注「训练集表观性能」；RSF-lite 的 1/2 基因边界不再报错。

用法：cd plugins/dsh-geobulk && python3 test/py_prognosis_regression_test.py
"""

from __future__ import annotations

import os
import sys
import tempfile
import traceback
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_prognosis as prog  # noqa: E402
import bg_survival as surv  # noqa: E402
from bg_common import EngineError  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-prog-reg-")
PASSED: list[str] = []
FAILED: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


# ---------------------------------------------------------------------------
# GEO-04：事件清洗
# ---------------------------------------------------------------------------

def test_geo04_event_cleaning():
    print("\n[GEO-04] 事件列逐值清洗")
    samples = [f"S{i}" for i in range(8)]
    info = pd.DataFrame({
        "os_time": [10, 20, 30, 40, 50, 60, 70, 80],
        "event": [1, 0, 1, np.nan, 1, "Dead", "Alive", np.nan],
    }, index=samples)
    ctx = {"sample_info": info, "X": np.zeros((8, 1)), "genes": ["G0"], "samples": samples}
    data = surv.clinical(ctx)
    frame = data["frame"]
    check("GEO-04a 未知/缺失事件被排除而不是填 0",
          len(frame) == 6, f"kept={len(frame)}")
    check("GEO-04b Dead 被识别为事件、Alive 为删失",
          list(frame["event"].values) == [1.0, 0.0, 1.0, 1.0, 1.0, 0.0],
          str(list(frame["event"].values)))
    check("GEO-04c 发过警告说明排除了多少条",
          any("无法识别的事件状态" in w for w in data["warnings"]), str(data["warnings"])[:80])
    check("GEO-04d 清洗后的 event 写回返回 frame",
          frame["event"].dtype.kind == "f" and set(np.unique(frame["event"])) <= {0.0, 1.0})
    check("GEO-04e 清洗后的 time 写回返回 frame 且为数值",
          np.issubdtype(frame["os_time"].dtype, np.number))

    # TCGA 风格编码
    for value, expected in (("0:LIVING", 0.0), ("1:DECEASED", 1.0), ("1:deceased", 1.0), ("0:alive", 0.0)):
        check(f"GEO-04f TCGA 风格 {value} → {int(expected)}",
              surv._parse_event_value(value) == expected, str(surv._parse_event_value(value)))
    # 未知文本 / 非 0/1 数字 / 空值
    for value in ("Unknown", "not reported", 2, 3.5, np.nan, None, ""):
        check(f"GEO-04g 非确定值 {value!r} → None（不可当删失）",
              surv._parse_event_value(value) is None, str(surv._parse_event_value(value)))

    # 单一 0/1 数字编码、字符串 0/1 都必须正确
    for value, expected in ((1, 1.0), (0, 0.0), ("1", 1.0), ("0", 0.0), ("Dead", 1.0), ("alive", 0.0)):
        check(f"GEO-04h 常规编码 {value!r} → {int(expected)}",
              surv._parse_event_value(value) == expected)

    # KM/Cox/预后 用同一批样本
    rng = np.random.default_rng(5)
    n = 40
    s = [f"P{i:02d}" for i in range(n)]
    time = rng.uniform(5, 60, n)
    status = rng.integers(0, 2, n)
    expr = rng.normal(7, 1, (n, 6))
    genes = [f"G{i}" for i in range(6)]
    cctx = {
        "sample_info": pd.DataFrame({"os_time": time, "os_event": status}, index=s),
        "X": expr, "genes": genes, "samples": s, "results": {}, "group_key": "os_event",
    }
    cox = surv.run_cox_uni(cctx, {"gene": "G0"}, os.path.join(OUTDIR, "g04-cox"))
    check("GEO-04i Cox 分析披露 n_samples_used / n_events",
          cox["metrics"].get("n_samples_used") == n and cox["metrics"].get("n_events") == int(status.sum()),
          f"{cox['metrics'].get('n_samples_used')}/{cox['metrics'].get('n_events')}")


# ---------------------------------------------------------------------------
# GEO-05：概率与删失
# ---------------------------------------------------------------------------

def _survival_fixture(n=80, seed=9, censor_early=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(0, 1, n)
    z = rng.normal(0, 1, n)
    base = rng.exponential(20, n) + 1
    hazard = np.exp(0.8 * x - 0.5 * z)
    t_event = base / hazard
    if censor_early:
        censor_t = rng.uniform(0.5, 8.0, n)   # 大量早期删失
    else:
        censor_t = rng.uniform(30, 120, n)
    obs_t = np.minimum(t_event, censor_t)
    event = (t_event <= censor_t).astype(int)
    frame = pd.DataFrame({"time": obs_t, "event": event, "X": x, "Z": z})
    return frame


def test_geo05_probabilities():
    print("\n[GEO-05] 校准 / 基线 / 列线图 / DCA")
    from lifelines import CoxPHFitter
    from lifelines import KaplanMeierFitter

    frame = _survival_fixture()
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(frame, duration_col="time", event_col="event")

    # 基线必须逐值等于拟合模型的 baseline_survival_
    timeline, baseline = prog._baseline_survival(cph, frame, "time")
    ref = cph.baseline_survival_
    check("GEO-05a 基线生存逐值等于 baseline_survival_",
          np.allclose(timeline, np.asarray(ref.index, float)) and
          np.allclose(baseline, np.asarray(ref.iloc[:, 0], float)),
          f"max|diff|={np.max(np.abs(baseline - np.asarray(ref.iloc[:, 0], float))):.2e}")

    # 校准曲线的 predicted 必须等于同一模型 predict_survival_function 的分箱均值
    horizons = [10.0, 20.0]
    calib = prog.calibration_curve(cph, frame, "time", "event", horizons)
    ok = True
    detail = ""
    for horizon in horizons:
        direct = np.asarray(cph.predict_survival_function(frame, times=[horizon]).values[0], float)
        sub = calib[calib["horizon"] == horizon]
        for _, row in sub.iterrows():
            mask = pd.qcut(pd.Series(direct).rank(method="first"), q=4, labels=False,
                           duplicates="drop") == row["bin"]
            expected = float(np.mean(direct[np.asarray(mask)]))
            if abs(expected - row["predicted"]) > 1e-9:
                ok = False
                detail = f"bin={row['bin']} 预期 {expected:.6f} 实际 {row['predicted']:.6f}"
    check("GEO-05b 校准概率来自同一已拟合模型（无重复施加协变量效应）", ok, detail)

    # 缺必需列 → 明确失败，不合成数值
    missing = frame.drop(columns=["Z"]).copy()
    raised = None
    try:
        prog._predict_survival(cph, missing, [10.0])
    except EngineError as error:
        raised = error
    check("GEO-05c 缺个体预测必需列 → prerequisite",
          raised is not None and raised.code == "prerequisite",
          f"{getattr(raised, 'message', '')[:60]}")
    raised = None
    try:
        prog.calibration_curve(cph, missing, "time", "event", [10.0])
    except EngineError as error:
        raised = error
    check("GEO-05d 校准缺列时不返回替代数值", raised is not None and raised.code == "prerequisite")
    raised = None
    try:
        prog.decision_curve(cph, missing, "time", "event", 10.0)
    except EngineError as error:
        raised = error
    check("GEO-05e DCA 缺列时不返回替代数值", raised is not None and raised.code == "prerequisite")

    # DCA：KM 版净获益 vs 旧「把时点前删失当阴性」的朴素做法必须不同
    censored = _survival_fixture(n=120, seed=17, censor_early=True)
    cph2 = CoxPHFitter(penalizer=0.01)
    cph2.fit(censored, duration_col="time", event_col="event")
    horizon = float(np.nanmedian(censored["time"]))
    dca = prog.decision_curve(cph2, censored, "time", "event", horizon)
    check("GEO-05f DCA 有结果且包含 KM 风险列",
          not dca.empty and {"net_benefit", "treat_all", "event_risk_treat"} <= set(dca.columns),
          f"rows={len(dca)}")

    # 手算 treat_all 对照：prevalence 用 KM 时点事件率，而不是「时点前是否观察到事件」
    km_all = float(1.0 - KaplanMeierFitter().fit(censored["time"], censored["event"]).predict(horizon))
    naive = float(((censored["time"] <= horizon) & (censored["event"] > 0)).mean())
    threshold = float(dca["threshold"].iloc[len(dca) // 2])
    row = dca.iloc[len(dca) // 2]
    odds = threshold / (1 - threshold)
    check("GEO-05g treat_all 用 KM 时点事件率",
          abs(float(row["treat_all"]) - (km_all - (1 - km_all) * odds)) < 1e-6,
          f"KM={km_all:.3f} naive={naive:.3f}")
    naive_nb = naive - (1 - naive) * odds
    check("GEO-05h 与旧「删失当阴性」算法有实质差别",
          abs(naive_nb - float(row["net_benefit"])) > 1e-3,
          f"KM-NB={float(row['net_benefit']):.4f} naive-NB={naive_nb:.4f}")

    # 缺基线属性 → 明确失败（不合成 exp(-0.05t)）
    class NoBaseline:
        pass

    raised = None
    try:
        prog._baseline_survival(NoBaseline(), frame, "time")
    except EngineError as error:
        raised = error
    check("GEO-05i 模型没有 baseline_survival_ → internal（不合成假曲线）",
          raised is not None and raised.code == "internal")


def _prognosis_ctx(n=60, n_genes=8, seed=31):
    rng = np.random.default_rng(seed)
    s = [f"P{i:02d}" for i in range(n)]
    x = rng.normal(0, 1, n)
    genes = [f"G{i}" for i in range(n_genes)]
    expr = rng.normal(7, 1, (n, n_genes))
    expr[:, 0] += 2.0 * x
    t_event = rng.exponential(20, n) / np.exp(0.9 * x)
    censor = rng.uniform(25, 120, n)
    obs = np.minimum(t_event, censor)
    event = (t_event <= censor).astype(int)
    info = pd.DataFrame({"os_time": obs, "os_event": event}, index=s)
    de = pd.DataFrame({
        "gene": genes, "logFC": np.r_[2.0, np.zeros(n_genes - 1)],
        "AveExpr": expr.mean(0), "t": np.r_[4.0, np.zeros(n_genes - 1)],
        "P.Value": np.r_[1e-5, np.full(n_genes - 1, 0.5)],
        "adj.P.Val": np.r_[1e-4, np.full(n_genes - 1, 0.6)],
        "B": np.zeros(n_genes),
    })
    de["regulation"] = np.where(de["adj.P.Val"] <= 0.05, "up", "ns")
    return {
        "samples": s, "genes": genes, "X": expr, "X_kind": "lognorm",
        "sample_info": info, "results": {}, "de": {"table": de},
    }


def test_geo06_prognosis():
    print("\n[GEO-06] 表观性能标注 / RSF-lite 边界")
    ctx = _prognosis_ctx()
    out = prog.run_prognosis(ctx, {}, os.path.join(OUTDIR, "g06-cox"))
    metrics = out["summary"]["metrics"]
    check("GEO-06a 明确标注训练集表观性能",
          "表观" in str(metrics.get("evaluation_scope") or "") or
          "apparent" in str(metrics.get("evaluation_scope") or ""),
          str(metrics.get("evaluation_scope"))[:60])
    check("GEO-06b 校准确性指标带上 evaluation_scope",
          "cindex_scope" in metrics)
    check("GEO-06c 基线生存表存在（来自 baseline_survival_）",
          any(t.get("id") == "baseline_survival" for t in out["tables"]),
          ",".join(t.get("id", "") for t in out["tables"]))
    check("GEO-06d 基线生存不是 exp(-0.05t) 假曲线",
          not _looks_like_fake_baseline(out))

    # 1 / 2 个基因时 RSF-lite 不应因抽样大小报错
    for count in (1, 2, 3):
        sub = _prognosis_ctx(n=60, n_genes=count, seed=41)
        raised = None
        try:
            prog.run_risk_rsf(sub, {"n_boot": 20}, os.path.join(OUTDIR, f"g06-rsf{count}"))
        except EngineError as error:
            raised = error
        check(f"GEO-06e {count} 个基因时集成风险评分可跑", raised is None,
              f"{getattr(raised, 'message', '')[:70]}")
    out_rsf = prog.run_risk_rsf(ctx, {"n_boot": 30}, os.path.join(OUTDIR, "g06-rsf"))
    rsf_notes = " ".join(out_rsf.get("notes") or [])
    check("GEO-06f 不再冒用 RSF 名称（明确写「不是随机生存森林」）",
          "不是随机生存森林" in rsf_notes or "随机生存森林" in rsf_notes, rsf_notes[:80])
    check("GEO-06g 模型类型记录为子集 bootstrap Cox 集成",
          out_rsf["summary"]["metrics"].get("model_kind") == "subset_bootstrap_cox_ensemble",
          str(out_rsf["summary"]["metrics"].get("model_kind")))


def _looks_like_fake_baseline(payload: dict) -> bool:
    for table_obj in payload.get("tables") or []:
        if table_obj.get("id") != "baseline_survival":
            continue
        rows = table_obj.get("rows") or []
        if not rows:
            return False
        times = np.array([row["time"] for row in rows], float)
        values = np.array([row["S0"] for row in rows], float)
        if np.allclose(values, np.exp(-0.05 * times), atol=1e-12):
            return True
    return False


def main() -> int:
    print("== dsh-geobulk 生存/预后回归（GEO-04/05/06）==")
    for fn in (test_geo04_event_cleaning, test_geo05_probabilities, test_geo06_prognosis):
        try:
            fn()
        except Exception:  # noqa: BLE001
            FAILED.append(fn.__name__)
            print(f"  FAIL {fn.__name__} 抛异常")
            traceback.print_exc()
    print(f"\n== 结果：{len(PASSED)} 通过 / {len(FAILED)} 失败 ==")
    if FAILED:
        for item in FAILED:
            print(f"  ✗ {item}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
