#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 差异分析回归测试（任务单 GEO-01 / GEO-02 / GEO-03）

这三个问题是「能跑通但结论错」的类型，原有测试全部通过也发现不了，所以这里
全部按「数值/方向/数据对应关系」断言，而不是只看 ok：

  GEO-01  自动推断的 case/control 方向：两种样本出现顺序都必须推断一致，
          且与显式指定时的 logFC 同号同值。
  GEO-02  DESeq2 结果不得被清空；logFC/p 必须与直接调用 DESeq2::results 一致；
          程序错误（exit 4）不得静默换成 Welch。
  GEO-03  配对必须真的进设计（强患者效应下配对与非配对结果不同）；数值协变量
          必须按 numeric 进入 model.matrix（不能变每个年龄一个哑变量）。

用法：cd plugins/dsh-geobulk && python3 test/py_de_regression_test.py
成功打印 PASS 并以 0 退出；失败打印 FAIL 并以 1 退出。
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_de  # noqa: E402
import rbridge  # noqa: E402
from bg_common import EngineError  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-de-reg-")
PASSED: list[str] = []
FAILED: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def build_ctx(X, genes, samples, info, *, kind="lognorm", counts=None, groups=None) -> dict:
    return {
        "dataset_id": "ds_de",
        "name": "synthetic",
        "genes": list(genes),
        "samples": list(samples),
        "X": np.asarray(X, dtype=np.float64),
        "X_kind": kind,
        "counts": None if counts is None else np.asarray(counts, dtype=np.float64),
        "sample_info": info,
        "group_key": "group",
        "groups": groups or [],
        "species": "human",
        "results": {},
        "modules": {},
    }


def run(ctx, params, tag):
    out = os.path.join(OUTDIR, tag)
    os.makedirs(out, exist_ok=True)
    return bg_de.run_de(ctx, params, out)


# ---------------------------------------------------------------------------
# GEO-01：自动比较方向
# ---------------------------------------------------------------------------

def test_geo01_direction():
    print("\n[GEO-01] 自动实验/对照方向")
    guess = bg_de._infer_pair("normal", "tumor")
    check("GEO-01a 关键字推断返回带字段名的 PairGuess",
          isinstance(guess, bg_de.PairGuess) and guess.case == "tumor" and guess.control == "normal",
          f"{guess}")

    rng = np.random.default_rng(1)
    n_genes = 40
    genes = [f"G{i}" for i in range(n_genes)]
    X = rng.normal(7.0, 0.8, (8, n_genes))
    X[:4, :10] += 3.0  # tumor 前 4 个样本里前 10 个基因更高
    labels = ["tumor"] * 4 + ["normal"] * 4
    dirs = {}
    for order_name, samples in (
        ("tumor_first", [f"S{i}" for i in range(8)]),
        ("normal_first", [f"N{i}" for i in range(8)]),
    ):
        info = pd.DataFrame({"group": labels}, index=samples)
        ctx = build_ctx(X, genes, samples, info)
        out = run(ctx, {"method": "welch", "group_col": "group"}, f"geo01-{order_name}")
        store = out["de_store"]
        dirs[order_name] = (store["case"], store["control"], out["summary"]["metrics"]["contrast"])
        top = store["table"].set_index("gene")
        check(f"GEO-01b {order_name} 自动推断 case=tumor",
              store["case"] == "tumor", f"case={store['case']} control={store['control']}")
        check(f"GEO-01c {order_name} tumor 中高表达基因 logFC 为正",
              float(top.loc["G0", "logFC"]) > 0, f"logFC={float(top.loc['G0', 'logFC']):.3f}")

    # 显式 vs 自动：必须同号同值
    samples = [f"S{i}" for i in range(8)]
    info = pd.DataFrame({"group": labels}, index=samples)
    ctx = build_ctx(X, genes, samples, info)
    explicit = run(ctx, {"method": "welch", "group_col": "group", "case": "tumor", "control": "normal"}, "geo01-explicit")
    auto = run(ctx, {"method": "welch", "group_col": "group"}, "geo01-auto")
    a = explicit["de_store"]["table"].set_index("gene")["logFC"]
    b = auto["de_store"]["table"].set_index("gene")["logFC"]
    diff = float(np.max(np.abs(a.values - b.values)))
    check("GEO-01d 显式与自动的 logFC 逐基因一致", diff < 1e-9, f"max|diff|={diff:.2e}")
    check("GEO-01e 显式对比方向写进 metrics",
          explicit["summary"]["metrics"]["contrast"] == "tumor-normal",
          str(explicit["summary"]["metrics"].get("contrast")))

    # treated/control、resistant/sensitive 也要一致
    for a_name, b_name in (("treated", "control"), ("resistant", "sensitive")):
        g1 = bg_de._infer_pair(a_name, b_name)
        g2 = bg_de._infer_pair(b_name, a_name)
        ok = g1 is not None and g2 is not None and g1.case == g2.case and g1.control == g2.control
        check(f"GEO-01f {a_name}/{b_name} 两种出现顺序推断一致", ok, f"{g1} vs {g2}")


# ---------------------------------------------------------------------------
# GEO-02：DESeq2 结果列名
# ---------------------------------------------------------------------------

def _counts_fixture(n_case=6, n_ctrl=6, n_genes=60, seed=7, de_genes=8):
    rng = np.random.default_rng(seed)
    base = rng.uniform(60, 400, n_genes)
    lib_case = rng.uniform(0.8, 1.2, n_case)
    lib_ctrl = rng.uniform(0.8, 1.2, n_ctrl)
    X = np.zeros((n_case + n_ctrl, n_genes))
    for i in range(n_case):
        mu = base * lib_case[i]
        mu[:de_genes] *= 3.0
        X[i] = rng.negative_binomial(2.0, 2.0 / (2.0 + mu))
    for j in range(n_ctrl):
        mu = base * lib_ctrl[j]
        X[n_case + j] = rng.negative_binomial(2.0, 2.0 / (2.0 + mu))
    genes = [f"G{i:03d}" for i in range(n_genes)]
    samples = [f"T{i}" for i in range(n_case)] + [f"C{j}" for j in range(n_ctrl)]
    labels = ["tumor"] * n_case + ["normal"] * n_ctrl
    return X, genes, samples, labels


def _direct_deseq2(X, genes, samples, labels, case="tumor", control="normal") -> pd.DataFrame:
    """直接调 DESeq2::results，作为官方基准。"""
    tmp = Path(tempfile.mkdtemp(prefix="deseq2-ref-"))
    mat = pd.DataFrame(X.T, index=genes, columns=samples)
    mat.to_csv(tmp / "counts.csv")
    pd.DataFrame({"group": labels}, index=samples).to_csv(tmp / "info.csv")
    script = tmp / "ref.R"
    script.write_text(
        """
args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
counts <- as.matrix(read.csv(file.path(d, "counts.csv"), row.names=1, check.names=FALSE))
info <- read.csv(file.path(d, "info.csv"), row.names=1, check.names=FALSE)
info$group <- factor(info$group, levels=c("%s", "%s"))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=info, design=~group)
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds, contrast=c("group", "%s", "%s"))
out <- data.frame(gene=rownames(res), logFC=res$log2FoldChange,
                  P.Value=res$pvalue, adj.P.Val=res$padj)
write.csv(out, file.path(d, "ref.csv"), row.names=FALSE)
""" % (control, case, case, control),
        encoding="utf-8",
    )
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(script), str(tmp)],
        capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"参考 DESeq2 调用失败：{proc.stderr[-1500:]}")
    return pd.read_csv(tmp / "ref.csv")


def test_geo02_deseq2():
    print("\n[GEO-02] DESeq2 结果列名 / 静默回退")
    if not (rbridge.r_available() and rbridge.has_pkg("DESeq2")):
        print("  SKIP 本机没有 R + DESeq2")
        return
    X, genes, samples, labels = _counts_fixture()
    info = pd.DataFrame({"group": labels}, index=samples)
    ctx = build_ctx(X, genes, samples, info, kind="counts", counts=X)
    out = run(ctx, {"method": "deseq2", "group_col": "group"}, "geo02-deseq2")
    frame = out["de_store"]["table"]
    check("GEO-02a DESeq2 结果表非空（修复前为 0 行）", len(frame) >= len(genes) * 0.9,
          f"rows={len(frame)} / {len(genes)}")
    check("GEO-02b actual_method=deseq2（没有静默回退）",
          out["summary"]["metrics"].get("actual_method") == "deseq2",
          str(out["summary"]["metrics"].get("actual_method")))
    check("GEO-02c lite=False", not out["de_store"]["lite"])

    ref = _direct_deseq2(X, genes, samples, labels)
    merged = frame.set_index("gene").join(ref.set_index("gene"), rsuffix="_ref")
    ok_rows = merged.dropna(subset=["logFC", "logFC_ref"])
    d_lfc = float(np.max(np.abs(ok_rows["logFC"] - ok_rows["logFC_ref"])))
    ok_p = ok_rows.dropna(subset=["P.Value", "P.Value_ref"])
    d_p = float(np.nanmax(np.abs(np.log10(np.clip(ok_p["P.Value"], 1e-300, 1)) -
                                 np.log10(np.clip(ok_p["P.Value_ref"], 1e-300, 1)))))
    check("GEO-02d logFC 与直接调用 DESeq2 一致", d_lfc < 1e-6, f"max|diff|={d_lfc:.2e}")
    check("GEO-02e p 值与直接调用 DESeq2 一致（log10 差 <0.01）", d_p < 1e-2, f"max dlog10p={d_p:.3e}")
    check("GEO-02f 方向：tumor 上调的基因 logFC > 0",
          float(frame.set_index("gene").loc["G000", "logFC"]) > 0,
          f"logFC={float(frame.set_index('gene').loc['G000', 'logFC']):.3f}")

    # 程序错误不得静默回退
    original = bg_de._r_diff
    try:
        bg_de._r_diff = lambda *a, **k: {"frame": None, "failed": True,
                                         "error": "R 计算内部错误：模拟", "code": "runtime"}
        raised = None
        try:
            run(ctx, {"method": "deseq2", "group_col": "group"}, "geo02-runtime")
        except EngineError as error:
            raised = error
        check("GEO-02g R 程序错误 → 明确失败，不换方法",
              raised is not None and raised.code == "runtime",
              f"code={getattr(raised, 'code', None)} msg={getattr(raised, 'message', '')[:60]}")
        bg_de._r_diff = original

        # 环境原因（无 R）仍允许回退，但必须标出 requested vs actual
        bg_de._r_diff = lambda *a, **k: {"frame": None, "failed": True,
                                         "error": "本机没有可用的 Rscript", "code": "no_r"}
        out2 = run(ctx, {"method": "deseq2", "group_col": "group"}, "geo02-nor")
        m = out2["summary"]["metrics"]
        check("GEO-02h 环境不可用可回退且标出 requested/actual",
              m.get("method_requested") == "deseq2" and m.get("actual_method") == "welch" and bool(m.get("lite")),
              f"requested={m.get('method_requested')} actual={m.get('actual_method')}")
    finally:
        bg_de._r_diff = original


# ---------------------------------------------------------------------------
# GEO-03：配对与协变量进设计
# ---------------------------------------------------------------------------

def _paired_fixture(n_pairs=8, n_genes=80, effect=2.0, seed=3, patient_sd=3.0):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    patient = rng.normal(0, patient_sd, (n_pairs, 1))
    effect_vec = np.zeros(n_genes)
    effect_vec[:10] = effect
    X, samples, labels, pairs = [], [], [], []
    for p in range(n_pairs):
        for arm, lab in ((0, "control"), (1, "case")):
            base = 7.0 + patient[p] + rng.normal(0, 0.3, n_genes)
            X.append(base + arm * effect_vec + rng.normal(0, 0.5, n_genes))
            samples.append(f"P{p:02d}{"T" if lab == "case" else "C"}")
            labels.append(lab)
            pairs.append(f"P{p:02d}")
    info = pd.DataFrame({"group": labels, "patient": pairs}, index=samples)
    return np.vstack(X), genes, samples, info


def test_geo03_design():
    print("\n[GEO-03] 配对与协变量进入设计矩阵")
    if not (rbridge.r_available() and rbridge.has_pkg("limma")):
        print("  SKIP 本机没有 R + limma")
        return
    X, genes, samples, info = _paired_fixture()
    ctx = build_ctx(X, genes, samples, info)

    paired = run(ctx, {"method": "limma", "group_col": "group", "paired": True, "pair_col": "patient"},
                 "geo03-paired")
    unpaired = run(ctx, {"method": "limma", "group_col": "group"}, "geo03-unpaired")
    pf = paired["de_store"]["table"].set_index("gene")
    uf = unpaired["de_store"]["table"].set_index("gene")
    check("GEO-03a 配对结果非空", len(pf) == len(genes), f"rows={len(pf)}")
    check("GEO-03b 配对 logFC ≈ 组内差（效应 2.0）",
          abs(float(pf.loc["G000", "logFC"]) - 2.0) < 0.6, f"logFC={float(pf.loc['G000', 'logFC']):.3f}")
    same_p = np.allclose(pf["P.Value"].values, uf["P.Value"].values, equal_nan=True)
    check("GEO-03c 配对 ≠ 非配对（患者效应进模型）", not same_p,
          f"paired p[G000]={float(pf.loc['G000', 'P.Value']):.2e} unpaired={float(uf.loc['G000', 'P.Value']):.2e}")
    check("GEO-03d 配对显著更强（强患者效应下）",
          float(pf.loc["G000", "P.Value"]) < float(uf.loc["G000", "P.Value"]))
    check("GEO-03e metrics 记录设计矩阵含 pair",
          "pair" in str(paired["summary"]["metrics"].get("design") or ""),
          str(paired["summary"]["metrics"].get("design")))
    check("GEO-03f metrics 记录配对数",
          int(paired["summary"]["metrics"].get("n_pairs") or 0) == 8,
          str(paired["summary"]["metrics"].get("n_pairs")))

    # 配对列缺失 → 必须报错，不允许退成独立样本
    bad_info = info.copy()
    bad_info.loc["P00C", "patient"] = np.nan
    ctx_bad = build_ctx(X, genes, samples, bad_info)
    raised = None
    try:
        run(ctx_bad, {"method": "limma", "group_col": "group", "paired": True, "pair_col": "patient"}, "geo03-missing")
    except EngineError as error:
        raised = error
    check("GEO-03g 配对 ID 缺失 → prerequisite", raised is not None and raised.code == "prerequisite",
          f"{getattr(raised, 'message', '')[:60]}")

    # 不成对（每组各 4 个不同患者）→ 也必须报错
    solo_samples = [f"S{i}" for i in range(8)]
    solo_info = pd.DataFrame({"group": ["case"] * 4 + ["control"] * 4,
                              "patient": [f"X{i}" for i in range(8)]}, index=solo_samples)
    ctx_solo = build_ctx(X[:8], genes, solo_samples, solo_info)
    raised = None
    try:
        run(ctx_solo, {"method": "limma", "group_col": "group", "paired": True, "pair_col": "patient"}, "geo03-solo")
    except EngineError as error:
        raised = error
    check("GEO-03h 没有共同配对 → prerequisite", raised is not None and raised.code == "prerequisite",
          f"{getattr(raised, 'message', '')[:60]}")

    # 连续协变量（age）：修复前会报 object 'cov_age' not found
    rng = np.random.default_rng(11)
    n = 16
    age = np.linspace(30, 75, n)
    groups = ["case"] * (n // 2) + ["control"] * (n - n // 2)
    G = 60
    g = [f"G{i:03d}" for i in range(G)]
    age_effect = np.zeros(G)
    age_effect[:15] = 0.05
    expr = 6.0 + np.outer(age, age_effect) + rng.normal(0, 0.4, (n, G))
    s = [f"A{i:02d}" for i in range(n)]
    ainfo = pd.DataFrame({"group": groups, "age": age, "sex": ["M" if i % 2 else "F" for i in range(n)]}, index=s)
    actx = build_ctx(expr, g, s, ainfo)
    raised = None
    adjusted = None
    try:
        adjusted = run(actx, {"method": "limma", "group_col": "group", "covars": ["age"]}, "geo03-age")
    except EngineError as error:
        raised = error
    check("GEO-03i 数值协变量 age 可跑（修复前 object 'cov_age' not found）",
          raised is None and adjusted is not None,
          f"{getattr(raised, 'message', '')[:80]}")
    if adjusted is not None:
        design = str(adjusted["summary"]["metrics"].get("design") or "")
        check("GEO-03j 设计里 cov_age 只有 1 列（numeric，不是每个年龄一个哑变量）",
              "cov_age" in design, design)
        # 与无协变量结果不同 ⇒ 协变量真的进了模型
        plain = run(actx, {"method": "limma", "group_col": "group"}, "geo03-plain")
        check("GEO-03k 加 age 后 p 值改变（协变量生效）",
              not np.allclose(adjusted["de_store"]["table"]["P.Value"].values,
                              plain["de_store"]["table"]["P.Value"].values, equal_nan=True))

    # 配对 + 患者内可变的连续协变量（如样本纯度）→ 应当可跑
    combo_info = info.copy()
    combo_info["purity"] = np.random.default_rng(21).uniform(0.2, 0.9, len(combo_info))
    ctx_combo = build_ctx(X, genes, samples, combo_info)
    raised = None
    combo = None
    try:
        combo = run(ctx_combo, {"method": "limma", "group_col": "group", "paired": True,
                                "pair_col": "patient", "covars": ["purity"]}, "geo03-combo")
    except EngineError as error:
        raised = error
    check("GEO-03l 配对 + 患者内可变协变量可跑", raised is None and combo is not None,
          f"{getattr(raised, 'message', '')[:80]}")

    # 配对 + 患者级不变协变量（age 每个患者只有一个值）→ 被 pair 吸收，不可识别，
    # 必须明确失败并点名共线项（这是正确的统计行为，不能悄悄丢弃）
    coll_info = info.copy()
    coll_info["age"] = np.repeat(np.linspace(40, 75, 8), 2)
    coll_info["sex"] = np.repeat(["M", "F"] * 4, 2)
    ctx_coll = build_ctx(X, genes, samples, coll_info)
    raised = None
    try:
        run(ctx_coll, {"method": "limma", "group_col": "group", "paired": True,
                       "pair_col": "patient", "covars": ["age", "sex"]}, "geo03-patlevel")
    except EngineError as error:
        raised = error
    check("GEO-03l2 配对 + 患者级协变量 → prerequisite（被 pair 吸收）",
          raised is not None and raised.code == "prerequisite",
          f"{getattr(raised, 'message', '')[:80]}")

    # 配对 + 患者内恒定的分类协变量（sex）→ 不可识别，必须明确失败并点名共线项
    coll_info = info.copy()
    coll_info["sex"] = ["M" if i % 2 else "F" for i in range(len(coll_info))]
    ctx_coll = build_ctx(X, genes, samples, coll_info)
    raised = None
    try:
        run(ctx_coll, {"method": "limma", "group_col": "group", "paired": True,
                       "pair_col": "patient", "covars": ["sex"]}, "geo03-collinear")
    except EngineError as error:
        raised = error
    check("GEO-03m 共线协变量 → prerequisite 且点名共线项",
          raised is not None and raised.code == "prerequisite" and "cov_sex" in str(getattr(raised, "message", "")),
          f"{getattr(raised, 'message', '')[:100]}")


def main() -> int:
    print("== dsh-geobulk 差异分析回归（GEO-01/02/03）==")
    test_geo01_direction()
    test_geo02_deseq2()
    test_geo03_design()
    print(f"\n== 结果：{len(PASSED)} 通过 / {len(FAILED)} 失败 ==")
    if FAILED:
        for item in FAILED:
            print(f"  ✗ {item}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
