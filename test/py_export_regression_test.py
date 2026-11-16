#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 导出回归测试（任务单 EXPORT-01 / EXPORT-02）

这两个问题是「看起来有导出、其实内容不成立」的类型，原有导出测试只检查文件存在、
非空、含 import，所以全部通过也发现不了。这里按真实执行记录逐条核对：

  EXPORT-01  Methods 必须由结构化执行记录生成：
             * method=limma / deseq2 / welch 分别跑，导出文字必须写实际方法；
             * normalize 的 none / cpm(logCPM) / quantile 分别跑，不得统一写
               "log-transformed and quantile-normalized"；
             * 没跑过生存分析不得出现 Cox / Kaplan-Meier / log-rank；
             * 没跑过富集不得出现 hypergeometric / GSEA；
             * 降级执行（requested=limma, actual=welch）必须写清 requested/actual；
             * 执行后改了分组，旧记录不得继续当作当前结论。

  EXPORT-02  导出的 script 必须能重放：
             * 脚本含真实代码（ast 顶层有函数定义、调用 bg_io.load_source /
               bg_registry.run_step），不是 import + 注释 + print；
             * 同时保存 manifest（输入路径与 sha256、样本/基因映射、按顺序的操作、
               每步实际参数、依赖版本、随机种子）；
             * 在独立空目录执行脚本，从同一输入恢复同一流程并逐项核对；
             * 输入被篡改必须失败，不得静默产出"已验证"；
             * 无法重放的产物（无源文件 / 含未记录参数的高级分析）只能标为
               记录 / 部分记录，不得标为可复现脚本。

用法：cd plugins/dsh-geobulk && python3 test/py_export_regression_test.py
成功打印 PASS 并以 0 退出；失败打印 FAIL 并以 1 退出。
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_io  # noqa: E402
import bg_state  # noqa: E402
import bg_registry  # noqa: E402
import bg_engine  # noqa: E402

OUT = Path(tempfile.mkdtemp(prefix="bg-export-reg-"))
PASSED: list[str] = []
FAILED: list[str] = []

_DS_COUNTER = {"n": 0}
_ENGINES: dict[str, dict] = {}


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def expect_text(name: str, text: str, must: tuple = (), must_not: tuple = ()) -> None:
    """同时断言必须出现与必须不出现的片段（大小写按原样）。"""
    lowered = text.lower()
    for needle in must:
        ok = needle in text or needle.lower() in lowered
        if not ok:
            check(f"{name} :: 缺少「{needle}」", False, _snippet(text))
            return
    found = [needle for needle in must_not if needle in text or needle.lower() in lowered]
    if found:
        check(f"{name} :: 不应出现「{found[0]}」", False, _snippet(text))
        return
    check(name, True, f"must={len(must)} must_not={len(must_not)}")


def _snippet(text: str, limit: int = 900) -> str:
    text = text.replace("\n", " ⏎ ")
    return text[:limit]


# ---------------------------------------------------------------------------
# fixtures：真实 CSV 输入 → bg_io 导入 → bg_state.Dataset（与 bg_engine 同一条链路）
# ---------------------------------------------------------------------------

def counts_frame(seed: int = 21, n_case: int = 5, n_ctrl: int = 5, n_genes: int = 120,
                 disp: float = 0.15, signal: float = 3.0):
    """DESeq2 可用的整数 counts（负二项，有真实离散度，否则 dispersion 拟合会失败）。"""
    rng = np.random.default_rng(seed)
    genes = [f"G{j}" for j in range(n_genes)]
    samples = [f"T{j}" for j in range(n_case)] + [f"N{j}" for j in range(n_ctrl)]
    base = np.exp(rng.uniform(np.log(30), np.log(4000), n_genes))
    lib = rng.uniform(0.7, 1.4, n_case + n_ctrl)
    mu = np.outer(lib, base)
    counts = rng.negative_binomial(1 / disp, 1 / (1 + mu * disp), size=(n_case + n_ctrl, n_genes))
    counts[:n_case, :10] = (counts[:n_case, :10] * signal).astype(int)
    frame = pd.DataFrame(counts.astype(int), index=samples, columns=genes)
    return frame, n_case, n_ctrl


def write_input(tag: str, **kwargs):
    folder = OUT / tag
    folder.mkdir(parents=True, exist_ok=True)
    frame, n_case, n_ctrl = counts_frame(**kwargs)
    path = folder / "expr.csv"
    frame.to_csv(path)
    return path, frame, n_case, n_ctrl


def dataset_from_csv(path: Path, workdir: str, n_case: int, n_ctrl: int, *,
                     clinical: bool = False, name: str = "geo-export") -> "bg_state.Dataset":
    loaded = bg_io.load_source(str(path), workdir=workdir, options={"data_format": "matrix_text"})
    samples = [str(s) for s in loaded["samples"]]
    info = pd.DataFrame({"group": ["tumor"] * n_case + ["normal"] * n_ctrl}, index=samples)
    if clinical:
        info["os_time"] = np.linspace(11, 11 + 9 * (len(samples) - 1) / 9, len(samples)).astype(float)
        info["os_event"] = [1, 1, 0, 1, 0, 1, 0, 1, 1, 0][: len(samples)]
    kind = loaded.get("value_kind") or (loaded.get("meta") or {}).get("value_kind") or "unknown"
    X = np.asarray(loaded["X"])
    counts = loaded.get("counts")
    if counts is None and str(kind) == "counts":
        counts = X
    _DS_COUNTER["n"] += 1
    return bg_state.Dataset(
        f"ds{_DS_COUNTER['n']}", name, X=X, genes=list(loaded["genes"]), samples=samples,
        counts=counts, sample_info=info, source_path=str(path),
        source_format=str(loaded.get("format") or "matrix_text"), species="human",
        value_kind=str(kind),
    )


def new_case(tag: str, *, chain: list, clinical: bool = False, **kwargs):
    """建一个干净的数据集 + 工作目录，并跑完 chain。"""
    path, frame, n_case, n_ctrl = write_input(tag, **kwargs)
    workdir = str(OUT / tag / "work")
    ds = dataset_from_csv(path, workdir, n_case, n_ctrl, clinical=clinical, name=tag)
    run_steps(ds, workdir, chain)
    return ds, workdir, path, frame


def engine_for(workdir: str) -> dict:
    key = str(workdir)
    if key not in _ENGINES:
        _ENGINES[key] = bg_engine.Worker(str(workdir)).engine_info()
    return _ENGINES[key]


def run_steps(ds, workdir: str, chain: list):
    engine = engine_for(workdir)
    outs = []
    for step_id, params in chain:
        outs.append(bg_registry.run_step(ds, step_id, dict(params), workdir, engine))
    return outs


def export_artifact(ds, workdir: str, what: str) -> tuple[dict, Path]:
    res = bg_registry.export_dataset(ds, what, None, workdir, engine_for(workdir))
    return res, Path(res["file"])


def export_text(ds, workdir: str, what: str = "methods") -> tuple[str, dict]:
    res, path = export_artifact(ds, workdir, what)
    return path.read_text(encoding="utf-8"), res


QC_FILTER = [("qc", {}), ("filter", {})]


# ---------------------------------------------------------------------------
# EXPORT-01：Methods 必须由真实执行记录生成
# ---------------------------------------------------------------------------

def test_01_no_fabricated_methods():
    """没有任何分析步骤时，不得声称做过标准化 / 差异 / 富集 / 生存。"""
    print("\n[EXPORT-01] 无执行记录时不得编造方法")
    path, frame, n_case, n_ctrl = write_input("no-steps")
    workdir = str(OUT / "no-steps" / "work")
    ds = dataset_from_csv(path, workdir, n_case, n_ctrl, name="no-steps")
    text, _ = export_text(ds, workdir)

    expect_text(
        "EXPORT-01a 空记录 Methods 不写任何未执行方法",
        text,
        must=("The dataset contained 10 samples and 120 genes",),
        must_not=(
            "moderated t-statistics", "quantile-normalized", "log-transformed",
            "Welch", "Cox", "Kaplan-Meier", "log-rank", "hypergeometric", "GSEA",
            "quantile", "limma::lmFit", "DESeq2",
        ),
    )
    expect_text(
        "EXPORT-01b 空记录 Methods 明确写未记录",
        text,
        must=("no analysis step was recorded",),
    )


def test_01_normalize_variants():
    """none / cpm / quantile 三种标准化分别导出，文字必须与真实执行一致。"""
    print("\n[EXPORT-01] 标准化方法逐条核对（none / cpm / quantile）")

    ds, workdir, _, _ = new_case(
        "norm-cpm",
        chain=QC_FILTER + [("normalize", {"method": "cpm", "log_transform": "auto"})],
    )
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01c cpm → 写 log2 counts per million",
        text,
        must=("log2 counts per million", "requested method: `cpm`"),
        must_not=("quantile",),
    )
    check("EXPORT-01c2 cpm 记录 value_kind counts→lognorm",
          ds.value_kind == "lognorm"
          and "`counts` → `lognorm`" in text,
          f"value_kind={ds.value_kind}")

    ds, workdir, _, _ = new_case(
        "norm-quantile",
        chain=QC_FILTER + [("normalize", {"method": "quantile", "log_transform": "yes"})],
    )
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01d quantile → 写 quantile-normalized",
        text,
        must=("quantile-normalized", "requested method: `quantile`"),
        must_not=("log2 counts per million",),
    )

    ds, workdir, _, _ = new_case(
        "norm-none",
        chain=QC_FILTER + [("normalize", {"method": "none", "log_transform": "no"})],
    )
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01e none → 不得再写 log-transformed and quantile-normalized",
        text,
        must=("No normalization was applied", "original `counts` scale"),
        must_not=("log-transformed and quantile-normalized", "log2 counts per million",
                  "quantile-normalized", "moderated t-statistics"),
    )
    check("EXPORT-01e2 none 保持 counts 尺度并提示 logFC 不可按 log2 解读",
          ds.value_kind == "counts" and "must not be interpreted as log2 fold changes" in text,
          f"value_kind={ds.value_kind}")


def test_01_de_methods():
    """limma / DESeq2 / welch 三种差异方法分别导出，文字必须写实际方法。"""
    print("\n[EXPORT-01] 差异方法逐条核对（limma / deseq2 / welch）")

    ds, workdir, _, _ = new_case(
        "de-limma",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"})],
    )
    text, _ = export_text(ds, workdir)
    store = ds.results["de"]
    de_metrics = {}
    for entry in reversed(ds.history):
        if entry.get("step") == "de" and entry.get("ok"):
            de_metrics = entry.get("metrics") or {}
            break
    expect_text(
        "EXPORT-01f limma → moderated t-statistics（R 官方）",
        text,
        must=("limma moderated t-statistics", "limma::lmFit", "Benjamini-Hochberg"),
        must_not=("Welch", "DESeq2"),
    )
    check("EXPORT-01f2 limma 对比方向写 case - control",
          "`tumor - normal`" in text and store["case"] == "tumor" and store["control"] == "normal",
          f"{store['case']} - {store['control']}")
    design = str(de_metrics.get("design") or "")
    check("EXPORT-01f3 limma 设计公式与执行记录逐字一致",
          bool(design) and f"Design: `{design}`" in text,
          f"design={design!r}")

    ds, workdir, _, _ = new_case(
        "de-deseq2",
        chain=QC_FILTER + [("normalize", {"method": "none", "log_transform": "no"}),
                           ("de", {"method": "deseq2"})],
    )
    text, _ = export_text(ds, workdir)
    store = ds.results["de"]
    expect_text(
        "EXPORT-01g DESeq2 → 负二项 Wald（不是 moderated t）",
        text,
        must=("DESeq2 negative binomial", "Wald"),
        must_not=("moderated t-statistics", "Welch"),
    )
    check("EXPORT-01g2 DESeq2 输入尺度写 counts",
          "Input matrix scale: `counts`" in text and store["method"] == "deseq2",
          f"method={store['method']}")

    ds, workdir, _, _ = new_case(
        "de-welch",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "welch"})],
    )
    text, _ = export_text(ds, workdir)
    store = ds.results["de"]
    expect_text(
        "EXPORT-01h welch → Welch t 检验（不是 moderated t）",
        text,
        must=("Welch's two-sample t-test", "tumor - normal"),
        must_not=("moderated t-statistics", "limma::lmFit", "DESeq2"),
    )
    check("EXPORT-01h2 welch 记录 actual_method=welch",
          store["method"] == "welch" and store["method_requested"] == "welch",
          f"{store['method_requested']}→{store['method']}")


def test_01_de_fallback():
    """requested=limma 但 R 不可用 → actual=welch，Methods 必须写清降级。"""
    print("\n[EXPORT-01] 降级执行（requested=limma → actual=welch）")
    import bg_de

    original = bg_de._r_diff

    def fake_r_diff(*args, **kwargs):
        return {"frame": None, "failed": True, "code": "no_r",
                "error": "Rscript 不可用（测试模拟环境缺失）"}

    bg_de._r_diff = fake_r_diff
    try:
        ds, workdir, _, _ = new_case(
            "de-fallback",
            chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"})],
        )
    finally:
        bg_de._r_diff = original

    store = ds.results["de"]
    check("EXPORT-01i 降级记录 requested=limma / actual=welch / lite=True",
          store["method_requested"] == "limma" and store["method"] == "welch" and store["lite"],
          f"{store['method_requested']}→{store['method']} lite={store['lite']}")

    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01j Methods 写 requested/actual 与回退原因，不写 moderated t",
        text,
        must=("Welch's two-sample t-test", "limma", "fallback"),
        must_not=("moderated t-statistics",),
    )


def test_01_survival_claims():
    """没跑生存不得出现 Cox/Kaplan；跑了但没拟合 Cox 也不得声称 Cox。"""
    print("\n[EXPORT-01] 生存分析：只说真实执行过的内容")
    ds, workdir, _, _ = new_case(
        "surv-never",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"})],
        clinical=True,
    )
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01k 未跑生存 → 无 Cox/Kaplan/log-rank",
        text,
        must=("limma moderated t-statistics",),
        must_not=("Cox", "Kaplan-Meier", "log-rank", "hazard ratio"),
    )

    ds, workdir, _, _ = new_case(
        "surv-km-only",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}),
                           ("survival", {"group_col": "group", "time_col": "os_time",
                                         "event_col": "os_event"})],
        clinical=True,
    )
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01l 只做 KM/log-rank → 不得写 Cox",
        text,
        must=("Kaplan-Meier", "log-rank", "`os_time`"),
        must_not=("Cox", "hazard ratio"),
    )

    ds, workdir, _, _ = new_case(
        "surv-cox",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}),
                           ("survival", {"gene": "G0", "time_col": "os_time",
                                         "event_col": "os_event"})],
        clinical=True,
    )
    metrics = None
    for entry in reversed(ds.history):
        if entry.get("step") == "survival":
            metrics = entry.get("metrics") or {}
            break
    check("EXPORT-01m 生存记录里确实有 Cox HR（新增强制记录）",
          bool(metrics) and metrics.get("HR") is not None,
          f"metrics={metrics}")
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01n 拟合了 Cox → Methods 写 Cox proportional hazards",
        text,
        must=("Kaplan-Meier", "log-rank", "Cox proportional hazards"),
    )


def test_01_stale_group_record():
    """执行后改了分组：旧 DE 记录不属于当前数据版本，不得继续作为结论。"""
    print("\n[EXPORT-01] 分组变化后旧记录失效")
    ds, workdir, _, _ = new_case(
        "stale-group",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"})],
    )
    before, _ = export_text(ds, workdir)
    check("EXPORT-01o 改动前 Methods 有肿瘤/正常对比",
          "tumor - normal" in before, "有对比")

    # 交换分组标签（数据集仍是 10 样本，但分组含义已变）
    labels = ds.sample_info["group"].astype(str).tolist()
    swapped = ["normal" if v == "tumor" else "tumor" for v in labels]
    ds.sample_info["group"] = swapped

    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01p 分组变了 → 不再把旧的 tumor - normal 当结论",
        text,
        must=("does not match the current sample groups",),
        must_not=("tumor - normal",),
    )
    _res, _path = export_artifact(ds, workdir, "script")
    check("EXPORT-01p2 记录过期时脚本不得标为可复现",
          _res.get("replayable") is False,
          f"replayable={_res.get('replayable')} reasons={_res.get('replay_reasons')}")


def test_01_report_uses_same_record():
    """report 的 Methods 段必须与 methods 导出一致（同一个结构化记录）。"""
    print("\n[EXPORT-01] report 与 methods 同源")
    ds, workdir, _, _ = new_case(
        "report-same",
        chain=QC_FILTER + [("normalize", {"method": "quantile", "log_transform": "yes"}),
                           ("de", {"method": "welch"})],
    )
    methods_text, _ = export_text(ds, workdir, "methods")
    report_text, _ = export_text(ds, workdir, "report")
    # report 里 Methods 正文是按 methods 导出逐行拼进去的：逐行核对（跳过带时间戳的行）
    core_lines = [line for line in methods_text.splitlines()
                  if line.startswith("- **") or line.startswith("- Differential expression used")]
    missing = [line[:70] for line in core_lines if line not in report_text]
    check("EXPORT-01q report 的 Methods 段与 methods 导出逐行一致",
          bool(core_lines) and not missing, f"missing={missing[:3]}")
    expect_text("EXPORT-01q2 report 同样写 Welch + quantile，不写 moderated t",
                report_text,
                must=("Welch's two-sample t-test", "quantile-normalized"),
                must_not=("moderated t-statistics",))


def test_01_analysis_must_not_leak_into_step_methods():
    """高级分析会覆盖 dataset.results，但 Methods 必须描述「步骤执行当时」的结果。"""
    print("\n[EXPORT-01] 高级分析不得改写步骤的 Methods 记录")
    ds, workdir, _, _ = new_case(
        "enrich-snapshot",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"}),
                           ("enrich", {"method": "ora", "databases": ["GO_BP"]})],
    )
    step_store = dict(ds.results.get("enrich") or {})
    bg_registry.run_analysis(ds, "gsea_hallmark", {}, workdir, engine_for(workdir))
    now_store = dict(ds.results.get("enrich") or {})
    check("EXPORT-01r 高级分析确实覆盖了 results['enrich']",
          step_store.get("method") != now_store.get("method")
          or step_store.get("databases") != now_store.get("databases"),
          f"step={step_store.get('method')}/{step_store.get('databases')} "
          f"analysis={now_store.get('method')}/{now_store.get('databases')}")
    text, _ = export_text(ds, workdir)
    enrich_line = next((line for line in text.splitlines() if line.startswith("- **enrich**")), "")
    expect_text(
        "EXPORT-01r2 enrich 行写步骤当时的 ORA 数据库，不写分析覆盖后的 GSEA",
        enrich_line,
        must=("hypergeometric", "GO_BP"),
        must_not=("HALLMARK", "Kolmogorov"),
    )


def test_01_legacy_record_fallback():
    """升级前的旧记录（无 metrics/快照）也不能编造：方法从 store 取，设计写未记录。"""
    print("\n[EXPORT-01] 旧记录（无结构化字段）的降级处理")
    ds, workdir, _, _ = new_case(
        "legacy-record",
        chain=QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"})],
    )
    for entry in ds.history:
        for key in ("metrics", "notes", "warnings", "data_version", "result_snapshot", "error"):
            entry.pop(key, None)
    text, _ = export_text(ds, workdir)
    expect_text(
        "EXPORT-01s 旧记录：实际方法从 store 取，设计公式标未记录",
        text,
        must=("limma moderated t-statistics", "tumor - normal", "未记录"),
        must_not=("Welch", "Design: `"),
    )


# ---------------------------------------------------------------------------
# EXPORT-02：可运行 / 可重放的导出脚本
# ---------------------------------------------------------------------------

REPLAY_CHAIN = QC_FILTER + [("normalize", {"method": "cpm"}), ("de", {"method": "limma"})]


def test_02_script_structure():
    """脚本必须是真实代码 + manifest，而不是 import + 注释 + print。"""
    print("\n[EXPORT-02] 脚本结构与 manifest")
    ds, workdir, csv_path, _ = new_case("script-structure", chain=REPLAY_CHAIN)
    res, script_path = export_artifact(ds, workdir, "script")
    text = script_path.read_text(encoding="utf-8")

    try:
        tree = ast.parse(text)
        parsed = True
    except SyntaxError as error:
        tree = None
        parsed = False
        check("EXPORT-02a 导出脚本可被 ast 解析", False, f"{error}")
    if parsed:
        funcs = [node.name for node in tree.body if isinstance(node, ast.FunctionDef)]
        check("EXPORT-02a 导出脚本可被 ast 解析且顶层含函数定义",
              True, f"top-level funcs={funcs[:6]}")

    expect_text(
        "EXPORT-02b 脚本含真实调用（导入 + 按记录重跑）",
        text,
        must=("bg_io.load_source(", "bg_registry.run_step(", "hashlib.sha256",
              "MANIFEST", "sys.exit(main())", "REPLAY STATUS"),
        must_not=("# expr = pd.read_csv('expr.csv', index_col=0)",),
    )
    check("EXPORT-02c 脚本含记录的实际差异方法参数",
          '"method": "limma"' in text and '"id": "normalize"' in text,
          "含 limma/normalize 操作参数")
    check("EXPORT-02d 脚本头标明可复现脚本",
          "可复现脚本" in text.split("MANIFEST", 1)[0], "头部标签")

    manifest_path = Path(res.get("manifest") or "")
    check("EXPORT-02e 同时写出 manifest 文件", bool(str(manifest_path)) and manifest_path.exists(),
          str(manifest_path))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}

    sha = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    check("EXPORT-02f manifest 记录输入路径与 sha256",
          manifest.get("input", {}).get("source_path") == str(csv_path)
          and manifest.get("input", {}).get("sha256") == sha,
          f"sha={str(manifest.get('input', {}).get('sha256'))[:12]}")
    check("EXPORT-02g manifest 记录样本/基因映射（数量 + 校验值 + 清单）",
          manifest.get("input", {}).get("genes_sha256")
          and manifest.get("input", {}).get("samples_sha256")
          and len(manifest.get("input", {}).get("genes") or []) == ds.n_genes
          and len(manifest.get("input", {}).get("samples") or []) == ds.n_samples,
          f"genes={len(manifest.get('input', {}).get('genes') or [])} samples={len(manifest.get('input', {}).get('samples') or [])}")
    ops = manifest.get("operations") or []
    check("EXPORT-02h manifest 按顺序记录有效操作与实际参数",
          [op.get("id") for op in ops] == list(ds.steps)
          and ops[-1].get("params", {}).get("method") == "limma",
          f"ops={[(op.get('id'), (op.get('params') or {}).get('method')) for op in ops]}")
    check("EXPORT-02i manifest 记录依赖版本与随机种子",
          manifest.get("versions", {}).get("python") == platform.python_version()
          and isinstance(manifest.get("versions", {}).get("packages"), dict)
          and manifest.get("seed") == 0,
          f"python={manifest.get('versions', {}).get('python')} seed={manifest.get('seed')}")
    expected = manifest.get("expected") or {}
    check("EXPORT-02j manifest 记录期望终态（样本/基因/尺度/矩阵校验值/对比方向）",
          expected.get("n_samples") == ds.n_samples
          and expected.get("n_genes") == ds.n_genes
          and expected.get("value_kind") == ds.value_kind
          and bool(expected.get("matrix_sha256"))
          and (expected.get("de") or {}).get("case") == "tumor"
          and (expected.get("de") or {}).get("control") == "normal",
          f"de={expected.get('de', {}).get('method')}/{expected.get('de', {}).get('case')}")
    check("EXPORT-02k 导出结果声明 replayable=True",
          res.get("replayable") is True and "可复现脚本" in str(res.get("label")),
          f"label={res.get('label')}")


def test_02_replay_isolated():
    """把脚本 + manifest 复制到独立空目录，从同一输入重放并逐项核对。"""
    print("\n[EXPORT-02] 独立空目录重放")
    ds, workdir, csv_path, _ = new_case("script-replay", chain=REPLAY_CHAIN)
    res, script_path = export_artifact(ds, workdir, "script")
    manifest_path = Path(res.get("manifest") or "")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    recorded_de = (manifest.get("expected") or {}).get("de") or {}

    fresh = Path(tempfile.mkdtemp(prefix="bg-replay-"))
    shutil.copy(script_path, fresh / script_path.name)
    shutil.copy(manifest_path, fresh / manifest_path.name)
    env = {**os.environ, "MPLBACKEND": "Agg"}
    started = time.time()
    proc = subprocess.run(
        [sys.executable, script_path.name],
        cwd=str(fresh), capture_output=True, text=True, timeout=900, env=env,
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    check("EXPORT-02l 重放脚本在独立目录执行成功", proc.returncode == 0,
          f"exit={proc.returncode} {output[-400:]}")
    check("EXPORT-02m 重放结束状态为 verified",
          "REPLAY STATUS: verified" in output, f"{output[-260:]}")

    report_path = fresh / "replay_out" / "replay_report.json"
    check("EXPORT-02n 重放写出核对报告", report_path.exists(), str(report_path))
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {}
    checks = report.get("checks") or []
    failed_checks = [item for item in checks if not item.get("ok")]
    check("EXPORT-02o 重放逐项核对全部通过",
          bool(checks) and not failed_checks,
          f"failed={[item.get('name') for item in failed_checks]}")
    check("EXPORT-02p 重放恢复样本/基因/尺度/矩阵",
          report.get("dataset", {}).get("n_samples") == ds.n_samples
          and report.get("dataset", {}).get("n_genes") == ds.n_genes
          and report.get("dataset", {}).get("value_kind") == ds.value_kind
          and report.get("dataset", {}).get("matrix_sha256") == manifest["expected"]["matrix_sha256"],
          f"{report.get('dataset')}")
    check("EXPORT-02q 重放恢复对比方向与实际方法",
          (report.get("de") or {}).get("method") == recorded_de.get("method")
          and (report.get("de") or {}).get("case") == "tumor"
          and (report.get("de") or {}).get("control") == "normal",
          f"{report.get('de', {}).get('method')} {report.get('de', {}).get('case')}-{report.get('de', {}).get('control')}")

    # 确定性步骤逐值比较：重放得到的 DE logFC 与记录值一致
    replay_rows = {row["gene"]: row for row in (report.get("de") or {}).get("rows") or []}
    recorded_rows = {row["gene"]: row for row in recorded_de.get("rows") or []}
    common = [gene for gene in recorded_rows if gene in replay_rows]
    worst = 0.0
    worst_gene = ""
    for gene in common:
        a = float(replay_rows[gene]["logFC"])
        b = float(recorded_rows[gene]["logFC"])
        diff = abs(a - b) / max(1.0, abs(b))
        if diff > worst:
            worst, worst_gene = diff, gene
    check("EXPORT-02r 重放 DE logFC 与记录逐值一致（相对误差 <1e-8）",
          len(common) >= 10 and worst < 1e-8,
          f"n={len(common)} max_rel={worst:.2e} gene={worst_gene}")
    print(f"  INFO 重放耗时 {time.time() - started:.1f}s")


def test_02_replay_tamper():
    """输入被改动后必须校验失败，不得静默产出"已验证"。"""
    print("\n[EXPORT-02] 输入篡改必须失败")
    ds, workdir, csv_path, _ = new_case("script-tamper", chain=REPLAY_CHAIN)
    res, script_path = export_artifact(ds, workdir, "script")
    manifest_path = Path(res.get("manifest") or "")

    fresh = Path(tempfile.mkdtemp(prefix="bg-replay-bad-"))
    shutil.copy(script_path, fresh / script_path.name)
    shutil.copy(manifest_path, fresh / manifest_path.name)
    lines = csv_path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[1] = str(int(float(cells[1])) + 12345)  # 改第一个样本第一个基因的值
    lines[1] = ",".join(cells)
    csv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    proc = subprocess.run(
        [sys.executable, script_path.name],
        cwd=str(fresh), capture_output=True, text=True, timeout=900,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    check("EXPORT-02s 输入 sha256 不符时脚本失败退出", proc.returncode != 0,
          f"exit={proc.returncode}")
    expect_text(
        "EXPORT-02t 失败信息指出输入校验不一致",
        output,
        must=("sha256",),
        must_not=("REPLAY STATUS: verified",),
    )


def test_02_label_not_replayable():
    """无源文件的产物只能标为记录，不得标为可复现脚本。"""
    print("\n[EXPORT-02] 无源文件 → 标记为记录/示例")
    path, frame, n_case, n_ctrl = write_input("no-source")
    workdir = str(OUT / "no-source" / "work")
    ds = dataset_from_csv(path, workdir, n_case, n_ctrl, name="no-source")
    run_steps(ds, workdir, QC_FILTER + [("normalize", {"method": "cpm"})])
    ds.source_path = None  # 例如 demo 数据集：没有可回放的原文件

    res, script_path = export_artifact(ds, workdir, "script")
    text = script_path.read_text(encoding="utf-8")
    head = text.split("MANIFEST", 1)[0]
    check("EXPORT-02u replayable=False 且标签为记录", res.get("replayable") is False,
          f"label={res.get('label')}")
    check("EXPORT-02v 脚本头明确写不可独立复现", "不可独立复现" in head, head[:200])
    try:
        ast.parse(text)
        ok, detail = True, ""
    except SyntaxError as error:
        ok, detail = False, str(error)
    check("EXPORT-02w 记录版脚本仍是合法 Python（只是不能独立重放）", ok, detail)


def test_02_label_partial_analysis():
    """有高级分析记录（参数未纳入重放）时只能标为部分记录。"""
    print("\n[EXPORT-02] 高级分析未纳入重放 → 部分记录")
    ds, workdir, _, _ = new_case("partial", chain=REPLAY_CHAIN)
    engine = engine_for(str(workdir))
    bg_registry.run_analysis(ds, "de_summary", {}, workdir, engine)
    check("EXPORT-02x 高级分析确实执行并产生记录",
          any(item.get("id") == "de_summary" and item.get("ok") for item in ds.analyses_run),
          f"analyses_run={[item.get('id') for item in ds.analyses_run]}")

    res, script_path = export_artifact(ds, workdir, "script")
    text = script_path.read_text(encoding="utf-8")
    head = text.split("MANIFEST", 1)[0]
    check("EXPORT-02y 部分记录时 replayable=False 且标签含部分",
          res.get("replayable") is False and "部分" in str(res.get("label")),
          f"label={res.get('label')} reasons={res.get('replay_reasons')}")
    check("EXPORT-02z 脚本头列出未纳入重放的项", "de_summary" in head, head[:300])


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------

TESTS = [
    test_01_no_fabricated_methods,
    test_01_normalize_variants,
    test_01_de_methods,
    test_01_de_fallback,
    test_01_survival_claims,
    test_01_stale_group_record,
    test_01_report_uses_same_record,
    test_01_analysis_must_not_leak_into_step_methods,
    test_01_legacy_record_fallback,
    test_02_script_structure,
    test_02_replay_isolated,
    test_02_replay_tamper,
    test_02_label_not_replayable,
    test_02_label_partial_analysis,
]


def main() -> int:
    started = time.time()
    print("=" * 78)
    print("dsh-geobulk 导出回归测试（EXPORT-01 / EXPORT-02）")
    print(f"工作目录：{OUT}")
    print("=" * 78)
    for test in TESTS:
        try:
            test()
        except Exception as error:  # noqa: BLE001
            FAILED.append(test.__name__)
            print(f"  FAIL {test.__name__} 未捕获异常：{type(error).__name__}: {error}")
            traceback.print_exc(limit=6)
    print("\n" + "=" * 78)
    print(f"总计 {len(PASSED) + len(FAILED)} checks：PASS {len(PASSED)} / FAIL {len(FAILED)}"
          f"（{time.time() - started:.1f}s）")
    if FAILED:
        print("失败项：")
        for item in FAILED:
            print(f"  - {item}")
        return 1
    print("全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
