#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 单细胞映射（bg_scmap.py）单元测试（CONTRACT §11）

不启动引擎、不依赖网络，直接构造 ctx 调 `bg_scmap.run_scmap` / 各分析：

  ctx = {
    "X": (20 样本 × 300 基因) numpy 矩阵, "genes"/"samples",
    "de": {"table": DataFrame(gene/logFC/AveExpr/t/P.Value/adj.P.Val/B/regulation)},
    "sc_ref": test/fixtures/demo_sc.h5ad,
    ...
  }

断言：
  1. 得分 ∈ [0, 1]（归一化列）；raw_score 有限
  2. 细胞类型表 6 行（与 sc_truth.json 一致）、细胞数合计 = 400
  3. 命中：用每个细胞类型的 truth marker 当 bulk 上调签名，top-1 命中率 ≥ 60%（契约 §11）
  4. 缺参考 → prerequisite；路径不存在 → not_found；非法方法 → bad_request；缺 de → prerequisite
  5. 三种打分方法都跑通；每个结果至少 1 图 + 1 表，图文件真存在
  6. STEPS/ANALYSES 注册符合契约（step=scmap/requires=de；分析 ≥6 且 needs 含 sc）
  7. 6 个 ANALYSIS 全部可跑，ctx 不被原地修改（契约 §2.3 只读）

用法：cd plugins/dsh-geobulk && python3 test/py_scmap_test.py
成功打印 PASS 并以 0 退出；失败打印 FAIL 并以 1 退出。
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import traceback
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/bg-mpl")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bg_scmap  # noqa: E402
from bg_common import EngineError  # noqa: E402

FIXTURES = ROOT / "test" / "fixtures"
SC_REF = FIXTURES / "demo_sc.h5ad"
TRUTH = json.loads((FIXTURES / "sc_truth.json").read_text(encoding="utf-8"))
OUTDIR = tempfile.mkdtemp(prefix="bg-scmap-test-")

PASSED: list[str] = []
FAILED: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def build_ctx(up_genes: list[str], *, n_samples: int = 20, n_genes: int = 300,
              ref: str | None | bool = True, with_de: bool = True) -> dict:
    """numpy 造 bulk 矩阵 + de 表 + sc_ref。up_genes 在 case 组上调 2 个 log2 单位。"""
    rng = np.random.default_rng(0)
    up = [str(g) for g in up_genes]
    rest = [f"GENE{i:04d}" for i in range(max(0, n_genes - len(up)))]
    genes = up + rest
    X = rng.normal(6.0, 1.0, (n_samples, n_genes)).astype(np.float64)
    labels = ["case"] * (n_samples // 2) + ["control"] * (n_samples - n_samples // 2)
    X[: n_samples // 2, : len(up)] += 2.0
    samples = [f"S{i:02d}" for i in range(n_samples)]
    lfc = np.concatenate([np.full(len(up), 2.0), np.zeros(len(rest))])
    pvals = np.full(n_genes, 1e-6)
    if not up:
        pvals = np.ones(n_genes)
    de = pd.DataFrame({
        "gene": genes,
        "logFC": lfc,
        "AveExpr": X.mean(axis=0),
        "t": np.zeros(n_genes),
        "P.Value": pvals,
        "adj.P.Val": pvals,
        "B": np.zeros(n_genes),
    })
    de["regulation"] = np.where((de["adj.P.Val"] <= 0.05) & (de["logFC"] > 0), "up",
                                np.where((de["adj.P.Val"] <= 0.05) & (de["logFC"] < 0), "down", "ns"))
    ref_path = os.path.abspath(str(SC_REF)) if ref is True else (ref if isinstance(ref, str) else None)
    return {
        "dataset_id": "ds_test",
        "name": "synthetic",
        "genes": genes,
        "samples": samples,
        "X": X.astype(np.float32),
        "X_kind": "lognorm",
        "counts": None,
        "sample_info": pd.DataFrame({"group": labels}, index=samples),
        "group_key": "group",
        "groups": ["control", "case"],
        "species": "human",
        "de": {"table": de} if with_de else None,
        "results": {},
        "modules": {},
        "rng": np.random.default_rng(0),
        "workdir": OUTDIR,
        "outdir": os.path.join(OUTDIR, "an"),
        "engine": {},
        "sc_ref": ref_path,
    }


def run(ctx, params, tag) -> dict:
    out = os.path.join(OUTDIR, f"an-{tag}")
    return bg_scmap.run_scmap(ctx, params, out)


def score_rows(res) -> list[dict]:
    for t in res.get("tables") or []:
        if t.get("id") == "sc_celltype_score":
            return t["rows"]
    return []


def require_figure_table(res, tag) -> None:
    tables = res.get("tables") or []
    figures = res.get("figures") or []
    check(f"{tag}: 至少 1 张表", len(tables) >= 1, f"{len(tables)} 表")
    check(f"{tag}: 至少 1 张图", len(figures) >= 1, f"{len(figures)} 图")
    ok_files = all(os.path.isfile(f.get("file", "")) and os.path.getsize(f["file"]) > 0 for f in figures)
    check(f"{tag}: 图文件存在且非空", ok_files,
          ", ".join(os.path.basename(f["file"]) for f in figures[:3]))


def main() -> int:
    print("=" * 74)
    print("bg_scmap 单细胞映射测试")
    print("=" * 74)
    print(f"参考 h5ad : {SC_REF}")
    check("fixture demo_sc.h5ad 存在", SC_REF.is_file())
    check("fixture sc_truth.json 解析出 6 个细胞类型", len(TRUTH.get("cell_types") or []) == 6,
          ", ".join(TRUTH.get("cell_types") or []))

    cell_types = list(TRUTH["cell_types"])
    markers = {k: list(v) for k, v in TRUTH["markers"].items()}
    name_map = {str(ct).lower(): str(ct) for ct in cell_types}

    # ------------------------------------------------------------------
    # 1) 基础运行：全部 marker 混合签名
    # ------------------------------------------------------------------
    print("\n[1] run_scmap 基础运行（de_up，mean）")
    all_markers = [g for ct in cell_types for g in markers.get(ct, [])]
    ctx = build_ctx(all_markers[:80])
    res = run(ctx, {"gene_set": "de_up", "score_method": "mean"}, "base")
    rows = score_rows(res)
    check("返回 summary.title", bool(res["summary"]["title"]), res["summary"]["title"])
    check("细胞类型表 6 行", len(rows) == 6, f"{len(rows)} 行")
    names = {str(r["cell_type"]) for r in rows}
    check("细胞类型名与 truth 一致", names == set(cell_types), ", ".join(sorted(names)))
    scores = [float(r["score"]) for r in rows]
    check("归一化得分 ∈ [0,1]", all(0.0 <= s <= 1.0 for s in scores), f"min={min(scores):.3f} max={max(scores):.3f}")
    raw = [float(r["raw_score"]) for r in rows]
    check("raw_score 全为有限值", all(np.isfinite(raw)), f"min={min(raw):.3f} max={max(raw):.3f}")
    check("最高得分 = 1.0（min-max 归一化上界）", abs(max(scores) - 1.0) < 1e-6, f"{max(scores):.4f}")
    total_cells = sum(int(r["n_cells"]) for r in rows)
    check("细胞数合计 = 参考 400", total_cells == int(TRUTH["n_cells"]), f"{total_cells} 细胞")
    require_figure_table(res, "基础")
    check("命中参考基因 = 签名基因（marker 全在参考里）",
          res["summary"]["metrics"]["signature_matched"] == res["summary"]["metrics"]["signature_genes"],
          f"{res['summary']['metrics']['signature_matched']}/{res['summary']['metrics']['signature_genes']}")

    # ------------------------------------------------------------------
    # 2) 命中率：每个细胞类型自己的 marker 当 bulk 上调签名
    # ------------------------------------------------------------------
    print("\n[2] 命中率：每类 truth marker → top-1 亚群")
    hits, details = 0, []
    for ct in cell_types:
        sig = markers.get(ct, [])[:12]
        out = run(build_ctx(sig), {"gene_set": "de_up", "score_method": "mean"}, f"hit-{ct.replace(' ', '_')}")
        top = score_rows(out)[0]
        got = str(top["cell_type"])
        hits += (got == ct)
        details.append(f"{ct}→{got}")
        if top["score"] is not None:
            check(f"命中 {ct}: 得分 ∈ [0,1]", 0.0 <= float(top["score"]) <= 1.0)
        require_figure_table(out, f"命中-{ct.replace(' ', '_')}")
    acc = hits / max(1, len(cell_types))
    check("top-1 命中率 ≥ 60%（契约 §11）", acc >= 0.60, f"{hits}/{len(cell_types)} = {acc:.0%}；" + "；".join(details))

    # ------------------------------------------------------------------
    # 3) 三种打分方法 + de_down / de_both / custom
    # ------------------------------------------------------------------
    print("\n[3] 三种打分方法 + 签名来源分支")
    for method in ("mean", "aucell", "addmodule"):
        out = run(build_ctx(all_markers[:60]), {"gene_set": "de_up", "score_method": method}, f"method-{method}")
        mrows = score_rows(out)
        mscores = [float(r["score"]) for r in mrows]
        check(f"score_method={method}: 6 行且得分 ∈ [0,1]",
              len(mrows) == 6 and all(0.0 <= s <= 1.0 for s in mscores),
              f"{len(mrows)} 行 max={max(mscores):.3f}")
        require_figure_table(out, f"method-{method}")
    # de_down：把 marker 放到 control 组上调 → 需要重新造 ctx，因此用 custom + genes 直接给列表
    out = run(build_ctx([]), {"gene_set": "custom", "genes": markers[cell_types[0]][:10],
                              "score_method": "mean"}, "custom")
    check("gene_set=custom：细胞类型表 6 行", len(score_rows(out)) == 6)
    check("gene_set=custom：命中基因数 = 10", out["summary"]["metrics"]["signature_matched"] == 10,
          str(out["summary"]["metrics"]["signature_matched"]))

    # ------------------------------------------------------------------
    # 4) 错误分支（缺 ref / 路径不存在 / 非法方法 / 缺 de）
    # ------------------------------------------------------------------
    print("\n[4] 错误分支")

    def expect_error(code: str, params, ctx_obj, tag) -> None:
        try:
            run(ctx_obj, params, f"err-{tag}")
        except EngineError as exc:
            check(f"{tag} → {code}", exc.code == code, f"实际 {exc.code}: {exc.message[:60]}")
        except Exception as exc:  # noqa: BLE001
            check(f"{tag} → {code}", False, f"未捕获 {type(exc).__name__}: {exc}")
        else:
            check(f"{tag} → {code}", False, "未抛异常")

    expect_error("prerequisite", {}, build_ctx(all_markers[:20], ref=False), "缺 sc_ref")
    expect_error("not_found", {"ref_path": os.path.join(OUTDIR, "nope.h5ad")},
                 build_ctx(all_markers[:20]), "参考路径不存在")
    expect_error("bad_request", {"score_method": "bogus"}, build_ctx(all_markers[:20]), "非法打分方法")
    expect_error("prerequisite", {"gene_set": "de_up"}, build_ctx(all_markers[:20], with_de=False), "缺 de")
    expect_error("prerequisite", {"gene_set": "custom", "genes": ["NOT_A_GENE_XYZ"]},
                 build_ctx(all_markers[:20]), "签名基因全未命中")

    # ------------------------------------------------------------------
    # 5) 注册表 & 只读性
    # ------------------------------------------------------------------
    print("\n[5] 注册表 & ctx 只读")
    step_ids = [s["id"] for s in bg_scmap.STEPS]
    check("STEPS 含 scmap", "scmap" in step_ids, str(step_ids))
    step = next((s for s in bg_scmap.STEPS if s["id"] == "scmap"), {})
    check("scmap.requires == ['de']", list(step.get("requires") or []) == ["de"], str(step.get("requires")))
    check("scmap 有 params 且含 score_method/ref_path",
          {"score_method", "ref_path"} <= {p.get("key") for p in step.get("params") or []})
    ana_ids = [a["id"] for a in bg_scmap.ANALYSES]
    check("ANALYSES ≥ 6", len(ana_ids) >= 6, str(ana_ids))
    bad_needs = [a["id"] for a in bg_scmap.ANALYSES if not any(str(n).startswith("sc") for n in a.get("needs") or [])]
    check("所有分析 needs 含 sc", not bad_needs, str(bad_needs))

    ctx_ro = build_ctx(all_markers[:40])
    x_before = ctx_ro["X"].copy()
    genes_before = list(ctx_ro["genes"])
    samples_before = list(ctx_ro["samples"])
    de_before = len(ctx_ro["de"]["table"])
    run(ctx_ro, {"gene_set": "de_up", "score_method": "mean"}, "readonly")
    check("ctx.X 未被原地修改", np.allclose(ctx_ro["X"], x_before))
    check("ctx.genes/samples/de 未被原地修改",
          ctx_ro["genes"] == genes_before and ctx_ro["samples"] == samples_before
          and len(ctx_ro["de"]["table"]) == de_before)

    # ------------------------------------------------------------------
    # 6) 6 个 ANALYSIS 全部可跑
    # ------------------------------------------------------------------
    print("\n[6] 全部分析可跑")
    for ana in bg_scmap.ANALYSES:
        try:
            payload = ana["run"](build_ctx(all_markers[:80]),
                                 {"genes": markers[cell_types[0]][:6]},
                                 os.path.join(OUTDIR, f"ana-{ana['id']}"))
            tables = payload.get("tables") or []
            figures = payload.get("figures") or []
            ok_fig = all(os.path.isfile(f.get("file", "")) for f in figures)
            check(f"分析 {ana['id']} 跑通（≥1 表/图，图可读）",
                  len(tables) >= 1 and len(figures) >= 1 and ok_fig,
                  f"{len(tables)} 表 {len(figures)} 图")
        except Exception as exc:  # noqa: BLE001
            check(f"分析 {ana['id']} 跑通", False, f"{type(exc).__name__}: {exc}")
            traceback.print_exc()

    # ------------------------------------------------------------------
    # 汇总
    # ------------------------------------------------------------------
    print("\n" + "=" * 74)
    total = len(PASSED) + len(FAILED)
    if FAILED:
        print(f"FAILED  {len(FAILED)}/{total} 检查未通过：")
        for name in FAILED:
            print(f"  - {name}")
        print("=" * 74)
        return 1
    print(f"PASS  {total}/{total} 检查全部通过")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
