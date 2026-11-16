#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — GSEA 回归测试（任务单 GEO-08）

  * 有 R + fgsea 时，结果必须来自 `fgsea::fgseaMultilevel`，并记录 actual_method /
    版本 / 置换策略 / 覆盖率；ES 必须与直接调用 fgsea 一致。
  * 自实现分支的 NES 必须是官方定义（按符号分别用同号置换 ES 均值归一化），
    单侧 p 值分母是同号置换次数，而不是全部置换次数。

用法：cd plugins/dsh-geobulk && python3 test/py_gsea_regression_test.py
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

import bg_enrich  # noqa: E402
import rbridge  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-gsea-reg-")
PASSED: list[str] = []
FAILED: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def _de_fixture(n_genes=400, seed=3):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    logfc = rng.normal(0, 1, n_genes)
    logfc[:30] += 3.0     # 前 30 个基因上调
    logfc[30:60] -= 3.0   # 接下来 30 个下调
    de = pd.DataFrame({
        "gene": genes, "logFC": logfc,
        "AveExpr": rng.normal(7, 1, n_genes),
        "t": logfc * 3, "P.Value": np.full(n_genes, 1e-4),
        "adj.P.Val": np.full(n_genes, 1e-3), "B": np.zeros(n_genes),
        "regulation": np.where(logfc > 0, "up", "down"),
    })
    sets = {
        "UP_SET": genes[:30],
        "DOWN_SET": genes[30:60],
        "RANDOM_SET": list(rng.choice(genes[60:], 25, replace=False)),
    }
    names = {k: k.title() for k in sets}
    return de, sets, names


def test_official_fgsea():
    print("\n[GEO-08] 官方 fgsea 分支")
    if not (rbridge.r_available() and rbridge.has_pkg("fgsea")):
        print("  SKIP 本机没有 R + fgsea")
        return
    de, sets, names = _de_fixture()
    res = bg_enrich.gsea(de, sets, names, "TEST", outdir=os.path.join(OUTDIR, "fgsea"))
    check("GEO-08a 有结果", res is not None and not res.empty, f"rows={0 if res is None else len(res)}")
    prov = getattr(res, "attrs", {}).get("provenance") or {}
    check("GEO-08b actual_method 记录为 fgsea",
          "fgsea" in str(prov.get("actual_method")), str(prov.get("actual_method")))
    check("GEO-08c 记录软件版本", bool(prov.get("software_version")), str(prov.get("software_version")))
    check("GEO-08d 记录置换策略", bool(prov.get("permutation")), str(prov.get("permutation"))[:60])
    check("GEO-08e 记录覆盖率", bool((prov.get("coverage") or {}).get("n_ranked")),
          str(prov.get("coverage")))

    # ES 与直接调用 fgsea 对照
    ref = _direct_fgsea(de, sets)
    merged = res.set_index("ID")[["ES", "NES", "pvalue"]].join(
        ref.set_index("pathway")[["ES", "NES", "pvalue"]], rsuffix="_ref", how="inner"
    )
    if merged.empty:
        check("GEO-08f ES 与直接调用 fgsea 一致", False, "没有可对照的集合")
    else:
        d_es = float(np.max(np.abs(merged["ES"] - merged["ES_ref"])))
        check("GEO-08f ES 与直接调用 fgsea 一致", d_es < 1e-8, f"max|dES|={d_es:.2e}")
        d_nes = float(np.max(np.abs(merged["NES"] - merged["NES_ref"])))
        check("GEO-08g NES 与直接调用 fgsea 一致（同一算法）", d_nes < 1e-6, f"max|dNES|={d_nes:.2e}")
    up = res.set_index("ID").loc["UP_SET"]
    down = res.set_index("ID").loc["DOWN_SET"]
    check("GEO-08h 上调集 ES>0、下调集 ES<0", float(up["ES"]) > 0 and float(down["ES"]) < 0,
          f"up={float(up['ES']):.3f} down={float(down['ES']):.3f}")


def _direct_fgsea(de, sets):
    tmp = Path(tempfile.mkdtemp(prefix="fgsea-ref-"))
    ranked = de[["gene", "logFC"]].copy()
    ranked.to_csv(tmp / "ranked.csv", index=False)
    rows = [{"term": k, "gene": g} for k, members in sets.items() for g in members]
    pd.DataFrame(rows).to_csv(tmp / "sets.csv", index=False)
    script = tmp / "ref.R"
    script.write_text(
        """
args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
r <- read.csv(file.path(d, "ranked.csv"), stringsAsFactors=FALSE)
s <- read.csv(file.path(d, "sets.csv"), stringsAsFactors=FALSE)
stats <- setNames(r$logFC, r$gene)
paths <- split(s$gene, s$term)
suppressPackageStartupMessages(library(fgsea))
set.seed(0)
res <- fgsea::fgseaMultilevel(pathways=paths, stats=stats, minSize=15, maxSize=500, nPermSimple=1000, nproc=1)
write.csv(data.frame(pathway=res$pathway, ES=res$ES, NES=res$NES, pvalue=res$pval),
          file.path(d, "ref.csv"), row.names=FALSE)
""",
        encoding="utf-8",
    )
    proc = subprocess.run(["Rscript", "--vanilla", str(script), str(tmp)],
                          capture_output=True, text=True, timeout=900)
    if proc.returncode != 0:
        raise RuntimeError(f"参考 fgsea 调用失败：{proc.stderr[-800:]}")
    return pd.read_csv(tmp / "ref.csv")


def test_python_fallback_statistic():
    print("\n[GEO-08] 自实现分支的 NES 定义")
    de, sets, names = _de_fixture()
    res = bg_enrich.gsea(de, sets, names, "TEST", use_r=False, nperm=400, seed=1)
    prov = getattr(res, "attrs", {}).get("provenance") or {}
    check("GEO-08i 自实现明确标注 actual_method=gsea_python",
          prov.get("actual_method") == "gsea_python", str(prov.get("actual_method")))
    check("GEO-08j 不再冒称官方 GSEA/NES",
          "自实现" in str(prov.get("permutation")) or "自实现" in str(prov.get("software")),
          str(prov.get("software")))

    # 重新计算一次同号归一化，验证实现与声明一致
    frame = de.dropna(subset=["logFC"]).copy()
    frame["gene"] = frame["gene"].astype(str).str.upper()
    frame = frame.drop_duplicates(subset=["gene"]).sort_values("logFC", ascending=False)
    ranked = frame["gene"].tolist()
    weights = frame["logFC"].values
    rng = np.random.default_rng(1)
    pool = np.array(ranked)
    for key in ("UP_SET", "DOWN_SET"):
        members = {str(g).upper() for g in sets[key]}
        size = len(members & set(ranked))
        es = bg_enrich._es_score(ranked, members, weights)
        null = np.array([
            bg_enrich._es_score(ranked, set(rng.choice(pool, size=size, replace=False).tolist()), weights)
            for _ in range(400)
        ])
        row = res.set_index("ID").loc[key]
        if es >= 0:
            same = null[null >= 0]
            expected_nes = es / float(np.mean(same)) if same.size and np.mean(same) > 0 else 0.0
            expected_p = (int(np.sum(same >= es)) + 1) / (same.size + 1)
        else:
            same = null[null < 0]
            expected_nes = es / float(-np.mean(same)) if same.size and np.mean(same) < 0 else 0.0
            expected_p = (int(np.sum(same <= es)) + 1) / (same.size + 1)
        check(f"GEO-08k {key} NES = ES/同号置换均值",
              abs(float(row["NES"]) - expected_nes) < 1e-9,
              f"{float(row['NES']):.4f} vs {expected_nes:.4f}")
        check(f"GEO-08l {key} p 值分母为同号置换次数",
              abs(float(row["pvalue"]) - expected_p) < 1e-9,
              f"{float(row['pvalue']):.4f} vs {expected_p:.4f}")


def main() -> int:
    print("== dsh-geobulk GSEA 回归（GEO-08）==")
    test_official_fgsea()
    test_python_fallback_statistic()
    print(f"\n== 结果：{len(PASSED)} 通过 / {len(FAILED)} 失败 ==")
    if FAILED:
        for item in FAILED:
            print(f"  ✗ {item}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
