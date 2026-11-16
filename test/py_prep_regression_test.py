#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 标准化回归测试（任务单 GEO-07）

  * `method=none,log_transform=no` 时矩阵不变，但 value_kind 必须保持真实尺度
    （修复前无论是否转换都标成 lognorm，下游会把原始计数差当成 logFC）。
  * `auto` 对已标准化（lognorm）输入默认不再改动（修复前会做中位数中心化）。
  * 分位数标准化必须让相同值得到相同结果（ties 均值），并与 preprocessCore 对照。

用法：cd plugins/dsh-geobulk && python3 test/py_prep_regression_test.py
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

import bg_prep  # noqa: E402

OUTDIR = tempfile.mkdtemp(prefix="bg-prep-reg-")
PASSED: list[str] = []
FAILED: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def run(ctx, params, tag):
    out = os.path.join(OUTDIR, tag)
    os.makedirs(out, exist_ok=True)
    return bg_prep.run_normalize(ctx, params, out)


def counts_ctx(seed=3, n=6, g=40):
    rng = np.random.default_rng(seed)
    X = rng.poisson(80, (n, g)).astype(float)
    return {
        "X": X, "genes": [f"G{i}" for i in range(g)],
        "samples": [f"S{i}" for i in range(n)], "X_kind": "counts",
        "sample_info": pd.DataFrame(index=[f"S{i}" for i in range(n)]),
    }


def test_value_kind():
    print("\n[GEO-07] 尺度标记 / auto 幂等")
    ctx = counts_ctx()
    out = run(ctx, {"method": "none", "log_transform": "no"}, "g07-none")
    replaced = out["replace"]
    check("GEO-07a none+no 矩阵不变", np.allclose(replaced["X"], ctx["X"], atol=1e-6))
    check("GEO-07b none+no 的 value_kind 仍是 counts（修复前写成 lognorm）",
          replaced["value_kind"] == "counts", replaced["value_kind"])
    check("GEO-07c 明确警告 logFC 不是 log2 尺度",
          any("不是 log2 尺度" in w for w in out.get("warnings") or []),
          str(out.get("warnings"))[:80])
    check("GEO-07d metrics 记录前后尺度",
          out["summary"]["metrics"].get("value_kind_before") == "counts" and
          out["summary"]["metrics"].get("value_kind_after") == "counts")

    out = run(ctx, {"method": "cpm"}, "g07-cpm")
    check("GEO-07e cpm → lognorm", out["replace"]["value_kind"] == "lognorm",
          out["replace"]["value_kind"])

    # 已 lognorm 的输入，auto 不得再改
    log_ctx = dict(ctx)
    log_ctx["X"] = np.log2(ctx["X"] + 1.0)
    log_ctx["X_kind"] = "lognorm"
    out = run(log_ctx, {"method": "auto"}, "g07-auto-lognorm")
    check("GEO-07f auto 对已 lognorm 输入默认不改动",
          np.allclose(out["replace"]["X"], log_ctx["X"], atol=1e-6),
          f"method={out['summary']['metrics'].get('method')}")
    check("GEO-07g auto 对已 lognorm 选择 method=none",
          out["summary"]["metrics"].get("method") == "none",
          str(out["summary"]["metrics"].get("method")))


def test_quantile_ties():
    print("\n[GEO-07] 分位数 ties")
    X = np.array([[0.0, 0.0, 0.0, 10.0],
                  [1.0, 2.0, 3.0, 4.0]], dtype=float)
    out = bg_prep._quantile_normalize(X)
    check("GEO-07h 相同的 0 标准化后仍相同",
          len(set(np.round(out[0][:3], 9))) == 1, str(out[0]))
    check("GEO-07i 不再出现 0.5/1/1.5 分裂", not np.allclose(sorted(out[0][:3]), [0.5, 1.0, 1.5]))

    # 与 preprocessCore 对照（如果有 R）
    try:
        ok_r = subprocess.run(
            ["Rscript", "--vanilla", "-e",
             'cat(requireNamespace("preprocessCore", quietly=TRUE))'],
            capture_output=True, text=True, timeout=60,
        ).stdout.strip() == "TRUE"
    except Exception:  # noqa: BLE001
        ok_r = False
    if not ok_r:
        print("  SKIP preprocessCore 对照（本机没有 R/preprocessCore）")
        return
    rng = np.random.default_rng(7)
    mat = rng.poisson(20, (5, 30)).astype(float)
    mat[0, :10] = 0.0          # 大量 ties
    mat[1, 5:15] = 0.0
    mat[2, :] = np.round(mat[2])
    ref = _preprocess_core_quantile(mat)
    ours = bg_prep._quantile_normalize(mat)
    diff = float(np.max(np.abs(ref - ours)))
    check("GEO-07j 与 preprocessCore::normalize.quantiles 一致", diff < 1e-8, f"max|diff|={diff:.2e}")


def _preprocess_core_quantile(mat: np.ndarray) -> np.ndarray:
    tmp = Path(tempfile.mkdtemp(prefix="qnorm-ref-"))
    # preprocessCore 按列（=样本）标准化，我们的矩阵是 样本×基因，所以要转置传入
    pd.DataFrame(mat.T).to_csv(tmp / "in.csv", index=False, header=False)
    script = tmp / "ref.R"
    script.write_text(
        """
args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
m <- as.matrix(read.csv(file.path(d, "in.csv"), header=FALSE, check.names=FALSE))
suppressPackageStartupMessages(library(preprocessCore))
out <- preprocessCore::normalize.quantiles(m)
write.table(out, file.path(d, "out.csv"), sep=",", row.names=FALSE, col.names=FALSE)
""",
        encoding="utf-8",
    )
    proc = subprocess.run(["Rscript", "--vanilla", str(script), str(tmp)],
                          capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        raise RuntimeError(f"preprocessCore 调用失败：{proc.stderr[-800:]}")
    return pd.read_csv(tmp / "out.csv", header=None).values.astype(float).T


def main() -> int:
    print("== dsh-geobulk 标准化回归（GEO-07）==")
    test_value_kind()
    test_quantile_ties()
    print(f"\n== 结果：{len(PASSED)} 通过 / {len(FAILED)} 失败 ==")
    if FAILED:
        for item in FAILED:
            print(f"  ✗ {item}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
