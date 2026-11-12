#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 测试 fixture 生成器（CONTRACT §11）

生成 12+ 种格式的 fixture 到临时目录（或指定目录），并写 `truth.json` 记录
每种 fixture 的期望维度 / 样本名 / 基因名 / 表达值类型 / 数值和，供
`test/py_import_test.py` 交叉校验。

用法：
    python3 test/make_fixtures.py            # 生成到临时目录并打印路径
    python3 test/make_fixtures.py /tmp/xxx   # 生成到指定目录

fixture 清单（核心 12 种）：
    csv（样本×基因）、tsv（基因×样本）、xlsx、gz、GEO series matrix、
    zip、tar.gz、rds、RData、matrix_text 转置、miRNA、circ
  额外（覆盖模块/边界）：
    lncRNA、fusion、rMATS splicing、无表头矩阵、TPM 值域、多样本表 join、
    大矩阵（inspect 内存）、坏文件（空/二进制/坏 zip）、多候选 zip、多对象 RData
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

R_CANDIDATES = [
    os.environ.get("DS_BG_RSCRIPT") or "",
    "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    "/usr/local/bin/Rscript",
    "/opt/homebrew/bin/Rscript",
]


def rscript_path() -> str | None:
    for candidate in R_CANDIDATES:
        if candidate and os.path.exists(candidate):
            return candidate
    return shutil.which("Rscript")


def _digest(values) -> str:
    arr = np.round(np.asarray(values, dtype=np.float64), 4).astype(np.float32)
    return hashlib.sha256(arr.tobytes()).hexdigest()[:16]


def _entry(path, X, genes, samples, fmt, value_kind="counts", **extra) -> dict:
    """把期望结果写进 truth。X 必须是 samples×genes。"""
    arr = np.asarray(X, dtype=np.float64)
    assert arr.shape == (len(samples), len(genes)), f"{path}: shape {arr.shape} != {(len(samples), len(genes))}"
    data = {
        "file": str(path),
        "format": fmt,
        "n_samples": int(arr.shape[0]),
        "n_genes": int(arr.shape[1]),
        "samples": [str(item) for item in samples],
        "genes": [str(item) for item in genes],
        "value_kind": value_kind,
        "sum": round(float(arr.sum()), 3),
        "mean": round(float(arr.mean()), 5),
        "max": round(float(arr.max()), 4),
        "digest": _digest(arr),
        "loader_options": {},
    }
    data.update(extra)
    return data


def _write_r_fixtures(outdir: Path, truth: dict) -> None:
    """用本机 R 生成 ExpressionSet 的 rds / RData（BioBase 已装）。"""
    rscript = rscript_path()
    if not rscript:
        truth["_r_available"] = False
        return
    try:
        subprocess.run([rscript, "--vanilla", "-e", "stopifnot(requireNamespace('Biobase', quietly=TRUE))"],
                       capture_output=True, text=True, timeout=120, check=True)
    except Exception:  # noqa: BLE001
        truth["_r_available"] = False
        return
    truth["_r_available"] = True
    rng = np.random.default_rng(2024)
    # 30 基因 × 8 样本 counts
    genes = ["ESR1", "TP53", "EGFR", "MYC", "BRCA1", "BRCA2", "CDH1", "VIM", "TWIST1", "SNAI1",
             "GAPDH", "ACTB", "B2M", "TBP", "HPRT1", "RPLP0", "GUSB", "TFRC", "LDHA", "PKM",
             "HK2", "PGK1", "ENO1", "ALDOA", "MKI67", "PCNA", "CCND1", "CDK4", "RB1", "PTEN"]
    samples = [f"P{i:02d}" for i in range(1, 9)]
    values = rng.integers(0, 500, size=(len(genes), len(samples))).astype(np.int64)
    tsv = outdir / "r_input.tsv"
    pd.DataFrame(values, index=genes, columns=samples).to_csv(tsv, sep="\t")
    rds = outdir / "exprset.rds"
    rdata = outdir / "eset.RData"
    two = outdir / "two_objects.RData"
    script = outdir / "_make_rds.R"
    script.write_text(
        """
suppressPackageStartupMessages({library(Biobase)})
mat <- as.matrix(read.table(commandArgs(trailingOnly=TRUE)[1], header=TRUE, row.names=1,
                            sep="\\t", check.names=FALSE))
groups <- rep(c("case", "control"), length.out = ncol(mat))
pdata <- data.frame(group = groups,
                    age = seq(45, 45 + 5 * (ncol(mat) - 1), by = 5),
                    sex = rep(c("M", "F"), length.out = ncol(mat)),
                    row.names = colnames(mat))
eset <- ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(pdata))
saveRDS(eset, commandArgs(trailingOnly=TRUE)[2])
save(eset, file = commandArgs(trailingOnly=TRUE)[3])
eset2 <- eset[1:12, 1:5]
save(eset, eset2, file = commandArgs(trailingOnly=TRUE)[4])
cat("r-fixtures-ok\\n")
""", encoding="utf-8")
    completed = subprocess.run([rscript, "--vanilla", str(script), str(tsv), str(rds), str(rdata), str(two)],
                               capture_output=True, text=True, timeout=300)
    if completed.returncode != 0:
        truth["_r_available"] = False
        truth["_r_error"] = (completed.stderr or "")[-500:]
        return
    truth["exprset.rds"] = _entry(rds, values.T, genes, samples, "rds", "counts",
                                  sample_info_columns=["group", "age", "sex"])
    truth["eset.RData"] = _entry(rdata, values.T, genes, samples, "rdata", "counts",
                                 sample_info_columns=["group", "age", "sex"])
    truth["two_objects.RData"] = {
        "file": str(two), "format": "rdata", "ambiguous": True,
        "candidates": ["eset", "eset2"],
    }


def build(outdir=None, verbose: bool = False) -> dict:
    """生成全部 fixture，返回 truth dict（含 `_dir`）。"""
    outdir = Path(outdir) if outdir else Path(tempfile.mkdtemp(prefix="bg-fixtures-"))
    outdir.mkdir(parents=True, exist_ok=True)
    truth: dict = {"_dir": str(outdir)}
    rng = np.random.default_rng(20241024)

    # ---------------------------------------------------------------- 1. csv 样本×基因
    samples1 = [f"TCGA-{tag}-01A" for tag in ("AA01", "AA02", "AA03", "BB01", "BB02", "BB03", "CC01", "CC02")]
    genes1 = [f"GENE{i}" for i in range(1, 31)]
    X1 = rng.integers(0, 400, size=(len(samples1), len(genes1)))
    path1 = outdir / "matrix_samples_by_genes.csv"
    pd.DataFrame(X1, index=samples1, columns=genes1).to_csv(path1)
    truth[path1.name] = _entry(path1, X1, genes1, samples1, "matrix_text", "counts")

    # ---------------------------------------------------------------- 2. tsv 基因×样本
    genes2 = ["TP53", "EGFR", "MYC", "BRCA1", "BRCA2", "CDH1", "VIM", "TWIST1", "SNAI1", "GAPDH",
              "ACTB", "B2M", "TBP", "HPRT1", "RPLP0", "GUSB", "TFRC", "LDHA", "PKM", "HK2",
              "PGK1", "ENO1", "ALDOA", "MKI67", "PCNA", "CCND1", "CDK4", "RB1", "PTEN", "AKT1",
              "MAPK1", "STAT3", "JUN", "FOS", "RELA", "NFKB1", "IL6", "TNF", "VEGFA", "MMP9"]
    samples2 = [f"Ctrl_{i}" for i in range(1, 5)] + [f"Case_{i}" for i in range(1, 4)]
    X2 = rng.integers(0, 800, size=(len(samples2), len(genes2)))
    path2 = outdir / "matrix_genes_by_samples.tsv"
    pd.DataFrame(X2.T, index=genes2, columns=samples2).to_csv(path2, sep="\t")
    truth[path2.name] = _entry(path2, X2, genes2, samples2, "matrix_text", "counts")

    # ---------------------------------------------------------------- 3. 空白分隔 txt
    genes3 = [f"PROBE_{i:03d}" for i in range(1, 21)]
    samples3 = [f"S{i}" for i in range(1, 6)]
    X3 = rng.integers(0, 200, size=(len(samples3), len(genes3))).astype(float)
    path3 = outdir / "matrix_space.txt"
    with open(path3, "w", encoding="utf-8") as handle:
        handle.write("gene " + " ".join(samples3) + "\n")
        for gene, row in zip(genes3, X3.T):
            handle.write(gene + " " + " ".join(str(int(value)) for value in row) + "\n")
    truth[path3.name] = _entry(path3, X3, genes3, samples3, "matrix_text", "counts")

    # ---------------------------------------------------------------- 4. gz
    genes4 = [f"G{i}" for i in range(1, 26)]
    samples4 = [f"GSM{300 + i}" for i in range(6)]
    X4 = rng.integers(0, 600, size=(len(samples4), len(genes4)))
    path4 = outdir / "matrix_gz.csv.gz"
    pd.DataFrame(X4, index=samples4, columns=genes4).to_csv(path4, compression="gzip")
    truth[path4.name] = _entry(path4, X4, genes4, samples4, "matrix_text", "counts")

    # ---------------------------------------------------------------- 5. GEO series matrix
    genes5 = [f"ILMN_{i:06d}" for i in range(1, 21)]
    gsm = [f"GSM{1000 + i}" for i in range(6)]
    titles = [f"Gastric cancer sample {i + 1}" for i in range(3)] + [f"Normal tissue {i + 1}" for i in range(3)]
    X5 = rng.normal(8.0, 1.2, size=(len(gsm), len(genes5)))
    values5 = np.round(X5, 3)
    path5 = outdir / "GSE99999_series_matrix.txt"
    geo_lines = [
        '!Series_title\t"Demo gastric cancer expression series"',
        '!Series_summary\t"Test series matrix fixture"',
        '!Sample_title\t' + "\t".join(f'"{title}"' for title in titles),
        '!Sample_geo_accession\t' + "\t".join(f'"{value}"' for value in gsm),
        '!Sample_source_name_ch1\t' + "\t".join(
            ['"gastric tumor"'] * 3 + ['"gastric normal"'] * 3),
        '!Sample_characteristics_ch1\t' + "\t".join(
            f'"age: {40 + i}"' for i in range(len(gsm))),
        '!Sample_characteristics_ch1\t' + "\t".join(
            ['"sex: M"', '"sex: F"', '"sex: M"', '"sex: F"', '"sex: M"', '"sex: F"']),
        "!series_matrix_table_begin",
        'ID_REF\t' + "\t".join(f'"{value}"' for value in gsm),
    ]
    for i, gene in enumerate(genes5):
        geo_lines.append(gene + "\t" + "\t".join(f"{value:.3f}" for value in values5[:, i]))
    geo_lines.append("!series_matrix_table_end")
    path5.write_text("\n".join(geo_lines) + "\n", encoding="utf-8")
    truth[path5.name] = _entry(path5, values5, genes5, gsm, "geo_series", "lognorm",
                               sample_info_columns=["title", "age", "sex", "source_name_ch1"])

    # ---------------------------------------------------------------- 6. zip
    genes6 = [f"GENE{i}" for i in range(1, 26)]
    samples6 = ["Ctrl_1", "Ctrl_2", "Ctrl_3", "Case_1", "Case_2", "Case_3"]
    X6 = rng.integers(0, 300, size=(len(samples6), len(genes6)))
    path6 = outdir / "bundle.zip"
    with zipfile.ZipFile(path6, "w") as archive:
        archive.writestr("data/matrix.csv", pd.DataFrame(X6, index=samples6, columns=genes6).to_csv())
        archive.writestr("README.txt", "demo bundle")
    truth[path6.name] = _entry(path6, X6, genes6, samples6, "archive", "counts")

    # ---------------------------------------------------------------- 7. 多候选 zip
    path_multizip = outdir / "bundle_multi.zip"
    genes7a = [f"GA{i}" for i in range(1, 11)]
    genes7b = [f"GB{i}" for i in range(1, 9)]
    samples7 = [f"S{i}" for i in range(1, 6)]
    X7a = rng.integers(0, 100, size=(len(samples7), len(genes7a)))
    X7b = rng.integers(0, 100, size=(len(samples7), len(genes7b)))
    with zipfile.ZipFile(path_multizip, "w") as archive:
        archive.writestr("first_matrix.csv", pd.DataFrame(X7a, index=samples7, columns=genes7a).to_csv())
        archive.writestr("second_matrix.csv", pd.DataFrame(X7b, index=samples7, columns=genes7b).to_csv())
    truth[path_multizip.name] = {
        "file": str(path_multizip), "format": "archive", "ambiguous": True,
        "member": "second_matrix.csv",
        "member_truth": _entry(path_multizip, X7b, genes7b, samples7, "matrix_text", "counts"),
    }

    # ---------------------------------------------------------------- 8. tar.gz
    genes8 = [f"TAR_GENE{i}" for i in range(1, 31)]
    samples8 = ["N1", "N2", "T1", "T2", "T3"]
    X8 = rng.integers(0, 500, size=(len(samples8), len(genes8)))
    inner = outdir / "_tar_inner"
    inner.mkdir(exist_ok=True)
    inner_file = inner / "matrix.tsv"
    pd.DataFrame(X8.T, index=genes8, columns=samples8).to_csv(inner_file, sep="\t")
    path8 = outdir / "bundle.tar.gz"
    with tarfile.open(path8, "w:gz") as archive:
        archive.add(inner_file, arcname="matrix.tsv")
    truth[path8.name] = _entry(path8, X8, genes8, samples8, "archive", "counts")

    # ---------------------------------------------------------------- 9. 无表头
    genes9 = [f"ENSG{i:011d}" for i in range(1, 13)]
    samples9 = [f"sample_{i}" for i in range(1, 7)]
    X9 = rng.integers(0, 90, size=(len(samples9), len(genes9)))
    path9 = outdir / "matrix_no_header.tsv"
    with open(path9, "w", encoding="utf-8") as handle:
        for gene, row in zip(genes9, X9.T):
            handle.write(gene + "\t" + "\t".join(str(int(value)) for value in row) + "\n")
    truth[path9.name] = _entry(path9, X9, genes9, samples9, "matrix_text", "counts",
                               auto_generated_names=True)

    # ---------------------------------------------------------------- 10. TPU 值域（TPM）
    genes10 = [f"TPM{i}" for i in range(1, 19)]
    samples10 = ["P1", "P2", "P3", "P4", "P5", "P6"]
    X10 = np.round(rng.gamma(2.0, 300.0, size=(len(samples10), len(genes10))), 3)
    path10 = outdir / "matrix_tpm.csv"
    pd.DataFrame(X10, index=samples10, columns=genes10).to_csv(path10)
    truth[path10.name] = _entry(path10, X10, genes10, samples10, "matrix_text", "tpm")

    # ---------------------------------------------------------------- 11. 强制转置（样本×基因 → 转置）
    genes11 = [f"TG{i}" for i in range(1, 21)]
    samples11 = [f"QB{i}" for i in range(1, 9)]
    X11 = rng.normal(6.0, 0.8, size=(len(samples11), len(genes11)))
    path11 = outdir / "matrix_transpose_forced.tsv"
    pd.DataFrame(X11, index=samples11, columns=genes11).to_csv(path11, sep="\t")
    # transpose=True ⇒ 把「行」当特征，「列」当样本
    truth[path11.name] = _entry(path11, X11.T, samples11, genes11, "matrix_text", "lognorm",
                                loader_options={"transpose": True},
                                note="显式 transpose=True：行=特征(8)，列=样本(20)")

    # ---------------------------------------------------------------- 12. xlsx（多 sheet）
    genes12 = [f"XG{i}" for i in range(1, 31)]
    samples12 = [f"S{i}" for i in range(1, 7)]
    X12 = rng.integers(0, 200, size=(len(samples12), len(genes12)))
    path12 = outdir / "expr_multi_sheet.xlsx"
    with pd.ExcelWriter(path12, engine="openpyxl") as writer:
        pd.DataFrame(X12.T, index=genes12, columns=samples12).to_excel(writer, sheet_name="expr")
        pd.DataFrame({"note": ["demo", "fixture"]}).to_excel(writer, sheet_name="notes", index=False)
    truth[path12.name] = _entry(path12, X12, genes12, samples12, "xlsx", "counts",
                                loader_options={"sheet": "expr"}, multiple_sheets=True)

    # ---------------------------------------------------------------- 13. 样本信息表 / 临床表
    sheet_path = outdir / "sample_sheet.csv"
    pd.DataFrame({
        "sample": samples1[:5],
        "group": ["case", "case", "case", "control", "control"],
        "patient": [f"PT{i}" for i in range(1, 6)],
    }).to_csv(sheet_path, index=False)
    clinical_path = outdir / "clinical.csv"
    pd.DataFrame({
        "GSM": gsm,
        "os_time": [12, 30, 8, 60, 45, 20],
        "os_event": [1, 0, 1, 0, 1, 0],
        "age": [40, 41, 42, 43, 44, 45],
    }).to_csv(clinical_path, index=False)
    truth["sample_sheet.csv"] = {"file": str(sheet_path), "n_rows": 5}
    truth["clinical.csv"] = {"file": str(clinical_path), "n_rows": 6}

    # ---------------------------------------------------------------- 14. 非编码 / 融合 / 剪接
    mirna_genes = ["hsa-miR-21-5p", "hsa-let-7a-5p", "hsa-miR-155-5p", "hsa-miR-145-5p",
                   "hsa-miR-203a-3p", "hsa-miR-34a-5p", "hsa-miR-125b-5p", "hsa-miR-9-5p",
                   "hsa-miR-200c-3p", "hsa-miR-141-3p", "hsa-miR-10b-5p", "hsa-miR-31-5p",
                   "hsa-mir-21", "hsa-mir-155"]
    mirna_samples = [f"M{i}" for i in range(1, 9)]
    X13 = rng.integers(0, 900, size=(len(mirna_samples), len(mirna_genes)))
    path13 = outdir / "mirna_matrix.tsv"
    pd.DataFrame(X13.T, index=mirna_genes, columns=mirna_samples).to_csv(path13, sep="\t")
    truth[path13.name] = _entry(path13, X13, mirna_genes, mirna_samples, "mirna", "counts")

    circ_genes = [f"chr{(i % 6) + 1}:{1000 * i}-{1000 * i + 250}:+" for i in range(1, 21)] + \
                 [f"hsa_circ_{i:06d}" for i in range(1, 6)]
    circ_samples = [f"CS{i}" for i in range(1, 7)]
    X14 = rng.integers(0, 150, size=(len(circ_samples), len(circ_genes)))
    path14 = outdir / "circ_matrix.tsv"
    pd.DataFrame(X14.T, index=circ_genes, columns=circ_samples).to_csv(path14, sep="\t")
    truth[path14.name] = _entry(path14, X14, circ_genes, circ_samples, "circ", "counts")

    lncrna_genes = ["LINC00152", "LINC00662", "PVT1", "HOTAIR", "H19", "MALAT1", "NEAT1",
                    "XIST", "KCNQ1OT1", "HULC", "CASC9", "LINC01234", "SOX2-OT", "CCAT1",
                    "GAS5", "MEG3", "LINC-PINT", "TUG1", "SNHG1", "DANCR"]
    lncrna_samples = [f"L{i}" for i in range(1, 7)]
    X15 = rng.integers(0, 700, size=(len(lncrna_samples), len(lncrna_genes)))
    path15 = outdir / "lncrna_matrix.tsv"
    pd.DataFrame(X15.T, index=lncrna_genes, columns=lncrna_samples).to_csv(path15, sep="\t")
    truth[path15.name] = _entry(path15, X15, lncrna_genes, lncrna_samples, "lncrna", "counts")

    path16 = outdir / "fusion_starfusion.tsv"
    pd.DataFrame({
        "#FusionName": ["BCR--ABL1", "EML4--ALK", "TPM3--NTRK1", "FGFR3--TACC3", "RET--CCDC6"],
        "JunctionReadCount": [12, 5, 0, 3, 9],
        "SpanningFragCount": [8, 3, 1, 2, 6],
        "LeftGene": ["BCR", "EML4", "TPM3", "FGFR3", "RET"],
        "RightGene": ["ABL1", "ALK", "NTRK1", "TACC3", "CCDC6"],
    }).to_csv(path16, sep="\t", index=False)
    truth[path16.name] = {"file": str(path16), "format": "fusion", "n_events": 5,
                          "expected_pairs": ["BCR--ABL1", "EML4--ALK"]}

    path17 = outdir / "SE.MATS.JC.txt"
    pd.DataFrame({
        "ID": [f"SE.{i}" for i in range(1, 7)],
        "GeneID": [f"G{i}" for i in range(1, 7)],
        "GeneSymbol": ["A", "B", "C", "D", "E", "F"],
        "IJC_SAMPLE_1": ["10,12,9", "3,4,2", "0,1,0", "5,6,7", "8,9,10", "2,3,1"],
        "SJC_SAMPLE_1": ["20,22,19", "4,5,3", "1,2,1", "15,16,17", "18,19,20", "12,13,11"],
        "IncLevel1": ["0.8,0.85,0.82", "0.3,0.35,0.31", "0.1,0.05,0.08",
                      "0.6,0.62,0.61", "0.9,0.88,0.91", "0.25,0.28,0.26"],
        "IJC_SAMPLE_2": ["2,3,1", "8,9,7", "5,6,4", "9,10,8", "1,2,1", "6,7,5"],
        "SJC_SAMPLE_2": ["20,22,19", "10,11,9", "6,7,5", "19,20,18", "11,12,10", "16,17,15"],
        "IncLevel2": ["0.1,0.12,0.11", "0.7,0.72,0.71", "0.5,0.52,0.49",
                      "0.3,0.32,0.31", "0.08,0.09,0.07", "0.35,0.38,0.36"],
    }).to_csv(path17, sep="\t", index=False)
    truth[path17.name] = {"file": str(path17), "format": "splicing", "n_events": 6,
                          "n_samples": 6, "event_type": "SE"}

    # ---------------------------------------------------------------- 15. 大矩阵（inspect 内存）
    genes_big = [f"BG{i:05d}" for i in range(1, 20001)]
    samples_big = [f"BS{i}" for i in range(1, 41)]
    rng_big = np.random.default_rng(7)
    X_big = rng_big.integers(0, 100, size=(len(samples_big), len(genes_big)))
    path_big = outdir / "matrix_large.csv"
    pd.DataFrame(X_big, index=samples_big, columns=genes_big).to_csv(path_big)
    truth[path_big.name] = {"file": str(path_big), "format": "matrix_text", "large": True,
                            "inspect_only": True, "n_full_samples": len(samples_big),
                            "n_full_genes": len(genes_big)}

    # ---------------------------------------------------------------- 16. 坏文件
    truth["bad_binary.csv"] = {"file": str(outdir / "bad_binary.csv"), "broken": True}
    (outdir / "bad_binary.csv").write_bytes(b"\x00\x01\x02\x03\xff\xfe" * 64)
    truth["empty.csv"] = {"file": str(outdir / "empty.csv"), "broken": True}
    (outdir / "empty.csv").write_text("", encoding="utf-8")
    truth["bad_archive.zip"] = {"file": str(outdir / "bad_archive.zip"), "broken": True}
    (outdir / "bad_archive.zip").write_bytes(b"this is not a zip file")

    # ---------------------------------------------------------------- 17. R fixtures
    _write_r_fixtures(outdir, truth)

    # 清理临时目录
    for temp_dir in (inner,):
        shutil.rmtree(temp_dir, ignore_errors=True)
    for temp_file in (outdir / "_make_rds.R", outdir / "r_input.tsv", outdir / "_make_rds.R"):
        try:
            Path(temp_file).unlink(missing_ok=True)
        except Exception:  # noqa: BLE001
            pass

    truth_path = outdir / "truth.json"
    truth_path.write_text(json.dumps(truth, ensure_ascii=False, indent=2), encoding="utf-8")
    truth["_truth_path"] = str(truth_path)
    if verbose:
        print(f"fixtures -> {outdir}")
        print(f"truth    -> {truth_path}")
        print(f"条目数：{len([k for k in truth if not k.startswith('_')])}")
    return truth


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else None
    build(outdir=target, verbose=True)
