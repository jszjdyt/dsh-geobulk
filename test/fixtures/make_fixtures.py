"""Generate deterministic fixtures for bulk-ncrna tests (scmap + ncrna).

Creates (all under this directory, idempotent unless force=True):
  demo_sc.h5ad        400 cells x 800 genes x 6 cell types (sparse CSR X)
  sc_truth.json       markers / shared DE genes / expected mapping (truth)
  sample_info.tsv     S1..S12 -> case/control
  rmats_SE.txt        rMATS SE.MATS.JC style table (3 vs 3 replicates)
  arriba_fusions.tsv  Arriba style fusion table
  star_fusion.tsv     STAR-Fusion style fusion table
  psi_matrix.tsv      SUPPA style PSI matrix (mixed event types)
  circ_expr.tsv       60 circRNAs x 12 samples (+ linear column for ratio)
  mirna_expr.tsv      40 miRNAs x 12 samples (hsa-miR-* / hsa-mir-* precursor)
  lncrna_expr.tsv     60 lncRNAs x 12 samples (LINC* / *-AS1)

All data are SYNTHETIC (no real patient data). Real gene symbols are used for
markers so that demo mapping is interpretable.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261112

SAMPLES = [f"S{i:02d}" for i in range(1, 13)]
GROUPS = ["case"] * 6 + ["control"] * 6

# 细胞类型对齐 bg_state.make_demo 的 bulk demo 基因（增殖/T/B/髓系/上皮/成纤维），
# 这样 bulk DE 上调→Proliferating、下调(免疫)→T/B 细胞的映射结果可解释。
CELL_TYPES = ["Proliferating", "T cell", "B cell", "Myeloid", "Epithelial", "Fibroblast"]

# 22 candidate markers per cell type (dedupe + pad to 20 programmatically).
_MARKER_CANDIDATES = {
    "Proliferating": ["MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB",
                      "PLK1", "BUB1", "BUB1B", "TYMS", "RRM2", "PCNA", "MCM2", "MCM3",
                      "MCM4", "MCM5", "MCM6", "MCM7", "CCNE1", "CCNE2", "CDK2",
                      "E2F1", "E2F3"],
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "IL7R", "GZMA", "GZMB",
               "PRF1", "NKG7", "GNLY", "KLRD1", "CD27", "CD28", "ICOS", "LTB",
               "LCK", "ZAP70", "TRAC", "CCR7", "LEF1"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD19", "IGHM", "IGHG1", "IGKC", "MZB1",
               "JCHAIN", "CD74", "HLA-DRA", "PAX5", "EBF1", "CD22", "CD37", "CD53",
               "CD40", "FCER2", "BANK1", "BLNK", "CD83", "CR2"],
    "Myeloid": ["CD14", "LYZ", "FCGR3A", "FCGR1A", "S100A8", "S100A9", "S100A12",
                "CD68", "CSF1R", "TLR4", "IL1B", "CCL2", "FCN1", "VCAN", "MNDA",
                "SERPINA1", "CEBPB", "CD86", "MSR1", "LILRB1", "CTSB", "CD163"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "CLDN4",
                   "CLDN7", "MUC1", "ELF3", "ESRP1", "GRHL2", "SFN", "TJP3",
                   "KRT15", "KRT5", "LAMB3", "LAMC2", "ITGB4", "S100A14", "KRT13",
                   "CD24"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1",
                   "COL6A2", "DCN", "LUM", "FBN1", "FN1", "SPARC", "TIMP1", "TIMP2",
                   "VIM", "THY1", "PDGFRA", "PDGFRB", "ACTA2", "POSTN", "TAGLN",
                   "MFAP5"],
}

# Other real symbols used to fill the gene universe (keeps the h5ad realistic).
# 前两行刻意覆盖 bg_state.make_demo 的 UP_INTERFERON / UP_REPAIR 基因，方便 demo 联动。
_EXTRA_REAL = [
    # demo UP_INTERFERON / UP_REPAIR（非 marker，但要在参考里存在）
    "STAT1", "STAT2", "IFNG", "CXCL9", "CXCL10", "CXCL11", "ISG15", "IFI44L",
    "IFIT1", "IFIT3", "OAS1", "OAS2", "MX1", "MX2", "BRCA1", "BRCA2", "RAD51",
    "ATM", "ATR", "CHEK1", "CHEK2", "PALB2", "FANCA", "FANCD2", "MLH1", "MSH2",
    "MSH6",
    "ACTB", "GAPDH", "B2M", "TMSB10", "TMSB4X", "HLA-A", "HLA-B", "HLA-C",
    "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-E", "HLA-F", "CD52", "CD48",
    "CD44", "CD47", "CD59", "CD63", "CD81", "CD9", "CD99", "CD164", "CD320",
    "RPLP0", "RPL13A", "RPS18", "RPS27", "RPL10", "RPS6", "EEF1A1", "EEF1B2",
    "EIF4A1", "EIF5B", "MALAT1", "NEAT1", "XIST", "TUG1", "HOTAIR", "PVT1",
    "MEG3", "NORAD", "DANCR", "H19", "KCNQ1OT1", "SNHG1", "SNHG5", "SNHG6",
    "UBB", "UBC", "UBA52", "UBA2", "PSMA1", "PSMB1", "PSMC1", "ATP5F1A",
    "ATP5F1B", "COX4I1", "COX5A", "NDUFA1", "NDUFB1", "SDHA", "SDHB",
    "LDHA", "LDHB", "PKM", "ENO1", "PGK1", "PPIA", "PPIB", "SRP14", "SSR4",
    "CALR", "CANX", "HSP90AA1", "HSPA1A", "HSPA8", "HSPB1", "HSPD1", "HSPE1",
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "IER2", "IER3", "DUSP1",
    "NFKBIA", "TNFAIP3", "ZFP36", "IER5", "ATF3", "KLF6", "KLF4", "MYC",
    "GADD45B", "BTG2", "SOCS3", "RHOB", "RGS1", "RGS2", "RGS10", "CXCR4",
    "CCR2", "CCR5", "CXCL10", "CXCL12", "CCL3", "CCL4", "CCL5", "CXCL9",
    "IL10", "IL6", "IL1RN", "TNF", "TGFB1", "TGFBR2", "IFNGR1", "IFNGR2",
    "STAT1", "STAT3", "STAT6", "IRF1", "IRF7", "IRF8", "IRF9", "JAK1",
    "JAK2", "TYK2", "MYD88", "TRAF6", "NFKB1", "RELA", "RELB", "NFATC1",
    "NFATC2", "RUNX1", "RUNX3", "GATA3", "TBX21", "RORC", "FOXP3", "PTPRC",
    "PTPN22", "PTPN11", "PTP4A1", "SELL", "SELPLG", "SELP", "SELE", "ITGAL",
    "ITGAM", "ITGB1", "ITGB2", "ITGA4", "ITGA5", "ITGAV", "ITGB5", "ICAM1",
    "VCAM1", "MADCAM1", "PECAM1", "VWF", "KDR", "FLT1", "TEK", "CDH5",
    "ENG", "ESAM", "TIE1", "ANGPT1", "ANGPT2", "PLVAP", "CLDN5", "AQP1",
    "MKI67", "TOP2A", "PCNA", "CCNB1", "CCND1", "CCNE1", "CDC20", "BUB1",
    "AURKA", "AURKB", "PLK1", "TP53", "CDKN1A", "CDKN1B", "CDKN2A", "RB1",
    "BCL2", "BCL2L1", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "APAF1",
    "CYCS", "DIABLO", "XIAP", "BIRC5", "MCL1", "BID", "BAD", "BCL2A1",
    "COL4A1", "COL4A2", "COL7A1", "COL8A1", "COL11A1", "COL12A1", "COL14A1",
    "COL15A1", "COL16A1", "COL17A1", "COL18A1", "ELN", "FBLN1", "FBLN2",
    "FBLN5", "LTBP1", "LTBP2", "EMILIN1", "MMP1", "MMP2", "MMP3", "MMP7",
    "MMP9", "MMP10", "MMP11", "MMP14", "TIMP3", "TIMP4", "SERPINE1",
    "SERPINE2", "PLAT", "PLAU", "PLAUR", "THBS1", "THBS2", "CTGF", "CYR61",
    "WISP1", "WISP2", "LOX", "LOXL1", "LOXL2", "TGFBI", "SERPINF1",
]

SYNTH_PREFIX = "GENE"


def _build_gene_panel():
    used = set()
    markers = {}
    for ct in CELL_TYPES:
        picked = []
        for g in _MARKER_CANDIDATES[ct]:
            if g not in used:
                picked.append(g)
                used.add(g)
            if len(picked) == 20:
                break
        assert len(picked) == 20, (ct, len(picked))
        markers[ct] = picked
    for g in _EXTRA_REAL:
        if g not in used:
            used.add(g)
    return markers, used


def build_demo_sc(path: str):
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(SEED)
    markers, _used = _build_gene_panel()

    all_real = []
    for ct in CELL_TYPES:
        all_real.extend(markers[ct])
    all_real.extend(_EXTRA_REAL)
    # dedupe keeping order
    real = []
    seen = set()
    for g in all_real:
        if g not in seen:
            real.append(g)
            seen.add(g)
    n_target = 800
    n_synth = n_target - len(real)
    assert n_synth > 0
    synth = [f"{SYNTH_PREFIX}{i:04d}" for i in range(1, n_synth + 1) if f"{SYNTH_PREFIX}{i:04d}" not in seen]
    genes = real + synth[:n_synth]
    assert len(genes) == n_target, (len(genes), n_target)
    gene_idx = {g: i for i, g in enumerate(genes)}

    n_cells = 400
    types = []
    per = n_cells // len(CELL_TYPES)
    for ct in CELL_TYPES:
        types.extend([ct] * per)
    types.extend([CELL_TYPES[-1]] * (n_cells - len(types)))
    types = np.array(types, dtype=object)
    rng.shuffle(types)

    # baseline expression: gamma, then marker boost inside matching cell type
    base = rng.gamma(shape=1.1, scale=0.8, size=(n_cells, len(genes)))
    marker_rows = {ct: np.where(types == ct)[0] for ct in CELL_TYPES}
    for ct in CELL_TYPES:
        rows = marker_rows[ct]
        cols = np.array([gene_idx[g] for g in markers[ct]])
        boost = rng.uniform(4.0, 8.0, size=(len(rows), 20))
        base[np.ix_(rows, cols)] *= boost
    # a few weak cross-type signals so the map is not trivially perfect
    for ct in CELL_TYPES:
        rows = marker_rows[ct]
        others = [g for c2 in CELL_TYPES if c2 != ct for g in markers[c2][:3]]
        cols = np.array([gene_idx[g] for g in others])
        base[np.ix_(rows, cols)] *= rng.uniform(1.0, 2.2, size=(len(rows), len(cols)))

    # demo 的 UP_INTERFERON 在 T/髓系中有中等信号（让 demo de_up 里干扰素基因也有映射）
    ifn = ["STAT1", "STAT2", "IFNG", "CXCL9", "CXCL10", "CXCL11", "ISG15",
           "IFI44L", "IFIT1", "IFIT3", "OAS1", "OAS2", "MX1", "MX2"]
    ifn = [g for g in ifn if g in gene_idx]
    for ct in ("T cell", "Myeloid", "Epithelial"):
        rows = marker_rows[ct]
        if rows.size and ifn:
            cols = np.array([gene_idx[g] for g in ifn])
            base[np.ix_(rows, cols)] *= rng.uniform(1.5, 2.4, size=(len(rows), len(cols)))

    # dropout (structural zeros) + log1p -> sparse CSR
    keep = rng.random((n_cells, len(genes))) < 0.55
    counts = np.where(keep, base * 4.0, 0.0)
    X = np.log1p(counts).astype(np.float32)
    X_sparse = sparse.csr_matrix(X)

    # 30 shared DE genes: 10 Proliferating (bulk de_up 方向) + 10 T + 10 B (de_down 方向)
    up_markers = list(markers["Proliferating"][:10])
    down_markers = list(markers["T cell"][:10]) + list(markers["B cell"][:10])
    shared = [{"gene": g, "cell_type": "Proliferating"} for g in up_markers] + [
        {"gene": g, "cell_type": ct}
        for g, ct in zip(down_markers, ["T cell"] * 10 + ["B cell"] * 10)
    ]
    assert len(shared) == 30

    obs = __import__("pandas").DataFrame(index=[f"cell{i:04d}" for i in range(n_cells)])
    obs["cell_type"] = [str(t) for t in types]
    cluster_map = {ct: str(i) for i, ct in enumerate(CELL_TYPES)}
    obs["clusters"] = [_cluster_label(cluster_map[str(t)]) for t in types]
    obs["leiden"] = obs["clusters"].to_numpy()
    obs["donor"] = [f"D{int(rng.integers(1, 4))}" for _ in range(n_cells)]
    obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)
    obs["total_counts"] = X.sum(axis=1).astype(np.float32)

    var = __import__("pandas").DataFrame(index=genes)
    var["gene_ids"] = [f"ENSG{i:011d}" for i in range(len(genes))]

    adata = ad.AnnData(X=X_sparse, obs=obs, var=var)
    adata.uns["fixture"] = "synthetic bulk-ncrna demo scRNA-seq (not real data)"
    adata.write_h5ad(path, compression="gzip")

    truth = {
        "synthetic": True,
        "n_cells": n_cells,
        "n_genes": len(genes),
        "cell_types": CELL_TYPES,
        "markers": markers,
        "shared_de_genes": [s["gene"] for s in shared],
        "shared_de_map": {s["gene"]: s["cell_type"] for s in shared},
        "de_up_genes": up_markers,
        "de_down_genes": down_markers,
        "expected_up_top": "Proliferating",
        "expected_down_top_in": ["T cell", "B cell"],
        "genes": genes,
    }
    with open(os.path.join(os.path.dirname(path), "sc_truth.json"), "w", encoding="utf-8") as f:
        json.dump(truth, f, ensure_ascii=False, indent=1)
    return truth


def _cluster_label(ct_label: str) -> str:
    return f"c{ct_label}"


# ---------------------------------------------------------------- ncrna fixtures

def _write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if v is None else str(v) for v in r) + "\n")


def build_sample_info(path):
    _write_tsv(path, ["sample", "group", "batch"], [
        [s, g, "B1" if i < 6 else "B2"] for i, (s, g) in enumerate(zip(SAMPLES, GROUPS))
    ])


def build_rmats(path):
    rng = np.random.default_rng(SEED + 1)
    header = ["ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
              "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "ID",
              "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
              "IncFormLen", "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2",
              "IncLevelDifference"]
    n = 120
    rows = []
    transcripts = ["NM_000001", "NM_000002", "NM_000003", "NM_000004", "NM_000005"]
    genes = ["CD44", "MKI67", "FN1", "COL1A1", "CDH1", "VIM", "ACTB", "TP53", "BCL2",
             "EGFR", "MET", "ALK", "ROS1", "RET", "FGFR1", "FGFR2", "FGFR3", "ERBB2",
             "NOTCH1", "NOTCH2", "MYC", "KRAS", "NRAS", "BRAF", "PIK3CA", "PTEN",
             "RB1", "CDKN2A", "MDM2", "CCND1"]
    for i in range(n):
        gene = genes[i % len(genes)]
        chrom = f"chr{int(rng.integers(1, 23))}"
        start = int(rng.integers(10000, 90000000))
        inc_len = 120
        skip_len = 250
        delta = rng.normal(0, 0.22)
        if i % 6 == 0:
            delta += 0.55   # true upregulated inclusion
        if i % 7 == 0:
            delta -= 0.55   # true skipped exon
        delta = float(np.clip(delta, -0.95, 0.95))
        psi1 = float(np.clip(0.5 + delta / 2 + rng.normal(0, 0.05), 0.02, 0.98))
        psi2 = float(np.clip(0.5 - delta / 2 + rng.normal(0, 0.05), 0.02, 0.98))
        n1 = rng.integers(12, 60, size=3)
        n2 = rng.integers(12, 60, size=3)
        ijc1 = n1.copy()
        sjc1 = np.maximum(1, (n1 * (1 - psi1) / max(psi1, 0.02)).astype(int))
        ijc2 = n2.copy()
        sjc2 = np.maximum(1, (n2 * (1 - psi2) / max(psi2, 0.02)).astype(int))
        inc1 = np.clip(psi1 + rng.normal(0, 0.04, size=3), 0.01, 0.99)
        inc2 = np.clip(psi2 + rng.normal(0, 0.04, size=3), 0.01, 0.99)
        row = [
            f"SE.MATS.JC.{i + 1}", f"ENSG{i:011d}", gene, chrom,
            rng.choice(["+", "-"]),
            start, start + inc_len - 1,
            start - 300, start - 1,
            start + inc_len, start + inc_len + 299,
            f"SE.MATS.JC.{i + 1}",
            ",".join(str(v) for v in ijc1), ",".join(str(v) for v in sjc1),
            ",".join(str(v) for v in ijc2), ",".join(str(v) for v in sjc2),
            inc_len, skip_len,
            f"{rng.uniform(0.001, 0.9):.4g}", f"{rng.uniform(0.001, 0.95):.4g}",
            ",".join(f"{v:.3f}" for v in inc1),
            ",".join(f"{v:.3f}" for v in inc2),
            f"{float(np.mean(inc1) - np.mean(inc2)):.3f}",
        ]
        rows.append(row)
    _write_tsv(path, header, rows)


def build_arriba(path):
    rng = np.random.default_rng(SEED + 2)
    header = ["#gene1", "gene2", "breakpoint1", "breakpoint2", "type", "confidence",
              "reading_frame", "direction1", "direction2", "junction_reads",
              "spanning_reads", "discordant_mates"]
    fusions = [
        ("EML4", "ALK"), ("TMPRSS2", "ERG"), ("BCR", "ABL1"), ("SLC34A2", "ROS1"),
        ("CD74", "ROS1"), ("FGFR3", "TACC3"), ("FGFR2", "BICC1"), ("RET", "CCDC6"),
        ("NTRK1", "TPM3"), ("NTRK3", "ETV6"), ("EWSR1", "FLI1"), ("PAX3", "FOXO1"),
        ("SS18", "SSX1"), ("SYT1", "SSX2"), ("HMGA2", "LPP"), ("MLL", "AF9"),
        ("KMT2A", "AFF1"), ("TCF3", "PBX1"), ("ETV6", "RUNX1"), ("RUNX1", "RUNX1T1"),
        ("PML", "RARA"), ("MYB", "NFIB"), ("PLAG1", "CTNNB1"), ("FUS", "DDIT3"),
        ("EWSR1", "ATF1"), ("EWSR1", "WT1"), ("NCOA4", "RET"), ("TRIM24", "BRAF"),
        ("BRAF", "KIAA1549"), ("EGFR", "SEPT14"), ("EGFR", "SLC44A2"),
        ("MET", "MST1R"), ("MET", "CAPZA2"), ("RET", "KIF5B"), ("ALK", "STRN"),
        ("ROS1", "GOPC"), ("ROS1", "SDC4"), ("NTRK2", "NTRK2"), ("BRAF", "FAM131B"),
        ("BRD4", "NUTM1"), ("MYC", "IGH"), ("BCL2", "IGH"), ("BCL6", "IGH"),
        ("CCND1", "IGH"), ("MAF", "IGH"), ("MALT1", "IGH"), ("MYB", "GATA1"),
        ("KMT2A", "MLLT3"), ("KMT2A", "MLLT1"), ("KMT2A", "MLLT10"),
        ("TCF3", "HLF"), ("ETV6", "ABL1"), ("ETV6", "NTRK3"), ("RARA", "PLZF"),
        ("JAZF1", "SUZ12"), ("PHF1", "JAZF1"), ("SFPQ", "TFE3"), ("TCF12", "NR4A3"),
        ("ACTB", "GLI1"), ("PTCH1", "GLI1"), ("LPP", "HMGA2"), ("HMGA2", "WIF1"),
        ("THADA", "IGF2BP3"), ("TFG", "NTRK1"), ("TPM4", "ALK"), ("CLTC", "ALK"),
        ("ATIC", "ALK"), ("KIF5B", "RET"), ("CCDC6", "RET"), ("PRKACA", "DNAJB1"),
        ("DNAJB1", "PRKACA"), ("ANKRD28", "ALK"), ("RANBP2", "ALK"), ("VCL", "ALK"),
        ("TFG", "ALK"), ("TPR", "MET"), ("TFE3", "ASPSCR1"), ("YWHAE", "NUTM2B"),
        ("C11orf95", "RELA"), ("MN1", "ETV6"), ("ABL1", "NUP214"), ("ABL1", "ETV6"),
    ]
    rows = []
    for i, (g1, g2) in enumerate(fusions):
        jr = int(rng.integers(2, 120))
        sr = int(rng.integers(0, 80))
        rows.append([
            g1, g2,
            f"chr{int(rng.integers(1, 23))}:{int(rng.integers(10000, 90000000))}:+",
            f"chr{int(rng.integers(1, 23))}:{int(rng.integers(10000, 90000000))}:+",
            "in-frame" if i % 3 != 1 else "out-of-frame",
            "high" if i % 4 != 0 else "low",
            "in-frame" if i % 3 != 1 else "out-of-frame",
            "+", "-", jr, sr, int(rng.integers(0, 40)),
        ])
    _write_tsv(path, header, rows)


def build_star_fusion(path):
    rng = np.random.default_rng(SEED + 3)
    header = ["#FusionName", "JunctionReadCount", "SpanningFragCount", "SpliceType",
              "LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint",
              "FusionJuncSpan", "RECURRENCE"]
    pairs = [("EML4", "ALK"), ("TMPRSS2", "ERG"), ("BCR", "ABL1"), ("SLC34A2", "ROS1"),
             ("FGFR3", "TACC3"), ("RET", "CCDC6"), ("NTRK1", "TPM3"), ("CD74", "ROS1"),
             ("KMT2A", "AFF1"), ("EWSR1", "FLI1"), ("PAX3", "FOXO1"), ("SS18", "SSX1"),
             ("MYC", "IGH"), ("BCL2", "IGH"), ("CCND1", "IGH"), ("ALK", "STRN"),
             ("ROS1", "GOPC"), ("BRAF", "KIAA1549"), ("MET", "CAPZA2"), ("ABL1", "NUP214")]
    rows = []
    for i, (g1, g2) in enumerate(pairs):
        jr = int(rng.integers(3, 90))
        sr = int(rng.integers(0, 60))
        rows.append([
            f"{g1}--{g2}", jr, sr, "ONLY_REF_SPLICE",
            f"{g1}^ENSG00000{i:06d}.1", f"{g2}^ENSG00000{i + 100:06d}.1",
            f"chr1:100000:+", f"chr1:200000:+", f"{jr} ({jr + sr})", "1",
        ])
    _write_tsv(path, header, rows)


def build_psi_matrix(path):
    rng = np.random.default_rng(SEED + 4)
    types = ["SE", "A3SS", "A5SS", "MXE", "RI"]
    genes = ["CD44", "FN1", "MKI67", "COL1A1", "CDH1", "VIM", "EGFR", "MET", "ALK",
             "ROS1", "RET", "FGFR3", "ERBB2", "NOTCH1", "MYC", "KRAS", "PTEN", "RB1",
             "MDM2", "CCND1"]
    rows = []
    n = 120
    for i in range(n):
        et = types[i % len(types)]
        gene = genes[i % len(genes)]
        chrom = f"chr{int(rng.integers(1, 23))}"
        start = int(rng.integers(10000, 90000000))
        ev = f"{et};{et}_{i + 1:04d};{gene};{chrom}:{start}-{start + 200}:+"
        vals = []
        shift = 0.28 if i % 5 == 0 else (-0.28 if i % 8 == 0 else 0.0)
        for j in range(12):
            base = 0.5 + shift / 2 if j < 6 else 0.5 - shift / 2
            v = float(np.clip(base + rng.normal(0, 0.07), 0.0, 1.0))
            if rng.random() < 0.02:
                vals.append("")
            else:
                vals.append(f"{v:.4f}")
        rows.append([ev] + vals)
    _write_tsv(path, ["event_id"] + SAMPLES, rows)


def build_circ(path):
    rng = np.random.default_rng(SEED + 5)
    header = ["circ_id"] + SAMPLES + ["linear"]
    rows = []
    n = 60
    genes = ["CD44", "FN1", "MKI67", "COL1A1", "CDH1", "VIM", "EGFR", "MET", "ALK",
             "ROS1", "RET", "FGFR3", "ERBB2", "NOTCH1", "MYC", "KRAS", "PTEN", "RB1",
             "MDM2", "CCND1", "HIPK3", "ZCCHC7", "SMARCA5", "LARP1", "DNMT1", "FAM120A"]
    for i in range(n):
        if i % 2 == 0:
            chrom = f"chr{int(rng.integers(1, 23))}"
            start = int(rng.integers(10000, 90000000))
            length = int(rng.integers(400, 4000))
            cid = f"{chrom}:{start}-{start + length}:{rng.choice(['+', '-'])}"
        else:
            cid = f"hsa_circ_{i:06d}"
        fold = 1.0
        if i % 4 == 0:
            fold = 2.6          # differential circRNAs case vs control
        elif i % 9 == 0:
            fold = 0.4
        vals = []
        for j in range(12):
            base = 40.0 * fold if j < 6 else 40.0
            vals.append(int(max(1, rng.poisson(base))))
        linear = int(max(10, rng.poisson(180)))
        rows.append([cid] + vals + [linear])
    _write_tsv(path, header, rows)


def build_mirna(path):
    rng = np.random.default_rng(SEED + 6)
    header = ["mirna"] + SAMPLES
    names = [
        "hsa-miR-21-5p", "hsa-miR-21-3p", "hsa-miR-155-5p", "hsa-miR-146a-5p",
        "hsa-miR-155-3p", "hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-18a-5p",
        "hsa-miR-19a-3p", "hsa-miR-19b-3p", "hsa-miR-20a-5p", "hsa-miR-22-3p",
        "hsa-miR-24-3p", "hsa-miR-25-3p", "hsa-miR-26a-5p", "hsa-miR-26b-5p",
        "hsa-miR-27a-3p", "hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p",
        "hsa-miR-30a-5p", "hsa-miR-30b-5p", "hsa-miR-31-5p", "hsa-miR-34a-5p",
        "hsa-miR-101-3p", "hsa-miR-124-3p", "hsa-miR-125a-5p", "hsa-miR-125b-5p",
        "hsa-miR-126-3p", "hsa-miR-126-5p", "hsa-miR-133a-3p", "hsa-miR-143-3p",
        "hsa-miR-145-5p", "hsa-miR-146b-5p", "hsa-miR-150-5p", "hsa-miR-181a-5p",
        "hsa-miR-182-5p", "hsa-miR-200a-3p", "hsa-miR-200b-3p", "hsa-miR-200c-3p",
        "hsa-miR-205-5p", "hsa-miR-210-3p", "hsa-miR-214-3p", "hsa-miR-221-3p",
        "hsa-miR-222-3p", "hsa-miR-223-3p", "hsa-miR-296-5p", "hsa-miR-335-5p",
        "hsa-miR-375-3p", "hsa-miR-451a", "hsa-miR-486-5p", "hsa-miR-493-5p",
        "hsa-miR-497-5p", "hsa-miR-503-5p", "hsa-miR-let-7a-5p", "hsa-let-7a-5p",
        "hsa-let-7b-5p", "hsa-let-7c-5p", "hsa-let-7d-5p", "hsa-let-7e-5p",
        "hsa-let-7f-5p", "hsa-let-7g-5p", "hsa-mir-21", "hsa-mir-155",
        "hsa-mir-146a", "hsa-mir-34a", "hsa-mir-200c", "hsa-mir-9-1",
    ]
    rows = []
    for i, name in enumerate(names):
        fold = 1.0
        if i % 5 == 0:
            fold = 2.4
        elif i % 11 == 0:
            fold = 0.45
        vals = []
        for j in range(12):
            mu = 120.0 * fold if j < 6 else 120.0
            vals.append(f"{max(0.0, rng.lognormal(mean=np.log(mu), sigma=0.35)):.3f}")
        rows.append([name] + vals)
    _write_tsv(path, header, rows)


def build_lncrna(path):
    rng = np.random.default_rng(SEED + 7)
    header = ["lncrna"] + SAMPLES
    names = [
        "MALAT1", "NEAT1", "XIST", "HOTAIR", "PVT1", "MEG3", "H19", "TUG1", "GAS5",
        "DANCR", "NORAD", "KCNQ1OT1", "HULC", "UCA1", "CCAT1", "CCAT2", "HOTTIP",
        "HOXA11-AS", "HOXD-AS1", "CDKN2B-AS1", "ANRIL", "BANCR", "LINC00152",
        "LINC00261", "LINC00312", "LINC00472", "LINC00511", "LINC00662",
        "LINC00673", "LINC00707", "LINC00958", "LINC01006", "LINC01116",
        "LINC01234", "LINC01314", "LINC01410", "LINC01503", "LINC01614",
        "LINC01705", "LINC01857", "LINC01963", "LINC02035", "LINC02159",
        "LINC02268", "LINC02381", "LINC02446", "LINC02532", "LINC02605",
        "LINC02731", "LINC02807", "SNHG1", "SNHG3", "SNHG5", "SNHG6", "SNHG7",
        "SNHG8", "SNHG12", "SNHG15", "SNHG16", "SNHG17", "SNHG20", "ZEB1-AS1",
        "ZEB2-AS1", "ZNFX1-AS1", "DLX6-AS1", "FOXD3-AS1", "PAX8-AS1", "TP53TG1",
        "PCAT1", "PCAT6", "PCAT7", "PRNCR1", "SPRY4-IT1", "TINCR", "HAR1A",
    ]
    rows = []
    for i, name in enumerate(names):
        fold = 1.0
        if i % 5 == 0:
            fold = 2.2
        elif i % 9 == 0:
            fold = 0.5
        vals = []
        for j in range(12):
            mu = 60.0 * fold if j < 6 else 60.0
            vals.append(f"{max(0.0, rng.lognormal(mean=np.log(mu), sigma=0.4)):.3f}")
        rows.append([name] + vals)
    _write_tsv(path, header, rows)


def ensure_fixtures(force: bool = False) -> dict:
    """Create all fixtures if missing. Returns truth dict for demo_sc."""
    truth_path = os.path.join(HERE, "sc_truth.json")
    h5ad = os.path.join(HERE, "demo_sc.h5ad")
    if force or not os.path.exists(h5ad) or not os.path.exists(truth_path):
        truth = build_demo_sc(h5ad)
    else:
        with open(truth_path, encoding="utf-8") as f:
            truth = json.load(f)
    makers = {
        "sample_info.tsv": build_sample_info,
        "rmats_SE.txt": build_rmats,
        "arriba_fusions.tsv": build_arriba,
        "star_fusion.tsv": build_star_fusion,
        "psi_matrix.tsv": build_psi_matrix,
        "circ_expr.tsv": build_circ,
        "mirna_expr.tsv": build_mirna,
        "lncrna_expr.tsv": build_lncrna,
    }
    for name, fn in makers.items():
        p = os.path.join(HERE, name)
        if force or not os.path.exists(p):
            fn(p)
    return truth


def path_for(name: str) -> str:
    return os.path.join(HERE, name)


if __name__ == "__main__":
    t = ensure_fixtures(force="--force" in __import__("sys").argv)
    print("fixtures ready in", HERE)
    print("demo_sc.h5ad:", t["n_cells"], "cells x", t["n_genes"], "genes,",
          len(t["cell_types"]), "cell types")
    for name in ["rmats_SE.txt", "arriba_fusions.tsv", "star_fusion.tsv",
                 "psi_matrix.tsv", "circ_expr.tsv", "mirna_expr.tsv",
                 "lncrna_expr.tsv", "sample_info.tsv"]:
        p = os.path.join(HERE, name)
        print(f"  {name}: {os.path.getsize(p)} bytes")
