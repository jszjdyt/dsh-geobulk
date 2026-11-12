#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 内置基因集 / 通路库（离线回退）+ Enrichr 在线库拉取与缓存。

契约：CONTRACT.md §8.2（bulk-de 负责基因集）。
本文件必须**自包含**：无网络、无 R、无 bg_common 时也能 import 并返回可用基因集。

对外 API
--------
* ``GO_BP`` / ``GO_CC`` / ``GO_MF`` / ``KEGG`` / ``HALLMARK`` / ``REACTOME``
      每个都是 ``{term_id: {"name": str, "genes": [SYMBOL, ...]}}``
* ``IMMUNE_28``：Charoentong 28 类免疫细胞签名，``{cell_type: [SYMBOL, ...]}``
      （供 bulk-clin 免疫模块复用，键名固定 IMMUNE_28）
* ``get_gene_sets(db, workdir=None)`` -> ``{term_id: {"name","genes"}}``
* ``get_gene_sets_with_source(db, workdir=None)`` -> ``(sets, source, warnings)``
* ``load_gmt(path)`` / ``to_gmt(sets)``：标准 GMT 读写
* ``fetch_enrichr_library(name, workdir=None, timeout=20)``：Enrichr 在线库（缓存）
* ``available_dbs()`` / ``builtin_stats()``

基因集内容为人工整理的常用条目（GO/KEGG/Hallmark/Reactome 精选子集，基因 symbol，
人类），用于无网络/无 R 时的确定性回退；有 R 时优先 clusterProfiler/msigdbr 官方库。
"""

from __future__ import annotations

import importlib.util
import os
import re
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# EngineError：优先用 bg_common 的（引擎按 .code/.message 识别）；独立导入时本地兜底
# ---------------------------------------------------------------------------
try:  # pragma: no cover - 正常运行时总走 bg_common
    from bg_common import EngineError  # type: ignore
except Exception:  # noqa: BLE001

    class EngineError(Exception):  # type: ignore[no-redef]
        """独立导入时的最小错误类（与 bg_common.EngineError 字段一致）。"""

        VALID_CODES = {
            "bad_request", "not_found", "prerequisite", "unsupported",
            "ambiguous", "internal", "timeout", "no_r", "no_python",
        }

        def __init__(self, code: str, message: str, **details):
            self.code = code if code in self.VALID_CODES else "internal"
            self.message = str(message)
            self.details = details or {}
            super().__init__(self.message)


ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
ENRICHR_GMT_URL = ENRICHR_BASE + "/geneSetLibrary?mode=text&libraryName=%s"
MIN_GENES_PER_SET = 8
_CACHE_DIRNAME = "genesets"


# ===========================================================================
# 1. 内置库（离线回退；key = 标准 ID，value = (name, "space separated symbols")）
# ===========================================================================

# ---- GO Biological Process（常用条目） -----------------------------------
_RAW_GO_BP = {
    "GO:0006954": ("Inflammatory response",
                   "IL1A IL1B IL6 CXCL8 TNF CXCL1 CXCL2 CXCL3 CCL2 CCL3 CCL4 CCL5 PTGS2 "
                   "NFKB1 RELA TLR4 MYD88 IL1R1 IL6R ICAM1 VCAM1 SELE"),
    "GO:0006955": ("Immune response",
                   "CD4 CD8A CD8B CD19 MS4A1 PTPRC HLA-A HLA-B HLA-DRA IL2 IFNG TNF CCR7 "
                   "CXCR4 CD28 CD80 CD86"),
    "GO:0045087": ("Innate immune response",
                   "TLR1 TLR2 TLR4 TLR5 TLR9 MYD88 TICAM1 IRF3 IRF7 IFNB1 IFNA1 STAT1 NFKB1 "
                   "NLRP3 IL1B IL18 C3 C1QA"),
    "GO:0002250": ("Adaptive immune response",
                   "CD4 CD8A CD19 MS4A1 CD40LG CD28 CTLA4 ICOS IL2 IL4 IL7R JAK1 JAK3 STAT5A "
                   "PTPRC HLA-DRA"),
    "GO:0006950": ("Response to stress",
                   "ATF3 JUN FOS HIF1A TP53 SOD1 SOD2 CAT HSPA1A HSPA1B HSP90AA1 DDIT3 GADD45A"),
    "GO:0007165": ("Signal transduction",
                   "GNB1 GNG2 GNAI1 GNAS ADCY1 PRKACA MAPK1 MAPK3 PIK3CA AKT1 PLCG1 EGFR INS "
                   "GRB2 SOS1 RASA1"),
    "GO:0007155": ("Cell adhesion",
                   "CDH1 CDH2 CTNNB1 CTNNA1 ITGB1 ITGA5 ITGAV ICAM1 VCAM1 PECAM1 NCAM1 L1CAM "
                   "JAM1 CD44 FN1"),
    "GO:0008283": ("Cell population proliferation",
                   "MKI67 PCNA CCND1 CCNE1 CDK1 CDK2 CDK4 CDK6 E2F1 MYC FOS JUN EGF EGFR IL2 IL2RA"),
    "GO:0006915": ("Apoptotic process",
                   "BAX BAK1 BAD BCL2 BCL2L1 MCL1 BID CASP3 CASP8 CASP9 APAF1 CYCS TP53 TNF FAS"),
    "GO:0042981": ("Regulation of apoptotic process",
                   "BCL2 BCL2L1 MCL1 BAX BAK1 BAD BID BIRC2 BIRC3 XIAP TP53 MDM2 AKT1 PIK3CA NFKB1"),
    "GO:0006260": ("DNA replication",
                   "MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 PCNA RFC1 RPA1 POLA1 POLD1 POLE CDK2 CDC6 CDT1 TOP2A"),
    "GO:0006281": ("DNA repair",
                   "RAD51 BRCA1 BRCA2 PALB2 ATM ATR CHEK1 CHEK2 XRCC5 XRCC6 PRKDC LIG4 MSH2 MSH6 "
                   "MLH1 ERCC1 XPA XPC PARP1 OGG1"),
    "GO:0006974": ("Cellular response to DNA damage stimulus",
                   "TP53 CDKN1A GADD45A ATM ATR CHEK1 CHEK2 H2AFX MDC1 BRCA1 RAD51 PARP1 XPC"),
    "GO:0007049": ("Cell cycle",
                   "CDK1 CDK2 CDK4 CDK6 CCNA2 CCNB1 CCNB2 CCND1 CCNE1 E2F1 RB1 CDKN1A CDKN1B "
                   "CDKN2A CDC25A CDC25C PLK1 AURKA BUB1 MAD2L1"),
    "GO:0051301": ("Cell division",
                   "CDK1 PLK1 AURKA AURKB BUB1 MAD2L1 NDC80 NUF2 KIF11 KIF23 ANLN ECT2 RACGAP1 "
                   "PRC1 TOP2A"),
    "GO:0000278": ("Mitotic cell cycle",
                   "CDK1 CCNB1 CCNB2 CDC20 CDC25A CDC25B CDC25C PLK1 AURKA BUB1 MAD2L1 KIF11 "
                   "KIF23 ECT2 INCENP"),
    "GO:0007059": ("Chromosome segregation",
                   "BUB1 BUB1B MAD2L1 NDC80 NUF2 SPC24 SPC25 SKA1 SKA2 CENPE CENPF KIF2C KIF11 "
                   "PLK1 AURKB INCENP SMC1A SMC3"),
    "GO:0006412": ("Translation",
                   "RPL3 RPL4 RPL5 RPL6 RPL7 RPL8 RPL10 RPS3 RPS5 RPS6 RPS8 RPS9 EEF1A1 EEF1B2 "
                   "EEF2 EIF4E"),
    "GO:0006413": ("Translational initiation",
                   "EIF1 EIF2S1 EIF2S2 EIF3A EIF3B EIF4A1 EIF4B EIF4E EIF4G1 EIF5 EIF5B PABPC1 "
                   "RPS3 RPS5"),
    "GO:0006397": ("mRNA processing",
                   "SNRNP70 SRSF1 SRSF2 SRSF3 HNRNPA1 HNRNPA2B1 HNRNPC PTBP1 U2AF1 U2AF2 PRPF8 "
                   "SF3B1 DDX5 DDX17"),
    "GO:0008380": ("RNA splicing",
                   "SF3B1 SRSF1 SRSF2 SRSF3 HNRNPA1 HNRNPA2B1 PTBP1 U2AF1 U2AF2 PRPF8 SNRNP70 "
                   "RBM5 RBMX"),
    "GO:0006355": ("Regulation of DNA-templated transcription",
                   "MYC JUN FOS EGR1 STAT1 STAT3 NFKB1 RELA TP53 SP1 YY1 ATF3 NR3C1 ESR1"),
    "GO:0032774": ("RNA biosynthetic process",
                   "POLR2A POLR2B POLR2C POLR2D POLR2E POLR2H TBP TAF1 TAF2 GTF2B GTF2F1 CDK7 CCNH"),
    "GO:0016070": ("RNA metabolic process",
                   "POLR2A POLR2B TBP TAF1 SRSF1 HNRNPA1 DDX5 ELAVL1 ZFP36 YTHDF1 YTHDF2 METTL3"),
    "GO:0007154": ("Cell communication",
                   "CD40LG IL2 IFNG TNF CCL2 CXCL8 EGFR ERBB2 GNB1 GNG2 ADCY1 GNAS PRKACA PLCG1 "
                   "JAK1 STAT1"),
    "GO:0007229": ("Integrin-mediated signaling pathway",
                   "ITGB1 ITGB2 ITGB3 ITGA2 ITGA4 ITGA5 ITGAV ITGAM ITGAX PTK2 PXN VCL TLN1 SRC "
                   "FYN ILK"),
    "GO:0030198": ("Extracellular matrix organization",
                   "COL1A1 COL1A2 COL3A1 COL4A1 COL5A1 COL6A1 FN1 FBN1 LAMA1 LAMB1 LAMC1 MMP2 "
                   "MMP9 TIMP1 SPARC THBS1 VTN"),
    "GO:0001525": ("Angiogenesis",
                   "VEGFA VEGFB VEGFC KDR FLT1 FLT4 ANGPT1 ANGPT2 TEK PDGFB PDGFRB TGFB1 FGF2 "
                   "HIF1A NRP1 PECAM1"),
    "GO:0001568": ("Blood vessel development",
                   "VEGFA KDR FLT1 TEK ANGPT1 ANGPT2 PDGFB PDGFRB TGFB1 NOTCH1 DLL4 JAG1 EFNB2 "
                   "PECAM1 CDH5"),
    "GO:0001816": ("Cytokine production",
                   "IL1B IL6 IL12A IL12B IL18 TNF IFNG CSF2 CCL2 CXCL8 TLR4 MYD88 NFKB1 RELA STAT1 IRF3"),
    "GO:0001817": ("Regulation of cytokine production",
                   "IL10 TGFB1 IFNG IL4 IL13 IL2 SOCS1 SOCS3 STAT1 STAT3 STAT5A NFKB1 RELA NLRP3 PTPN22"),
    "GO:0006935": ("Chemotaxis",
                   "CXCL12 CXCR4 CCL2 CCR2 CCL3 CCR1 CCL5 CCR5 CXCL8 CXCR1 CXCL1 CXCR2 S1PR1 "
                   "S1PR3 LTB4R"),
    "GO:0030593": ("Neutrophil chemotaxis",
                   "CXCL1 CXCL2 CXCL3 CXCL5 CXCL6 CXCL8 CXCR1 CXCR2 CCL3 CCL4 S100A8 S100A9 LTB4R"),
    "GO:0007166": ("Cell surface receptor signaling pathway",
                   "EGFR ERBB2 ERBB3 MET FGFR1 FGFR2 IGF1R INSR PDGFRA PDGFRB TGFBR1 TGFBR2 "
                   "TNFRSF1A IL6R IFNGR1 JAK1 STAT1"),
    "GO:0006468": ("Protein phosphorylation",
                   "CDK1 CDK2 MAPK1 MAPK3 MAPK8 MAPK14 AKT1 SRC FYN LYN JAK1 JAK2 STAT1 STAT3 "
                   "EGFR PRKACA GSK3B"),
    "GO:0006508": ("Proteolysis",
                   "MMP2 MMP9 MMP14 CTSB CTSD CTSS CTSK CASP3 CASP8 CASP9 PSMB8 PSMB9 LGMN "
                   "ADAM17 ADAM10"),
    "GO:0006511": ("Ubiquitin-dependent protein catabolic process",
                   "PSMA1 PSMA2 PSMA3 PSMA4 PSMB1 PSMB2 PSMB5 PSMB8 PSMB9 PSMC1 PSMC2 PSMD1 "
                   "UBA1 UBE2C UBB UBC"),
    "GO:0016192": ("Vesicle-mediated transport",
                   "RAB5A RAB7A RAB11A RAB27A STX4 STX7 VAMP2 VAMP3 SNAP23 SNAP25 NSF SYT1 CFL1 ARF1"),
    "GO:0006897": ("Endocytosis",
                   "CLTC CLTA AP2A1 AP2B1 EPS15 EPN1 DNM1 DNM2 RAB5A ARF6 CAV1 CAV2 ITSN1 SNAP91"),
    "GO:0006886": ("Intracellular protein transport",
                   "SEC61A1 SEC61B SEC62 SEC63 SRP54 SRP9 COPA COPB1 COPB2 ARF1 RAB1A RAB2A "
                   "KPNA2 KPNB1 NUP98 RAN"),
    "GO:0006096": ("Glycolytic process",
                   "HK1 HK2 GPI PFKL PFKM PFKP ALDOA TPI1 GAPDH PGK1 PGAM1 ENO1 PKM LDHA PDK1"),
    "GO:0006631": ("Fatty acid metabolic process",
                   "ACACA ACACB FASN SCD ELOVL1 ELOVL6 FADS1 FADS2 CPT1A CPT2 ACADM ACADL HADHA HADHB"),
    "GO:0006629": ("Lipid metabolic process",
                   "PLIN2 PLIN3 DGAT1 DGAT2 LIPE PNPLA2 AGPAT2 LPIN1 SREBF1 SREBF2 HMGCR LDLR "
                   "ABCA1 APOE"),
    "GO:0006695": ("Cholesterol biosynthetic process",
                   "HMGCS1 HMGCR MVK MVD IDI1 SQLE LSS DHCR7 DHCR24 CYP51A1 MSMO1 NSDHL FDFT1 "
                   "INSIG1 SREBF2"),
    "GO:0006694": ("Steroid biosynthetic process",
                   "CYP11A1 CYP17A1 CYP19A1 CYP21A2 STS HSD3B1 HSD3B2 HSD17B1 SRD5A1 AKR1C3 DHCR7"),
    "GO:0055085": ("Transmembrane transport",
                   "SLC2A1 SLC2A3 SLC7A5 SLC7A11 SLC16A1 SLC16A3 SLC25A5 SLC25A6 SLC4A1 SLC9A1 "
                   "SLC12A2 ABCB1 ABCG2 ATP1A1 ATP1B1"),
    "GO:0006811": ("Ion transport",
                   "ATP1A1 ATP1B1 SLC9A1 SLC12A2 SCN1A KCNQ1 KCNA1 CLCN3 TRPM7 TRPV1 SLC8A1 "
                   "CACNA1C HCN2"),
    "GO:0006810": ("Transport",
                   "SLC2A1 SLC7A5 SLC25A5 ATP1A1 ATP1B1 ABCB1 ABCG2 CLCN3 KCNQ1 SCN1A TRPM7 "
                   "RAB5A RAB7A SEC61A1"),
    "GO:0007160": ("Cell-matrix adhesion",
                   "ITGB1 ITGA5 ITGAV ITGA2 ITGA6 ITGB4 LAMB1 LAMC1 FN1 COL1A1 PTK2 PXN TLN1 "
                   "VCL ILK"),
    "GO:0098609": ("Cell-cell adhesion",
                   "CDH1 CDH2 CTNNA1 CTNNB1 CTNND1 JUP DSC1 DSG1 OCLN TJP1 CLDN1 CLDN3 ICAM1 "
                   "VCAM1 PECAM1"),
    "GO:0016477": ("Cell migration",
                   "RHOA RAC1 CDC42 PXN PTK2 ITGB1 MMP2 MMP9 CXCR4 CCR7 S1PR1 VEGFA MYH9 VIM"),
    "GO:0009615": ("Response to virus",
                   "IFIT1 IFIT2 IFIT3 ISG15 MX1 MX2 OAS1 OAS2 DDX58 IFIH1 IRF3 IRF7 STAT1 STAT2 IFNB1"),
    "GO:0051607": ("Defense response to virus",
                   "IFIT1 IFIT3 ISG15 MX1 MX2 OAS1 OAS2 OAS3 DDX58 IFIH1 IFNB1 IFNL1 STAT1 STAT2 IRF7"),
    "GO:0002224": ("Toll-like receptor signaling pathway",
                   "TLR1 TLR2 TLR3 TLR4 TLR5 TLR7 TLR8 TLR9 MYD88 TICAM1 TIRAP IRAK1 IRAK4 TRAF6 "
                   "IRF3 IFNB1 NFKB1"),
    "GO:0007173": ("Epidermal growth factor receptor signaling pathway",
                   "EGFR ERBB2 ERBB3 ERBB4 HBEGF EGF TGFA AREG EREG GRB2 SOS1 HRAS MAP2K1 MAPK1 "
                   "STAT3 PIK3CA"),
    "GO:0006457": ("Protein folding",
                   "HSPA1A HSPA1B HSPA5 HSPA8 HSP90AA1 HSP90AB1 HSPD1 HSPE1 DNAJB1 DNAJC3 CCT5 "
                   "CALR CANX PPIA"),
    "GO:0034976": ("Response to endoplasmic reticulum stress",
                   "HSPA5 ATF4 ATF6 DDIT3 XBP1 ERN1 EIF2AK3 DNAJB9 SEL1L EDEM1 HERPUD1 PDIA3 PDIA4"),
    "GO:0030968": ("Endoplasmic reticulum unfolded protein response",
                   "HSPA5 ATF4 ATF6 DDIT3 XBP1 ERN1 EIF2AK3 DNAJB9 HERPUD1 SEL1L PDIA4 HSP90B1"),
    "GO:0016567": ("Protein ubiquitination",
                   "UBA1 UBE2C UBE2D1 UBE3A MDM2 SKP2 FBXW7 CBL ITCH NEDD4 HERC1 RNF2 PRKN TRIM21"),
    "GO:0006486": ("Protein glycosylation",
                   "ALG1 ALG2 ALG5 ALG6 ALG8 B4GALT1 MGAT1 MGAT2 MGAT4A FUT8 ST6GAL1 UGGT1 DDOST RPN1"),
    "GO:0006091": ("Generation of precursor metabolites and energy",
                   "NDUFA1 NDUFB1 SDHA SDHB UQCRC1 UQCRC2 COX4I1 ATP5F1A ATP5F1B PFKM HK2 GAPDH "
                   "PDHA1 IDH2"),
    "GO:0006119": ("Oxidative phosphorylation",
                   "NDUFA1 NDUFA2 NDUFB1 NDUFB8 SDHA SDHB SDHC UQCRC1 UQCRC2 UQCRFS1 COX4I1 COX5A "
                   "COX6A1 ATP5F1A ATP5F1B ATP5F1C"),
    "GO:0006936": ("Muscle contraction",
                   "ACTA1 ACTC1 MYH1 MYH2 MYH7 MYL2 MYL3 TPM1 TNNT2 TNNC1 DES DMD TTN RYR1"),
    "GO:0007010": ("Cytoskeleton organization",
                   "ACTB ACTG1 TUBA1A TUBB MYH9 RHOA RAC1 CDC42 CFL1 DSTN VIM LMNA SPTAN1"),
    "GO:0030036": ("Actin cytoskeleton organization",
                   "ACTB ACTG1 ACTR2 ACTR3 ARPC1B ARPC2 CRK CFL1 DSTN LIMK1 SSH1 ROCK1 RHOA"),
    "GO:0007015": ("Actin filament organization",
                   "ACTB ACTG1 MYH9 MYL6 MYL12A CFL1 DSTN TPM1 TPM4 PFN1 GSN CAPZA1 CAPZB"),
    "GO:0032502": ("Developmental process",
                   "BMP4 BMP7 WNT3A SHH FGF8 GLI1 CTNNB1 SMAD1 SMAD4 NOTCH1 DLL1 HES1 PAX6 SOX2"),
    "GO:0048869": ("Cellular developmental process",
                   "CD34 KIT RUNX1 GATA2 SPI1 CEBPA MYOD1 MYOG PAX7 SOX9 RUNX2 PPARG CEBPB"),
    "GO:0030154": ("Cell differentiation",
                   "RUNX2 PPARG CEBPA MYOD1 MYOG GATA1 GATA2 SPI1 PAX7 SOX9 BMP2 SMAD4"),
    "GO:0007399": ("Nervous system development",
                   "NGF BDNF NTF3 NTRK1 NTRK2 SOX2 PAX6 NEUROD1 ASCL1 DCX TUBB3 MAP2 GFAP MBP"),
    "GO:0001649": ("Osteoblast differentiation",
                   "RUNX2 SP7 ALPL COL1A1 IBSP BGLAP SPP1 MSX2 DLX5 SOX9 BMP2"),
    "GO:0042692": ("Muscle cell differentiation",
                   "MYOD1 MYOG MYF5 MYF6 PAX3 PAX7 MEF2C MEF2D DES MYH1 CKM"),
    "GO:0009725": ("Response to hormone",
                   "ESR1 AR PGR NR3C1 THRB RARA RARB VDR PPARG STAT5A STAT5B CEBPB FOXO1"),
    "GO:0010033": ("Response to organic substance",
                   "CYP1A1 CYP1B1 AHR NQO1 NFE2L2 HMOX1 GCLC GCLM TXNRD1 SLC7A11"),
    "GO:0042221": ("Response to chemical",
                   "NFE2L2 HMOX1 NQO1 GCLC GCLM TXNRD1 CYP1A1 AHR HSPA1A SOD1 CAT GPX1"),
    "GO:0009410": ("Response to xenobiotic stimulus",
                   "AHR CYP1A1 CYP1B1 CYP2E1 CYP3A4 NQO1 NFE2L2 GSTA1 GSTP1 UGT1A1"),
    "GO:0015980": ("Energy derivation by oxidation of organic compounds",
                   "SDHA IDH2 MDH2 CS ACLY PDHA1 PDHB DLAT DLST OGDH SUCLA2"),
    "GO:0006163": ("Purine nucleotide metabolic process",
                   "HPRT1 APRT ADA AMPD1 IMPDH1 IMPDH2 GMPS ATIC ADSL ADSS"),
    "GO:0009259": ("Ribonucleotide metabolic process",
                   "RRM1 RRM2 NME1 NME2 HPRT1 APRT IMPDH1 IMPDH2 UMPS CAD"),
    "GO:0045055": ("Regulated exocytosis",
                   "STX1A SNAP25 VAMP2 SYT1 RAB3A RAB27A UNC13A CPLX1 NSF"),
    "GO:0002443": ("Leukocyte mediated immunity",
                   "CD8A GZMB PRF1 GNLY NKG7 CD4 IFNG IL2 CD40LG FASLG"),
    "GO:0043312": ("Neutrophil degranulation",
                   "ELANE MPO PRTN3 CTSG LTF CAMP DEFA3 DEFA4 AZU1 BPI CEACAM8 MMP9 S100A8 S100A9"),
    "GO:0002274": ("Myeloid leukocyte activation",
                   "CD14 TLR4 MYD88 LY96 ITGAM ITGB2 FCGR1A FCGR2A SYK TYROBP SPI1 CEBPB"),
    "GO:0065007": ("Biological regulation",
                   "STAT1 STAT3 NFKB1 RELA TP53 MYC JUN FOS SP1 YY1 ESR1 AR NR3C1"),
    "GO:0050789": ("Regulation of biological process",
                   "STAT1 STAT3 SOCS1 SOCS3 PIAS1 PTPN11 CBL ITCH NEDD4 FBXW7 SKP2 MDM2"),
    "GO:0050896": ("Response to stimulus",
                   "FOS JUN EGR1 ATF3 HSPA1A HSPA5 DDIT3 NFE2L2 HMOX1 SOD1 CAT TP53 NFKB1"),
    "GO:0023052": ("Signaling",
                   "GNB1 GNG2 GNAI1 GNAS ADCY1 PRKACA PLCG1 EGFR GRB2 SOS1 HRAS MAP2K1 MAPK1 "
                   "PIK3CA AKT1 JAK1 STAT1"),
    "GO:0009628": ("Response to abiotic stimulus",
                   "HSPA1A HSPA5 HIF1A NFE2L2 HMOX1 SOD1 SOD2 CAT GPX1 DDIT3 ATF3 JUN"),
    "GO:0009636": ("Response to toxic substance",
                   "NFE2L2 HMOX1 NQO1 GCLC GCLM TXNRD1 CYP1A1 CYP1B1 AHR MT1A MT2A SOD1 CAT"),
    "GO:0006796": ("Phosphate-containing compound metabolic process",
                   "ATP5F1A ATP5F1B ATP1A1 ATP1B1 PRKACA PRKCD PIK3CA PIK3CB PLCG1 PTEN INPP5D"),
    "GO:0008284": ("Positive regulation of cell population proliferation",
                   "MYC CCND1 CCND2 CDK4 CDK6 E2F1 EGF EGFR IL2 IL2RA JUN FOS STAT3 PIK3CA"),
    "GO:0043066": ("Negative regulation of apoptotic process",
                   "BCL2 BCL2L1 MCL1 BCL2A1 BIRC2 BIRC3 BIRC5 XIAP AKT1 PIK3CA IL10 SOCS3"),
    "GO:0030335": ("Positive regulation of cell migration",
                   "CXCR4 CCR7 S1PR1 VEGFA HGF MET RHOA RAC1 CDC42 PXN PTK2 ITGB1"),
    "GO:0009611": ("Response to wounding",
                   "F3 PLAT PLAU SERPINE1 TGFB1 PDGFB VEGFA CXCL8 CCL2 IL6 COL1A1 FN1"),
    "GO:0043488": ("Regulation of mRNA stability",
                   "ZFP36 ZFP36L1 ZFP36L2 ELAVL1 HNRNPA1 HNRNPD MBNL1 PUM1 PUM2 IGF2BP1 IGF2BP2 "
                   "YTHDF1 YTHDF2"),
}

# ---- GO Cellular Component（常用条目） ------------------------------------
_RAW_GO_CC = {
    "GO:0005576": ("Extracellular region",
                   "COL1A1 COL1A2 COL3A1 FN1 LAMA1 LAMB1 LAMC1 SPARC THBS1 VTN MMP2 MMP9 TIMP1 "
                   "APOE C3 CFB SERPINA1 ALB"),
    "GO:0005615": ("Extracellular space",
                   "IL1B IL6 CXCL8 CCL2 CCL3 CCL4 CCL5 TNF VEGFA FGF2 IGF1 EGF HBEGF PDGFB "
                   "SERPINE1 FN1 MMP9"),
    "GO:0031012": ("Extracellular matrix",
                   "COL1A1 COL1A2 COL3A1 COL4A1 COL5A1 COL6A1 FN1 FBN1 LAMA1 LAMB1 LAMC1 LAMC2 "
                   "HAS2 VCAN DCN LUM BGN POSTN"),
    "GO:0005925": ("Focal adhesion",
                   "PTK2 PXN VCL TLN1 ILK PARVA ACTN1 ZYX ITGB1 ITGA5 ITGA6 LAMB1 LAMC1 FN1 SRC"),
    "GO:0005886": ("Plasma membrane",
                   "CD4 CD8A EGFR ERBB2 ITGB1 ITGA5 CXCR4 CCR7 HLA-A PTPRC ICAM1 CDH1 TNFRSF1A "
                   "IL6R ATP1A1 SLC2A1"),
    "GO:0005887": ("Integral component of plasma membrane",
                   "EGFR ERBB2 ITGB1 ITGA5 CXCR4 CCR7 CD4 CD8A HLA-A HLA-B SLC2A1 SLC7A5 ATP1A1 "
                   "TNFRSF1A IL6R"),
    "GO:0070062": ("Extracellular exosome",
                   "CD9 CD63 CD81 TSG101 PDCD6IP SDCBP HSPA8 GAPDH ENO1 PKM ACTB ACTG1 ALDOA PGK1"),
    "GO:0009986": ("Cell surface",
                   "CD4 CD8A CD19 MS4A1 CD14 TLR4 ITGAM PECAM1 ICAM1 HLA-A CD44 ITGB1 CD40 CD28"),
    "GO:0005737": ("Cytoplasm",
                   "GAPDH ACTB ACTG1 RPL3 RPL5 RPS3 RPS6 EEF1A1 EEF2 HSPA8 TUBA1A TUBB VIM YWHAZ"),
    "GO:0005829": ("Cytosol",
                   "GAPDH PGK1 ENO1 PKM LDHA ALDOA TPI1 HK1 PFKM GPI SOD1 CAT PRKACA GSK3B"),
    "GO:0005634": ("Nucleus",
                   "MYC JUN FOS TP53 RB1 E2F1 STAT1 STAT3 NFKB1 RELA SP1 YY1 ESR1 AR LMNA"),
    "GO:0005654": ("Nucleoplasm",
                   "POLR2A TBP TAF1 CDK7 CCNH SRSF1 SRSF2 HNRNPA1 TP53 RB1 E2F1 STAT1 STAT3"),
    "GO:0000785": ("Chromatin",
                   "H2AFX H2AFZ H3F3A H3F3B HMGB1 HMGB2 HMGN1 HMGN2 TOP1 TOP2A EZH2 SUZ12"),
    "GO:0005694": ("Chromosome",
                   "H2AFX H3F3A TOP2A SMC1A SMC3 RAD21 STAG1 STAG2 CTCF CENPA CENPB INCENP"),
    "GO:0000775": ("Chromosome, centromeric region",
                   "CENPA CENPB CENPC CENPE CENPF CENPH NDC80 NUF2 SPC24 SPC25 INCENP AURKB"),
    "GO:0000776": ("Kinetochore",
                   "NDC80 NUF2 SPC24 SPC25 BUB1 BUB1B BUB3 MAD2L1 CENPE CENPF ZW10 ZWINT KNTC1"),
    "GO:0005874": ("Microtubule",
                   "TUBA1A TUBA1B TUBB TUBB2A TUBB3 MAP1B MAP2 MAP4 MAPT STMN1 DCTN1 DYNC1H1 "
                   "KIF11 KIF23"),
    "GO:0005819": ("Spindle",
                   "TUBA1A TUBB KIF11 KIF23 KIF2C PLK1 AURKA AURKB TPX2 NUMA1 TACC3 CKAP5 INCENP"),
    "GO:0000922": ("Spindle pole",
                   "CEP192 PLK1 AURKA TPX2 NUMA1 KIF11 TUBG1 TUBGCP2 TUBGCP3 CDK5RAP2"),
    "GO:0005739": ("Mitochondrion",
                   "SDHA SDHB UQCRC1 UQCRC2 COX4I1 COX5A ATP5F1A ATP5F1B TOMM20 TOMM40 TIMM23 "
                   "VDAC1 VDAC2 HSPD1"),
    "GO:0005743": ("Mitochondrial inner membrane",
                   "SDHA SDHB SDHC SDHD UQCRC1 UQCRC2 UQCRFS1 COX4I1 COX5A ATP5F1A ATP5F1B TOMM40 "
                   "TIMM23 SLC25A5"),
    "GO:0005741": ("Mitochondrial outer membrane",
                   "TOMM20 TOMM22 TOMM40 TOMM70 VDAC1 VDAC2 VDAC3 BCL2 BAX BAK1 MFN1 MFN2"),
    "GO:0005759": ("Mitochondrial matrix",
                   "HSPD1 HSPE1 CS IDH2 MDH2 PDHA1 PDHB DLAT DLST OGDH SUCLA2 SUCLG1 ACO2"),
    "GO:0005783": ("Endoplasmic reticulum",
                   "HSPA5 CALR CANX PDIA3 PDIA4 SEC61A1 SEC61B ATF6 ERN1 RPN1 RPN2 DDOST UGGT1"),
    "GO:0005789": ("Endoplasmic reticulum membrane",
                   "SEC61A1 SEC61B SEC62 SEC63 ATF6 ERN1 RPN1 RPN2 DDOST UGGT1 HMGCR SREBF2"),
    "GO:0005794": ("Golgi apparatus",
                   "GOLGA2 GOLGB1 GOLGA4 B4GALT1 MGAT1 MGAT2 ST6GAL1 FUT8 ARF1 COPB1 COPB2 GOLPH3"),
    "GO:0005768": ("Endosome",
                   "RAB5A RAB7A RAB11A EEA1 CLTC AP2A1 AP2B1 EHD1 SNX1 SNX2 VPS35 VPS26A"),
    "GO:0005764": ("Lysosome",
                   "CTSB CTSD CTSS CTSK LAMP1 LAMP2 NPC1 NPC2 HEXA HEXB GBA GLA SMPD1"),
    "GO:0005765": ("Lysosomal membrane",
                   "LAMP1 LAMP2 CD63 CD68 NPC1 NPC2 SLC11A1 CLCN7 OSTM1 TMEM175 ATP6V0A1"),
    "GO:0005777": ("Peroxisome",
                   "PEX1 PEX3 PEX5 PEX6 PEX7 PEX10 PEX12 PEX13 PEX14 PEX19 CAT ACOX1 HSD17B4"),
    "GO:0005840": ("Ribosome",
                   "RPL3 RPL4 RPL5 RPL6 RPL7 RPL8 RPS3 RPS5 RPS6 RPS8 RPS9 RPLP0 RPSA"),
    "GO:0005911": ("Cell-cell junction",
                   "CDH1 CTNNA1 CTNNB1 CTNND1 JUP TJP1 OCLN CLDN1 CLDN3 GJA1 GJB1 DSC2 DSG2"),
    "GO:0030054": ("Cell junction",
                   "CDH1 CTNNA1 CTNNB1 JUP TJP1 OCLN CLDN1 CLDN3 GJA1 GJB1 PECAM1 ITGB1 ITGA5"),
    "GO:0009897": ("External side of plasma membrane",
                   "HLA-A HLA-B HLA-C HLA-DRA HLA-DRB1 CD4 CD8A CD19 MS4A1 CD14 TLR4 CD40 FAS FASLG"),
    "GO:0045121": ("Membrane raft",
                   "CAV1 CAV2 FLOT1 FLOT2 PROM1 CD44 CD9 CD63 CD81 ITGB1 TLR4 LYN HCK"),
    "GO:0016020": ("Membrane",
                   "EGFR ERBB2 ITGB1 ITGA5 CXCR4 CCR7 CD4 CD8A HLA-A SLC2A1 SLC7A5 ATP1A1 TOMM20"),
    "GO:0098796": ("Membrane protein complex",
                   "ATP5F1A ATP5F1B COX4I1 COX5A NDUFA1 NDUFB1 SDHA SDHB UQCRC1 UQCRC2 TOMM20 "
                   "SEC61A1 SEC61B"),
    "GO:0005747": ("Mitochondrial respiratory chain complex I",
                   "NDUFA1 NDUFA2 NDUFA3 NDUFA4 NDUFA5 NDUFA6 NDUFB1 NDUFB2 NDUFB3 NDUFB8 NDUFS1 "
                   "NDUFS2 NDUFS3"),
    "GO:0005751": ("Mitochondrial respiratory chain complex IV",
                   "COX4I1 COX4I2 COX5A COX5B COX6A1 COX6B1 COX6C COX7A2 COX7B COX7C COX8A"),
    "GO:0005753": ("Mitochondrial proton-transporting ATP synthase complex",
                   "ATP5F1A ATP5F1B ATP5F1C ATP5F1D ATP5F1E ATP5PB ATP5PD ATP5PF ATP5MC1 "
                   "ATP5MC2 ATP5MC3"),
    "GO:0000786": ("Nucleosome",
                   "H2AFX H2AFZ H3F3A H3F3B H4C1 H4C2 H4C3 H2BC1 H2BC3 H1F0 H1FX"),
    "GO:0005730": ("Nucleolus",
                   "NCL NPM1 FBL DKC1 GNL3 NPM3 PES1 BOP1 RPL3 RPL5 RPS3 SIRT7"),
    "GO:0016604": ("Nuclear body",
                   "PML SP100 DAXX ATRX TP53 MDM2 CREBBP EP300 NCL NPM1 SIRT1"),
    "GO:0016607": ("Nuclear speck",
                   "SRSF1 SRSF2 SRSF3 SRSF5 SRSF7 HNRNPA1 HNRNPC U2AF1 U2AF2 PRPF8 SON SRRM1"),
    "GO:0005681": ("Spliceosomal complex",
                   "SNRNP70 SNRPB2 SNRPB SNRPD1 SNRPD2 SNRPD3 SF3B1 SF3A1 PRPF8 PRPF19 U2AF1 U2AF2"),
    "GO:0005839": ("Proteasome core complex",
                   "PSMA1 PSMA2 PSMA3 PSMA4 PSMA5 PSMA6 PSMA7 PSMB1 PSMB2 PSMB3 PSMB4 PSMB5 PSMB6 "
                   "PSMB7 PSMB8 PSMB9 PSMB10"),
    "GO:0000502": ("Proteasome complex",
                   "PSMA1 PSMB1 PSMC1 PSMC2 PSMC3 PSMD1 PSMD2 PSMD3 PSMD4 ADRM1 UCHL5"),
    "GO:0030136": ("Clathrin-coated vesicle",
                   "CLTC CLTA CLTB AP2A1 AP2B1 AP1B1 EPS15 EPN1 DNM2 SYNJ1"),
    "GO:0031982": ("Vesicle",
                   "RAB5A RAB7A RAB11A RAB27A VAMP2 VAMP3 STX4 CD9 CD63 CD81 SNAP23 SYT1"),
    "GO:0009898": ("Cytoplasmic side of plasma membrane",
                   "PXN PTK2 SRC LYN FYN HCK ITGB1 ITGB2 TLN1 VCL PTPN6"),
    "GO:0005604": ("Basement membrane",
                   "LAMA1 LAMB1 LAMC1 LAMA5 LAMB2 COL4A1 COL4A2 COL4A3 COL4A4 COL4A5 COL4A6 NID1 "
                   "NID2 HSPG2"),
    "GO:0031966": ("Mitochondrial membrane",
                   "TOMM20 TOMM40 VDAC1 VDAC2 SDHA UQCRC1 COX4I1 ATP5F1A BCL2 BAX MFN2 OPA1"),
    "GO:0005770": ("Late endosome",
                   "RAB7A RAB9A LAMP1 CD63 NPC1 NPC2 VPS35 VPS26A SNX1"),
    "GO:0005769": ("Early endosome",
                   "RAB5A EEA1 RAB4A RAB11A SNX1 SNX2 SNX5 VPS35 CLTC AP2A1"),
    "GO:0031234": ("Extrinsic component of cytoplasmic side of plasma membrane",
                   "SRC FYN LYN HCK ITGB1 PXN PTK2 TLN1 VCL PTPN6"),
    "GO:0031252": ("Cell leading edge",
                   "RHOA RAC1 CDC42 PXN PTK2 VCL ACTN1 MYH9 ITGB1 VIM"),
    "GO:0015629": ("Actin cytoskeleton",
                   "ACTB ACTG1 ACTN1 MYH9 MYL6 CFL1 DSTN TPM1 TPM4 PFN1 GSN VCL TLN1"),
    "GO:0031965": ("Nuclear membrane",
                   "LMNA LMNB1 LMNB2 NUP98 NUP153 NUP214 NUP62 NDC1 TPR EMD LBR"),
    "GO:0005746": ("Mitochondrial respiratory chain",
                   "NDUFA1 NDUFB1 SDHA SDHB SDHC SDHD UQCRC1 UQCRC2 UQCRFS1 COX4I1 COX5A ATP5F1A "
                   "ATP5F1B"),
}

# ---- GO Molecular Function（常用条目） ------------------------------------
_RAW_GO_MF = {
    "GO:0004672": ("Protein kinase activity",
                   "CDK1 CDK2 CDK4 CDK6 MAPK1 MAPK3 MAPK8 MAPK14 AKT1 SRC FYN LYN JAK1 JAK2 EGFR "
                   "PRKACA GSK3B AURKA PLK1"),
    "GO:0004674": ("Protein serine/threonine kinase activity",
                   "CDK1 CDK2 CDK4 CDK6 MAPK1 MAPK3 MAPK8 MAPK14 AKT1 PRKACA GSK3B MTOR AURKA "
                   "PLK1 CHEK1"),
    "GO:0004713": ("Protein tyrosine kinase activity",
                   "SRC FYN LYN HCK JAK1 JAK2 JAK3 EGFR ERBB2 KDR FLT1 MET ABL1 SYK BTK"),
    "GO:0004714": ("Transmembrane receptor protein tyrosine kinase activity",
                   "EGFR ERBB2 ERBB3 ERBB4 KDR FLT1 FLT4 MET FGFR1 IGF1R INSR PDGFRA PDGFRB NTRK1"),
    "GO:0004930": ("G protein-coupled receptor activity",
                   "CXCR4 CCR7 CXCR1 CXCR2 CCR2 CCR5 ADRB2 ADRB1 AGTR1 S1PR1 S1PR3 GPR183 C5AR1 FPR1"),
    "GO:0004888": ("Transmembrane signaling receptor activity",
                   "EGFR ERBB2 TNFRSF1A IL6R IFNGR1 TLR4 TLR2 NOTCH1 TGFBR1 TGFBR2 KDR FLT1 CXCR4 CCR7"),
    "GO:0005102": ("Signaling receptor binding",
                   "IL6 IL1B IFNG TNF VEGFA EGF FGF2 PDGFB TGFB1 CXCL8 CXCL12 CCL2 CSF2 IL2"),
    "GO:0005515": ("Protein binding",
                   "TP53 ACTB GAPDH MYC JUN FOS EGFR ERBB2 SRC PIK3CA AKT1 STAT1 STAT3 NFKB1 RELA"),
    "GO:0046983": ("Protein dimerization activity",
                   "JUN FOS MYC MAX STAT1 STAT3 NFKB1 RELA BCL2 BAX EGFR ERBB2"),
    "GO:0008092": ("Cytoskeletal protein binding",
                   "ACTN1 VCL TLN1 PXN PTK2 MYH9 CFL1 DSTN TPM1 PFN1 GSN CAPZA1"),
    "GO:0003779": ("Actin binding",
                   "ACTB ACTG1 ACTN1 MYH9 MYL6 CFL1 DSTN TPM1 TPM4 PFN1 GSN CAPZA1 VCL"),
    "GO:0008017": ("Microtubule binding",
                   "TUBA1A TUBB MAP1B MAP2 MAP4 MAPT STMN1 DCTN1 DYNC1H1 KIF11 KIF23 KIF2C"),
    "GO:0008134": ("Transcription factor binding",
                   "CREBBP EP300 MED1 MED12 NCOA1 NCOA2 NCOR1 NCOR2 TAF1 TBP SP1 YY1"),
    "GO:0003700": ("DNA-binding transcription factor activity",
                   "MYC JUN FOS EGR1 STAT1 STAT3 NFKB1 RELA TP53 SP1 YY1 ATF3 NR3C1 ESR1"),
    "GO:0043565": ("Sequence-specific DNA binding",
                   "MYC JUN FOS EGR1 STAT1 STAT3 NFKB1 RELA TP53 SP1 YY1 ATF3 GATA1"),
    "GO:0003677": ("DNA binding",
                   "TP53 MYC JUN FOS E2F1 RB1 STAT1 STAT3 NFKB1 RELA CTCF H2AFX TOP2A"),
    "GO:0003682": ("Chromatin binding",
                   "EZH2 SUZ12 H2AFX H3F3A TOP2A CTCF HMGB1 HMGB2 EP300 CREBBP KDM1A"),
    "GO:0003723": ("RNA binding",
                   "SRSF1 SRSF2 SRSF3 HNRNPA1 HNRNPA2B1 HNRNPC PTBP1 ELAVL1 DDX5 DDX17 EIF4E YTHDF1"),
    "GO:0019843": ("rRNA binding",
                   "RPL3 RPL4 RPL5 RPL6 RPS3 RPS5 RPS6 RPS9 RPLP0"),
    "GO:0003735": ("Structural constituent of ribosome",
                   "RPL3 RPL4 RPL5 RPL6 RPL7 RPL8 RPS3 RPS5 RPS6 RPS7 RPS8 RPS9 RPLP0 RPSA"),
    "GO:0005198": ("Structural molecule activity",
                   "ACTB ACTG1 TUBA1A TUBB KRT8 KRT18 KRT19 VIM COL1A1 COL1A2 FN1 LMNA"),
    "GO:0016740": ("Transferase activity",
                   "ACACA ACACB MGAT1 MGAT2 FUT8 ST6GAL1 DNMT1 DNMT3A PRMT1 SETD7 GSK3B"),
    "GO:0016772": ("Transferase activity, transferring phosphorus-containing groups",
                   "CDK1 CDK2 PIK3CA PIK3CB PRKACA AKT1 MTOR MAP2K1 MAPK1 EGFR JAK1"),
    "GO:0016301": ("Kinase activity",
                   "CDK1 CDK2 CDK4 MAPK1 MAPK3 AKT1 PRKACA PIK3CA MTOR EGFR JAK1 PLK1 AURKA GSK3B"),
    "GO:0019901": ("Protein kinase binding",
                   "CDKN1A CDKN1B CDKN2A CCNB1 CCNE1 CCND1 SKP2 PIK3CA AKT1 IRS1 GRB2 SHC1"),
    "GO:0019904": ("Protein domain specific binding",
                   "GRB2 SHC1 SOS1 PIK3CA PTK2 PXN VCL YWHAZ YWHAQ YWHAB"),
    "GO:0008233": ("Peptidase activity",
                   "MMP2 MMP9 MMP14 CTSB CTSD CTSS CTSK CASP3 CASP8 CASP9 PSMB8 PSMB9 LGMN ADAM17"),
    "GO:0008234": ("Cysteine-type peptidase activity",
                   "CTSB CTSS CTSK CASP3 CASP7 CASP8 CASP9 CASP10 LGMN USP1 USP7"),
    "GO:0008237": ("Metallopeptidase activity",
                   "MMP2 MMP9 MMP14 MMP1 MMP3 MMP13 ADAM10 ADAM17 ACE2 LTA4H"),
    "GO:0070003": ("Threonine-type peptidase activity",
                   "PSMA1 PSMA2 PSMA3 PSMA4 PSMA5 PSMA6 PSMA7 PSMB1 PSMB2 PSMB5"),
    "GO:0004252": ("Serine-type endopeptidase activity",
                   "PRTN3 ELANE CTSG F2 F10 PLG PLAT PLAU TMPRSS2 KLK3"),
    "GO:0016491": ("Oxidoreductase activity",
                   "IDH2 MDH2 SDHA SDHB G6PD HMOX1 NQO1 CYP1A1 CYP1B1 DHCR7 HSD17B4"),
    "GO:0003824": ("Catalytic activity",
                   "HK1 HK2 GPI PFKM ALDOA GAPDH PGK1 ENO1 PKM LDHA ACACA ACACB"),
    "GO:0009055": ("Electron transfer activity",
                   "SDHB SDHD UQCRFS1 UQCRC1 UQCRC2 COX4I1 COX5A NDUFA1 NDUFS1 CYCS"),
    "GO:0015075": ("Ion transmembrane transporter activity",
                   "ATP1A1 ATP1B1 SLC9A1 SLC12A2 SLC8A1 SCN1A KCNQ1 CLCN3 TRPM7 SLC25A5"),
    "GO:0022857": ("Transmembrane transporter activity",
                   "SLC2A1 SLC2A3 SLC7A5 SLC7A11 SLC16A1 SLC16A3 SLC25A5 SLC4A1 ABCB1 ABCG2 ATP1A1"),
    "GO:0022804": ("Active transmembrane transporter activity",
                   "ATP1A1 ATP1B1 ABCB1 ABCG2 SLC2A1 SLC7A5 SLC25A5 SLC12A2"),
    "GO:0005215": ("Transporter activity",
                   "SLC2A1 SLC2A3 SLC7A5 SLC16A1 SLC25A5 ABCB1 ABCG2 ATP1A1 CLCN3 SCN1A"),
    "GO:0008324": ("Cation transmembrane transporter activity",
                   "ATP1A1 ATP1B1 SLC9A1 SLC8A1 SCN1A KCNQ1 TRPM7 TRPV1 CACNA1C"),
    "GO:0005261": ("Cation channel activity",
                   "SCN1A KCNQ1 KCNA1 TRPM7 TRPV1 CACNA1C HCN2 CLIC1"),
    "GO:0005085": ("Guanyl-nucleotide exchange factor activity",
                   "SOS1 SOS2 GRB2 PREX1 TIAM1 VAV1 DOCK2 ARHGEF1"),
    "GO:0005088": ("Ras guanyl-nucleotide exchange factor activity",
                   "SOS1 SOS2 GRB2 SHC1 PTPN11 RASA1 RASGRP1 RASGRF1"),
    "GO:0003924": ("GTPase activity",
                   "HRAS KRAS NRAS RAB5A RAB7A RAB11A RHOA RAC1 CDC42 GNAI1 GNAS"),
    "GO:0005525": ("GTP binding",
                   "HRAS KRAS NRAS RAB5A RAB7A RAB11A RHOA RAC1 CDC42 GNAI1 GNAS TUBA1A TUBB"),
    "GO:0005524": ("ATP binding",
                   "CDK1 CDK2 CDK4 MAPK1 AKT1 MTOR EGFR ABCB1 ABCG2 ATP1A1 ATP5F1A HSPA1A"),
    "GO:0046872": ("Metal ion binding",
                   "MT1A MT2A SOD1 CAT CA2 S100A8 S100A9 MMP2 MMP9"),
    "GO:0005509": ("Calcium ion binding",
                   "S100A8 S100A9 S100A12 CALR CANX CALM1 CALM2 CALM3 S100B ANXA1 ANXA2"),
    "GO:0046914": ("Transition metal ion binding",
                   "MT1A MT2A SOD1 CAT S100A8 S100A9 MMP2 MMP9 CYP1A1"),
    "GO:0005506": ("Iron ion binding",
                   "FTH1 FTL HMOX1 CYP1A1 CYP1B1 CYP2E1 SDHB ACO1 IREB2"),
    "GO:0051082": ("Unfolded protein binding",
                   "HSPA1A HSPA1B HSPA5 HSPA8 HSP90AA1 HSP90AB1 HSPD1 HSPE1 DNAJB1 CCT5"),
    "GO:0003755": ("Peptidyl-prolyl cis-trans isomerase activity",
                   "PPIA PPIB PPIC PPID FKBP1A FKBP4 FKBP5"),
    "GO:0008289": ("Lipid binding",
                   "APOE APOA1 FABP3 FABP4 FABP5 PLIN2 PLIN3 SREBF1 SREBF2 CETP"),
    "GO:0005543": ("Phospholipid binding",
                   "ANXA1 ANXA2 ANXA5 S100A8 S100A9 PEBP1 PLIN2 PLIN3"),
    "GO:0042393": ("Histone binding",
                   "EZH2 SUZ12 HMGN1 HMGN2 H2AFX H3F3A EP300 CREBBP KDM1A"),
    "GO:0019899": ("Enzyme binding",
                   "CDKN1A CDKN1B CCNB1 CCNE1 CCND1 SKP2 PIK3CA AKT1 IRS1 GRB2 SHC1 TP53"),
    "GO:0004867": ("Serine-type endopeptidase inhibitor activity",
                   "SERPINA1 SERPINA3 SERPINB5 SERPINE1 SERPINE2 SERPINF1 SERPING1 TIMP1 TIMP2"),
    "GO:0004866": ("Endopeptidase inhibitor activity",
                   "SERPINA1 SERPINE1 TIMP1 TIMP2 TIMP3 CSTB CSTA BIRC2 BIRC3 XIAP"),
    "GO:0030234": ("Enzyme regulator activity",
                   "SERPINA1 SERPINE1 TIMP1 TIMP2 BIRC2 BIRC3 XIAP SOCS1 SOCS3 PTPN11"),
    "GO:0060589": ("Nucleoside-triphosphatase regulator activity",
                   "RASA1 SOS1 GRB2 SHC1 PTPN11 RACGAP1 ARHGAP1 TBC1D1"),
    "GO:0005096": ("GTPase activator activity",
                   "RASA1 RASA2 ARHGAP1 ARHGAP5 TBC1D1 TBC1D4 RACGAP1"),
    "GO:0004872": ("Signaling receptor activity",
                   "EGFR ERBB2 TNFRSF1A IL6R IFNGR1 TLR4 TLR2 NOTCH1 TGFBR1 TGFBR2 KDR FLT1"),
    "GO:0004896": ("Cytokine receptor activity",
                   "IL6R IL6ST IL2RA IL2RB IL4R IL7R IL10RA IL10RB IFNGR1 IFNGR2 CSF2RA"),
    "GO:0005125": ("Cytokine activity",
                   "IL1B IL6 IL2 IFNG TNF IL4 IL10 IL12A IL12B IL18 CSF2 CXCL8 CCL2 LIF"),
    "GO:0005126": ("Cytokine receptor binding",
                   "IL6 IL2 IFNG TNF IL4 IL10 IL12A IL18 CSF2 CXCL8 CCL2 LIF"),
    "GO:0008009": ("Chemokine activity",
                   "CXCL1 CXCL2 CXCL3 CXCL5 CXCL6 CXCL8 CXCL9 CXCL10 CXCL11 CXCL12 CXCL13 CCL2 "
                   "CCL3 CCL4 CCL5 CCL19 CCL21"),
    "GO:0042379": ("Chemokine receptor binding",
                   "CCR1 CCR2 CCR5 CCR7 CXCR1 CXCR2 CXCR3 CXCR4 CXCR5 ACKR1 ACKR2"),
    "GO:0005178": ("Integrin binding",
                   "FN1 COL1A1 LAMB1 LAMC1 ICAM1 VCAM1 THBS1 FBN1 VTN SPP1"),
    "GO:0005179": ("Hormone activity",
                   "INS IGF1 IGF2 EGF FGF2 VEGFA PDGFB TGFB1 LEP ADIPOQ GH1 POMC"),
    "GO:0008083": ("Growth factor activity",
                   "EGF FGF2 FGF7 IGF1 IGF2 PDGFB VEGFA VEGFB HGF TGFB1 BMP2 CSF2"),
    "GO:0070851": ("Growth factor receptor binding",
                   "EGFR ERBB2 ERBB3 KDR FLT1 MET FGFR1 IGF1R PDGFRA PDGFRB NTRK1"),
    "GO:0004721": ("Phosphoprotein phosphatase activity",
                   "PTPRC PTPN6 PTPN11 PTEN PPP1CA PPP1CB PPP2CA PPP3CA DUSP1 DUSP6"),
    "GO:0016791": ("Phosphatase activity",
                   "PTPRC PTEN PPP1CA PPP2CA PPP3CA DUSP1 DUSP6 ACP1 ALPL"),
    "GO:0004553": ("Hydrolase activity, hydrolyzing O-glycosyl compounds",
                   "GBA GLA HEXA HEXB MAN2B1 MANBA GLB1 GALC IDUA"),
    "GO:0016787": ("Hydrolase activity",
                   "MMP2 MMP9 CTSB CTSD CTSS CASP3 CASP8 GBA HEXA ACHE APEH"),
    "GO:0008168": ("Methyltransferase activity",
                   "DNMT1 DNMT3A DNMT3B PRMT1 PRMT5 SETD7 EZH2 DOT1L METTL3 GNMT"),
    "GO:0016746": ("Acyltransferase activity",
                   "ACACA ACACB ACAT1 ACAT2 DGAT1 DGAT2 AGPAT2 LCAT SOAT1"),
    "GO:0008270": ("Zinc ion binding",
                   "MT1A MT2A S100A8 S100A9 MMP2 MMP9 CA2"),
}

# ---- KEGG 通路（常用条目，hsa ID） ----------------------------------------
_RAW_KEGG = {
    "hsa04110": ("Cell cycle",
                 "CCNA2 CCNB1 CCNB2 CCND1 CCND2 CCND3 CCNE1 CDK1 CDK2 CDK4 CDK6 CDC25A CDC25B "
                 "CDC25C CDKN1A CDKN1B CDKN2A E2F1 RB1 TP53 MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 CDC6 "
                 "CDT1 SKP2 PLK1 AURKA AURKB BUB1 MAD2L1 CDC20 ESPL1 PTTG1 MDM2 GSK3B ATM ATR "
                 "CHEK1 CHEK2"),
    "hsa04115": ("p53 signaling pathway",
                 "TP53 MDM2 CDKN1A BAX BCL2 BBC3 PMAIP1 GADD45A GADD45B GADD45G ATM ATR CHEK1 "
                 "CHEK2 FAS TNFRSF10A TNFRSF10B ZMAT3 SFN SERPINE1 PTGS2 CCNG1 RRM2B TIGAR "
                 "TP53I3 BID APAF1 CASP3 CASP8 CASP9 CYCS PRKDC XPC DDB2"),
    "hsa04210": ("Apoptosis",
                 "BAX BAK1 BAD BCL2 BCL2L1 MCL1 BID BIRC2 BIRC3 APAF1 CYCS CASP3 CASP6 CASP7 "
                 "CASP8 CASP9 CASP10 FADD FAS FASLG TNF TNFRSF1A TRADD RIPK1 NFKB1 IKBKB MAP3K5 "
                 "MAPK8 MAPK9 JUN PIK3CA AKT1"),
    "hsa04310": ("Wnt signaling pathway",
                 "WNT2 WNT3A WNT4 WNT5A WNT7A WNT10B FZD1 FZD2 FZD4 FZD7 LRP5 LRP6 DVL1 DVL2 "
                 "DVL3 APC AXIN1 AXIN2 CTNNB1 TCF7 TCF7L1 TCF7L2 LEF1 GSK3B CSNK1A1 BTRC SFRP1 "
                 "DKK1 MYC CCND1"),
    "hsa04350": ("TGF-beta signaling pathway",
                 "TGFB1 TGFB2 TGFB3 TGFBR1 TGFBR2 TGFBR3 SMAD1 SMAD2 SMAD3 SMAD4 SMAD5 SMAD6 "
                 "SMAD7 SMAD9 ACVR1 ACVR2A ACVR2B BMP2 BMP4 BMP6 BMPR1A BMPR1B BMPR2 ID1 ID2 "
                 "ID3 ID4 JUN JUNB FOS SKIL SERPINE1 INHBA"),
    "hsa04330": ("Notch signaling pathway",
                 "NOTCH1 NOTCH2 NOTCH3 NOTCH4 DLL1 DLL3 DLL4 JAG1 JAG2 HES1 HES2 HES4 HES5 HES6 "
                 "HES7 HEY1 HEY2 HEYL LFNG RFNG MFNG PSEN1 PSEN2 NCSTN APH1A APH1B PSENEN MAML1 "
                 "MAML2 MAML3 RBPJ DTX1 NUMB NUMBL"),
    "hsa04340": ("Hedgehog signaling pathway",
                 "SHH IHH DHH PTCH1 PTCH2 SMO GLI1 GLI2 GLI3 SUFU BOC CDON GAS1 HHIP WNT2B "
                 "WNT5A BMP2 BMP4 PRKACA GSK3B CSNK1A1"),
    "hsa04360": ("Axon guidance",
                 "NTN1 DCC NEO1 SLIT1 SLIT2 SLIT3 ROBO1 ROBO2 ROBO3 SEMA3A SEMA3B SEMA3C "
                 "SEMA3F PLXNA1 PLXNA2 PLXNA3 NRP1 NRP2 EPHA2 EPHA4 EFNA1 EFNB2 RHOA RAC1 CDC42 "
                 "ABLIM1 LIMK1 CFL1 GSK3B"),
    "hsa04010": ("MAPK signaling pathway",
                 "HRAS KRAS NRAS RAF1 BRAF MAP2K1 MAP2K2 MAP2K3 MAP2K4 MAP2K6 MAPK1 MAPK3 MAPK7 "
                 "MAPK8 MAPK9 MAPK10 MAPK11 MAPK12 MAPK13 MAPK14 MAPKAPK2 ELK1 JUN FOS MYC "
                 "RPS6KA1 MKNK1 DUSP1 DUSP6 SPRY1 SPRY2 PTPN11 SOS1 GRB2 SHC1 EGFR"),
    "hsa04012": ("ErbB signaling pathway",
                 "EGFR ERBB2 ERBB3 ERBB4 HBEGF EGF TGFA AREG EREG NRG1 NRG2 GRB2 SOS1 HRAS "
                 "MAP2K1 MAPK1 STAT3 PIK3CA AKT1 PLCG1 SRC CBL EGF"),
    "hsa04014": ("Ras signaling pathway",
                 "HRAS KRAS NRAS RAF1 BRAF MAP2K1 MAPK1 RASA1 RASA2 SOS1 SOS2 GRB2 SHC1 PTPN11 "
                 "EGFR ERBB2 FGFR1 IGF1R INSR KDR FLT1 MET PIK3CA AKT1"),
    "hsa04015": ("Rap1 signaling pathway",
                 "RAP1A RAP1B RAPGEF1 RAPGEF3 RAPGEF4 TERF2? PLCE1 PIK3CA AKT1 MAPK1 MAPK3 "
                 "RAC1 CDC42 RHOA ITGB1 ITGA5 TLN1 PXN SRC FYN"),
    "hsa04020": ("Calcium signaling pathway",
                 "CACNA1C CACNA1D CACNA1A ATP2A1 ATP2A2 ATP2B1 SLC8A1 RYR1 RYR2 RYR3 ITPR1 "
                 "ITPR2 ITPR3 PRKACA PRKCA CAMK2A CAMK2B CALM1 CALM2 CALM3 PLCB1 PLCG1 GNAS"),
    "hsa04022": ("cGMP-PKG signaling pathway",
                 "GUCA1A GUCY1A1 GUCY1B1 PRKG1 PRKG2 PDE5A PDE1A PDE1C ADRB1 ADRB2 ATP2B1 "
                 "SLC8A1 PLCB1 PRKACA CREB1"),
    "hsa04024": ("cAMP signaling pathway",
                 "ADCY1 ADCY2 ADCY3 ADCY5 ADCY6 ADRB1 ADRB2 GNAS GNAI1 PRKACA PRKACB PRKACG "
                 "CREB1 CREM ATF1 PDE3A PDE4A PDE4B RAPGEF3 EPAC? PIK3CA AKT1 MAPK1"),
    "hsa04062": ("Chemokine signaling pathway",
                 "CCL2 CCL3 CCL4 CCL5 CCL7 CCL8 CCL11 CCL13 CCL17 CCL19 CCL20 CCL21 CCL22 "
                 "CXCL1 CXCL2 CXCL3 CXCL5 CXCL8 CXCL9 CXCL10 CXCL11 CXCL12 CXCL13 CXCL16 "
                 "CCR1 CCR2 CCR5 CCR7 CXCR1 CXCR2 CXCR3 CXCR4 CXCR5 CXCR6 ACKR1 PIK3CA AKT1 "
                 "MAPK1 MAPK3 PRKACA PLCB1"),
    "hsa04064": ("NF-kappa B signaling pathway",
                 "TNF TNFRSF1A TRADD RIPK1 TRAF2 TRAF6 IKBKB IKBKG NFKB1 NFKB2 RELA RELB NFKBIA "
                 "NFKBIB BIRC2 BIRC3 MAP3K7 TLR4 MYD88 IRAK1 IRAK4 CD40 CD40LG LTA LTBR IL1B "
                 "IL1R1 CCL19 CXCL12"),
    "hsa04066": ("HIF-1 signaling pathway",
                 "HIF1A EPAS1 ARNT VEGFA VEGFB KDR FLT1 SLC2A1 SLC2A3 HK1 HK2 PGK1 ENO1 PKM "
                 "LDHA PDK1 PDK3 EGLN1 EGLN2 EGLN3 VHL EP300 CREBBP STAT3 MTOR"),
    "hsa04068": ("FoxO signaling pathway",
                 "FOXO1 FOXO3 FOXO4 FOXO6 AKT1 AKT2 AKT3 PIK3CA PTEN SIRT1 SIRT2 SKP2 MDM2 "
                 "CDKN1A CDKN1B GADD45A SOD1 SOD2 CAT TNFSF10 BCL2L11 IGF1 IGF1R INSR"),
    "hsa04120": ("Ubiquitin mediated proteolysis",
                 "UBA1 UBE2A UBE2B UBE2C UBE2D1 UBE2D3 UBE2E1 UBE2G1 UBE2I UBE2L3 UBE2M UBE2N "
                 "UBE2V1 UBE3A UBE4B MDM2 SKP2 FBXW7 FBXO7 CUL1 CUL2 CUL3 CUL4A CUL4B RBX1 "
                 "ANAPC1 ANAPC2 CDC16 CDC23 CDC27"),
    "hsa04140": ("Autophagy - animal",
                 "MAP1LC3A MAP1LC3B MAP1LC3C GABARAP GABARAPL1 ATG3 ATG4A ATG4B ATG5 ATG7 ATG9A "
                 "ATG10 ATG12 ATG13 ATG14 ATG16L1 BECN1 PIK3C3 PIK3R4 UVRAG ULK1 ULK2 RB1CC1 "
                 "SQSTM1 NBR1 OPTN RAB7A STX17"),
    "hsa04141": ("Protein processing in endoplasmic reticulum",
                 "HSPA5 HSPA1A HSPA1B HSPA8 HSP90B1 CALR CANX PDIA3 PDIA4 PDIA6 SEC61A1 SEC61B "
                 "SEC61G SEC62 SEC63 SRP54 SRP9 SRP14 UGGT1 UGGT2 EDEM1 EDEM2 EDEM3 SEL1L "
                 "HERPUD1 DNAJB9 DNAJC3 ATF4 DDIT3 XBP1 ERN1 ATF6 EIF2AK3"),
    "hsa04142": ("Lysosome",
                 "CTSB CTSD CTSS CTSK CTSA CTSZ LAMP1 LAMP2 NPC1 NPC2 HEXA HEXB GBA GLA GALC GLB1 "
                 "MAN2B1 MANBA ARSA ARSB IDUA SGSH NAGLU GUSB SMPD1 ASAH1"),
    "hsa04144": ("Endocytosis",
                 "CLTC CLTA CLTB AP2A1 AP2A2 AP2B1 AP2M1 EPS15 EPN1 DNM1 DNM2 RAB5A RAB4A ARF6 "
                 "CAV1 CAV2 CAV3 ITSN1 SNAP91 PICALM EHD1 EHD2 ACAP1"),
    "hsa04145": ("Phagosome",
                 "ITGAM ITGB2 ITGAV ITGB3 ITGA5 FCGR1A FCGR2A FCGR3A TLR2 TLR4 TLR6 TUBB TUBA1A "
                 "ACTB ACTG1 RAB5A RAB7A RAB11A LAMP1 CTSS CTSK CTSC MRC1 MARCO CD14 CD36 SCARB1"),
    "hsa04146": ("Peroxisome",
                 "PEX1 PEX3 PEX5 PEX6 PEX7 PEX10 PEX12 PEX13 PEX14 PEX19 CAT ACOX1 ACOX2 HSD17B4 "
                 "SCP2 SLC25A17 ABCD1 ABCD2 ABCD3 SOD1"),
    "hsa04150": ("mTOR signaling pathway",
                 "MTOR RPTOR RICTOR MLST8 AKT1 AKT2 TSC1 TSC2 RHEB RPS6KB1 RPS6KB2 EIF4EBP1 "
                 "RPS6 PIK3CA PIK3CB PTEN PDPK1 PRKAA1 PRKAA2 TSC1 IGF1 INS"),
    "hsa04151": ("PI3K-Akt signaling pathway",
                 "PIK3CA PIK3CB PIK3CD PIK3CG PIK3R1 PIK3R2 PIK3R3 AKT1 AKT2 AKT3 PTEN PDPK1 "
                 "MTOR GSK3B FOXO1 FOXO3 CDKN1A BCL2 BCL2L1 MCL1 BAD CASP9 EGFR ERBB2 IGF1R "
                 "INSR KDR FLT1 MET ITGB1 ITGA5 FN1 COL1A1 VEGFA"),
    "hsa04152": ("AMPK signaling pathway",
                 "PRKAA1 PRKAA2 PRKAB1 PRKAB2 PRKAG1 PRKAG2 PRKAG3 STK11 CAB39 STRADA LKB1? "
                 "CAMKK2 AKT1 MTOR TSC2 RPTOR PPARGC1A CPT1A ACACA ACACB HMGCR SLC2A4 FOXO3 "
                 "CDKN1A"),
    "hsa04216": ("Ferroptosis",
                 "SLC7A11 SLC3A2 GPX4 ACSL4 LPCAT3 ALOX15 ALOX12 NFE2L2 HMOX1 FTH1 FTL TFRC "
                 "SLC40A1 GCLM GCLC SLC7A5"),
    "hsa04217": ("Necroptosis",
                 "TNF TNFRSF1A RIPK1 RIPK3 MLKL ZBP1 TRADD FADD CASP8 CYLD TRAF2 TRAF5 NFKB1 "
                 "BIRC2 BIRC3 MAPK8 MAPK9 PELI1"),
    "hsa04218": ("Cellular senescence",
                 "TP53 CDKN1A CDKN2A RB1 E2F1 MDM2 ATM ATR CHEK1 CHEK2 MAPK14 MAPK1 IL6 CXCL8 "
                 "CCL2 SERPINE1 IGFBP3 IGFBP7 CDKN2B GLB1 LMNA"),
    "hsa04370": ("VEGF signaling pathway",
                 "VEGFA VEGFB VEGFC KDR FLT1 FLT4 NRP1 PIK3CA AKT1 PLCG1 PRKCA MAPK1 MAPK3 "
                 "SRC FYN PTK2 PXN CASP9 BCL2 NOS3 PTGS2"),
    "hsa04380": ("Osteoclast differentiation",
                 "TNFRSF11A TNFSF11 CSF1 CSF1R SPI1 NFKB1 NFATC1 RELA JUN FOS MAPK14 MAPK1 "
                 "MAPK3 PIK3CA AKT1 SRC CTSK MMP9 ACP5 TRAP? ATP6V0D2"),
    "hsa04510": ("Focal adhesion",
                 "ITGB1 ITGB2 ITGB3 ITGB4 ITGB5 ITGB6 ITGB7 ITGB8 ITGA1 ITGA2 ITGA3 ITGA4 ITGA5 "
                 "ITGA6 ITGAV ITGAX PTK2 PXN VCL TLN1 ILK PARVA ACTN1 ZYX SRC FYN LAMB1 LAMC1 "
                 "FN1 COL1A1 COL4A1 THBS1 VEGFA"),
    "hsa04512": ("ECM-receptor interaction",
                 "COL1A1 COL1A2 COL2A1 COL3A1 COL4A1 COL4A2 COL5A1 COL6A1 COL6A2 COL6A3 COL9A1 "
                 "COL11A1 FN1 FBN1 LAMA1 LAMA2 LAMA4 LAMB1 LAMB2 LAMC1 LAMC2 THBS1 THBS2 VTN "
                 "SPP1 ITGB1 ITGA5 ITGAV CD44 HSPG2"),
    "hsa04514": ("Cell adhesion molecules",
                 "CDH1 CDH2 CDH3 CDH4 CDH5 CDH6 CDH11 CDH15 CTNNA1 CTNNB1 CTNND1 JUP ITGB1 "
                 "ITGB2 ITGA4 ITGAL ITGAM ICAM1 ICAM2 ICAM3 VCAM1 PECAM1 NCAM1 NCAM2 L1CAM "
                 "MADCAM1 SELL SELP SELE CD2 CD28 CD40 CD80 CD86 CTLA4 PDCD1 CD274"),
    "hsa04520": ("Adherens junction",
                 "CDH1 CDH2 CDH3 CTNNA1 CTNNA2 CTNNB1 CTNND1 JUP PTPN6 PTPN11 SRC FYN YES1 "
                 "IQGAP1 RHOA RAC1 CDC42 ACTN1 ACTN2 VCL EGFR ERBB2 MET"),
    "hsa04530": ("Tight junction",
                 "TJP1 TJP2 TJP3 OCLN CLDN1 CLDN2 CLDN3 CLDN4 CLDN5 CLDN7 CLDN8 CLDN10 CLDN11 "
                 "JAM1 JAM2 JAM3 CDH1 CTNNA1 CTNNB1 MYH9 MYL6 PRKACA PRKCA SRC"),
    "hsa04540": ("Gap junction",
                 "GJA1 GJA3 GJA4 GJA5 GJA8 GJB1 GJB2 GJB3 GJB4 GJB5 GJC1 GJD2 TUBB TUBA1A "
                 "ADCY1 ADCY3 PRKACA PRKACB MAPK1 MAPK3 SRC"),
    "hsa04550": ("Signaling pathways regulating pluripotency of stem cells",
                 "POU5F1 SOX2 NANOG KLF4 MYC LIN28A LIN28B WNT3A WNT5A FZD1 LRP5 LRP6 CTNNB1 "
                 "BMP4 SMAD1 SMAD4 JAK1 STAT3 LIF IL6ST PIK3CA AKT1 MAPK1 FGF2 FGFR1"),
    "hsa04610": ("Complement and coagulation cascades",
                 "C1QA C1QB C1QC C1R C1S C2 C3 C4A C4B C5 C5AR1 C6 C7 C8A C8B C9 CFB CFD CFH CFI "
                 "CFP CR1 CR2 SERPING1 CD46 CD55 CD59 MASP1 MASP2 F2 F3 F5 F7 F8 F9 F10 F11 F12 "
                 "F13A1 FGA FGB FGG PLG PLAT PLAU SERPINE1"),
    "hsa04612": ("Antigen processing and presentation",
                 "HLA-A HLA-B HLA-C HLA-E HLA-F HLA-DRA HLA-DRB1 HLA-DPA1 HLA-DPB1 HLA-DQA1 "
                 "HLA-DQB1 B2M TAP1 TAP2 TAPBP PSMB8 PSMB9 PSMB10 PSMA2 PSMA3 PSMA4 PSMA5 PSMB1 "
                 "PSMB2 CD74 CIITA CTSB CTSD LGMN IFNG IFNGR1 JAK1 JAK2 STAT1"),
    "hsa04620": ("Toll-like receptor signaling pathway",
                 "TLR1 TLR2 TLR3 TLR4 TLR5 TLR6 TLR7 TLR8 TLR9 MYD88 TICAM1 TIRAP IRAK1 IRAK2 "
                 "IRAK4 TRAF3 TRAF6 TANK IKBKB IKBKE NFKB1 RELA MAP3K7 MAPK8 MAPK14 PIK3CA AKT1 "
                 "IRF3 IRF5 IRF7 IFNB1 IFNA1 CXCL8 IL6 TNF CD80 CD86"),
    "hsa04621": ("NOD-like receptor signaling pathway",
                  "NOD1 NOD2 NLRC3 NLRP1 NLRP3 NLRP6 NLRP12 NLRC4 PYCARD CASP1 IL1B IL18 CARD6 "
                  "CARD9 RIPK2 MAP3K7 IKBKB NFKB1 RELA MAPK14 MAPK1 CXCL8 IL6 TNF"),
    "hsa04622": ("RIG-I-like receptor signaling pathway",
                  "DDX58 IFIH1 DHX58 MAVS TBK1 IKBKE IRF3 IRF7 IFNB1 IFNA1 NFKB1 RELA TRAF3 "
                  "TRAF6 TICAM1 AZI2 TMEM173 ISG15 MX1 OAS1"),
    "hsa04623": ("Cytosolic DNA-sensing pathway",
                  "CGAS MB21D1? TMEM173 STING1 TBK1 IKBKE IRF3 IRF7 IFNB1 NFKB1 RELA ZBP1 "
                  "AIM2 PYCARD CASP1 IL1B IL18 POLR3A POLR3B"),
    "hsa04630": ("JAK-STAT signaling pathway",
                  "JAK1 JAK2 JAK3 TYK2 STAT1 STAT2 STAT3 STAT4 STAT5A STAT5B STAT6 SOCS1 SOCS2 "
                  "SOCS3 SOCS4 SOCS5 PTPN2 PTPN6 PIAS1 PIAS3 IL2 IL2RA IL4 IL4R IL6 IL6ST "
                  "IL7 IL7R IL10 IL10RA IFNG IFNGR1 CSF2 CSF2RB EPO EPOR GH1 GHR PRL PRLR"),
    "hsa04640": ("Hematopoietic cell lineage",
                  "CD34 KIT FLT3 CSF1R CSF2RA CSF3R IL3RA IL7R IL2RA MS4A1 CD19 CD22 CD33 CD38 "
                  "CD44 CD14 ITGAM ITGAX PTPRC GATA1 GATA2 SPI1 RUNX1 CEBPA MPO ELANE LYZ"),
    "hsa04650": ("Natural killer cell mediated cytotoxicity",
                  "NKG7 GNLY GZMA GZMB GZMH GZMK GZMM PRF1 KLRD1 KLRF1 KLRK1 NCR1 NCR3 FCGR3A "
                  "NCAM1 SH2D1A CD160 IFNG TNF FASLG ICAM1 HLA-A HLA-B HLA-C KIR2DL1 KIR2DL2 "
                  "KIR2DL3 KIR3DL1 KIR3DL2 PIK3CA AKT1 MAPK1"),
    "hsa04658": ("Th1 and Th2 cell differentiation",
                  "CD4 IFNG IL12A IL12B IL12RB1 IL12RB2 STAT1 STAT4 TBX21 IL2 IL2RA IL4 IL4R "
                  "IL13 STAT6 GATA3 JAK1 JAK2 RUNX3 IFNG IFNG? T-bet?"),
    "hsa04659": ("Th17 cell differentiation",
                  "CD4 IL17A IL17F IL17RA IL23A IL23R IL6 IL6R IL6ST IL21 IL21R IL22 TGFB1 "
                  "TGFBR1 TGFBR2 STAT3 RORC JAK2 JAK1 SOCS3 RUNX1 BATF IRF4"),
    "hsa04660": ("T cell receptor signaling pathway",
                  "CD3D CD3E CD3G CD3D CD247 CD4 CD8A CD28 CD80 CD86 CTLA4 LCK FYN ZAP70 LAT "
                  "PLCG1 PIK3CA AKT1 GRB2 SOS1 HRAS MAP2K1 MAPK1 JUN FOS NFKB1 NFATC1 NFATC2 "
                  "PPP3CA IL2 IL2RA"),
    "hsa04662": ("B cell receptor signaling pathway",
                  "CD19 CD79A CD79B MS4A1 BLNK SYK BTK PLCG2 PIK3CD PIK3CA VAV1 GRB2 SOS1 HRAS "
                  "MAP2K1 MAPK1 JUN NFATC1 NFKB1 RELA PTPN6 LYN HCK CARD11 BCL10 MALT1"),
    "hsa04664": ("Fc epsilon RI signaling pathway",
                  "FCER1A MS4A2 FCER1G SYK LYN FYN PLCG1 PIK3CA AKT1 MAP2K1 MAPK1 MAPK3 MAPK8 "
                  "MAPK14 PRKCA TNF IL4 IL5 IL13"),
    "hsa04666": ("Fc gamma R-mediated phagocytosis",
                  "FCGR1A FCGR2A FCGR2B FCGR3A FCGR3B FCER1G SYK LYN HCK PLCG1 PIK3CA AKT1 "
                  "MAPK1 MAPK3 RAC1 CDC42 ARPC1B ACTR2 ACTR3 WASL LIMK1 CFL1"),
    "hsa04668": ("TNF signaling pathway",
                  "TNF TNFRSF1A TNFRSF1B TRADD TRAF2 RIPK1 IKBKB NFKB1 RELA NFKBIA MAP3K7 "
                  "MAP2K3 MAP2K4 MAPK8 MAPK9 MAPK14 JUN FOS ATF2 AP1? CCL2 CXCL1 CXCL8 IL6 "
                  "ICAM1 VCAM1 BIRC2 BIRC3"),
    "hsa04670": ("Leukocyte transendothelial migration",
                  "ITGAM ITGB2 ITGAL ITGA4 ITGB1 ICAM1 VCAM1 PECAM1 CD99 JAM1 JAM2 JAM3 CLDN1 "
                  "CLDN3 CLDN5 OCLN TJP1 MMP2 MMP9 CXCL12 CXCR4 RHOA RAC1 CDC42 PXN PTK2 VCL"),
    "hsa04910": ("Insulin signaling pathway",
                  "INS INSR IRS1 IRS2 PIK3CA AKT1 AKT2 SLC2A4 RPS6KB1 MTOR GSK3B FOXO1 PDPK1 "
                  "PRKAA1 TSC2 RHEB SHC1 GRB2 SOS1 HRAS MAP2K1 MAPK1 PTPN11"),
    "hsa04915": ("Estrogen signaling pathway",
                  "ESR1 ESR2 SHC1 GRB2 SOS1 HRAS MAP2K1 MAPK1 PIK3CA AKT1 JUN FOS MMP2 MMP9 "
                  "GPER1 PGR SP1 CTSD PRLR"),
    "hsa04919": ("Thyroid hormone signaling pathway",
                  "THRA THRB THR? MED1 MED12 NCOA1 NCOA2 NCOR1 NCOR2 ESR1 ESR2 PIK3CA AKT1 "
                  "MTOR MAPK1 MAPK3 SLC2A1 SLC2A4 ATP1A1"),
    "hsa04920": ("Adipocytokine signaling pathway",
                  "LEP LEPR ADIPOQ ADIPOR1 ADIPOR2 TNF TNFRSF1A IL6 IL6R JAK1 JAK2 STAT3 SOCS3 "
                  "PTPN11 PRKAA1 PRKAA2 CPT1A ACACA ACACB SLC2A4 PPARGC1A"),
    "hsa04931": ("Insulin resistance",
                  "INS INSR IRS1 IRS2 PIK3CA AKT1 PTPN11 SOCS1 SOCS3 JAK1 JAK2 STAT3 TNF "
                  "TNFRSF1A IL6 IL6R SLC2A4 G6PC PCK1 FOXO1 PPARGC1A PRKAA1 PRKAA2"),
    "hsa04932": ("Non-alcoholic fatty liver disease",
                  "TNF TNFRSF1A IL1B IL6 CXCL8 CCL2 TGFB1 COL1A1 COL1A2 BCL2 BAX CASP3 CASP8 "
                  "CASP9 APAF1 CYCS IRS1 IRS2 PIK3CA AKT1 PRKAA1 PRKAA2 SREBF1 HMGCR"),
    "hsa04933": ("AGE-RAGE signaling pathway in diabetic complications",
                  "AGER AGT REN ACE AGTR1 EDN1 CXCL8 CCL2 ICAM1 VCAM1 SELE TGFB1 COL1A1 "
                  "COL3A1 FN1 MMP2 STAT3 JAK2 PIK3CA AKT1 MAPK1 MAPK3 NFKB1 RELA"),
    "hsa04060": ("Cytokine-cytokine receptor interaction",
                  "IL1A IL1B IL1R1 IL2 IL2RA IL2RB IL4 IL4R IL6 IL6R IL6ST IL7 IL7R IL10 IL10RA "
                  "IL10RB IL12A IL12B IL12RB1 IL13 IL13RA1 IL15 IL15RA IL17A IL17RA IL18 IL18R1 "
                  "IL21 IL21R IL23A IL23R IFNG IFNGR1 IFNGR2 TNF TNFRSF1A TNFRSF1B CXCL8 CXCR1 "
                  "CXCR2 CCL2 CCR2 CXCL12 CXCR4 CSF1 CSF1R CSF2 CSF2RA CSF3 CSF3R VEGFA KDR"),
    "hsa04080": ("Neuroactive ligand-receptor interaction",
                  "ADRB1 ADRB2 ADRB3 CHRM1 CHRM2 CHRM3 DRD1 DRD2 DRD3 DRD4 DRD5 GABRA1 GABRB2 "
                  "GABRG2 GRIN1 GRIN2A GRIN2B HTR1A HTR2A HTR2C SLC6A4 OPRM1 OPRD1 OPRK1 "
                  "NPY1R NPY2R AGTR1 BDKRB1 BDKRB2 GLP1R"),
    "hsa04114": ("Oocyte meiosis",
                  "CDK1 CDK2 CCNB1 CCNB2 CCNA2 CDC20 CDC25A CDC25B CDC25C BUB1 MAD2L1 PLK1 "
                  "AURKA AURKB MOS? MAPK1 MAPK3 MOS SMC1A SMC3 RAD21 STAG1 STAG2 REC8"),
    "hsa03010": ("Ribosome",
                  "RPL3 RPL4 RPL5 RPL6 RPL7 RPL8 RPL9 RPL10 RPL11 RPL12 RPL13 RPL13A RPL14 "
                  "RPL15 RPL17 RPL18 RPL18A RPL19 RPL21 RPL22 RPL23 RPL23A RPL24 RPL26 RPL27 "
                  "RPL27A RPL28 RPL29 RPL30 RPL31 RPL32 RPL34 RPL35 RPL36 RPL37 RPL38 RPL39 "
                  "RPS3 RPS4X RPS5 RPS6 RPS7 RPS8 RPS9 RPS10 RPS11 RPS12 RPS13 RPS14 RPS15 "
                  "RPS15A RPS16 RPS17 RPS18 RPS19 RPS20 RPS21 RPS23 RPS24 RPS25 RPS26 RPS27 "
                  "RPS27A RPS28 RPS29 RPLP0 RPLP1 RPLP2 RPSA"),
    "hsa03040": ("Spliceosome",
                  "SNRNP70 SNRPB SNRPB2 SNRPD1 SNRPD2 SNRPD3 SNRPE SNRPF SNRPG SNRNP200 "
                  "SF3A1 SF3A2 SF3A3 SF3B1 SF3B2 SF3B3 SF3B4 SF3B5 SF3B6 PRPF3 PRPF4 PRPF6 "
                  "PRPF8 PRPF19 PRPF31 U2AF1 U2AF2 SRSF1 SRSF2 SRSF3 HNRNPA1"),
    "hsa03050": ("Proteasome",
                  "PSMA1 PSMA2 PSMA3 PSMA4 PSMA5 PSMA6 PSMA7 PSMB1 PSMB2 PSMB3 PSMB4 PSMB5 "
                  "PSMB6 PSMB7 PSMB8 PSMB9 PSMB10 PSMC1 PSMC2 PSMC3 PSMC4 PSMC5 PSMC6 PSMD1 "
                  "PSMD2 PSMD3 PSMD4 PSMD6 PSMD7 PSMD8 PSMD11 PSMD12 PSMD13 PSME1 PSME2"),
    "hsa03030": ("DNA replication",
                  "MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 PCNA RFC1 RFC2 RFC3 RFC4 RFC5 RPA1 RPA2 RPA3 "
                  "POLA1 POLA2 POLD1 POLD2 POLD3 POLE POLE2 POLE3 POLE4 PRIM1 PRIM2 LIG1 FEN1 "
                  "RNASEH1 RNASEH2A"),
    "hsa03410": ("Base excision repair",
                  "OGG1 MUTYH NTHL1 APEX1 XRCC1 POLB PARP1 PARP2 LIG1 LIG3 FEN1 PCNA POLD1 "
                  "POLE MBD4 TDG UNG SMUG1 NEIL1 NEIL2 NEIL3 APLF"),
    "hsa03420": ("Nucleotide excision repair",
                  "XPA XPC XPD? ERCC2 ERCC3 ERCC4 ERCC5 ERCC6 ERCC8 DDB1 DDB2 RAD23A RAD23B "
                  "CETN2 XAB2 RPA1 RPA2 RPA3 PCNA RFC1 LIG1 POLH POLR2A"),
    "hsa03430": ("Mismatch repair",
                  "MSH2 MSH3 MSH6 MLH1 MLH3 PMS1 PMS2 EXO1 RPA1 RPA2 RPA3 PCNA RFC1 LIG1 "
                  "POLD1 POLD2 POLD3 RPA1"),
    "hsa03440": ("Homologous recombination",
                  "RAD51 RAD51B RAD51C RAD51D RAD52 RAD54B RAD54L BRCA1 BRCA2 PALB2 BARD1 XRCC2 "
                  "XRCC3 NBN MRE11 RAD50 RPA1 RPA2 RPA3 BLM"),
    "hsa03450": ("Non-homologous end-joining",
                  "XRCC5 XRCC6 PRKDC LIG4 XRCC4 NHEJ1 DCLRE1C ARTEMIS? POLL POLM TDP1 APLF"),
    "hsa03460": ("Fanconi anemia pathway",
                  "FANCA FANCB FANCC FANCD2 FANCE FANCF FANCG FANCI FANCL FANCM BRCA1 BRCA2 "
                  "PALB2 RAD51 RAD51C BRIP1 ERCC4 SLX4"),
    "hsa00010": ("Glycolysis / Gluconeogenesis",
                  "HK1 HK2 HK3 GPI PFKL PFKM PFKP ALDOA ALDOB ALDOC TPI1 GAPDH PGK1 PGAM1 PGAM2 "
                  "ENO1 ENO2 ENO3 PKM PKLR LDHA LDHB PCK1 PCK2 FBP1 FBP2 GCK ALDH1A1 ALDH2"),
    "hsa00020": ("Citrate cycle (TCA cycle)",
                  "CS ACO1 ACO2 IDH1 IDH2 IDH3A IDH3B OGDH DLST SUCLA2 SUCLG1 SDHA SDHB SDHC "
                  "SDHD FH MDH1 MDH2 PC PCK1 PCK2 DLD DLAT PDHA1 PDHB"),
    "hsa00030": ("Pentose phosphate pathway",
                  "G6PD PGLS PGD RPE RPIA TKT TKTL1 PRPS1 PRPS2 RBKS RPE"),
    "hsa00071": ("Fatty acid degradation",
                  "ACADL ACADM ACADS ACADVL HADHA HADHB ECHS1 ECI1 ECI2 ACAA1 ACAA2 ACAT1 CPT1A "
                  "CPT1B CPT2 SLC25A20 ACAD8 ACAD9 ACAD10"),
    "hsa00190": ("Oxidative phosphorylation",
                  "NDUFA1 NDUFA2 NDUFA3 NDUFA4 NDUFA5 NDUFA6 NDUFA7 NDUFA8 NDUFA9 NDUFA10 "
                  "NDUFB1 NDUFB2 NDUFB3 NDUFB4 NDUFB5 NDUFB6 NDUFB7 NDUFB8 NDUFB9 NDUFB10 "
                  "NDUFS1 NDUFS2 NDUFS3 NDUFS4 NDUFS5 NDUFS6 NDUFS7 NDUFS8 NDUFV1 NDUFV2 SDHA "
                  "SDHB SDHC SDHD UQCRC1 UQCRC2 UQCRFS1 UQCRB UQCRQ UQCRH COX4I1 COX4I2 COX5A "
                  "COX5B COX6A1 COX6B1 COX6C COX7A2 COX7B COX7C COX8A ATP5F1A ATP5F1B ATP5F1C "
                  "ATP5F1D ATP5F1E ATP5PB ATP5PD ATP5PF ATP5MC1 ATP5MC2 ATP5MC3"),
    "hsa00230": ("Purine metabolism",
                  "HPRT1 APRT ADA AMPD1 AMPD2 AMPD3 IMPDH1 IMPDH2 GMPS ATIC ADSL ADSS ADK AK1 "
                  "AK2 GUK1 NME1 NME2 NME3 NME4 RRM1 RRM2 PRPS1 PRPS2 PPAT PAICS GART PFAS "
                  "PDE1A PDE4A PDE5A"),
    "hsa00240": ("Pyrimidine metabolism",
                  "UMPS CAD DHODH CTPS1 CTPS2 NME1 NME2 RRM1 RRM2 DPYD DPYS UPP1 UPP2 TYMP "
                  "TK1 TK2 UCK1 UCK2 CMPK1 NME1"),
    "hsa00500": ("Starch and sucrose metabolism",
                  "GYS1 GYS2 GYG1 PYGL PYGM AGL GAA MGAM SI TREH HK1 GCK UGP2 UGDH"),
    "hsa00510": ("N-Glycan biosynthesis",
                  "ALG1 ALG2 ALG3 ALG5 ALG6 ALG8 ALG9 ALG11 ALG12 ALG13 ALG14 DOLK DOLPP1 "
                  "DPAGT1 GPT1? MGAT1 MGAT2 MGAT4A MGAT4B MGAT5 B4GALT1 FUT8"),
    "hsa00520": ("Amino sugar and nucleotide sugar metabolism",
                  "GFPT1 GFPT2 GNPDA1 GNPDA2 UAP1 UGP2 HK1 GPI PGM1 PGM3 NAGK NAGLU CHIT1 "
                  "AMDHD2 GNE"),
    "hsa00561": ("Glycerolipid metabolism",
                  "GPAM GPAT2 AGPAT1 AGPAT2 AGPAT3 AGPAT4 AGPAT5 LPIN1 LPIN2 LPIN3 DGAT1 DGAT2 "
                  "DGKA DGKB DGKG PLPP1 PLPP2 PLPP3 LPL LIPE PNPLA2"),
    "hsa00564": ("Glycerophospholipid metabolism",
                  "GPAM AGPAT1 AGPAT2 LPIN1 CHKA CHKB PCYT1A PCYT1B CEPT1 CHPT1 LCAT LPCAT1 "
                  "LPCAT2 LPCAT3 PLA2G4A PLA2G6 PLD1 PLD2 PTDSS1 PISD"),
    "hsa00590": ("Arachidonic acid metabolism",
                  "PLA2G4A PLA2G2A PLA2G6 PTGS1 PTGS2 ALOX5 ALOX12 ALOX15 ALOX5AP LTA4H "
                  "LTB4R CYP2C8 CYP2C9 CYP4A11 EPHX2 PTGES PTGIS TBXAS1"),
    "hsa00562": ("Inositol phosphate metabolism",
                  "PIK3CA PIK3CB PIK3CD PIK3CG PTEN INPP1 INPP4A INPP4B INPP5A INPP5B INPP5D "
                  "INPP5E IPMK IPPK IMPA1 IMPA2 ITPK1"),
    "hsa00600": ("Sphingolipid metabolism",
                  "SPTLC1 SPTLC2 KDSR CERS2 CERS4 CERS5 CERS6 DEGS1 DEGS2 SMPD1 SMPD2 SMPD3 "
                  "SMPD4 ASAH1 ACER1 ACER2 SGMS1 SGMS2 UGCG GBA GLA GALC HEXA HEXB"),
    "hsa00620": ("Pyruvate metabolism",
                  "PKM PKLR PDHA1 PDHB PDHK1? PDK1 PDK2 PDK3 PDK4 DLD DLAT ACLY PC PCK1 PCK2 "
                  "LDHA LDHB MDH1 MDH2 ME1 ME2 ME3"),
    "hsa00640": ("Propanoate metabolism",
                  "ACSS2 ACSS3 ALDH2 ALDH3A2 ACACA ACACB PCC? PCCA PCCB MUT MCEE SUCLA2 SUCLG1 "
                  "LACTB2 HIBADH ABAT"),
    "hsa00760": ("Nicotinate and nicotinamide metabolism",
                  "NAMPT NMNAT1 NMNAT2 NMNAT3 NADSYN1 QPRT NAPRT NNT IDH1 IDH2 NADK NNMT PARP1 "
                  "PARP2 SIRT1 SIRT2 SIRT3"),
    "hsa00860": ("Porphyrin and chlorophyll metabolism",
                  "ALAS1 ALAS2 ALAD HMBS UROS UROD UROD? CPOX PPOX FECH UROD UGT1A1 BLVRA BLVRB"),
    "hsa00970": ("Aminoacyl-tRNA biosynthesis",
                  "AARS1 AARS2 DARS1 DARS2 EARS1? EPRS1 FARS1? FARSA FARSB GARS1 HARS1 IARS1 "
                  "KARS1 LARS1 MARS1 NARS1 QARS1 RARS1 SARS1 TARS1 VARS1 WARS1 YARS1"),
    "hsa00100": ("Steroid biosynthesis",
                  "SQLE LSS DHCR7 DHCR24 MSMO1 NSDHL EBP CYP51A1 FDFT1 HMGCR HMGCS1 MVK MVD IDI1 "
                  "FADS1 SC5D? SC5D LBR TM7SF2"),
    "hsa00140": ("Steroid hormone biosynthesis",
                  "CYP11A1 CYP17A1 CYP19A1 CYP21A2 CYP11B1 CYP11B2 HSD3B1 HSD3B2 HSD17B1 HSD17B2 "
                  "HSD17B3 HSD11B1 HSD11B2 SRD5A1 SRD5A2 STS UGT1A1"),
    "hsa00900": ("Terpenoid backbone biosynthesis",
                  "HMGCS1 HMGCR MVK MVD IDI1 IDI2 FDPS FDFT1 GGPS1 PDSS1 PDSS2 DHDDS NUS1 "
                  "COQ2 COQ3"),
    "hsa04976": ("Bile secretion",
                  "ABCB11 ABCB4 ABCG5 ABCG8 SLC10A1 SLC10A2 SLC22A1 SLC22A3 SLC51A SLC51B "
                  "SLCO1A2 SLCO1B1 SLCO1B3 NR1H4 NR0B2 HMGCR CYP7A1 CYP8B1 CYP27A1"),
    "hsa04975": ("Fat digestion and absorption",
                  "LPL LIPASE? PNLIP PNLIPRP2 CLPS CEL2? CEL FABP1 FABP2 FABP4 FABP5 CD36 SLC27A1 "
                  "SLC27A2 SLC27A4 APOA1 APOA2 APOB APOC2 APOE MTTP"),
    "hsa04978": ("Mineral absorption",
                  "SLC34A1 SLC34A2 SLC34A3 SLC8A1 TRPV5 TRPV6 SLC30A1 SLC30A2 SLC30A3 SLC39A1 "
                  "SLC39A4 SLC40A1 FTH1 FTL SLC11A1 SLC11A2 ATP7A ATP7B"),
    "hsa04979": ("Cholesterol metabolism",
                  "HMGCR HMGCS1 LDLR LRP1 APOA1 APOA2 APOB APOC1 APOC2 APOC3 APOE ABCA1 ABCG1 "
                  "NPC1 NPC2 SOAT1 SCARB1 PCSK9 CETP LCAT"),
    "hsa04973": ("Carbohydrate digestion and absorption",
                  "AMY2A AMY2B SI MGAM TREH SLC5A1 SLC2A2 SLC2A5 HK1 GCK GPI PFKM ALDOA GAPDH"),
    "hsa00051": ("Fructose and mannose metabolism",
                  "HK1 HK2 HK3 KHK ALDOA ALDOB ALDOC TPI1 PFKL PFKM PFKP FBP1 FBP2 GAPDH GPI MPI "
                  "FBP1"),
    "hsa00052": ("Galactose metabolism",
                  "GALK1 GALK2 GALT GALE UGP2 B4GALT1 B4GALT2 GLA GLB1 GALM PGM1 LCT"),
    "hsa00512": ("Mucin type O-glycan biosynthesis",
                  "GALNT1 GALNT2 GALNT3 GALNT4 GALNT6 GALNT7 GALNT10 GALNT11 GALNT12 C1GALT1 "
                  "C1GALT1C1 B3GNT6 B3GNT7 B3GALNT2"),
    "hsa00531": ("Glycosaminoglycan degradation",
                  "GUSB IDUA IDS HGSNAT NAGLU SGSH GLB1 GALNS ARSB HYAL1 HYAL2 HEXA HEXB"),
    "hsa00601": ("Glycosphingolipid biosynthesis - lacto and neolacto series",
                  "B3GNT2 B3GNT3 B3GNT4 B3GNT5 B4GALT1 B4GALT2 B4GALT3 B4GALT4 FUT2 FUT3 FUT9 "
                  "ST3GAL4 ST3GAL6 ST6GAL1"),
    "hsa00603": ("Glycosphingolipid biosynthesis - globo and isoglobo series",
                  "B3GALNT1 B3GALT1 B3GALT2 B3GALT4 A4GALT GLA GLB1 HEXA HEXB GLA ST3GAL1 "
                  "ST3GAL2"),
    "hsa00604": ("Glycosphingolipid biosynthesis - ganglio series",
                  "ST3GAL5 ST6GALNAC1 ST6GALNAC2 ST6GALNAC3 ST6GALNAC4 ST8SIA1 ST8SIA2 B4GALNT1 "
                  "B3GALT4 GM2A HEXA HEXB GLB1 GALC"),
    "hsa04922": ("Glucagon signaling pathway",
                  "GCG GCGR ADCY1 PRKACA CREB1 PCK1 PCK2 G6PC FBP1 FOXO1 PPARGC1A CPT1A ACACA "
                  "ACACB SLC2A4 GYS1 PYGL"),
    "hsa04921": ("Oxytocin signaling pathway",
                  "OXT OXTR ADCY1 PRKACA PLCB1 ITPR1 RYR1 CAMK2A CALM1 GNAS RHOA ROCK1 MYL9 "
                  "MYH9 PRKCA MAPK1"),
}

# ---- MSigDB Hallmark（常用 40 条） -----------------------------------------
_RAW_HALLMARK = {
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": ("TNFA signaling via NFKB",
        "NFKB1 NFKB2 RELA RELB IKBKB NFKBIA NFKBIB TNFAIP3 TNF BIRC3 CCL2 CXCL1 CXCL2 CXCL3 "
        "CXCL8 IL6 JUN JUNB FOS ATF3 ZFP36 ICAM1 VCAM1 SOD2 PTGS2"),
    "HALLMARK_INFLAMMATORY_RESPONSE": ("Inflammatory response",
        "IL1A IL1B IL6 CXCL8 CCL2 CCL3 CCL4 CCL5 TNF PTGS2 NFKB1 RELA TLR2 TLR4 MYD88 CXCL1 "
        "CXCL10 CXCL9 CCL19 CCL20 IL18"),
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": ("IL6 JAK STAT3 signaling",
        "IL6 IL6R IL6ST JAK1 JAK2 JAK3 TYK2 STAT1 STAT3 STAT5A STAT5B SOCS1 SOCS3 PTPN11 IL10 "
        "IL10RA IL10RB CEBPD MYC BCL2L1 MCL1 VEGFA OSMR"),
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": ("Interferon gamma response",
        "IFNG IFNGR1 IFNGR2 STAT1 STAT2 IRF1 IRF3 IRF7 JAK1 JAK2 CXCL9 CXCL10 CXCL11 CCL5 IDO1 "
        "CD274 PDCD1LG2 HLA-A HLA-B HLA-DRA B2M TAP1 ISG15 MX1 OAS1 GBP1"),
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": ("Interferon alpha response",
        "IFNA1 IFNA2 IFNB1 IFNAR1 IFNAR2 STAT1 STAT2 IRF7 IRF9 JAK1 TYK2 ISG15 MX1 MX2 OAS1 "
        "OAS2 OAS3 IFIT1 IFIT2 IFIT3 BST2 RSAD2"),
    "HALLMARK_COMPLEMENT": ("Complement",
        "C1QA C1QB C1QC C1R C1S C2 C3 C4A C4B C5 CFB CFD CFH CFI CFP CR1 CR2 CD46 CD55 CD59 "
        "SERPING1 MASP1 MASP2 VTN"),
    "HALLMARK_ALLOGRAFT_REJECTION": ("Allograft rejection",
        "HLA-A HLA-B HLA-C HLA-DRA HLA-DRB1 HLA-DPA1 HLA-DPB1 HLA-DQA1 B2M CD28 CD40LG CD80 "
        "CD86 CTLA4 ICOS IFNG IL2 IL2RA FASLG GZMB PRF1"),
    "HALLMARK_IL2_STAT5_SIGNALING": ("IL2 STAT5 signaling",
        "IL2 IL2RA IL2RB IL2RG JAK1 JAK3 STAT5A STAT5B SOCS1 SOCS3 PTPN11 FYN LCK BCL2 BCL2L1 "
        "MYC CCND1 IL4 IL4R IL13"),
    "HALLMARK_TGF_BETA_SIGNALING": ("TGF beta signaling",
        "TGFB1 TGFB2 TGFB3 TGFBR1 TGFBR2 TGFBR3 SMAD1 SMAD2 SMAD3 SMAD4 SMAD5 SMAD6 SMAD7 "
        "ACVR1 ACVR2A BMP2 BMP4 BMPR1A BMPR2 ID1 ID2 ID3 ID4 SERPINE1"),
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": ("WNT beta catenin signaling",
        "WNT2 WNT3A WNT4 WNT5A WNT6 WNT7A WNT10B FZD1 FZD2 FZD4 FZD7 LRP5 LRP6 DVL1 DVL2 DVL3 "
        "APC AXIN1 AXIN2 CTNNB1 TCF7 TCF7L2 LEF1 GSK3B MYC CCND1"),
    "HALLMARK_NOTCH_SIGNALING": ("Notch signaling",
        "NOTCH1 NOTCH2 NOTCH3 NOTCH4 DLL1 DLL3 DLL4 JAG1 JAG2 HES1 HES2 HES4 HES5 HES6 HES7 "
        "HEY1 HEY2 HEYL LFNG RFNG MFNG PSEN1 PSEN2 NCSTN APH1A RBPJ"),
    "HALLMARK_HEDGEHOG_SIGNALING": ("Hedgehog signaling",
        "SHH IHH DHH PTCH1 PTCH2 SMO GLI1 GLI2 GLI3 SUFU BOC CDON GAS1 HHIP WNT2B WNT5A BMP2 BMP4"),
    "HALLMARK_APOPTOSIS": ("Apoptosis",
        "BAX BAK1 BAD BID BCL2 BCL2L1 MCL1 BIRC2 BIRC3 XIAP APAF1 CYCS CASP3 CASP6 CASP7 CASP8 "
        "CASP9 CASP10 FADD FAS FASLG TNF TNFRSF1A TP53"),
    "HALLMARK_P53_PATHWAY": ("p53 pathway",
        "TP53 MDM2 CDKN1A BAX BBC3 PMAIP1 GADD45A GADD45B ATM ATR CHEK1 CHEK2 FAS TNFRSF10A "
        "TNFRSF10B ZMAT3 SERPINE1 PTGS2 CCNG1 RRM2B TIGAR"),
    "HALLMARK_DNA_REPAIR": ("DNA repair",
        "RAD51 BRCA1 BRCA2 PALB2 ATM ATR CHEK1 CHEK2 XRCC5 XRCC6 PRKDC LIG4 MSH2 MSH6 MLH1 "
        "PMS2 ERCC1 XPA XPC PARP1 OGG1 FEN1 PCNA"),
    "HALLMARK_G2M_CHECKPOINT": ("G2M checkpoint",
        "CDK1 CCNB1 CCNB2 CCNA2 CDC20 CDC25A CDC25B CDC25C PLK1 AURKA AURKB BUB1 BUB1B MAD2L1 "
        "NDC80 NUF2 TPX2 TOP2A MKI67 KIF11 KIF23 KIF2C CENPE CENPF ANLN ECT2 RACGAP1 BIRC5 UBE2C"),
    "HALLMARK_E2F_TARGETS": ("E2F targets",
        "E2F1 E2F2 E2F3 E2F4 E2F5 E2F6 E2F7 E2F8 RB1 RBL1 RBL2 TFDP1 MCM2 MCM3 MCM4 MCM5 MCM6 "
        "MCM7 CDC6 CDT1 PCNA TYMS DHFR CCNA2 CCNE1 CCNE2 CDK2 MKI67 TK1 SKP2"),
    "HALLMARK_MYC_TARGETS_V1": ("MYC targets V1",
        "MYC MAX NPM1 NCL RPL3 RPL5 RPS6 ODC1 LDHA NME1 NME2 PES1 DKC1 GNL3 SLC7A5 SRM TYMS "
        "PA2G4 EIF4E TFDP1 HSPD1 PTMA NPM3 RRP1"),
    "HALLMARK_MYC_TARGETS_V2": ("MYC targets V2",
        "MYC MAX MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 NCL NPM1 PES1 DKC1 GNL3 RRP1 EIF4E NME1 ODC1"),
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": ("Oxidative phosphorylation",
        "NDUFA1 NDUFA2 NDUFA3 NDUFA4 NDUFA5 NDUFB1 NDUFB2 NDUFB3 NDUFB5 NDUFB8 NDUFS1 NDUFS2 "
        "NDUFS3 NDUFS4 NDUFS7 NDUFV1 NDUFV2 SDHA SDHB SDHC SDHD UQCRC1 UQCRC2 UQCRFS1 COX4I1 "
        "COX5A COX6A1 COX6B1 COX7A2 COX7B ATP5F1A ATP5F1B ATP5F1C ATP5F1D ATP5F1E"),
    "HALLMARK_GLYCOLYSIS": ("Glycolysis",
        "HK1 HK2 GPI PFKL PFKM PFKP ALDOA ALDOB ALDOC TPI1 GAPDH PGK1 PGAM1 ENO1 ENO2 PKM LDHA "
        "LDHB SLC2A1 SLC2A3 PDK1 PFKFB3"),
    "HALLMARK_FATTY_ACID_METABOLISM": ("Fatty acid metabolism",
        "ACACA ACACB FASN SCD ELOVL1 ELOVL6 FADS1 FADS2 CPT1A CPT1B CPT2 ACADM ACADL ACADS "
        "ACADVL HADHA HADHB SLC25A20 PPARA PPARD PPARGC1A CD36 LPL FABP3 FABP4"),
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS": ("Cholesterol homeostasis",
        "HMGCR HMGCS1 MVK MVD IDI1 SQLE LSS DHCR7 DHCR24 CYP51A1 MSMO1 NSDHL EBP FDFT1 INSIG1 "
        "SREBF1 SREBF2 LDLR LRP1 APOE ABCA1 ABCG1 PCSK9"),
    "HALLMARK_ADIPOGENESIS": ("Adipogenesis",
        "PPARG CEBPA CEBPB CEBPD ADIPOQ LEP LPL FABP4 FABP5 CD36 SLC2A4 PLIN1 PLIN2 PLIN4 "
        "PPARGC1A SREBF1 FASN ACACA ACACB DGAT1 DGAT2 LIPE PNPLA2"),
    "HALLMARK_ANGIOGENESIS": ("Angiogenesis",
        "VEGFA VEGFB VEGFC KDR FLT1 FLT4 ANGPT1 ANGPT2 TEK PDGFB PDGFRB TGFB1 FGF2 HIF1A NRP1 "
        "PECAM1 CDH5 EFNB2 NOTCH1 DLL4 JAG1"),
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": ("Epithelial mesenchymal transition",
        "CDH2 VIM FN1 SNAI1 SNAI2 TWIST1 TWIST2 ZEB1 ZEB2 TGFB1 TGFBR1 COL1A1 COL1A2 COL3A1 "
        "COL5A1 LAMA1 LAMB1 LAMC1 MMP2 MMP9 TIMP1 SPARC THBS1 ITGB1 ITGA5"),
    "HALLMARK_APICAL_JUNCTION": ("Apical junction",
        "CDH1 CDH2 CTNNA1 CTNNB1 CTNND1 JUP OCLN TJP1 CLDN1 CLDN3 CLDN4 CLDN7 JAM1 GJA1 GJB1 "
        "ITGB1 ITGA5 PECAM1 VCL ACTN1 MYH9"),
    "HALLMARK_APICAL_SURFACE": ("Apical surface",
        "AFP ANPEP ATP1A1 ATP1B1 CD24 CD46 CD55 CEACAM6 CLDN1 EGFR EPCAM FXYD3 GGT1 ITGB4 "
        "MUC1 MUC13"),
    "HALLMARK_KRAS_SIGNALING_UP": ("KRAS signaling up",
        "ETV4 ETV5 MMP1 MMP2 MMP7 MMP9 CXCL8 CXCL1 CCL2 IL6 FOS JUN MYC CCND1 SPRY2 SPRY4 "
        "DUSP4 DUSP6 EGR1 EPHA2 ITGA2 ITGB1"),
    "HALLMARK_KRAS_SIGNALING_DN": ("KRAS signaling down",
        "CDKN1A CDKN2A FOXO1 FOXO3 PPARG CEBPA ADIPOQ FABP4 LPL LIPE PLIN1 UCP1 BMP4 ID2"),
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING": ("PI3K AKT MTOR signaling",
        "PIK3CA PIK3CB PIK3CD PIK3R1 PIK3R2 AKT1 AKT2 AKT3 PTEN PDPK1 MTOR RPTOR RICTOR TSC1 "
        "TSC2 RHEB RPS6KB1 RPS6KB2 EIF4EBP1 RPS6 FOXO1 GSK3B"),
    "HALLMARK_MTORC1_SIGNALING": ("MTORC1 signaling",
        "MTOR RPTOR RHEB TSC1 TSC2 RPS6KB1 RPS6KB2 EIF4EBP1 RPS6 EIF4E SLC7A5 SLC2A1 HK2 PGK1 "
        "ENO1 PKM LDHA PSAT1 SHMT2 MTHFD2"),
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": ("Unfolded protein response",
        "HSPA5 ATF4 ATF6 DDIT3 XBP1 ERN1 EIF2AK3 DNAJB9 SEL1L EDEM1 HERPUD1 PDIA3 PDIA4 PDIA6 "
        "HSP90B1 SEC61A1 SEC62 SEC63 ASNS ATF3 TRIB3"),
    "HALLMARK_PROTEIN_SECRETION": ("Protein secretion",
        "SEC61A1 SEC61B SEC61G SEC62 SEC63 SRP9 SRP14 SRP54 COPA COPB1 COPB2 ARF1 RAB1A RAB2A "
        "RAB5A RAB7A VAMP2 STX4 SNAP23 YKT6 BET1 GOSR1"),
    "HALLMARK_PEROXISOME": ("Peroxisome",
        "PEX1 PEX3 PEX5 PEX6 PEX7 PEX10 PEX12 PEX13 PEX14 PEX19 CAT ACOX1 ACOX2 HSD17B4 SCP2 "
        "SLC25A17 ABCD1 ABCD2 ABCD3"),
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY": ("Reactive oxygen species pathway",
        "SOD1 SOD2 CAT GPX1 GPX2 GPX3 GPX4 TXNRD1 TXN TXN2 NQO1 HMOX1 GCLC GCLM G6PD"),
    "HALLMARK_UV_RESPONSE_UP": ("UV response up",
        "CDKN1A GADD45A DDIT3 ATF3 JUN FOS EGR1 MDM2 TP53 BBC3 POLH XPC ERCC1 PCNA FEN1 RAD51"),
    "HALLMARK_UV_RESPONSE_DN": ("UV response down",
        "BCL2 BCL2L1 MCL1 CCND1 CCNE1 CDK1 CDK2 CDK4 CDK6 E2F1 E2F2 MYC RB1 ATM ATR"),
    "HALLMARK_HYPOXIA": ("Hypoxia",
        "HIF1A EPAS1 ARNT VEGFA VEGFB SLC2A1 SLC2A3 HK1 HK2 PGK1 ENO1 PKM LDHA PDK1 PDK3 "
        "EGLN1 EGLN3 VHL BHLHE40 BNIP3 NDRG1"),
    "HALLMARK_XENOBIOTIC_METABOLISM": ("Xenobiotic metabolism",
        "AHR CYP1A1 CYP1B1 CYP2E1 CYP3A4 CYP2C8 CYP2C9 CYP4A11 NQO1 NFE2L2 GSTA1 GSTP1 GSTT1 "
        "UGT1A1 UGT1A6 ADH1B ALDH1A1 EPHX1"),
    "HALLMARK_BILE_ACID_METABOLISM": ("Bile acid metabolism",
        "ABCB11 ABCB4 ABCG5 ABCG8 SLC10A1 SLC10A2 SLC51A SLC51B NR1H4 NR0B2 CYP7A1 CYP8B1 "
        "CYP27A1 CYP39A1 AKR1D1 AKR1C1 AKR1C4"),
    "HALLMARK_HEME_METABOLISM": ("Heme metabolism",
        "ALAS2 ALAS1 ALAD HMBS UROS UROD CPOX PPOX FECH BLVRA BLVRB HBG1 HBG2 HBA1 HBA2 HBB AHSP"),
    "HALLMARK_COAGULATION": ("Coagulation",
        "F2 F3 F5 F7 F8 F9 F10 F11 F12 F13A1 FGA FGB FGG PLG PLAT PLAU SERPINE1 VWF PF4 PPBP THBS1"),
    "HALLMARK_PANCREAS_BETA_CELLS": ("Pancreas beta cells",
        "INS GCK PDX1 NKX2-2 NEUROD1 PAX6 MAFA IAPP SLC2A2 KCNJ11 ABCC8 PCSK1 PCSK2 G6PC2 UCN3"),
    "HALLMARK_SPERMATOGENESIS": ("Spermatogenesis",
        "DAZL DDX4 TDRD1 TDRD9 PIWIL1 PIWIL2 PIWIL4 MAEL SYCP1 SYCP3 STRA8 DMC1 RAD51"),
    "HALLMARK_MYOGENESIS": ("Myogenesis",
        "MYOD1 MYOG MYF5 MYF6 PAX3 PAX7 MEF2A MEF2C MEF2D DES MYH1 MYH2 MYL1 CKM ACTN2 TNNT2 "
        "TPM1 MYH7"),
    "HALLMARK_ESTROGEN_RESPONSE_EARLY": ("Estrogen response early",
        "ESR1 PGR TFF1 TFF3 GREB1 CCND1 MYC FOS JUN SP1 PTGS2"),
    "HALLMARK_ESTROGEN_RESPONSE_LATE": ("Estrogen response late",
        "ESR1 PGR TFF1 TFF3 GREB1 CCND1 MYC JUN SP1 PTGS2 CA12"),
    "HALLMARK_ANDROGEN_RESPONSE": ("Androgen response",
        "AR KLK2 KLK3 FKBP5 TMPRSS2 NKX3-1 UGT2B15 SRD5A2"),
}

# ---- Reactome 通路（常用 25+ 条，R-HSA ID） --------------------------------
_RAW_REACTOME = {
    "R-HSA-168256": ("Immune System",
        "CD4 CD8A CD19 MS4A1 PTPRC HLA-A HLA-DRA IL2 IFNG TNF TLR4 MYD88"),
    "R-HSA-168249": ("Innate Immune System",
        "TLR1 TLR2 TLR4 TLR5 TLR9 MYD88 TICAM1 IRF3 IRF7 IFNB1 NLRP3 IL1B C3 C1QA"),
    "R-HSA-1280215": ("Cytokine Signaling in Immune system",
        "IL6 IL6R JAK1 JAK2 STAT1 STAT3 SOCS1 SOCS3 TNF IFNG IFNGR1 IL10 IL10RA CCL2 CXCL8"),
    "R-HSA-1280218": ("Adaptive Immune System",
        "CD4 CD8A CD19 MS4A1 CD40LG CD28 CTLA4 ICOS IL2 IL7R JAK3 STAT5A PTPRC"),
    "R-HSA-162582": ("Signal Transduction",
        "GNB1 GNG2 GNAI1 GNAS ADCY1 PRKACA PLCG1 EGFR GRB2 SOS1 HRAS MAP2K1 MAPK1 PIK3CA AKT1"),
    "R-HSA-9006934": ("Signaling by Receptor Tyrosine Kinases",
        "EGFR ERBB2 ERBB3 ERBB4 KDR FLT1 MET FGFR1 IGF1R INSR PDGFRA PDGFRB NTRK1 GRB2 SOS1"),
    "R-HSA-194315": ("Signaling by Rho GTPases",
        "RHOA RAC1 CDC42 ROCK1 ROCK2 PAK1 PAK2 LIMK1 CFL1 ARPC1B WASL IQGAP1"),
    "R-HSA-422475": ("Axon guidance",
        "NTN1 DCC SLIT2 ROBO1 SEMA3A NRP1 PLXNA1 EPHA4 EFNB2 RHOA RAC1 CDC42 LIMK1 CFL1 GSK3B"),
    "R-HSA-1266738": ("Developmental Biology",
        "BMP4 BMP7 SHH WNT3A FGF8 GLI1 CTNNB1 SMAD1 SMAD4 NOTCH1 DLL1 HES1 PAX6 SOX2"),
    "R-HSA-8953897": ("Cellular responses to stress",
        "HSPA1A HSPA5 DDIT3 ATF4 ATF3 TP53 CDKN1A NFE2L2 HMOX1 SOD1 CAT HIF1A GADD45A"),
    "R-HSA-392499": ("Metabolism of proteins",
        "RPL3 RPL5 RPS3 RPS6 EEF1A1 EEF2 EIF4E HSPA8 PSMA1 PSMB1 UBA1 UBE2C"),
    "R-HSA-1430728": ("Metabolism",
        "HK1 GPI PFKM ALDOA GAPDH PKM LDHA CS IDH2 MDH2 ACACA FASN HMGCR"),
    "R-HSA-74160": ("Gene expression (Transcription)",
        "POLR2A TBP TAF1 CDK7 CCNH MYC JUN FOS STAT1 STAT3 NFKB1 RELA TP53"),
    "R-HSA-73857": ("RNA Polymerase II Transcription",
        "POLR2A TBP TAF1 GTF2B GTF2F1 CDK7 CCNH MED1 MED12 MYC JUN FOS"),
    "R-HSA-212436": ("Generic Transcription Pathway",
        "TBP TAF1 TAF2 GTF2B GTF2F1 MED1 MED12 NCOA1 NCOA2 NCOR1 SP1 YY1"),
    "R-HSA-72172": ("mRNA Splicing",
        "SNRNP70 SRSF1 SRSF2 SRSF3 HNRNPA1 PTBP1 U2AF1 U2AF2 PRPF8 SF3B1"),
    "R-HSA-72163": ("mRNA Splicing - Major Pathway",
        "SNRNP70 SF3B1 SF3A1 U2AF1 U2AF2 PRPF8 PRPF19 SNRPD1 SNRPD2"),
    "R-HSA-5653656": ("Vesicle-mediated transport",
        "RAB5A RAB7A RAB11A RAB27A STX4 VAMP2 VAMP3 SNAP23 NSF CLTC AP2A1 COPB1"),
    "R-HSA-199991": ("Membrane Trafficking",
        "RAB5A RAB7A RAB11A RAB27A RAB1A RAB2A STX4 VAMP2 VAMP3 SNAP23 CLTC AP2A1 ARF1 COPB1"),
    "R-HSA-597592": ("Post-translational protein modification",
        "UBA1 UBE2C MDM2 SKP2 FBXW7 CBL ITCH NEDD4 PSMA1 PSMB1 SUMO1 UBC"),
    "R-HSA-109582": ("Hemostasis",
        "F2 F3 F5 F7 F8 F9 F10 F11 F12 F13A1 FGA FGB FGG PLG VWF PF4 PPBP"),
    "R-HSA-1643685": ("Disease",
        "TP53 EGFR ERBB2 KRAS BRAF PIK3CA PTEN RB1 CDKN2A BRCA1 BRCA2 APC CTNNB1 MYC"),
    "R-HSA-9711097": ("Cellular response to stimuli",
        "FOS JUN EGR1 ATF3 HSPA1A NFE2L2 HMOX1 SOD1 CAT TP53 NFKB1"),
    "R-HSA-2262752": ("Cellular responses to stimuli",
        "FOS JUN EGR1 ATF3 HSPA1A HSPA5 DDIT3 NFE2L2 HMOX1 TP53 CDKN1A GADD45A"),
    "R-HSA-9013149": ("RHO GTPase Effectors",
        "RHOA RAC1 CDC42 ROCK1 ROCK2 PAK1 PAK2 LIMK1 CFL1 WASL ARPC1B IQGAP1"),
    "R-HSA-69306": ("DNA Replication",
        "MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 PCNA RFC1 RPA1 POLA1 POLD1 POLE CDC6 CDT1"),
    "R-HSA-73894": ("DNA Repair",
        "RAD51 BRCA1 BRCA2 ATM ATR CHEK1 CHEK2 XRCC5 XRCC6 PRKDC LIG4 MSH2 MSH6 MLH1 ERCC1 XPA "
        "XPC PARP1 OGG1"),
    "R-HSA-69481": ("G2/M Checkpoints",
        "CDK1 CCNB1 CCNB2 CDC25A CDC25B CDC25C PLK1 AURKA AURKB BUB1 MAD2L1 CHEK1 CHEK2 WEE1"),
    "R-HSA-69278": ("Cell Cycle, Mitotic",
        "CDK1 CCNB1 CCNA2 CDC20 PLK1 AURKA AURKB BUB1 MAD2L1 NDC80 NUF2 TPX2 TOP2A KIF11"),
    "R-HSA-1640170": ("Cell Cycle",
        "CDK1 CDK2 CDK4 CDK6 CCNA2 CCNB1 CCND1 CCNE1 E2F1 RB1 CDKN1A CDKN2A CDC25A PLK1"),
    "R-HSA-453279": ("Mitotic G1-G1/S phases",
        "CDK2 CCND1 CCND2 CDK4 CDK6 E2F1 RB1 CDKN1A CDKN1B CDKN2A MCM2 MCM3 MCM4 CDC6 CDT1"),
    "R-HSA-1474244": ("Extracellular matrix organization",
        "COL1A1 COL1A2 COL3A1 COL4A1 COL5A1 COL6A1 FN1 FBN1 LAMA1 LAMB1 LAMC1 MMP2 MMP9 TIMP1 "
        "SPARC THBS1 VTN"),
    "R-HSA-1474290": ("Collagen formation",
        "COL1A1 COL1A2 COL3A1 COL4A1 COL5A1 COL6A1 COL7A1 COL17A1 COL18A1 P4HA1 P4HB PLOD1 PLOD2"),
    "R-HSA-5357801": ("Programmed Cell Death",
        "BAX BAK1 BCL2 BCL2L1 MCL1 BID CASP3 CASP8 CASP9 APAF1 CYCS FADD FAS TNF"),
    "R-HSA-109606": ("Intrinsic Pathway for Apoptosis",
        "BAX BAK1 BAD BID BCL2 BCL2L1 MCL1 APAF1 CYCS CASP9 CASP3"),
    "R-HSA-169911": ("Regulation of Apoptosis",
        "BCL2 BCL2L1 MCL1 BIRC2 BIRC3 XIAP BAD BID BAX BAK1"),
    "R-HSA-70171": ("Glycolysis",
        "HK1 HK2 GPI PFKL PFKM PFKP ALDOA TPI1 GAPDH PGK1 PGAM1 ENO1 PKM LDHA"),
    "R-HSA-71403": ("Citric acid cycle (TCA cycle)",
        "CS ACO1 ACO2 IDH1 IDH2 IDH3A OGDH SUCLA2 SUCLG1 SDHA SDHB SDHC SDHD FH MDH1 MDH2"),
    "R-HSA-611105": ("Respiratory electron transport",
        "NDUFA1 NDUFB1 SDHA SDHB SDHC SDHD UQCRC1 UQCRC2 UQCRFS1 COX4I1 COX5A ATP5F1A ATP5F1B"),
    "R-HSA-8978868": ("Fatty acid metabolism",
        "ACACA ACACB FASN SCD CPT1A CPT1B CPT2 ACADM ACADL HADHA HADHB"),
    "R-HSA-556833": ("Metabolism of lipids",
        "PLIN2 DGAT1 DGAT2 LIPE PNPLA2 AGPAT2 LPIN1 SREBF1 HMGCR LDLR ABCA1 APOE"),
    "R-HSA-8957322": ("Metabolism of steroids",
        "HMGCR HMGCS1 MVK MVD IDI1 SQLE LSS DHCR7 CYP51A1 MSMO1"),
    "R-HSA-191273": ("Cholesterol biosynthesis",
        "HMGCS1 HMGCR MVK MVD IDI1 SQLE LSS DHCR7 DHCR24 CYP51A1 MSMO1 NSDHL FDFT1"),
    "R-HSA-211859": ("Biological oxidations",
        "CYP1A1 CYP1B1 CYP2E1 CYP3A4 NQO1 HMOX1 GSTA1 GSTP1 UGT1A1 EPHX1"),
    "R-HSA-156580": ("Phase II - Conjugation of compounds",
        "GSTA1 GSTP1 GSTT1 GSTM1 UGT1A1 UGT1A6 SULT1A1 NAT1 NAT2"),
    "R-HSA-6798695": ("Neutrophil degranulation",
        "ELANE MPO PRTN3 CTSG LTF CAMP DEFA3 DEFA4 AZU1 BPI CEACAM8 MMP9 S100A8 S100A9"),
    "R-HSA-2029480": ("Fcgamma receptor (FCGR) dependent phagocytosis",
        "FCGR1A FCGR2A FCGR3A FCER1G SYK LYN HCK PLCG1 PIK3CA AKT1 RAC1 CDC42 ARPC1B"),
    "R-HSA-166658": ("Complement cascade",
        "C1QA C1QB C1QC C1R C1S C2 C3 C4A C4B C5 CFB CFD CFH CFI CFP CD46 CD55 CD59"),
    "R-HSA-977606": ("Regulation of Complement cascade",
        "CD46 CD55 CD59 CFH CFI C4BPA C4BPB CR1 SERPING1"),
    "R-HSA-173623": ("Classical antibody-mediated complement activation",
        "C1QA C1QB C1QC C1R C1S C2 C3 C4A C4B C5"),
    "R-HSA-6805567": ("Keratinization",
        "KRT1 KRT5 KRT10 KRT14 KRT6A KRT6B KRT16 KRT17 FLG LOR IVL"),
    "R-HSA-6809371": ("Formation of the cornified envelope",
        "FLG LOR IVL KRT1 KRT10 TGM1 TGM3 SPRR1A SPRR2A"),
    "R-HSA-9006931": ("Signaling by Nuclear Receptors",
        "ESR1 ESR2 AR PGR NR3C1 THRA THRB RARA RARB RXRA VDR PPARG"),
}

# ---- Charoentong 28 类免疫细胞签名（IMMUNE_28，供 bulk-clin 复用） ----------
# 说明：Charoentong et al. 2017 的 28 类面板；为保持 28 个条目，略去与 CD56bright/
# CD56dim 重复的泛 NK 及与 iDC/pDC 重复的 aDC，保留全部 CD4/CD8 亚群与 Th1/Th2/Th17/Tfh/Treg。
_RAW_IMMUNE_28 = {
    "Activated B cell":
        "CD19 MS4A1 TNFRSF17 CD79A CD79B POU2AF1 SPIB BLNK COCH CR2",
    "Activated CD4 T cell":
        "CD4 IL2RA TNFRSF4 TNFRSF18 CTLA4 ICOS CD28 IL2 IFNG",
    "Activated CD8 T cell":
        "CD8A CD8B TNFRSF9 GZMB PRF1 IFNG NKG7 CD27 CD69",
    "CD56bright natural killer cell":
        "NCAM1 KLRF1 GZMK GZMA PRF1 NKG7 KLRK1 CCL5 XCL1",
    "CD56dim natural killer cell":
        "FCGR3A NKG7 GNLY GZMB PRF1 KLRD1 NCR1 SPON2",
    "Central memory CD4 T cell":
        "CD4 IL7R CCR7 SELL TCF7 CD27 LEF1 CD3D CD28",
    "Central memory CD8 T cell":
        "CD8A IL7R CCR7 SELL TCF7 CD27 LEF1 CD3D CD28",
    "Effector memory CD4 T cell":
        "CD4 IL7R CCR5 GZMA IFNG KLRG1 GZMK CX3CR1",
    "Effector memory CD8 T cell":
        "CD8A GZMB PRF1 IFNG CX3CR1 KLRG1 GZMA ZEB2",
    "Eosinophil":
        "PRG2 EPX RNASE3 CLC IL5RA CCR3 SIGLEC8 RNASE2",
    "Gamma delta T cell":
        "TRDC TRGC1 TRGC2 CD3D CD3E GZMB IFNG TRAT1",
    "Immature B cell":
        "CD19 MS4A1 CD79A CD79B IGLL1 VPREB1 VPREB3 CD22",
    "Immature dendritic cell":
        "CD1C CD1E CLEC10A FCER1A CD83 CD86 HLA-DRA CD1A",
    "Macrophage":
        "CD68 CD163 MRC1 MSR1 MARCO CSF1R C1QA C1QB",
    "Mast cell":
        "TPSAB1 TPSB2 CPA3 MS4A2 KIT HDC CMA1 SIGLEC6",
    "Memory B cell":
        "CD27 MS4A1 CD19 CD40 CR2 TNFRSF17 CD79A CD79B",
    "Monocyte":
        "CD14 LYZ FCN1 VCAN S100A8 S100A9 CSF1R FCGR3A",
    "Natural killer T cell":
        "CD3D CD3E CD3G NKG7 GNLY KLRD1 PRF1 KLRK1",
    "Neutrophil":
        "FCGR3B CXCR2 CSF3R S100A8 S100A9 ELANE MPO CXCR1",
    "Naive B cell":
        "MS4A1 CD19 TCL1A IGHM IGHD FCER2 CD22 CD79A",
    "Naive CD4 T cell":
        "CD4 CCR7 SELL IL7R LEF1 TCF7 CD3D CD28",
    "Naive CD8 T cell":
        "CD8A CCR7 SELL IL7R LEF1 TCF7 CD8B CD3D",
    "Plasmacytoid dendritic cell":
        "IL3RA CLEC4C LILRA4 IRF7 GZMB TCF4 IRF8 GZMA",
    "Regulatory T cell":
        "FOXP3 IL2RA CTLA4 TIGIT IKZF2 TNFRSF18 IL10 CCR4",
    "T follicular helper cell":
        "CXCR5 ICOS PDCD1 BCL6 IL21 MAF CD40LG IL4",
    "Type 1 T helper cell":
        "IFNG TBX21 CXCR3 IL12RB2 STAT4 IL2 TNF CCR5",
    "Type 2 T helper cell":
        "GATA3 IL4 IL5 IL13 CCR4 STAT6 IL4R PTGDR2",
    "Type 17 T helper cell":
        "IL17A IL17F RORC CCR6 IL23R STAT3 IL1R1 CCL20",
}

# ===========================================================================
# 2. 规范化构造（去重、大写、剔除占位符）
# ===========================================================================

_GENE_RE = re.compile(r"^[A-Z0-9][A-Z0-9._-]*$")


def _clean_genes(text) -> list[str]:
    """把空格分隔的基因串清成大写、去重、保序的 symbol 列表。"""
    out: list[str] = []
    seen: set[str] = set()
    for token in str(text).replace(",", " ").replace(";", " ").split():
        raw = token.strip()
        if not raw or raw.endswith("?"):
            continue
        gene = raw.upper()
        if not _GENE_RE.match(gene):
            continue
        if gene in seen:
            continue
        seen.add(gene)
        out.append(gene)
    return out


def _build(raw: dict, min_genes: int = MIN_GENES_PER_SET) -> dict:
    """把 {id: (name, "genes")} 规范化为 {id: {"name","genes"}}（丢弃过小条目）。"""
    sets: dict[str, dict] = {}
    for key, value in raw.items():
        name, genes = value
        cleaned = _clean_genes(genes)
        if len(cleaned) < int(min_genes):
            continue
        sets[str(key)] = {"name": str(name), "genes": cleaned}
    return sets


def _build_immune(raw: dict, min_genes: int = MIN_GENES_PER_SET) -> dict:
    """免疫签名：{cell_type: [symbol, ...]}（纯 dict，供 bulk-clin 直接迭代）。"""
    out: dict[str, list[str]] = {}
    for name, genes in raw.items():
        cleaned = _clean_genes(genes)
        if len(cleaned) < int(min_genes):
            continue
        out[str(name)] = cleaned
    return out


GO_BP = _build(_RAW_GO_BP)
GO_CC = _build(_RAW_GO_CC)
GO_MF = _build(_RAW_GO_MF)
KEGG = _build(_RAW_KEGG)
HALLMARK = _build(_RAW_HALLMARK)
REACTOME = _build(_RAW_REACTOME)
IMMUNE_28 = _build_immune(_RAW_IMMUNE_28)

BUILTIN = {
    "GO_BP": GO_BP,
    "GO_CC": GO_CC,
    "GO_MF": GO_MF,
    "KEGG": KEGG,
    "HALLMARK": HALLMARK,
    "REACTOME": REACTOME,
    "IMMUNE_28": {name: {"name": name, "genes": list(genes)} for name, genes in IMMUNE_28.items()},
}

BUILTIN_LABELS = {
    "GO_BP": "GO 生物过程",
    "GO_CC": "GO 细胞组分",
    "GO_MF": "GO 分子功能",
    "KEGG": "KEGG 通路",
    "HALLMARK": "MSigDB Hallmark",
    "REACTOME": "Reactome 通路",
    "IMMUNE_28": "Charoentong 免疫细胞签名",
}

# 常用别名（大小写不敏感）
DB_ALIASES = {
    "GO": "GO_BP",
    "BP": "GO_BP",
    "GOBP": "GO_BP",
    "GO:BP": "GO_BP",
    "GO_BIOLOGICAL_PROCESS": "GO_BP",
    "GO_BP": "GO_BP",
    "CC": "GO_CC",
    "GOCC": "GO_CC",
    "GO:CC": "GO_CC",
    "GO_CELLULAR_COMPONENT": "GO_CC",
    "GO_CC": "GO_CC",
    "MF": "GO_MF",
    "GOMF": "GO_MF",
    "GO:MF": "GO_MF",
    "GO_MOLECULAR_FUNCTION": "GO_MF",
    "GO_MF": "GO_MF",
    "KEGG": "KEGG",
    "KEGG_PATHWAY": "KEGG",
    "HALLMARK": "HALLMARK",
    "MSIGDB_HALLMARK": "HALLMARK",
    "MSIGDB": "HALLMARK",
    "H": "HALLMARK",
    "REACTOME": "REACTOME",
    "REACTOME_PATHWAY": "REACTOME",
    "IMMUNE": "IMMUNE_28",
    "IMMUNE_28": "IMMUNE_28",
    "IMMUNE28": "IMMUNE_28",
}

# 在线库名 → 内置回退（网络失败时用）
_FALLBACK_RULES = [
    ("HALLMARK", "HALLMARK"),
    ("MSIGDB", "HALLMARK"),
    ("REACTOME", "REACTOME"),
    ("KEGG", "KEGG"),
    ("GO_BIOLOGICAL_PROCESS", "GO_BP"),
    ("GOBP", "GO_BP"),
    ("GO_CELLULAR_COMPONENT", "GO_CC"),
    ("GOCC", "GO_CC"),
    ("GO_MOLECULAR_FUNCTION", "GO_MF"),
    ("GOMF", "GO_MF"),
]


def resolve_db(db) -> str | None:
    """把 db 名解析成内置 key（未知返回 None）。"""
    text = str(db or "").strip()
    if not text:
        return None
    upper = text.upper().replace(" ", "_")
    if upper in DB_ALIASES:
        return DB_ALIASES[upper]
    if upper in BUILTIN:
        return upper
    return None


def _fallback_db(db) -> str | None:
    upper = str(db or "").strip().upper().replace(" ", "_")
    for needle, key in _FALLBACK_RULES:
        if needle in upper:
            return key
    return None


def available_dbs() -> list[str]:
    """内置可选库列表（用于 UI 下拉/校验）。"""
    return list(BUILTIN.keys())


def builtin_stats() -> dict:
    """内置库条目数 / 每库最小基因数（测试与自检用）。"""
    out = {}
    for key, sets in BUILTIN.items():
        sizes = [len(item["genes"]) for item in sets.values()]
        out[key] = {
            "n_sets": len(sets),
            "min_genes": int(min(sizes)) if sizes else 0,
            "max_genes": int(max(sizes)) if sizes else 0,
        }
    return out


def _copy_sets(sets: dict) -> dict:
    """浅拷贝嵌套结构，避免调用方修改内置库。"""
    return {
        str(key): {"name": str(val.get("name") or key), "genes": list(val.get("genes") or [])}
        for key, val in (sets or {}).items()
    }


# ===========================================================================
# 3. GMT 读写
# ===========================================================================

def _parse_gmt_text(text: str) -> dict:
    sets: dict[str, dict] = {}
    for line in str(text or "").splitlines():
        line = line.rstrip("\r\n")
        if not line.strip() or line.startswith("#"):
            continue
        parts = [part.strip() for part in line.split("\t")]
        parts = [part for part in parts if part]
        if len(parts) < 3:
            continue
        term = parts[0]
        desc = parts[1]
        genes = _clean_genes(" ".join(parts[2:]))
        if len(genes) < 3:
            continue
        sets[term] = {"name": desc or term, "genes": genes}
    return sets


def load_gmt(path) -> dict:
    """读取标准 GMT → {term: {"name","genes"}}。文件不存在/不可解析抛 EngineError。"""
    file_path = Path(str(path or "")).expanduser()
    if not str(path or "").strip():
        raise EngineError("bad_request", "需要提供 GMT 文件路径")
    if not file_path.exists() or not file_path.is_file():
        raise EngineError("bad_request", "GMT 文件不存在或不是文件：%s" % path)
    try:
        text = file_path.read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        raise EngineError("internal", "GMT 文件读取失败：%s" % error)
    sets = _parse_gmt_text(text)
    if not sets:
        raise EngineError("bad_request", "GMT 文件里没有可用基因集（需 ≥3 个基因/行）：%s" % path)
    return sets


def to_gmt(sets) -> str:
    """把 {term: {"name","genes"}} 或 {term: [genes]} 写成 GMT 文本。"""
    lines = []
    for term, value in (sets or {}).items():
        if isinstance(value, dict):
            name = str(value.get("name") or term)
            genes = list(value.get("genes") or [])
        elif isinstance(value, (list, tuple, set)):
            name = str(term)
            genes = list(value)
        else:
            continue
        cleaned = _clean_genes(" ".join(str(g) for g in genes))
        if not cleaned:
            continue
        lines.append("\t".join([str(term), name] + cleaned))
    return "\n".join(lines) + ("\n" if lines else "")


# ===========================================================================
# 4. Enrichr 在线库（GMT 文本接口，20s 超时，缓存到 <workdir>/genesets/）
# ===========================================================================

def _safe_name(name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name or "")).strip("._-")
    return (clean[:120] or "library")


def genesets_cache_dir(workdir=None) -> Path:
    """在线库缓存目录：<workdir>/genesets/（workdir 为空时用系统临时目录）。"""
    base = str(workdir or "").strip()
    root = Path(base).expanduser().resolve() if base else Path(tempfile.gettempdir()) / "dsh-geobulk-cache"
    path = root / _CACHE_DIRNAME
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError:
        path = Path(tempfile.gettempdir()) / "dsh-geobulk-cache" / _CACHE_DIRNAME
        path.mkdir(parents=True, exist_ok=True)
    return path


def enrichr_cache_path(name: str, workdir=None) -> str:
    return str(genesets_cache_dir(workdir) / ("%s.gmt" % _safe_name(name)))


# macOS 系统 CA（Python.framework 的 ssl 默认路径可能不含它 → 证书校验失败时用它重试）
_SYSTEM_CA_CANDIDATES = [
    "/etc/ssl/cert.pem",
    "/usr/local/etc/openssl/cert.pem",
    "/opt/homebrew/etc/openssl/cert.pem",
]


def _ssl_contexts() -> list:
    """返回可尝试的 SSL 上下文：默认 → 系统 CA 包（纯标准库，不依赖 certifi）。"""
    import ssl

    contexts = [ssl.create_default_context()]
    for cafile in _SYSTEM_CA_CANDIDATES:
        if not os.path.exists(cafile):
            continue
        try:
            contexts.append(ssl.create_default_context(cafile=cafile))
        except Exception:  # noqa: BLE001 - CA 文件异常时跳过
            continue
    return contexts


def _http_get(url: str, timeout: float) -> str:
    request = urllib.request.Request(
        url, headers={"User-Agent": "dsh-geobulk/1.0 (+bulk-transcriptome-workbench)"})
    contexts = _ssl_contexts()
    last_error: Exception | None = None
    for index, context in enumerate(contexts):
        try:
            with urllib.request.urlopen(request, timeout=timeout, context=context) as response:
                return response.read().decode("utf-8", "replace")
        except Exception as error:  # noqa: BLE001 - SSL/HTTP/超时统一重试
            last_error = error
            if index + 1 < len(contexts):
                continue
            raise last_error  # type: ignore[misc]
    raise last_error or RuntimeError("HTTP GET failed")


def fetch_enrichr_library(name: str, workdir=None, timeout: int = 20) -> dict:
    """从 Enrichr 拉取 GMT 基因集库并缓存；失败抛 EngineError("unsupported")。

    返回 {term: {"name","genes"}}；命中缓存不访问网络。
    """
    library = str(name or "").strip()
    if not library:
        raise EngineError("bad_request", "在线基因集库名称不能为空")
    cache = Path(enrichr_cache_path(library, workdir))
    if cache.exists() and cache.stat().st_size > 0:
        try:
            cached = _parse_gmt_text(cache.read_text(encoding="utf-8", errors="replace"))
            if cached:
                return cached
        except OSError:
            pass
    url = ENRICHR_GMT_URL % urllib.parse.quote(library, safe="")
    try:
        text = _http_get(url, float(timeout))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as error:
        raise EngineError(
            "unsupported",
            "在线基因集库「%s」拉取失败（%s）；已使用内置库回退" % (library, error),
        )
    sets = _parse_gmt_text(text)
    if not sets:
        raise EngineError(
            "unsupported",
            "在线基因集库「%s」返回内容无法解析（可能库名不存在）" % library,
        )
    try:
        cache.write_text(to_gmt(sets), encoding="utf-8")
    except OSError:
        pass
    return sets


# ===========================================================================
# 5. 统一入口
# ===========================================================================

def get_gene_sets_with_source(db, workdir=None) -> tuple[dict, str, list[str]]:
    """取基因集，并返回 (sets, source, warnings)。

    * 内置库（GO_BP/GO_CC/GO_MF/KEGG/HALLMARK/REACTOME/IMMUNE_28）→ 直接返回；
    * 以 .gmt 结尾的路径 → load_gmt；
    * 其它名称 → 当作 Enrichr 库名在线拉取；失败则按名称回退到对应内置库并给 warning，
      无对应内置库时抛 EngineError("unsupported")。
    """
    warnings: list[str] = []
    text = str(db or "").strip()
    if not text:
        raise EngineError("bad_request", "需要指定基因集库（%s）" % "、".join(available_dbs()))
    # 1) 内置库 / 别名
    key = resolve_db(text)
    if key:
        return _copy_sets(BUILTIN[key]), "builtin", warnings
    # 2) 本地 GMT 文件
    if text.lower().endswith(".gmt"):
        return load_gmt(text), "gmt", warnings
    # 3) 在线 Enrichr
    try:
        sets = fetch_enrichr_library(text, workdir=workdir)
        return sets, "enrichr", warnings
    except EngineError as error:
        fallback = _fallback_db(text)
        if fallback:
            warnings.append("%s；已回退到内置库 %s" % (error.message, fallback))
            return _copy_sets(BUILTIN[fallback]), "builtin-fallback", warnings
        raise EngineError(
            "unsupported",
            "%s；可选内置库：%s（或提供本地 .gmt 文件）" % (error.message, "、".join(available_dbs())),
        )


def get_gene_sets(db, workdir=None) -> dict:
    """取 {term: {"name","genes"}}（内置库 / .gmt 文件 / Enrichr 在线库）。"""
    sets, _source, _warnings = get_gene_sets_with_source(db, workdir=workdir)
    return sets


__all__ = [
    "GO_BP", "GO_CC", "GO_MF", "KEGG", "HALLMARK", "REACTOME", "IMMUNE_28",
    "BUILTIN", "BUILTIN_LABELS", "DB_ALIASES",
    "available_dbs", "builtin_stats", "resolve_db",
    "load_gmt", "to_gmt", "get_gene_sets", "get_gene_sets_with_source",
    "fetch_enrichr_library", "enrichr_cache_path", "genesets_cache_dir",
    "EngineError", "MIN_GENES_PER_SET",
]

