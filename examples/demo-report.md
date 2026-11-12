# Bulk 转录组分析报告：demo-bulk

## 数据概览

- 样本数：**36**；基因数：**4000**
- 数据来源：`内置示例`（格式 demo）
- 值类型：lognorm；物种：human
- 分组列：`group`；分组：case, control
- 已完成步骤：qc, filter, normalize, batch, pca, de, enrich, immune, survival, prognosis, scmap
- 引擎：dsh-geobulk 工作台（Python 3.14.2, R 4.4.2)

## Methods


Bulk transcriptome analysis was performed with the dsh-geobulk 工作台（Python 3.14.2, R 4.4.2).
The dataset contained 36 samples and 4000 genes.

## Data processing
- **qc**: Sample-level quality control assessed library size, detected genes and outlier samples. Parameters: `{"group_col": "group"}`.
- **filter**: Low-expression genes were removed. Parameters: `{"group_col": "group"}`.
- **normalize**: Expression values were log-transformed and quantile-normalized. Parameters: `{"group_col": "group"}`.
- **batch**: Batch effects were corrected when multiple batches were present. Parameters: `{"group_col": "group"}`.
- **pca**: Principal component analysis was used for sample-level dimensionality reduction. Parameters: `{"group_col": "group"}`.
- **de**: Differential expression analysis was performed. Parameters: `{"group_col": "group", "case": "case", "control": "control", "method": "limma"}`.
- **enrich**: Functional enrichment analysis was performed on the differential gene list. Parameters: `{"db": ["GO_BP", "KEGG", "HALLMARK"], "method": "ora", "group_col": "group"}`.
- **immune**: Immune cell infiltration was estimated with deconvolution and signature-based methods. Parameters: `{"methods": ["cibersort", "estimate", "mcpcounter", "ssgsea", "epic"], "group_col": "group"}`.
- **survival**: Kaplan-Meier curves and Cox proportional hazards models were used for survival analysis. Parameters: `{"time_col": "os_time", "event_col": "os_event", "group_col": "group"}`.
- **prognosis**: A risk model was built and evaluated with time-dependent ROC and C-index. Parameters: `{"model": "cox", "group_col": "group"}`.
- **scmap**: Bulk signatures were projected onto a single-cell reference dataset to identify expressing cell subsets. Parameters: `{"ref_path": "test/fixtures/demo_sc.h5ad", "gene_set": "de_up", "group_col": "group"}`.

## Statistics
Differential expression used moderated t-statistics with Benjamini-Hochberg FDR control. Enrichment used hypergeometric tests (ORA) or weighted Kolmogorov-Smirnov statistics (GSEA) against public gene sets. Survival comparisons used the log-rank test and Cox proportional hazards models. All tests were two-sided and P values below 0.05 were considered significant unless otherwise stated.

_Generated at 2026-11-12 21:02:55._

## 表达矩阵（前 10 基因 × 前 10 样本）

| gene | C01 | C02 | C03 | C04 | C05 | C06 | C07 | C08 | C09 | C10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8A | 7.97 | 8.49 | 8.50 | 8.97 | 9.59 | 7.88 | 9.35 | 9.14 | 8.80 | 9.10 |
| CD8B | 4.11 | 4.51 | 5.00 | 6.59 | 5.50 | 5.18 | 6.15 | 5.46 | 4.63 | 4.47 |
| GZMA | 5.32 | 6.32 | 7.26 | 7.69 | 6.60 | 5.96 | 6.73 | 6.99 | 6.96 | 7.34 |
| GZMB | 4.15 | 5.01 | 5.13 | 5.67 | 5.49 | 5.34 | 5.45 | 4.26 | 5.39 | 4.34 |
| GZMK | 2.90 | 1.58 | 3.31 | 4.31 | 4.07 | 3.24 | 4.40 | 5.23 | 3.91 | 3.86 |
| GZMH | 9.69 | 8.80 | 8.88 | 9.12 | 8.75 | 9.38 | 8.79 | 9.15 | 8.33 | 9.31 |
| PRF1 | 5.67 | 6.57 | 6.27 | 7.79 | 6.94 | 6.59 | 7.02 | 6.79 | 6.36 | 6.27 |
| NKG7 | 6.38 | 5.73 | 6.60 | 7.06 | 6.42 | 6.00 | 7.78 | 7.19 | 5.52 | 6.75 |
| GNLY | 8.47 | 7.98 | 8.67 | 9.35 | 8.25 | 8.27 | 9.25 | 8.28 | 7.70 | 8.06 |
| KLRD1 | 7.67 | 8.35 | 8.32 | 9.64 | 8.99 | 8.93 | 8.85 | 9.37 | 8.49 | 9.29 |

## 差异表达（前 30）

| gene | logFC | AveExpr | P.Value | adj.P.Val | regulation |
| --- | --- | --- | --- | --- | --- |
| IGHM | -2.401 | 8.414 | 8.34e-24 | 3.33e-20 | down |
| CKAP4 | 2.034 | 8.118 | 6.46e-23 | 1.29e-19 | up |
| E2F1 | 2.142 | 8.823 | 1.19e-22 | 1.42e-19 | up |
| NT5C2 | 1.985 | 9.258 | 1.49e-22 | 1.42e-19 | up |
| LOC126860867 | -2.014 | 7.236 | 1.77e-22 | 1.42e-19 | down |
| LOC127826254 | 2.029 | 6.997 | 2.65e-22 | 1.77e-19 | up |
| PILRA | 2.069 | 7.873 | 3.84e-22 | 2.11e-19 | up |
| HLA-DOA | 2.12 | 8.184 | 4.23e-22 | 2.11e-19 | up |
| KRT19 | -2.015 | 7.547 | 6.31e-22 | 2.81e-19 | down |
| S100A7A | 2.065 | 8.154 | 8.85e-22 | 3.54e-19 | up |
| BUB1 | 2.023 | 8.659 | 1.11e-21 | 4.04e-19 | up |
| IGLL5 | 2.084 | 7.145 | 1.22e-21 | 4.06e-19 | up |
| IL7R | -2.386 | 8.191 | 1.48e-21 | 4.57e-19 | down |
| HYKK | 2.033 | 9.075 | 1.89e-21 | 5.39e-19 | up |
| KLRG1 | 2.063 | 7.63 | 2.10e-21 | 5.59e-19 | up |
| COL1A1 | -2.216 | 7.748 | 2.95e-21 | 6.72e-19 | down |
| CRP | 2.166 | 9.367 | 3.30e-21 | 6.72e-19 | up |
| CPA3 | -2.463 | 7.019 | 3.32e-21 | 6.72e-19 | down |
| LIG3 | 2.195 | 6.646 | 3.36e-21 | 6.72e-19 | up |
| MAPT | 2.125 | 8.016 | 3.36e-21 | 6.72e-19 | up |
| LOC127829456 | -1.825 | 9.481 | 5.48e-21 | 1.04e-18 | down |
| TYMS | 1.9 | 7.437 | 8.02e-21 | 1.46e-18 | up |
| CCNE1 | 1.941 | 7.327 | 8.54e-21 | 1.48e-18 | up |
| CD300LF | 2.175 | 7.259 | 1.04e-20 | 1.73e-18 | up |
| JCHAIN | -1.911 | 8.042 | 1.22e-20 | 1.95e-18 | down |
| MCM4 | 2.173 | 9.194 | 1.30e-20 | 1.96e-18 | up |
| GNLY | -2.121 | 9.277 | 1.32e-20 | 1.96e-18 | down |
| CCT4 | -1.893 | 10.11 | 1.52e-20 | 2.17e-18 | down |
| IGHV4-55 | 1.879 | 9.71 | 1.75e-20 | 2.42e-18 | up |
| GNRH1 | 2.013 | 6.225 | 2.06e-20 | 2.75e-18 | up |

## 功能富集（前 20）

| ID | Description | GeneRatio | pvalue | p.adjust | Count |
| --- | --- | --- | --- | --- | --- |
| hsa04110 | Cell cycle | 18/143 | 4.31e-19 | 8.19e-18 | 18 |
| HALLMARK_E2F_TARGETS | E2F targets | 15/143 | 4.32e-18 | 5.61e-17 | 15 |
| HALLMARK_G2M_CHECKPOINT | G2M checkpoint | 10/143 | 2.67e-14 | 1.73e-13 | 10 |
| GO:0006260 | DNA replication | 10/143 | 6.52e-13 | 1.04e-11 | 10 |
| GO:0051607 | Defense response to virus | 9/143 | 1.49e-11 | 1.20e-10 | 9 |
| GO:0009615 | Response to virus | 9/143 | 4.71e-11 | 2.51e-10 | 9 |
| GO:0007049 | Cell cycle | 9/143 | 1.28e-10 | 5.11e-10 | 9 |
| hsa04114 | Oocyte meiosis | 8/143 | 3.31e-10 | 3.15e-09 | 8 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | Interferon gamma response | 10/143 | 1.12e-09 | 4.60e-09 | 10 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | Interferon alpha response | 9/143 | 1.42e-09 | 4.60e-09 | 9 |
| HALLMARK_MYC_TARGETS_V2 | MYC targets V2 | 7/143 | 2.19e-09 | 5.27e-09 | 7 |
| HALLMARK_DNA_REPAIR | DNA repair | 8/143 | 2.43e-09 | 5.27e-09 | 8 |
| GO:0051301 | Cell division | 6/143 | 1.88e-09 | 6.03e-09 | 6 |
| GO:0000278 | Mitotic cell cycle | 6/143 | 1.28e-08 | 3.41e-08 | 6 |
| GO:0002250 | Adaptive immune response | 7/116 | 4.32e-09 | 4.75e-08 | 7 |
| GO:0006281 | DNA repair | 7/143 | 4.40e-08 | 1.01e-07 | 7 |
| hsa03030 | DNA replication | 8/143 | 7.10e-08 | 4.50e-07 | 8 |
| GO:0002443 | Leukocyte mediated immunity | 5/116 | 3.78e-07 | 2.08e-06 | 5 |
| HALLMARK_UV_RESPONSE_DN | UV response down | 6/143 | 1.46e-06 | 2.70e-06 | 6 |
| GO:0008283 | Cell population proliferation | 6/143 | 1.46e-06 | 2.91e-06 | 6 |

## 免疫浸润（前 10 样本）

| sample | B_cells_naive_CIBERSORT | B_cells_memory_CIBERSORT | Plasma_cells_CIBERSORT | T_cells_CD8_CIBERSORT | T_cells_CD4_naive_CIBERSORT | T_cells_CD4_memory_resting_CIBERSORT | T_cells_CD4_memory_activated_CIBERSORT |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C01 | 0.01521 | 0 | 0.05374 | 0 | 0 | 0 | 0 |
| C02 | 0.009275 | 0 | 0.05685 | 0 | 0 | 0.1893 | 0.0177 |
| C03 | 0.04258 | 0 | 0.03485 | 0 | 0 | 0.1494 | 0 |
| C04 | 0 | 0 | 0.09009 | 0.003642 | 0 | 0.004166 | 0 |
| C05 | 0 | 0 | 0.05688 | 0.0598 | 0 | 0.01718 | 0.01291 |
| C06 | 0 | 0 | 0 | 0 | 0 | 0.01777 | 0.005897 |
| C07 | 0 | 0 | 0.06208 | 0.2008 | 0 | 0.005695 | 0 |
| C08 | 0 | 0 | 0.06724 | 0.1668 | 0 | 0 | 0 |
| C09 | 0.04184 | 0 | 0.02441 | 0 | 0 | 9.68e-04 | 0 |
| C10 | 0.004487 | 0 | 0.04617 | 0.1112 | 0 | 0.1024 | 0 |

## 分析记录

- 2026-11-12T13:02:38.559Z `qc` ✓ 36 个样本：文库大小中位数 520790，检出基因中位数 4000
- 2026-11-12T13:02:38.654Z `filter` ✓ 基因：4000 → 4000（去掉 0 个低表达基因，阈值 0.0 × 1 样本）
- 2026-11-12T13:02:38.882Z `normalize` ✓ 标准化方法：cpm（含 log2 转换）
- 2026-11-12T13:02:39.065Z `batch` ✓ combat：2 个批次（batch），批次解释度 R² 0.159 → 0.000
- 2026-11-12T13:02:39.432Z `pca` ✓ PCA：20 个主成分，PC1 解释 20.9%，PC1+PC2 共 24.2%
- 2026-11-12T13:02:40.939Z `de` ✓ limma（R 官方）：case vs control（18 vs 18 样本）
- 2026-11-12T13:02:41.170Z `enrich` ✓ ORA：GO_BP+KEGG+HALLMARK，共 77 条，显著 63 条（p.adjust≤0.05）
- 2026-11-12T13:02:54.480Z `immune` ✓ R/IOBR成功 5 个免疫浸润方法：cibersort、epic、estimate、mcpcounter、ssgsea
- 2026-11-12T13:02:54.567Z `survival` ✓ 分组（group）log-rank P=0.0823
- 2026-11-12T13:02:55.070Z `prognosis` ✓ Cox 比例风险模型（lifelines）：10 个基因，C-index=0.747
- 2026-11-12T13:02:55.352Z `scmap` ✓ 参考：demo_sc.h5ad（400 细胞 × 800 基因，6 个细胞类型；细胞类型列 cell_type，来源 params.ref_path）

_报告生成时间：2026-11-12 21:02:55_