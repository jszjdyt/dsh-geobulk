# Changelog

## 0.1.0 (2026-11-12)

首个公开版本。

### 工作流
- 12 步流程：`qc → filter → normalize → batch → pca → de → enrich → immune → survival → prognosis → scmap → ncrna`
- 70 项高级分析（基本/差异/通路/免疫/临床/预后/单细胞/非编码），按前置条件自动可用性判断
- 确定性 self-audit、下一步建议、4 套预设流程（standard / geo / quick / clinical）

### 数据
- GEO 检索：E-utilities + 192 条中英词典，中文查询自动翻译并逐级放宽；`matrix/` 与 `suppl/` 双目录解析；下载目录可指定
- 导入：matrix_text / xlsx / geo_series / archive / rds / rdata / circ / mirna / lncrna / fusion / splicing，样本表与临床表自动 join
- 导出：CSV 表、summary JSON、protocol、Methods、Markdown 报告、可复现脚本

### 方法
- 差异表达：limma（R 官方）+ DESeq2（counts）+ Welch 回退
- 富集：ORA 超几何 + GSEA；内置 GO/KEGG/Hallmark/Reactome + Enrichr 缓存
- 免疫浸润：IOBR（CIBERSORT/EPIC/MCP-counter/xCell/quanTIseq/ESTIMATE/ssGSEA/TIMER/IPS）+ Python 回退
- 生存与预后：KM/log-rank/Cox/cutoff/时间 ROC/亚组；Cox/LASSO-Cox/RSF-lite、C-index、校准、DCA、列线图数据
- 单细胞映射：h5ad 参考（含从 dsh-singlecell 拉取）、mean/AUCell/AddModuleScore
- 非编码：circRNA/miRNA/lncRNA/fusion/可变剪接

### 工程
- Host 路由 + Agent 工具同一套实现；纯函数层 + 组件渲染烟测可离线跑
- CI：语法检查 + 不依赖 R 的单测
