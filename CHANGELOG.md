# Changelog

## 0.1.2 (2026-11-16)

按《DSH_GEO_Singlecell_审查与优化任务单》修复结果正确性问题（每条都带回归测试）。

### 差异表达
- **GEO-01** 自动推断的 case/control 方向反转：`_infer_pair` 改为带字段名的 `PairGuess`，调用方按字段取用；两种样本顺序推断一致，与显式指定的 logFC 同号同值
- **GEO-02** DESeq2 结果列名先统一再过滤（修复前整表被清成 0 行）；R 程序错误用退出码 4 区分，不再静默回退 Welch；metrics 记录 requested/actual_method
- **GEO-03** 设计矩阵改为显式 data.frame（数值协变量保持 numeric）、`model.matrix(..., data=)`；配对以 `pair` 固定效应进 limma/DESeq2，新增 `pair_col` 参数，缺配对标 ID／不成对／共线变量都明确失败

### 标准化与富集
- **GEO-07** `value_kind` 反映真实尺度（`none`+不取对数不再标成 lognorm）；`auto` 对已是 log 尺度的输入默认不再改动；分位数标准化按平均秩插值，与 `preprocessCore::normalize.quantiles` 逐值一致（ties 不再被拆开）
- **GEO-08** GSEA 优先调用 R `fgsea::fgseaMultilevel` 并记录 actual_method／版本／置换策略／覆盖率；自实现的 NES 改为官方定义（按符号分别用同号置换 ES 均值归一化），单侧 p 值分母改为同号置换次数，并在结果里标明是自实现

### 临床生存 / 预后
- **GEO-04** 事件列逐值清洗（Dead/Alive、`0:LIVING`/`1:DECEASED`、混合编码、非 0/1 数字、缺失）；未知状态排除并记录，不再当删失；清洗后的 time/event 写回 frame，KM/Cox/预后/ROC 共用同一口径
- **GEO-05** 基线生存直接取自 `CoxPHFitter.baseline_survival_`；校准/列线图/DCA 的个体概率统一走 `predict_survival_function`（修复前把协变量效应施加两遍）；缺必需列时明确失败，不再合成 `exp(-0.05t)`；DCA 改为 KM 版净获益（不再把时点前删失当阴性）
- **GEO-06** 结果明确标注「训练集表观性能」；「RSF-lite」更名为「随机特征子集 + 患者 bootstrap 的 Cox 集成」（不再冒用 RSF 名义，修正 1/2 基因时抽样报错），校准/DCA 用对集成风险重拟合的单变量 Cox 并写明

### 工程
- 新增 `runtime` 错误码，区分「程序错误」与「数据/依赖不满足」
- 新增回归测试：`py_de_regression_test.py`、`py_prognosis_regression_test.py`、`py_prep_regression_test.py`、`py_gsea_regression_test.py`

## 0.1.1 (2026-11-16)

- 修复 README 里姊妹项目链接（`dsh-singlecell` → `singlecell-full-workflow`）
- 新增 512×512 PNG 图标，用于插件目录（`dsh-plugin` 生态要求位图）
- 新增 `LAUNCH.md` 发布文案与 `dsh-geobulk-0.1.1.tgz` npm 包形态的 Release 资产
- README 增加「从 Release tarball 安装」方式（插件市场同款）

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
