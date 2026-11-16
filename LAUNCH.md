# 发布文案（LAUNCH）

> 用途：发到 V2EX / 知乎 / 小红书 / 生信群 / 朋友圈 时直接复制对应段落。
> 链接统一用：https://github.com/jszjdyt/dsh-geobulk

## 一句话

**把「搜 GEO → 清洗 → 差异 → 富集 → 免疫浸润 → 生存预后 → 映射回单细胞」塞进了 DeepSeek Harness 的对话界面里，点着用，也能一句话让 Agent 跑完。**

## 30 秒版本（V2EX / 推文 / 群）

> 做了个 DeepSeek Harness 插件：**dsh-geobulk**，把 bulk 转录组分析做成可视化工作台。
>
> - GEO 检索**中文可搜**（胃癌 顺铂 耐药 表达谱），选目录下载、一键导入清洗
> - 自己的 counts/TPM 矩阵 + 样本表/临床表直接拖进来，走同一条流水线
> - 清洗（QC/标准化/ComBat 去批次）→ limma 差异 → GO/KEGG/Hallmark 富集
> - 免疫浸润走 IOBR 官方（CIBERSORT/EPIC/MCP-counter/xCell/quanTIseq/ESTIMATE/ssGSEA/TIMER）
> - KM/Cox 生存、Cox/LASSO-Cox 预后、时间 ROC、校准、DCA
> - **能把 bulk 的差异签名映射到你单细胞数据集的细胞亚群**（可直接从 dsh-singlecell 拉当前数据集）
> - circRNA/miRNA/lncRNA/fusion/可变剪接也能读
>
> 安装：`plugin_manager action=install_bundle target=<release tgz 链接>`（见 README）
> 仓库：https://github.com/jszjdyt/dsh-geobulk · MIT

## 90 秒版本（知乎 / 公众号 / 生信技能树投稿）

**为什么做这个**：GEO 上找一个能用的表达谱，通常要在网页检索、FTP 下载、R 里改脚本、拼样本表、查基因名注释之间来回折腾几个小时；做完差异，还得另开一套代码做免疫浸润和生存；最后想看看"这个基因到底在哪些细胞里表达"，又要去单细胞数据里重新找。

**dsh-geobulk 把这些串成一条流水线，放在聊天界面旁边**：

1. **检索**：中文关键词直接搜（内置 192 条中英词典，中文自动翻英文；结果太少会自动放宽关键词），结果卡片显示样本数、平台、是否有临床信息。
2. **导入**：GEO series matrix、counts/TPM/FPKM 矩阵（csv/tsv/xlsx/gz）、zip/tar、R 的 rds/RData 都行；样本表和临床表按样本名自动 join。
3. **清洗**：样本质控 → 低表达过滤 → log/分位数/TMM 标准化 → ComBat 去批次（出校正前后 PCA 对比）→ 样本聚类。
4. **分析**：limma（R 官方）差异、GO/KEGG/Hallmark/Reactome 富集（ORA+GSEA）、IOBR 免疫浸润（多种反卷积方法 + 组间差异）、KM/Cox 生存、Cox/LASSO-Cox 预后（风险图、时间 ROC、C-index、校准、DCA）。
5. **回单细胞**：把差异基因/签名打到单细胞参考的细胞类型上，出"细胞类型×得分"和"目标分子在各亚群表达"的表和图；参考 h5ad 可以直接从 `dsh-singlecell` 工作台拉当前数据集。
6. **非编码**：circRNA（差异/背靠背比值/圈图）、miRNA（靶基因/ceRNA）、lncRNA（分类/共表达网络）、fusion（STAR-Fusion/Arriba 解析 + 癌基因/激酶过滤）、可变剪接（rMATS/PSI 的 dPSI）。
7. **出结果**：一键导出表达表/差异表/富集表/免疫表 + 英文 Methods + Markdown 报告 + 可复现脚本。

**工程上值得一提的**：页面按钮和 Agent 工具走**同一套代码**（改一处两边同步）；Python 计算引擎 + 可选 R 桥，R 不在就自动回退到内置实现并诚实标注 `lite`；内置示例数据可离线跑通全流程；有确定性自检（0–100 分）告诉你流程哪里不靠谱。

**适合谁**：做肿瘤/免疫/生信分析、需要快速从 GEO 挖数据出图出表的研究生和博后；已经有单细胞数据、想用 bulk 做验证的人；不想每次都重写 R 脚本的人。

**不适合谁**：需要从 FASTQ 比对开始（本插件是矩阵级入口）；需要临床级结论（预后模型必须外部验证）。

## 小红书 / 短图文版

标题：**把 GEO 数据分析和单细胞打通了，全程点鼠标**

正文：
做 bulk 转录组最烦的几件事：GEO 找数据、下载解压、改成能分析的矩阵、做完差异还要做免疫浸润和生存……
我把它做成了一个工作台插件，装在 DeepSeek Harness 里：

🔍 中文搜 GEO（比如"胃癌 顺铂 耐药"），选目录下载
🧹 质控/标准化/去批次 一键跑，出前后对比图
📊 差异基因 → 火山图/热图 → GO/KEGG 富集
🛡️ 免疫浸润（CIBERSORT 等 8 种方法）→ 组间差异
⏳ KM 生存曲线 + Cox 森林图 + 预后风险评分
🧬 最后把差异基因映射到单细胞亚群：哪些细胞在表达，一目了然
🧩 circRNA / miRNA / lncRNA / 融合基因 / 可变剪接 也能读

开源免费（MIT），仓库：github.com/jszjdyt/dsh-geobulk

## 一句话发群里

> 分享个自己写的开源工具：DeepSeek Harness 的 bulk 转录组工作台插件，GEO 中文检索下载 → 清洗 → 差异/富集/免疫浸润/生存预后 → 映射回单细胞，全流程可视化，也能直接让 Agent 跑。MIT，https://github.com/jszjdyt/dsh-geobulk

## 配图建议（仓库 examples/ 里现成）

| 用途 | 文件 |
|---|---|
| 首图（差异结果） | `examples/de-volcano.png` |
| 免疫浸润 | `examples/immune-stack.png` |
| 预后风险 | `examples/risk-score.png` |
| 批次校正 | `examples/batch-pca.png` |
| 单细胞映射 | `examples/scmap-violin.png` |
| 富集 | `examples/enrich-bar.png` |

## 发布后看数据

- 仓库：Insights → Traffic（14 天 clone / view），只有仓库主人能看
- Release 资产下载数：Release 页每个 asset 右侧显示（`dsh-geobulk-*.tgz`）
- 插件目录：上传 `data/plugins/jszjdyt__dsh-geobulk.yml` 后，站点会定期刷新 stars/downloads
- npm（如果之后发布）：`https://api.npmjs.org/downloads/point/last-week/dsh-geobulk`
