# 安装 dsh-geobulk

## 方式一：本机安装（推荐）

本仓库就是一个插件包（克隆到任意目录即可）。

1. 确认 Python 可用：

```bash
python3 -c "import numpy, pandas, scipy, sklearn, matplotlib, anndata, lifelines, statsmodels, openpyxl"
```

缺包时（一般不需要，本机已装）：

```bash
python3 -m pip install -r requirements.txt
```

2. 确认 R（可选但强烈建议，官方实现更快更准）：

```bash
Rscript -e 'for (p in c("limma","edgeR","DESeq2","clusterProfiler","fgsea","msigdbr","GSVA","WGCNA","preprocessCore","timeROC","rms","glmnet","survminer","IOBR","org.Hs.eg.db")) cat(p, requireNamespace(p, quietly=TRUE), "\n")'
```

3. 在 Harness 里安装（需要 Full access 权限）：

- 设置 → 插件 → 安装本地插件包 → 选择克隆目录；
- 或让 Agent 调 `plugin_manager`：`action=install_bundle, target=<克隆目录的绝对路径>`。

4. 重载页面：侧边栏出现「Bulk 工作台」，对话里出现 `geobulk` 工具。

## 方式二：手动加入 profile（高级）

在 `~/.dsh/profiles/<profile>/package.json` 的 `dependencies` 加：

```json
"dsh-geobulk": "link:/克隆目录的绝对路径"
```

并在 `dsh.profile.bundles` 里追加 `"dsh-geobulk"`，然后重载宿主。
（一般不推荐手改 profile，用上面的 `plugin_manager install_bundle` 更稳。）

## 首次使用

1. 打开侧边栏「Bulk 工作台」。
2. 顶栏点「载入示例数据」→ 立刻有 36 样本 × 4000 基因、带分组/批次/临床列的数据。
3. 左栏从上到下点「运行」：质控 → 过滤 → 标准化 → 降维 → 差异 → 富集 → 免疫 → 生存 → 预后 → 单细胞映射。
4. 想用真实数据：顶栏「GEO 检索」搜关键词（如 `胃癌 顺铂 耐药`）→ 详情 → 填下载目录 → 下载 → 导入并清洗；或「导入数据」选自己的矩阵 + 样本信息表。

## 常见问题

**Q：引擎徽章显示「Python 未就绪」？**
点徽章看错误原因。常见是 Python 路径找不到：在插件配置里指定 `python`（绝对路径），或在环境变量 `DSH_GEOBULK_PYTHON` 指定。日志在 `<临时目录>/dsh-geobulk/engine.log`。

**Q：GEO 检索很慢或失败？**
NCBI E-utilities 需要能访问外网。失败会给出中文原因；可稍后重试，或用 `open` 直接导入已下载的 series matrix。

**Q：免疫浸润只有 2–3 个方法出结果？**
IOBR 的部分方法需要足够基因重叠（xCell ≥520、CIBERSORT 需要 LM22 基因）。数据基因少时会自动跳过并写 warning，其余方法照常返回。

**Q：没有 R 能用吗？**
能。所有 R 实现都有 Python 回退（limma→Welch+BH、clusterProfiler→超几何+内置 GMT、GSVA→ssGSEA、IOBR→内置签名打分），结果会标注 `lite`。

**Q：怎么和单细胞工作台联动？**
先打开 `dsh-singlecell` 并载入数据；在 Bulk 工作台「单细胞」tab 点「从单细胞工作台拉取」，或直接填本机 h5ad 路径。插件会调用单细胞工作台导出当前数据集后读取。

**Q：数据/图存在哪？**
工作目录默认系统临时目录下的 `dsh-geobulk/`（图表在 `analyses/`，导出在 `exports/`，GEO 下载在 `downloads/GEO`）。可在插件配置 `workdir` 指定固定目录；GEO 下载目录在检索弹层里可每次指定。
