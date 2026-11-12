# dsh-geobulk 内部契约（GEO / Bulk 转录组工作台）

> 本文件是并行开发者的**唯一事实来源**。实现与本契约冲突时以契约为准；确需改协议，先改本文件再改代码，并在团队任务里标注。
> 参考实现在 `../dsh-singlecell/`（同作者、同宿主、同形态），**不要修改它**。

## 0. 目标

在 DeepSeek Harness Web GUI 里提供一个**可视化、可操作、Agent 可驱动**的 Bulk 转录组工作台：

1. **GEO 数据检索**：中文界面内搜索 GEO Series（E-utilities），中文化字段、可筛选（物种/平台/样本数/是否有补充矩阵/是否有临床信息），点选后**指定本地目录下载**。
2. **自有测序数据导入**：表达矩阵（counts/TPM/FPKM，csv/tsv/txt/xlsx/gz）、样本信息表、临床表、GEO series matrix、zip/tar 包、h5ad（仅用于单细胞映射）。
3. **清洗**：样本质控 → 过滤 → 标准化/log 转换 → 去批次（ComBat / limma removeBatchEffect）→ 降维聚类。
4. **单数据集分析 + 多数据集合并去批次**。
5. **出图与数据**：样本×基因表达表、差异表达可视化、功能富集、免疫浸润（尽可能全）、临床生存、高级预后。
6. **单细胞映射**：把 bulk 的分子/签名映射到单细胞数据集的细胞亚群上（可直接从 dsh-singlecell 工作台拉当前数据集）。
7. **非编码与结构变异**：circRNA、miRNA、lncRNA、fusion gene、可变剪接（rMATS/SUPPA 等）。
8. UI 风格：清新淡雅 + 高科技感（浅色、留白、细边框、微渐变、等宽数字），交互逻辑对齐单细胞插件。

**硬性要求**：Python 未就绪时 UI 必须给可读原因而不是白屏；单个步骤/分析失败不影响其它步骤与整页；所有重计算有超时与进度；无网络时 GEO/在线富集给可读降级。

---

## 1. 目录与写权限（严格隔离，禁止越界改文件）

```

  CONTRACT.md                     # Lead
  README.md / INSTALL.md          # Lead
  package.json / cordis.patch.yml / icon.svg / locale/*.json   # Lead
  index.js                        # Lead  Host 入口（/geobulk 前缀路由 + geobulk 工具 + prompt section）
  lib/impl.mjs                    # Lead  唯一操作实现（HTTP + 工具共用）+ 路由分发 + 工具 + prompt
  lib/engine.mjs                  # Lead  Python worker 管理（探测/spawn/队列/超时/重启）
  lib/bg_engine.py                # Lead  Python worker 主循环（协议 §2 + 命令分发）
  lib/bg_registry.py              # Lead  模块汇总 + ctx 构造 + 步骤/分析分发 + 结果整形
  lib/bg_state.py                 # Lead  Dataset 对象、样本表、demo 数据、history
  lib/bg_common.py                # Lead  EngineError / 结果构造器 / 绘图工具 / 表格工具
  lib/rbridge.py                  # Lead  通用 Rscript 桥（bundle 导出 + 调用 + 结果读回）
  lib/bg_io.py                    # bulk-core   全格式导入
  lib/bg_geo.py                   # bulk-core   GEO 检索/详情/下载
  lib/bg_prep.py                  # bulk-core   步骤 qc/filter/normalize/batch/pca + 基础分析
  lib/bg_de.py      + bg_de.R     # bulk-de     步骤 de + 差异类分析（limma/DESeq2 官方 R）
  lib/bg_enrich.py  + bg_enrich.R # bulk-de     步骤 enrich + 通路类分析（clusterProfiler/fgsea/GSVA）
  lib/bg_genesets.py              # bulk-de     内置基因集/通路回退 + 在线库缓存
  lib/bg_immune.py  + bg_immune.R # bulk-clin   步骤 immune + 免疫浸润全方法（IOBR 优先）
  lib/bg_survival.py              # bulk-clin   步骤 survival + 临床生存分析（lifelines）
  lib/bg_prognosis.py             # bulk-clin   步骤 prognosis + 高级预后
  lib/bg_scmap.py   + bg_scmap.R  # bulk-ncrna  步骤 scmap + 单细胞映射
  lib/bg_ncrna.py   + bg_ncrna.R  # bulk-ncrna  步骤 ncrna + circ/miRNA/lncRNA/fusion/AS
  client.js                       # bulk-ui     浏览器 UI
  test/                           # 各自作用域，见 §11
  examples/ figures/              # Lead
```

- 只写自己作用域的文件。需要别人改：发消息 + 写进团队任务。
- **所有文件写之前先 `read`**；改别人的文件一律禁止。

---

## 2. Python worker 协议（`lib/bg_engine.py`）

### 2.1 启动与传输

- Host 以 `python3 -u lib/bg_engine.py --workdir <abs_dir>` 启动。
- stdin 一行一个 JSON 请求；stdout 一行一个 JSON 响应（UTF-8，`ensure_ascii=False`，`flush=True`）。**stdout 只允许协议行**，日志一律 stderr（Host 收进 `<workdir>/engine.log`）。
- 请求 `{"id":"<str>","cmd":"<name>","args":{...}}`；成功 `{"id":..,"ok":true,"result":{...}}`；失败 `{"id":..,"ok":false,"error":{"code":"bad_request|not_found|prerequisite|unsupported|ambiguous|internal|timeout","message":"中文可读原因","trace":"可选","candidates":[...]?}}`。
- 进度（可选）：`{"id":..,"event":"progress","message":"...","pct":0-100}`。
- 单条命令异常只影响该条；进程只在 `shutdown` 或 stdin EOF 退出。启动时不得向 stdout 打印任何东西。

### 2.2 命令表（字段名逐字一致）

**`ping`** `{}` →
```json
{"python":"3.14.2","platform":"darwin","mode":"builtin","workdir":"/abs",
 "packages":{"numpy":"2.5.1","pandas":"2.3.3","scipy":"1.17.1","sklearn":"1.9.0","matplotlib":"3.11.0","anndata":"0.12.19","lifelines":"0.30.3","statsmodels":"0.14.6","openpyxl":"3.1.5","requests":"2.34.2"},
 "r":{"available":true,"version":"4.4.2","script":"/Library/Frameworks/R.framework/Resources/bin/Rscript",
      "packages":{"limma":true,"edgeR":true,"DESeq2":true,"clusterProfiler":true,"GSVA":true,"IOBR":true,"fgsea":true,"msigdbr":true,"survival":true,"timeROC":true,"rms":true,"glmnet":true,"survminer":true,"WGCNA":true,"ComplexHeatmap":true,"pheatmap":true,"preprocessCore":true,"org.Hs.eg.db":true,"org.Mm.eg.db":true}},
 "capabilities":{"formats":[...],"steps":[...],"analyses":[...],"geo":true,"online":true}}
```

**`open`** args `{"path"?:str,"url"?:str,"name"?:str,"data_format"?:"auto|matrix_text|geo_series|xlsx|archive|h5ad|rds|rdata|mtx|circ|mirna|lncrna|fusion|splicing","member"?:str,"object"?:str,"transpose"?:bool,"delimiter"?:str,"sheet"?:str,"gene_col"?:str,"sample_sheet"?:str,"clinical"?:str,"species"?:"human|mouse|auto","demo"?:bool,"dest_dir"?:str}` → `dataset`（§2.4）。
- 同一路径重复 open 返回已存在数据集。`data_format` 省略=自动识别。
- `demo:true` 载入内置示例（§2.7）。

**`list`** `{}` → `{"datasets":[<dataset 摘要>]}`

**`close`** `{"dataset":id}` → `{"closed":true}`

**`preview`** `{"dataset":id}` →
```json
{"genes":["首 50 个基因"],"sample_columns":[{"key":"group","kind":"categorical|numeric|time|event","n_unique":2}],
 "n_samples":48,"n_genes":20530,"steps":["qc","filter"],"current":"filter","modules":["circ"]}
```

**`samples`** `{"dataset":id,"op":"list|update|join|derive|set_group|auto_group|remove","samples"?,"groups"?,"column"?,"group_regex"?,"path"?,"key"?,"name"?,"value"?,"source"?}` → `dataset`（摘要）+ `{"sample_info":[...]}`。
- `list` 返回样本表；`join` 把 CSV/TSV/XLSX 的临床/分组表按样本名（自动匹配行名或指定 `key` 列）join 进 sample_info；`update` 改样本显示名/大组；`auto_group` 用文件名/样本名正则（第 1 捕获组=大组）；`derive` 用正则从样本名派生新列；`set_group` 设为当前分组列；`remove` 剔除样本。

**`run`** `{"dataset":id,"step":str,"params":{}}` → `{"step":..,"summary":{"title","bullets":[],"metrics":{}},"steps":[...],"current":..,"dataset":<摘要>,"tables":[...]?,"figures":[...]?}`。
- 缺前置 → `prerequisite`。重跑某步清空其下游结果（`history` 保留）。

**`plot`** `{"dataset":id,"kind":str,"params":{}}` → `{"file":"/abs/*.png","width":int,"height":int,"bytes":int,"kind":str}`
`kind` ∈ `qc_box` `qc_density` `pca_scatter` `sample_corr` `sample_dendro` `volcano` `ma_plot` `de_heatmap` `gene_box` `enrich_bar` `enrich_dot` `immune_stack` `immune_corr` `km_curve` `forest` `risk_plot` `roc_time` `batch_pca` `circos` `splicing_heat`。

**`analyses`** `{"dataset":id}` → `{"analyses":[{...注册项...,"available":bool,"reason":str?}],"groups":["基本分析",...]}`（§2.5）

**`analysis`** `{"dataset":id,"id":str,"params":{}}` → 分析结果（§2.5 的 result 信封）+ `{"analysis_id":..,"ms":int}`

**`attach`** `{"dataset":id,"module":"circ|mirna|lncrna|fusion|splicing","path":str,"data_format"?:str,...}` → `{"dataset":..,"module":{...}}`（非编码模块数据挂载）

**`export`** `{"dataset":id,"what":"expr_csv|de_csv|enrich_csv|immune_csv|sample_info_csv|summary_json|protocol|methods|report|script|h5ad?","path"?:str}` → `{"file":"/abs/path","bytes":int,"what":str}`

**`audit`** `{"dataset":id}` → `{"score":0-100,"level":"ok|warn|fail","checks":[{"id","level","title","detail","hint"}]}`（确定性，不调模型）

**`suggest`** `{"dataset":id}` → `{"suggestions":[{"id","title","detail","step"?, "params"?}]}`（确定性）

**`recipe`** `{"id"?:str}` → 无 id：`{"recipes":[{id,title,detail,steps:[{step,params}]}]}`；有 id：该 recipe。预设 ≥3：`standard`（qc→filter→normalize→pca→de→enrich→immune）、`geo`（qc→filter→normalize→batch→pca→de）、`quick`（qc→filter→normalize→de）。

**`inspect`** `{"path"|"url",...}` → `{"format","n_samples"?,"n_genes"?,"candidates":[...]?,"notes":[],"size":int}`（不加载整表）

**`geo_search`** `{"term":str,"species"?:"human|mouse|all","study_type"?:"expression profiling by array|high throughput sequencing|all","min_samples"?:int,"has_clinical"?:bool,"limit"?:int,"page"?:int,"sort"?:"relevance|date|samples"}` →
```json
{"total":297,"page":0,"items":[{"gse":"GSE13861","title":"...","title_zh":"(机器/词典中文化，可为空)",
 "summary":"...","species":"Homo sapiens","n_samples":196,"platform":"GPL6884",
 "study_type":"expression profiling by array","release":"2010-01-01",
 "has_matrix":true,"has_clinical":true,"disease":"gastric cancer","tags":["胃癌","表达谱"]}]}
```
- 数据源：NCBI E-utilities（esearch db=gds + esummary + efetch gds brief）。必须容错：网络失败 → `EngineError("unsupported","NCBI 暂时不可达：...")`，带重试（≤2 次，退避 1s/2s）。
- 中文化：内置疾病/组织/技术词典（中英对照 ≥120 条，见 `bg_geo.py` 的 `ZH_DICT`），把 title/summary/disease 里的常见词替换；`title_zh` 允许为空字符串但要尽量生成。

**`geo_detail`** `{"gse":"GSE13861"}` →
```json
{"gse":"GSE13861","title":"...","summary":"...","design":"...","species":"...","n_samples":196,
 "platforms":["GPL6884"],"pubmed":"...","release":"...","contacts":"...",
 "clinical":{"from_series_matrix":true,"fields":[{"key":"age:ch1","n":180,"example":"62"}]},
 "files":[{"name":"GSE13861_series_matrix.txt.gz","size":12345,"url":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13861/suppl/GSE13861_series_matrix.txt.gz","kind":"series_matrix|raw|supplementary|other"}],
 "suggested":"GSE13861_series_matrix.txt.gz"}
```
- 实现：`esummary` + GEO FTP suppl 目录 HTML 解析（`https://ftp.ncbi.nlm.nih.gov/geo/series/<前缀>nnn/<GSE>/suppl/`，前缀规则：GSE 后数字去掉末 3 位 + `nnn`，如 GSE13861 → GSE13nnn）。

**`geo_download`** `{"gse":str,"files":[str],"dest_dir":str,"overwrite"?:bool}` → `{"dir":"/abs","files":[{"name","path","bytes","ok"}],"ms":int}`
- `dest_dir` 必须存在或可创建；越权路径一律拒绝（只允许绝对路径且不是 `/`；用户显式给出）。
- 下载用 `urllib`（不依赖 requests）+ Range 断点不必实现；单文件上限 8GB；超时 120s；进度事件走 `progress`。

**`sc_refs`** `{"path"?:str}` → 无 path：探测本机可能的单细胞 h5ad（`<workspace>/plugins/dsh-singlecell/exports/*.h5ad`、`~/.dsh/**` 不扫）→ `{"refs":[{"path","name","bytes"}]}`；有 path：`{"path","n_obs","n_vars","obs_columns","cell_types"?}`（anndata 只读头）。

**`shutdown`** `{}` → `{"bye":true}` 后退出。

### 2.3 `ctx` 对象（传给所有步骤/分析模块）

```python
ctx = {
  "dataset_id": "ds1", "name": "GSE13861",
  "genes": ["A1BG", ...],            # list[str]，长度 = X.shape[1]
  "samples": ["GSM1", ...],          # list[str]，长度 = X.shape[0]
  "X": np.ndarray float32,           # (n_samples, n_genes) 当前矩阵；已标准化则为 log2 尺度
  "X_kind": "counts|lognorm|tpm|fpkm|unknown",
  "counts": np.ndarray|None,         # 原始 counts（若有）
  "sample_info": pd.DataFrame,       # index = samples；列含 group/batch/临床列
  "group_key": "group"|None,         # 当前分组列
  "groups": ["case","control"]|None, # group_key 的取值（有序）
  "species": "human|mouse|unknown",
  "de": dict|None,                   # 最近一次 de 结果（见 §2.5）
  "results": {},                     # 跨分析缓存：{"enrich":..,"immune":..,"survival":..}
  "modules": {},                     # 非编码模块数据 {"circ":{"X","info","values",...}}
  "rng": np.random.Generator,        # 固定 seed=0
  "workdir": str, "outdir": str,     # outdir = <workdir>/analyses/<dataset>/<analysis>
  "engine": {"mode":"builtin","python":"..","r":{...}},
}
```
- **只读**：模块不得原地修改 ctx 的 X/sample_info（需要改就复制）。写回数据集由 registry 层根据返回值里的 `sample_columns`/`obs_updates` 完成。
- 大数据上限：bulk 场景样本数通常 ≤ 2000、基因 ≤ 60000；若 n_samples > 5000 或 n_genes > 40000，模块内部自行下采样并在 warnings 说明。

### 2.4 `dataset` 对象

```json
{"id":"ds1","name":"GSE13861","path":"/abs/x.txt.gz","source_format":"geo_series",
 "n_samples":196,"n_genes":20530,"species":"human",
 "sample_columns":[{"key":"group","kind":"categorical","n_unique":2},{"key":"os_time","kind":"time"}],
 "steps":["qc","filter","normalize"],"current":"normalize",
 "modules":[{"kind":"circ","n":312,"path":"..."}],
 "warnings":["..."],"notes":["..."],
 "history":[{"at":"ISO","step":"de","params":{},"ok":true,"ms":812,"summary":"...","engine":{"python":"3.14.2","r":"4.4.2"}}],
 "analyses_run":[{"at":"ISO","id":"volcano","ok":true,"ms":120}]}
```
`history` 最多 50 条；返回给 UI 时截断到最近 20 条。

### 2.5 步骤 / 分析模块契约

每个模块文件（`bg_prep.py` 等）导出一个或两个列表：

```python
STEPS = [{
  "id": "qc", "title": "样本质控", "icon": "qc",
  "requires": [],                       # 前置 step id 列表
  "params": [                           # UI 自动生成表单
    {"key":"min_detect_genes","label":"最少检出基因数","type":"number","default":0,"min":0,"max":30000,"step":100},
    {"key":"group_col","label":"分组列","type":"column","default":null},
    {"key":"method","label":"方法","type":"select","default":"limma","options":[{"value":"limma","label":"limma"}]},
    {"key":"genes","label":"基因列表","type":"genes","default":[]},
    {"key":"file","label":"文件路径","type":"file","default":""},
    {"key":"paired","label":"配对","type":"boolean","default":false},
  ],
  "run": callable,                      # run(ctx, params, outdir) -> result
}]

ANALYSES = [{
  "id": "volcano_custom", "title": "自定义火山图", "group": "差异分析",
  "method": "limma 结果 + matplotlib 绘制",     # 诚实写实现方式；lite/代理必须标注
  "needs": ["de"],                              # 语义前置：de / groups / clinical / modules:circ / sc / genes
  "params": [...],
  "run": callable,
}]
```

**result 信封（步骤与分析共用，键名固定）**：
```python
{"summary": {"title": "差异表达", "bullets": ["limma：case vs control，FDR<0.05 共 512 个"], "metrics": {...}},
 "tables": [{"id":"de","title":"差异基因","columns":[{"key":"gene","label":"基因","type":"string"},{"key":"logFC","label":"logFC","type":"number"}],
             "rows":[{"gene":"CDH1","logFC":2.31}], "note":"可选"}],
 "figures":[{"id":"volcano","title":"火山图","file":"/abs/x.png","width":720,"height":480}],
 "sample_columns": [{"key":"risk_group","kind":"categorical"}],   # 可选：要写回 sample_info 的新列
 "sample_values": {"risk_group": ["high","low",...]},            # 与 ctx.samples 等长
 "data_out": [{"kind":"de_csv","path":"/abs/de.csv","label":"差异基因表"}],
 "warnings": [], "notes": []}
```
- 行数上限：表格写文件完整，返回给 JSON 的 `rows` 每表 ≤ 500 行（多出时 `"truncated":true` 并写 `data_out`）。
- 图：PNG，白底，dpi≈140，**图上文字一律英文**（避免缺 CJK 字体），存 `outdir`，返回绝对路径。
- 缺前置 → 抛 `EngineError("prerequisite","中文原因")`，不得编造结果。
- 确定性：同数据同参数同结果（固定 seed）。

### 2.6 步骤清单（顺序即左栏卡片顺序）

| step | 标题 | requires | 主要参数（默认值） | 输出 |
|---|---|---|---|---|
| `qc` | 样本质控 | — | `group_col=null, min_detect_genes=0` | 每个样本 文库大小/检出基因数/中位数/离群标记；箱线图+密度图 |
| `filter` | 过滤 | qc | `min_expr=1.0, min_samples_frac=0.2, remove_outliers=false` | 过滤后 n_genes/n_samples；保留基因表 |
| `normalize` | 标准化 | filter | `method=auto(log2+quantile|median|tmm|none), log_transform=auto` | 前后分布图；X_kind=lognorm |
| `batch` | 去批次 | normalize | `batch_key=null, method=combat|removeBatchEffect|none, group_key` | 校正前后 PCA 对比图；批次混杂指标 |
| `pca` | 降维聚类 | normalize | `n_pcs=20, color_by=group, cluster_method=hclust, k=auto` | 碎石图、样本散点（带标签）、相关性热图、层次聚类树 |
| `de` | 差异表达 | normalize | `group_col, case, control, method=limma|deseq2|welch, covars=[], lfc=0, padj=0.05, paired=false` | 差异表 + 火山/MA/热图 |
| `enrich` | 功能富集 | de | `db=GO_BP,GO_CC,GO_MF,KEGG,HALLMARK,REACTOME, method=ora|gsea, pvalue=0.05, qvalue=0.2, top_n=20` | 富集表 + barplot/dotplot |
| `immune` | 免疫浸润 | normalize | `methods=[cibersort,epic,mcpcounter,xcell,quantiseq,estimate,ssgsea,timer], group_col` | 细胞比例表 + 堆叠图 + 组间差异 |
| `survival` | 临床生存 | normalize | `time_col, event_col, endpoint=OS, gene=null, cutoff=median|optimal|quartile, covariates=[]` | KM + logrank + Cox HR 森林图 |
| `prognosis` | 高级预后 | survival | `model=cox|lasso|rsf, genes=[], time_col, event_col, nfold=10` | 风险评分表 + 风险图 + 时间ROC + C-index |
| `scmap` | 单细胞映射 | de | `ref_path, gene_set=de_up|de_down|custom, genes=[], cell_type_col, score_method=mean|aucell|addmodule` | 细胞类型打分表 + 亚群表达图 |
| `ncrna` | 非编码/融合/剪接 | normalize | `kind=circ|mirna|lncrna|fusion|splicing, module_path|null` | 各模态差异/统计表 + 图 |

### 2.7 demo 数据（`bg_state.py` 负责生成）

- 2 组（`case` 18 例 / `control` 18 例）× 3000 基因 × 1 个批次列（2 批）× 临床列（`os_time`,`os_event`,`age`,`sex`,`stage`）× 患者 ID。
- 至少 120 个基因在两组间真实差异（logFC 1.5~2.5）；免疫相关基因（CD8A/GZMB/...)与 8 类免疫细胞签名基因内嵌，保证免疫浸润有信号；生存时间与一个"高危签名"正相关（保证 KM 有差异）。
- 生成 < 2 秒；`test/make_demo.py` 复用同一生成器（或直接 `import bg_state` 调用）。
- 另有个小 h5ad demo（`test/fixtures/demo_sc.h5ad`，400 细胞 × 800 基因 × 6 细胞类型，含 bulk 的差异基因）供 scmap 演示；没有时 scmap 给 `prerequisite` 提示。

### 2.8 引擎自检（`audit` 固定检查项）

`qc_pass`、`normalized`（未标准化就做 DE）、`batch_confounded`（分组与批次完全混杂）、`group_imbalanced`（某组 <3 样本）、`de_missing`、`enrich_missing`、`immune_missing`、`survival_missing`（有临床列但没跑）、`sc_missing`、`modules_missing`、`low_samples`（< 6 样本）。`score = clamp(100 - 12*Warn - 30*Fail, 0, 100)`。

---

## 3. 模块实现约定（各 owner 必读）

### 3.1 通用
- 依赖只能来自环境已有包（§10），**禁止 pip / conda 安装**、禁止改用户目录。
- 所有函数对输入做防御：`n_samples < 3`、某组样本为 0、基因不存在、全 NaN 列 → 抛 `EngineError("prerequisite", 中文)` 或降级并写 warning，不得 500 崩掉。
- 图统一走 `bg_common.pyplot_fig()` 上下文与 `bg_common.save_fig(fig, outdir, name, ...)`（设置 Agg、字体、dpi、白底、返回 {file,width,height}）。
- 表格统一走 `bg_common.table(id, title, columns, rows, note=None)`，内部做 500 行截断与 NaN→None。
- 颜色/字体由 `bg_common` 统一，避免每模块各写一套。

### 3.2 R 桥（`lib/rbridge.py`，Lead 实现，模块只调用）
```python
from rbridge import r_available, r_probe, run_r
r = run_r("bg_de.R", {"expr": X, "genes": genes, "samples": samples, "info": df, "params": {...}}, outdir, timeout=600)
# 返回：{"ok":True,"result":{...},"tables":{name: DataFrame},"figures":[{"id","title","file"}],"warnings":[],"log":"..."}
```
- 实现细节：写 bundle 到 `<outdir>/r_bundle/`（`expr.tsv` 样本×基因、`info.csv`、`params.json`），调 `Rscript --vanilla lib/<script>.R --in <bundle> --out <outdir>/r_out --params <json>`。
- 退出码：0 成功；3 缺前置（→ `EngineError("prerequisite")`）；其它 → 抛 `EngineError("internal", stderr 尾 2KB)`。
- R 脚本必须写 `<out>/result.json`（同 §2.5 信封，字段可少不可错）+ 图 PNG（英文）。
- R 不可用时模块必须**自动回退**到 Python 实现（limma→Welch+ BH；clusterProfiler→内置超几何；GSVA→ssGSEA 自实现；IOBR→内置签名打分），并在 `notes` 写明用了回退。

### 3.3 模块注册
- 模块文件顶部**只能 import 标准库 + numpy/pandas/scipy/sklearn/matplotlib**，以及 `bg_common`/`rbridge`/`bg_registry` 里的**纯函数**（不要 import bg_state，避免循环）。
- 模块内不要写 `if __name__ == "__main__"`，不要读全局状态。

---

## 4. HTTP API（Host，前缀 `/geobulk`，仅回环地址）

响应：成功 `{"ok":true, ...}`；失败 `{"ok":false,"code":"...","error":"中文"}`，HTTP 400/404/500。
请求可带 `"actor":"user"`（UI 默认）/ `"agent"`；所有操作写活动日志。

| method path | body | 说明 |
|---|---|---|
| `GET /health` | — | `{"engine":{"ready","mode","python","r","packages"},"geo":true,"workdir"}` |
| `GET /state` | — | `{"active","datasets":[摘要],"engine":{},"suggestions":[],"recipes":[],"activity":[最近50],"sc":{"available":bool,"datasets":[]}}` |
| `POST /geo/search` | `{term,species,study_type,min_samples,has_clinical,limit,page}` | GEO 检索（中文化） |
| `POST /geo/detail` | `{gse}` | 系列详情 + 文件列表 |
| `POST /geo/download` | `{gse,files,dest_dir,open_after}` | 异步任务（返回 `{job}`）；完成后可 `open_after` 自动导入 |
| `GET /jobs?id=` | — | `{"job":{id,kind,status:"running|done|failed",pct,message,result?,error?}}` |
| `POST /dataset/open` | `{path,url,name,data_format,member,object,transpose,sheet,gene_col,sample_sheet,clinical,species}` | 打开 |
| `POST /dataset/demo` | `{}` | 内置示例 |
| `POST /dataset/upload` | raw bytes + `x-dsh-filename` header | 落盘 workdir 再打开 |
| `POST /dataset/select` / `POST /dataset/close` | `{dataset}` | 切换/关闭 |
| `POST /inspect` | `{path,url,...}` | 只探测 |
| `GET /formats` | — | 支持格式清单 |
| `POST /preview` | `{dataset}` | 基因/列清单 |
| `POST /samples` | `{dataset,op,...}` | 样本表管理（§2.2） |
| `POST /run` | `{dataset,step,params}` | 跑一步 |
| `POST /plot` | `{dataset,kind,params}` | 出图 |
| `GET /analyses?dataset=` | — | 分析目录 |
| `POST /analysis` | `{dataset,id,params}` | 跑分析 |
| `POST /attach` | `{dataset,module,path,...}` | 挂载非编码模块数据 |
| `POST /export` | `{dataset,what,path}` | 导出 |
| `GET /figure?id=` | — | 回 PNG（id 只允许 `[A-Za-z0-9._-]+\.png`） |
| `GET /activity` | — | 活动日志 |
| `POST /audit` / `POST /suggest` / `POST /recipe` | `{dataset[,id]}` | 自检/建议/预设流程 |
| `POST /reset` | `{}` | 清空 |
| `GET /sc/refs` | — | 单细胞参考列表（本机 h5ad + 从 dsh-singlecell 拉取） |
| `POST /sc/attach` | `{dataset,ref_path}` | 把 h5ad 挂到当前 bulk 数据集供 scmap 用（校验存在性；若来自单细胞工作台则先请求它导出） |

**单细胞工作台联动**：`webServer.getOrigin()` 拿到宿主 origin，请求 `GET <origin>/singlecell/state`；可导出时调 `POST <origin>/singlecell/export {what:"h5ad"}` 得到路径。singlecell 不可用 → 仅本地 h5ad，`sc.available=false` 不报错。

**异步任务**：下载与单细胞映射等长任务走 `/jobs`；`POST` 返回 `{"ok":true,"job":{"id","kind","status"}}`，前端轮询。工具调用等待到完成（带上限）。

---

## 5. Agent 工具 `geobulk`

单一工具，`action` 枚举 + 扁平参数（**参数名不能叫 `format`** → 用 `data_format`）。参数严格按 JSON Schema（`type` 只能 string/number/boolean/object/array，enum 显式列出）。

```
action: status | search | detail | download | open | demo | inspect | preview | samples |
        run | plot | analyses | analysis | attach | export | audit | recipe | sc_refs | sc_attach | reset
```
关键参数：`term, gse, dest_dir, path, url, data_format, dataset, step, params, analysis, kind, what, op, module, ref_path, actor`。
返回给模型：`summary.bullets` + 表行数/前 20 行 + 图绝对路径 + 建议下一步。**图和表必须能给模型看到结论**。

工具描述里必须写清：什么时候用 `search`（GEO 检索）、`open`（自有数据）、`run`（清洗/分析）、`analysis`（高级分析目录）。

---

## 6. 浏览器端 `client.js`（bulk-ui）

### 6.1 形态（照抄 `../dsh-singlecell/client.js` 的骨架）
- 整文件 `;(function(){'use strict' ... })()`；顶部**纯函数层**（不碰 DOM、不 require react），`typeof window==='undefined'` 时挂 `globalThis.__BG_TEST__` 供 node 测试。
- 浏览器层 `window.__ModuleLoader__.load({id:'dsh-geobulk', factory(require){...}})`，**只 `require('react')`**。
- `inject:['slots','locale']`；`locale.register('dsh-geobulk',{zh,en})`；`slots.inject('main', ...)` 注册整页（key=`geobulk`），`slots.inject('sidebar.panellist', ...)` 注册侧边栏入口（id=`geobulk`，order=4）。
- 数据全走同源 `/geobulk/*`。
- 错误边界：单点异常不能打黑整页。

### 6.2 布局（整页，height:100%，min-height 72vh）
```
┌ 顶栏：「GEO · Bulk 转录组工作台」·引擎徽章(Python/R)·数据集 chips·[GEO 检索] [导入数据] [刷新]
├─ 左栏 268px：数据卡片 + 12 张流程卡片（序号/标题/状态/参数折叠/运行按钮/结果摘要）
├─ 中栏 flex：Tab 组
└─ 右栏 300px：Agent 面板（状态摘要 / 建议下一步 + 复制指令 / Agent 活动时间线）
```
**中栏 Tab**：`概览 / 表达矩阵 / 差异 / 富集 / 免疫 / 生存 / 预后 / 单细胞 / 非编码 / 记录 / 日志`
- 概览：数据集摘要卡（样本数/基因数/分组/批次/引擎）、PCA 图、样本相关性图、快速表达箱线图、AI 建议。
- 表达矩阵：可搜索基因 → 样本×基因表（虚拟滚动或分页，前 200 行）、基因箱线图、热图。
- 差异：DE 表（可排序）+ 火山图 + MA + 热图 + 导出。
- 富集：库选择（GO/KEGG/Hallmark/Reactome）+ ORA/GSEA + 富集表 + barplot/dotplot。
- 免疫：方法多选 + 细胞比例堆叠图 + 结果表 + 组间差异 + 相关性热图。
- 生存：时间/事件列选择 + 基因（可选）+ cutoff + KM 曲线 + Cox 森林图 + 分组统计表。
- 预后：模型选择（cox/lasso/rsf）+ 风险评分表 + 风险图 + 时间 ROC + C-index 卡。
- 单细胞：参考 h5ad 选择（本机 / 从单细胞工作台拉取）+ 签名来源 + 细胞类型打分表 + 亚群表达图。
- 非编码：模块挂载（circ/miRNA/lncRNA/fusion/splicing）+ 各自分析按钮 + 结果表/图。
- 记录：history 列表 + 自检（audit）+ 导出按钮组（expr/de/immune/summary/protocol/methods/report/script）。
- 日志：活动列表。

**GEO 检索弹层**（重点）：
- 搜索框（中文提示：如「胃癌 顺铂 耐药 表达谱」）、物种/研究类型/最小样本数/是否有临床 筛选、结果卡片列表（GSE 号、中英文标题、样本数、平台、年份、标签、是否有矩阵/临床）。
- 点卡片 → 详情面板：摘要、设计、平台、样本数、文件列表（勾选）、**下载目录输入框（默认 `~/.dsh/geobulk/GEO`，可点"选择目录"）** + 下载按钮 + 进度条；下载完成 → "导入并开始清洗"。
- 一键"直接下载并导入 series matrix"。

**导入弹层**：本机路径（文件/目录）、URL、文件选择、格式下拉（自动/各格式）、高级项（sheet/gene_col/transpose/member/object/species）+ 样本信息表与临床表路径；`ambiguous` → 候选列表点选重试。

### 6.3 样式（清新淡雅 + 高科技）
- 类名一律 `bg-` 前缀；**颜色只用 `--dsw-alias-*` 主题变量**（派生色用 `color-mix(in srgb, var(--dsw-alias-brand-primary) 14%, transparent)`），数据色（火山图点、免疫细胞调色板）可放 JS 常量。
- 卡片：底 `--dsw-alias-bg-layer-1`，边框 `--dsw-alias-border-l1`，圆角 12px，轻投影；页面底 `--dsw-alias-bg-base`。
- 高科技感：顶栏一条 1px 渐变线；数值用 `font-variant-numeric: tabular-nums`；关键指标用大号细体字 + 小写标签；状态点用 6px 圆点 + 微弱光晕；图区留白足；过渡 120–180ms ease。
- 字号 12/13px 为主；不使用 emoji 当图标（用内联 SVG 或文字）。
- 可见文案全部走 `T('key')`，`locale/zh.json` + `locale/en.json` 两份齐全。

### 6.4 client.js 纯函数层（可测）
`STEPS` 元数据、`TAB_KEYS`、`fmtNumber/fmtPct/fmtP(value)`、`colorForGroup`、`volcanoPoints()`、`pcaScreen()`（缩放/平移坐标变换）、`buildSearchPayload()`、`suggestionInstruction()`（建议→中文指令）、`apiClient(fetchImpl)`（method/url/body 正确）、`jobPolling()`、`jobDone()`。测试 `test/client-logic.test.mjs` 覆盖。

---

## 7. 数据格式与处理规则（`bg_io.py`）

### 7.1 必读格式
| data_format | 识别 | 规则 |
|---|---|---|
| `matrix_text` | `.csv/.tsv/.txt`（可 `.gz`） | 分隔符嗅探 `, \t ; 空白`；首列=基因名（或首行=样本名）；自动判方向（行名像基因 → genes×samples 转置；行名像样本 → samples×genes）；`transpose` 强制；`gene_col` 指定基因列；数值列非数值行丢弃并记 warning |
| `xlsx` | `.xlsx/.xls` | openpyxl 读第一个 sheet（或 `sheet` 指定）；同 matrix_text 规则；多 sheet 时返回 candidates |
| `geo_series` | 文本含 `!series_matrix_table_begin` | 取 table 段，首列 ID_REF=基因，其余列=样本 → 转置；`!Sample_*` 行变 sample_info 列（去 `!Sample_` 前缀、值去引号、同名去重）；`!Sample_title` 可用作样本名 |
| `archive` | `.zip/.tar/.tar.gz/.tgz` | 解压到 `<workdir>/unpack/<uuid>`（防穿越、跳符号链接、总大小上限 20GB）；找到系列矩阵/csv/10x 目录；多候选 → `ambiguous`+candidates，`member` 指定 |
| `rds`/`rdata` | 大小写不敏感 | R 桥读取：`ExpressionSet`/`SummarizedExperiment`/`Seurat`/matrix/data.frame → matrix + pdata |
| `h5ad` | `.h5ad` | 仅用于单细胞参考（`sc_refs`/`scmap`）；不当作 bulk 数据集打开 |
| `circ` | `.tsv/.csv/.txt` | circRNA 表达矩阵或 CIRCexplorer/CIRI 结果表：识别 `circ_id/junction/gene` 列 + 样本列 / 或 `circRNA` + counts |
| `mirna` | `.tsv/.csv/.txt` | miRNA 表达矩阵：miRBase 名（hsa-miR-…）或 mature 名；precursor→mature 映射表内置 |
| `lncrna` | 矩阵 | lncRNA 表达矩阵：ENSG/GENCODE id 或 symbol（`LINC…`/`…-AS1` 启发式识别） |
| `fusion` | `.tsv/.csv` | STAR-Fusion / Arriba / SOAPfuse 输出：识别常见列（`#FusionName`/`FusionName`/`gene1`+`gene2`/`JunctionReadCount`/`SpanningFragCount`） |
| `splicing` | `.tsv/.csv` | rMATS（`SE.MATS.JC.txt` 等，含 `IJC_SAMPLE_1/SJC_SAMPLE_1/IncLevel1`）或 SUPPA PSI 矩阵或通用 PSI 表 |

- 返回 `{"format","X","genes","samples","sample_info","notes","warnings","meta":{...}}`。
- **表达值类型判定**：整数且无负值 → `counts`；有小数且 max<50 → `lognorm`；max>100 且全正 → `tpm/fpkm`。写进 `dataset.notes`。
- 样本信息表 join：按行名或 `key` 列大小写不敏感匹配；匹配不上的样本 → warning。
- 自有数据：`open` 支持 `sample_sheet` / `clinical` 两个额外路径，导入后自动 join。

### 7.2 GEO 下载
- `dest_dir` 由用户指定；默认 `<workdir>/downloads/GEO`（UI 显示可改）。
- 文件名取 URL basename；已存在且 `overwrite=false` 跳过。
- 下载后自动 `inspect` 并给出"建议导入格式/建议系列矩阵"。

---

## 8. 步骤与分析详细要求（按 owner）

### 8.1 bulk-core（`bg_io.py` / `bg_geo.py` / `bg_prep.py`）
**步骤**：`qc`、`filter`、`normalize`、`batch`、`pca` 全实现。
- qc：每样本 `lib_size`（列和）、`detect_genes`（>0 计数）、`median_expr`、`pct_top_gene`、离群标记（`|z|>3` 任一指标）；图 `qc_box`、`qc_density`。
- filter：`min_expr` + `min_samples_frac` 过滤基因；`remove_outliers` 时剔除 qc 离群样本并写回 dataset（samples 同步）。
- normalize：`method=auto` → 若 counts：`log2(CPM+1)`；若 tpm/fpkm：`log2(x+1)`；若已 lognorm：不动（notes 说明）；`quantile` 用 `preprocessCore::normalize.quantiles`（R）或 numpy 分位数对齐（Python 回退）；`tmm` 用 edgeR::calcNormFactors（R）或 TMM 自实现回退。
- batch：`combat` → R `sva` 不可用，用 **`limma::removeBatchEffect`（R）或自实现 ComBat（Python）**；`removeBatchEffect` 同名方法走 R；输出校正前后 PCA 对比图 + 批次/分组混杂指标（如批次与分组 ARI/Cramer's V）。
- pca：SVD/PCA（sklearn），碎石图 + 样本散点在 PC1/PC2（按 group/batch 着色）+ 样本相关性热图 + 层次聚类树（scipy）。

**基础分析（≥8）**：`sample_correlation`、`sample_dendro`、`gene_boxplot`、`gene_heatmap`、`expr_distribution`、`batch_diagnose`、`outlier_detect`、`pca_samples`。

**GEO**：`geo_search` / `geo_detail` / `geo_download` 按 §2.2；中文词典 ≥120 条（疾病/组织/处理/技术/物种）；`geo_detail` 解析 suppl 目录 HTML。

### 8.2 bulk-de（`bg_de.py` / `bg_enrich.py` / `bg_genesets.py`）
**步骤** `de`：
- `method=limma`（默认，R `limma::lmFit + eBayes + topTable`，支持 `covars` 协变量与 `paired`）→ 回退 Welch t + BH（标注 lite）。
- `method=deseq2`：R `DESeq2`（需 counts；否则 prerequisite）。
- `method=welch`：纯 Python。
- 输出列统一：`gene, logFC, AveExpr, t, P.Value, adj.P.Val, B, regulation(up/down/ns)`；`data_out` 写 `de_all.csv` + `de_sig.csv`。
- 图：`volcano`、`ma_plot`、`de_heatmap`（top 50 差异基因 z-score 热图）。

**步骤** `enrich`：
- ORA：R `clusterProfiler::enrichGO/enrichKEGG/enrichPathway`（msigdbr/org.Hs.eg.db）；回退纯 Python 超几何 + 内置 GMT。
- GSEA：R `fgsea::fgsea` 或 `clusterProfiler::GSEA`；回退加权 KS + 置换。
- 输出：`ID, Description, GeneRatio/BgRatio, pvalue, p.adjust, Count, geneID`；图 `enrich_bar`、`enrich_dot`。
- 在线库：Enrichr API（`https://maayanlab.cloud/Enrichr/`）可作为额外库来源，超时 20s，失败降级；库缓存到 `<workdir>/genesets/`。

**分析（≥12）**：`ora_go`、`ora_kegg`、`ora_msigdb`、`gsea_go`、`gsea_kegg`、`gsea_hallmark`、`gsva_pathway`、`ssgsea_pathway`、`pathway_heatmap`、`de_multi`、`de_venn`、`de_summary`、`wgcna_module`、`tf_activity`。

### 8.3 bulk-clin（`bg_immune.py` / `bg_survival.py` / `bg_prognosis.py`）
**步骤** `immune`：**R `IOBR::deconvo_tme` 优先**，方法 `cibersort, epic, mcpcounter, xcell, quantiseq, estimate, ssgsea, timer`（每个方法单独 tryCatch，失败则跳过并记 warning，至少 2 个成功）；Python 回退：内置 ssGSEA（28 类免疫细胞，Charoentong 基因集）+ ESTIMATE 三分数（内置 141/141 基因集可简化）+ MCP-counter 式 marker 均值。输出：样本×细胞比例表、方法间相关性、组间差异（Mann-Whitney + FDR）、堆叠图、相关性热图。
**步骤** `survival`：KM + logrank（lifelines）+ Cox 单/多因素 + 最佳 cutoff 扫描 + 时间 ROC。临床列自动识别（`time|os|dfs|rfs|pfs|survival`、`event|status|dead|death`）。
**步骤** `prognosis`：Cox/LASSO-Cox（lifelines CoxPHFitter penalizer 或 R glmnet）/随机生存森林（sksurv 不可用 → 用 Cox + 自助法 或 sklearn RandomForestSurvival 自实现，诚实标注）；输出风险评分、风险分组 KM、时间 ROC（`timeROC` R 或 lifelines）、C-index、校准曲线、DCA、森林图、列线图数据表。
**分析（≥14）**：免疫 8 方法各一项 + `immune_all`、`immune_correlation`、`immune_group_diff`、`tme_cluster`、`km_curve`、`logrank_multi`、`cox_uni`、`cox_multi`、`cutoff_scan`、`time_roc`、`subgroup_forest`、`risk_score_cox`、`risk_score_lasso`、`risk_heatmap`、`cindex`、`calibration`、`dca`、`nomogram`、`clinical_table`。

### 8.4 bulk-ncrna（`bg_scmap.py` / `bg_ncrna.py`）
**步骤** `scmap`：
- 参考来源：`ref_path`（h5ad，anndata 读）；或 Host 已 attach 的单细胞数据集。
- 用 `gene_set`（bulk DE 上调/下调/自定义基因列表）在单细胞里做打分：`mean`（亚群平均表达 z）、`aucell`（rank-AUC）、`addmodule`（对照基因集校正）；输出：细胞类型/簇 × 得分表、目标分子在各亚群的平均表达/表达比例表、`sc_gene_violin`（可选）、与 bulk 样本表达的相关热图（celltype × bulk sample）。
- 不做细胞类型注释篡改：只读 reference 的 `cell_type`/`clusters` obs 列。
**步骤** `ncrna`：按 `kind` 分发到 5 个模态（见 §8.5）。
**分析（≥12）**：`sc_signature_score`、`sc_de_map`、`sc_celltype_prop`、`sc_correlation`、`sc_marker_consistency`、`sc_gene_violin`、`circ_de`、`circ_ratio`、`mirna_de`、`mirna_target`、`lncrna_de`、`lncrna_class`、`lncrna_network`、`cerna_network`、`fusion_summary`、`fusion_filter`、`splicing_dpsi`、`splicing_volcano`、`splicing_event`。

### 8.5 非编码模态规则
- **circRNA**：ID 形如 `chr1:12345-12789:+` / `hsa_circ_0001234`。差异分析同 bulk DE（用 modular 函数复用 `bg_de` 的统计核心，避免重复实现）；`circ_ratio` 计算 back-spliced/total 比例（若表里有 linear 列）；`circos` 图按染色体位置（解析 ID）画圈图（matplotlib 极坐标）。
- **miRNA**：名称规范化（`hsa-miR-21-5p`）；`mirna_target` 用 **相关性（miRNA-mRNA 负相关）+ 内置种子序列表**（若给 FASTA/序列则种子匹配，否则只做相关性并标注 lite）；输出 miRNA-target 表 + ceRNA 三元组候选。
- **lncRNA**：`lncrna_class` 按 ID 前缀/长度/biotype 注释（若用户给 GTF 则精确）；`lncrna_network` 共表达（|r|>0.6, FDR<0.05）+ 网络图。
- **fusion**：识别列 → 标准化表（gene1/gene2/junction_reads/spanning_reads/breakpoint/frame/oncogene）；`fusion_summary` 条形图（top 融合对）+ 样本×融合对矩阵；`fusion_filter` 按癌基因列表（内置 ≥120 个）+ 激酶列表（≥90）过滤。
- **splicing**：rMATS 表解析（`IJC_SAMPLE_1` 等）+ `IncLevel1/2` → 组间 dPSI + Welch/t 检验 + BH；通用 PSI 表 → 直接 dPSI；输出事件表（SE/A5SS/A3SS/MXE/RI）+ 火山图 + 事件类型条形图 + top 事件热图。

---

## 9. 与 dsh-singlecell 的联动（Lead 实现 HTTP 侧，bulk-ncrna 使用）
- `GET /geobulk/sc/refs` 先探测本机常见 h5ad（`<workspace>/plugins/dsh-singlecell/exports/*.h5ad`、`<cwd>/*.h5ad` 一层），再尝试 `GET <origin>/singlecell/state`；返回 `{"available","datasets":[{"id","name","path","n_obs"}],"local":[{"path","name","bytes"}]}`。
- `POST /geobulk/sc/attach {dataset, sc_dataset?|ref_path?}`：若给 `sc_dataset` → 调 singlecell `POST /singlecell/export {dataset, what:"h5ad", path:<workdir>/sc-ref/<id>.h5ad}`；否则用 `ref_path`。把路径写进 bulk dataset 的 `sc_ref` 字段。
- singlecell 未安装/未运行：`/sc/refs` 返回 `available:false` + 本地 h5ad 列表，不报 500。
- 不修改 singlecell 的任何文件。

---

## 10. 环境事实（所有人必读）
- 开发环境：macOS；插件目录 `plugins/dsh-geobulk`（相对仓库根）；Python 3.14 / R 4.4 的实测环境见 §10。
- Python `python3` = `/Library/Frameworks/Python.framework/Versions/3.14/bin/python3`（3.14.2）：numpy 2.5.1 / pandas 2.3.3 / scipy 1.17.1 / sklearn 1.9.0 / matplotlib 3.11.0 / anndata 0.12.19 / h5py 3.16 / lifelines 0.30.3 / statsmodels 0.14.6 / openpyxl 3.1.5 / requests 2.34.2 / seaborn / networkx / PIL / bs4 / lxml。**无 gseapy、无 umap-learn、无 sva、无 torch、无 pyarrow。**
- `~/.matplotlib` 不可写：必须 `MPLBACKEND=Agg` + `MPLCONFIGDIR` 指到 workdir（Host 已设，模块里再由 `bg_common` 兜底）。
- R `/Library/Frameworks/R.framework/Resources/bin/Rscript` 4.4.2：limma/edgeR/DESeq2/survival/clusterProfiler/DOSE/fgsea/msigdbr/GSVA/WGCNA/preprocessCore/ComplexHeatmap/pheatmap/timeROC/rms/glmnet/survminer/org.Hs.eg.db/org.Mm.eg.db/**IOBR 0.99.0**/GSEABase/Biobase/AnnotationDbi 可用；sva/estimate/immunedeconv/ReactomePA/forestplot 不可用。
- 网络：NCBI E-utilities/FTP、Enrichr、KEGG REST 均可达（已验证）。
- 宿主 Electron 进程 PATH 窄：探测 Python/R 用候选列表逐个试。
- **禁止** pip/conda/CRAN 安装、禁止改 profile、禁止 pnpm、禁止写 `$DSH_HOME`。
- 插件安装由 Lead 用 `plugin_manager install_bundle` 完成；开发者只写包内文件并跑自己的测试。
- 提交前自检：`python3 -m py_compile lib/*.py`；`node --check client.js`；`node -e "import('./lib/impl.mjs')"`。

---

## 11. 测试要求（每 owner 自己跑通再交）
- **bulk-core**：`test/py_import_test.py` + `test/make_fixtures.py`：生成 12 种 fixture（csv/tsv/xlsx/gz/geo series matrix/zip/tar/rds/RData/matrix_text 转置/miRNA/circ）并断言 `genes/samples/X` 与 truth 一致；`inspect` 不加载整表；坏文件给可读错误不崩；步骤 qc/filter/normalize/batch/pca 在 demo 上跑通且数值合理（过滤后基因数下降、标准化后列中位数接近、ComBat 前后批次分离度下降）。`test/py_geo_test.py`：mock HTTP 返回固定 JSON/HTML，断言 search/detail/download 解析正确（不依赖真实网络）。
- **bulk-de**：`test/py_de_test.py`：demo 注入的 120 个差异基因 top200 命中率 ≥ 80%（limma）；方法/协变量/配对分支可用；前置缺失给 `prerequisite`。`test/py_enrich_test.py`：富集结果 p.adjust 单调、基因集覆盖、GSEA ES∈[-1,1]、回退实现可用。
- **bulk-clin**：`test/py_immune_test.py`：≥2 个方法成功、比例和为 1±0.1（相对法）、组间差异有 p 值。`test/py_survival_test.py`：KM p 值域 [0,1]、HR/C-index 合理、cutoff 扫描返回单调统计、缺列给 prerequisite。
- **bulk-ncrna**：`test/py_scmap_test.py`：用 `test/fixtures/demo_sc.h5ad`（自建）跑通打分，得分 ∈ [0,1]，目标基因命中亚群正确率 ≥ 60%。`test/py_ncrna_test.py`：5 种模态各造小 fixture 跑通，rMATS 解析正确、fusion 列识别正确、circ ID 解析正确。
- **bulk-ui**：`test/client-logic.test.mjs`：纯函数单测（STEPS/格式化/坐标变换/火山图点/搜索 payload/API 客户端/建议指令/任务轮询）。
- 所有测试**不需要 node/npm 依赖**（python 用标准库 + 已装包；node 用内置 `node:test` 或纯脚本 `process.exit`）。

---

## 12. 交付验收（Lead 负责最终集成）
1. `plugin_manager install_bundle` 成功，`/geobulk/health` 200 且 `engine.ready=true`。
2. demo 全流程：qc→filter→normalize→batch→pca→de→enrich→immune→survival→prognosis→scmap→ncrna 每步成功且有 summary/图。
3. 真实 GEO 检索一次 + 下载一个 series matrix 成功（或网络受限时用 mock 证明链路）。
4. 自有数据：用工作区里真实 bulk 表达表（若无可造）导入 → 清洗 → DE。
5. UI：侧边栏出现「Bulk 工作台」，整页可打开，主要 tab 不白屏，GEO 检索/导入/流程卡/分析可交互。
6. 单细胞联动：从工作台拉取或指定 h5ad 出细胞类型打分。

---

## §14 实现期 API 增补（Lead 维护，与代码同步）

以下为实现时新增/明确的接口，**追加不改既有字段**：

1. `GET /geobulk/state` 额外返回：`steps`（12 个步骤目录，含 `params` 与 `available`）、`engine.modules`（各模块加载状态）、`sc`（单细胞工作台联动状态）。
2. `GET /geobulk/health` 的 `engine` 额外含 `r`（`{available,version,script,packages}`）、`modules`。
3. 新增 `POST /geobulk/recipe/run` `{dataset,recipe}`：Host 顺序执行预设流程并返回每步结果（工具里也走这个）。
4. 新增 `GET /geobulk/jobs`（列表）与 `GET /geobulk/jobs?id=`（单任务）；`POST /geobulk/geo/download` 立即返回 `{job}`，前端轮询。
5. 单细胞联动路由固定为 `GET /geobulk/sc/refs`、`POST /geobulk/sc/attach {dataset, ref_path?|sc_dataset?}`；`sc/attach` 会调 dsh-singlecell 的 `POST /singlecell/export {dataset,what:"h5ad",path}` 拿到 h5ad 后写入引擎 `sc_set`。
6. 步骤/分析返回信封允许的顶层附加键（registry 消费后**不发给前端**）：`"replace"`（回写 X/genes/samples/counts/value_kind）、`"de_store"`/`"store"`/`"result_store"`（写入 `dataset.results`）、`"sample_columns"`+`"sample_values"`（写回 sample_info）。
7. `ctx` 额外键：`sc_ref`（单细胞参考 h5ad 路径）、`notes`、`warnings`；`ctx["modules"][kind] = {"kind","X","genes","samples","table","info","path","format","n_features"}`。
8. 工具 `geobulk` 的 action 集合最终为：`status/search/detail/download/open/demo/inspect/preview/samples/run/plot/analyses/analysis/attach/export/audit/suggest/recipe/sc_refs/sc_attach/formats/reset`。
9. 引擎新增命令：`select`（切换 active 数据集）、`reset`（清空）、`steps`（步骤目录）、`formats`（格式清单）、`sc_refs {path?}`（本机 h5ad 列表或单个 h5ad 头信息）、`sc_set {dataset,ref_path}`。
10. 图路由 `GET /geobulk/figure?id=<basename>.png`：按 `analyses/**`、`figures/**`、`exports/**` 递归查找同名 PNG（只需传文件名）。

11. GEO 中文检索：`geo_search` 会把中文关键词用 `ZH_DICT` 反查成英文（返回 `query_zh/query_en/translated`），并在结果过少时**逐级放宽**关键词（去掉 expression profiling 类泛词 → 前 3 词 → 前 2 词，返回 `relaxed/used_term/relax_steps`）；完全无结果时给 `not_found` 可读错误。
12. `save_fig` 在 axes 已有标题时不再追加 suptitle（避免图上出现两行标题）。
