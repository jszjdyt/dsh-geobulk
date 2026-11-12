/**
 * dsh-geobulk — 浏览器端 UI（Client 半边）
 *
 * 形态：window.__ModuleLoader__.load({ id:'dsh-geobulk', factory(require){...} })，
 * 只 require('react')，不碰任何 @deepseek-ai Client 包（宿主版本一变就会炸）。
 *
 * 数据面：全部走同源相对路径 /geobulk/*（见 CONTRACT.md §4 / §6）。
 *
 * 文件分两层：
 *   1) 纯函数层（本段之前到 BG_TEST）：STEPS 元数据 / 数值与时间格式化 /
 *      分组取色 / 火山图点 / PCA 视口变换 / GEO 搜索 payload / 建议→中文指令 /
 *      可注入 fetch 的 API 客户端 / 异步任务轮询状态机 / 服务端响应归一化。
 *      不碰 DOM、不 require react；node 环境（typeof window === 'undefined'）
 *      挂到 globalThis.__BG_TEST__ 供 test/client-logic.test.mjs 直接 import。
 *   2) 浏览器层（factory 内）：React 组件 + 样式 + 槽注册。
 *
 * 约定：类名一律 bg- 前缀；颜色只用 --dsw-alias-* 主题变量（见 CSS 段），
 * 数据色（分组/细胞/火山图点）用 JS 常量 PALETTE / CONTINUOUS_STOPS。
 */
;(function () {
  'use strict'

  var NS = 'dsh-geobulk'
  var API = '/geobulk'
  var POLL_MS = 1500
  var POLL_TIMEOUT_MS = 20 * 60 * 1000
  var DEFAULT_API_TIMEOUT_MS = 60000
  var HEAVY_API_TIMEOUT_MS = 300000
  var MAX_HISTORY_ROWS = 20
  var MAX_TABLE_RENDER_ROWS = 200
  var DEFAULT_GEO_DIR = '~/.dsh/geobulk/GEO'

  /** 分组/细胞类型调色板（数据色，不参与主题）。 */
  var PALETTE = [
    '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948',
    '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac', '#1f77b4', '#ff7f0e',
    '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf',
  ]

  /** 连续值色带（数据色）：低→高。 */
  var CONTINUOUS_STOPS = ['#3b4cc0', '#6a9fdb', '#9ec6e8', '#e8d4a8', '#f0a35e', '#c0392b']

  /** 火山图三色（数据色）。 */
  var VOLCANO_COLORS = { up: '#d1495b', down: '#2e6f95', ns: '#b8bcc4' }

  /** 自检等级样式（由主题色派生）。 */
  var AUDIT_LEVEL_STYLE = {
    ok: { tone: 'ok', labelKey: 'auditOk' },
    warn: { tone: 'warn', labelKey: 'auditWarn' },
    fail: { tone: 'fail', labelKey: 'auditFail' },
  }

  /** 结果表格中可排序/展示的常见列（DE 表）。 */
  var DE_COLUMNS = ['gene', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B', 'regulation']

  /** 富集数据库选项（与 CONTRACT §2.6 enrich 的 db 取值一致）。 */
  var ENRICH_DBS = [
    { value: 'GO_BP', label: 'GO BP' },
    { value: 'GO_CC', label: 'GO CC' },
    { value: 'GO_MF', label: 'GO MF' },
    { value: 'KEGG', label: 'KEGG' },
    { value: 'HALLMARK', label: 'Hallmark' },
    { value: 'REACTOME', label: 'Reactome' },
  ]

  /** 免疫浸润方法（与 CONTRACT §2.6 immune 的 methods 取值一致）。 */
  var IMMUNE_METHODS = [
    { value: 'cibersort', label: 'CIBERSORT' },
    { value: 'epic', label: 'EPIC' },
    { value: 'mcpcounter', label: 'MCP-counter' },
    { value: 'xcell', label: 'xCell' },
    { value: 'quantiseq', label: 'quanTIseq' },
    { value: 'estimate', label: 'ESTIMATE' },
    { value: 'ssgsea', label: 'ssGSEA' },
    { value: 'timer', label: 'TIMER' },
  ]

  /** 非编码模态（与 CONTRACT §2.6 ncrna 的 kind 取值一致）。 */
  var NCRNA_MODULES = [
    { value: 'circ', labelKey: 'moduleCirc' },
    { value: 'mirna', labelKey: 'moduleMirna' },
    { value: 'lncrna', labelKey: 'moduleLncrna' },
    { value: 'fusion', labelKey: 'moduleFusion' },
    { value: 'splicing', labelKey: 'moduleSplicing' },
  ]

  /** 记录 tab 的导出按钮组（what 与 CONTRACT §4 / §2.2 export 一致）。 */
  var EXPORT_KINDS = [
    { id: 'expr_csv', labelKey: 'exportExprCsv', hint: 'CSV' },
    { id: 'de_csv', labelKey: 'exportDeCsv', hint: 'CSV' },
    { id: 'enrich_csv', labelKey: 'exportEnrichCsv', hint: 'CSV' },
    { id: 'immune_csv', labelKey: 'exportImmuneCsv', hint: 'CSV' },
    { id: 'sample_info_csv', labelKey: 'exportSampleCsv', hint: 'CSV' },
    { id: 'summary_json', labelKey: 'exportSummary', hint: 'JSON' },
    { id: 'protocol', labelKey: 'exportProtocol', hint: 'JSON' },
    { id: 'methods', labelKey: 'exportMethods', hint: 'Markdown' },
    { id: 'report', labelKey: 'exportReport', hint: 'Markdown' },
    { id: 'script', labelKey: 'exportScript', hint: 'Python' },
  ]

  /** 记录 tab 上方快速导出（只放最常用的三个）。 */
  var QUICK_EXPORTS = ['expr_csv', 'de_csv', 'summary_json']

  // ==========================================================================
  // 流程步骤：顺序即左栏卡片顺序；requires 与 CONTRACT §2.6 逐字一致。
  // params 用通用渲染（number/select/boolean/column/columns/text/genes/multiselect）。
  // ==========================================================================
  var STEPS = [
    {
      id: 'qc',
      title: '样本质控',
      icon: 'qc',
      requires: [],
      params: [
        { key: 'group_col', label: '分组列', type: 'column', default: null },
        { key: 'min_detect_genes', label: '最少检出基因数', type: 'number', default: 0, min: 0, max: 30000, step: 100 },
      ],
    },
    {
      id: 'filter',
      title: '过滤',
      icon: 'filter',
      requires: ['qc'],
      params: [
        { key: 'min_expr', label: '最小表达量', type: 'number', default: 1.0, min: 0, step: 0.5 },
        { key: 'min_samples_frac', label: '最少样本比例', type: 'number', default: 0.2, min: 0, max: 1, step: 0.05 },
        { key: 'remove_outliers', label: '剔除离群样本', type: 'boolean', default: false },
      ],
    },
    {
      id: 'normalize',
      title: '标准化',
      icon: 'normalize',
      requires: ['filter'],
      params: [
        {
          key: 'method',
          label: '标准化方法',
          type: 'select',
          default: 'auto',
          options: [
            { value: 'auto', label: '自动（counts→log2 CPM，TPM→log2）' },
            { value: 'quantile', label: '分位数归一化' },
            { value: 'median', label: '中位数标准化' },
            { value: 'tmm', label: 'TMM' },
            { value: 'none', label: '不标准化' },
          ],
        },
        {
          key: 'log_transform',
          label: 'log 转换',
          type: 'select',
          default: 'auto',
          options: [
            { value: 'auto', label: '自动' },
            { value: 'true', label: '是' },
            { value: 'false', label: '否' },
          ],
        },
      ],
    },
    {
      id: 'batch',
      title: '去批次',
      icon: 'batch',
      requires: ['normalize'],
      params: [
        { key: 'batch_key', label: '批次列', type: 'column', default: null },
        {
          key: 'method',
          label: '校正方法',
          type: 'select',
          default: 'combat',
          options: [
            { value: 'combat', label: 'ComBat / removeBatchEffect' },
            { value: 'removeBatchEffect', label: 'limma removeBatchEffect' },
            { value: 'none', label: '不校正' },
          ],
        },
        { key: 'group_key', label: '分组列（保护变量）', type: 'column', default: null },
      ],
    },
    {
      id: 'pca',
      title: '降维聚类',
      icon: 'pca',
      requires: ['normalize'],
      params: [
        { key: 'n_pcs', label: '主成分数', type: 'number', default: 20, min: 2, max: 200, step: 1 },
        { key: 'color_by', label: '着色列', type: 'column', default: 'group' },
        {
          key: 'cluster_method',
          label: '聚类方法',
          type: 'select',
          default: 'hclust',
          options: [
            { value: 'hclust', label: '层次聚类 hclust' },
            { value: 'kmeans', label: 'K-means' },
            { value: 'none', label: '不聚类' },
          ],
        },
        { key: 'k', label: '聚类数 k（auto=自动）', type: 'number', default: null, min: 2, max: 20, step: 1 },
      ],
    },
    {
      id: 'de',
      title: '差异表达',
      icon: 'de',
      requires: ['normalize'],
      params: [
        { key: 'group_col', label: '分组列', type: 'column', default: 'group' },
        { key: 'case', label: '实验组', type: 'text', default: '' },
        { key: 'control', label: '对照组', type: 'text', default: '' },
        {
          key: 'method',
          label: '方法',
          type: 'select',
          default: 'limma',
          options: [
            { value: 'limma', label: 'limma（R 优先）' },
            { value: 'deseq2', label: 'DESeq2（需 counts）' },
            { value: 'welch', label: 'Welch t 检验' },
          ],
        },
        { key: 'covars', label: '协变量（逗号分隔）', type: 'columns', default: [] },
        { key: 'lfc', label: 'logFC 阈值', type: 'number', default: 0, min: 0, step: 0.1 },
        { key: 'padj', label: 'FDR 阈值', type: 'number', default: 0.05, min: 0, max: 1, step: 0.01 },
        { key: 'paired', label: '配对设计', type: 'boolean', default: false },
      ],
    },
    {
      id: 'enrich',
      title: '功能富集',
      icon: 'enrich',
      requires: ['de'],
      params: [
        { key: 'db', label: '数据库', type: 'multiselect', default: ['GO_BP', 'KEGG'], options: ENRICH_DBS },
        {
          key: 'method',
          label: '方法',
          type: 'select',
          default: 'ora',
          options: [
            { value: 'ora', label: 'ORA 过表达分析' },
            { value: 'gsea', label: 'GSEA 排序富集' },
          ],
        },
        { key: 'pvalue', label: 'P 值阈值', type: 'number', default: 0.05, min: 0, max: 1, step: 0.01 },
        { key: 'qvalue', label: 'q 值阈值', type: 'number', default: 0.2, min: 0, max: 1, step: 0.01 },
        { key: 'top_n', label: '展示条数', type: 'number', default: 20, min: 1, max: 100, step: 1 },
      ],
    },
    {
      id: 'immune',
      title: '免疫浸润',
      icon: 'immune',
      requires: ['normalize'],
      params: [
        { key: 'methods', label: '方法', type: 'multiselect', default: ['cibersort', 'epic', 'mcpcounter', 'xcell', 'quantiseq', 'estimate', 'ssgsea', 'timer'], options: IMMUNE_METHODS },
        { key: 'group_col', label: '分组列', type: 'column', default: 'group' },
      ],
    },
    {
      id: 'survival',
      title: '临床生存',
      icon: 'survival',
      requires: ['normalize'],
      params: [
        { key: 'time_col', label: '时间列', type: 'column', default: null },
        { key: 'event_col', label: '事件列', type: 'column', default: null },
        {
          key: 'endpoint',
          label: '终点',
          type: 'select',
          default: 'OS',
          options: [
            { value: 'OS', label: 'OS 总生存' },
            { value: 'DFS', label: 'DFS 无病生存' },
            { value: 'RFS', label: 'RFS 无复发生存' },
            { value: 'PFS', label: 'PFS 无进展生存' },
          ],
        },
        { key: 'gene', label: '目标基因（可选）', type: 'text', default: '' },
        {
          key: 'cutoff',
          label: 'cutoff 选择',
          type: 'select',
          default: 'median',
          options: [
            { value: 'median', label: '中位数' },
            { value: 'optimal', label: '最佳 cutoff' },
            { value: 'quartile', label: '四分位' },
          ],
        },
        { key: 'covariates', label: '协变量（逗号分隔）', type: 'columns', default: [] },
      ],
    },
    {
      id: 'prognosis',
      title: '高级预后',
      icon: 'prognosis',
      requires: ['survival'],
      params: [
        {
          key: 'model',
          label: '模型',
          type: 'select',
          default: 'cox',
          options: [
            { value: 'cox', label: 'Cox 比例风险' },
            { value: 'lasso', label: 'LASSO-Cox' },
            { value: 'rsf', label: '随机生存森林' },
          ],
        },
        { key: 'genes', label: '基因清单（留空取 DE top）', type: 'genes', default: [] },
        { key: 'time_col', label: '时间列', type: 'column', default: null },
        { key: 'event_col', label: '事件列', type: 'column', default: null },
        { key: 'nfold', label: '交叉验证折数', type: 'number', default: 10, min: 2, max: 20, step: 1 },
      ],
    },
    {
      id: 'scmap',
      title: '单细胞映射',
      icon: 'scmap',
      requires: ['de'],
      params: [
        { key: 'ref_path', label: 'h5ad 参考路径', type: 'text', default: '' },
        {
          key: 'gene_set',
          label: '签名来源',
          type: 'select',
          default: 'de_up',
          options: [
            { value: 'de_up', label: 'DE 上调基因' },
            { value: 'de_down', label: 'DE 下调基因' },
            { value: 'custom', label: '自定义基因' },
          ],
        },
        { key: 'genes', label: '自定义基因', type: 'genes', default: [] },
        { key: 'cell_type_col', label: '细胞类型列', type: 'text', default: 'cell_type' },
        {
          key: 'score_method',
          label: '打分方法',
          type: 'select',
          default: 'mean',
          options: [
            { value: 'mean', label: '均值 z-score' },
            { value: 'aucell', label: 'AUCell rank-AUC' },
            { value: 'addmodule', label: 'AddModule 校正' },
          ],
        },
      ],
    },
    {
      id: 'ncrna',
      title: '非编码/融合/剪接',
      icon: 'ncrna',
      requires: ['normalize'],
      params: [
        { key: 'kind', label: '模态', type: 'select', default: 'circ', options: NCRNA_MODULES.map(function (item) { return { value: item.value, label: item.value.toUpperCase() } }) },
        { key: 'module_path', label: '模块文件路径', type: 'text', default: '' },
      ],
    },
  ]

  var STEP_BY_ID = {}
  for (var stepIndex = 0; stepIndex < STEPS.length; stepIndex++) STEP_BY_ID[STEPS[stepIndex].id] = STEPS[stepIndex]

  /** 左栏序号/标题查询。 */
  function stepById(id) {
    return STEP_BY_ID[id] || null
  }

  /** 中栏 tab：顺序与 CONTRACT §6.2 一致。 */
  var TAB_KEYS = ['overview', 'expr', 'de', 'enrich', 'immune', 'survival', 'prognosis', 'sc', 'ncrna', 'records', 'log']

  var TAB_META = {
    overview: { titleKey: 'tabOverview' },
    expr: { titleKey: 'tabExpr' },
    de: { titleKey: 'tabDe' },
    enrich: { titleKey: 'tabEnrich' },
    immune: { titleKey: 'tabImmune' },
    survival: { titleKey: 'tabSurvival' },
    prognosis: { titleKey: 'tabPrognosis' },
    sc: { titleKey: 'tabSc' },
    ncrna: { titleKey: 'tabNcrna' },
    records: { titleKey: 'tabRecords' },
    log: { titleKey: 'tabLog' },
  }

  // ==========================================================================
  // 文案字典（zh 完整；en 与 zh 同 key 集）。apply() 里 locale.register 注册。
  // 参数 label 等来自服务端 schema 的文案属于数据，不走字典。
  // ==========================================================================
  var DICT_ZH = {
    nav: 'Bulk 工作台',
    appTitle: 'GEO · Bulk 转录组工作台',
    appSubtitle: 'GEO 检索 · 表达矩阵清洗 · 差异/富集/免疫/生存 · 单细胞映射',
    refresh: '刷新',
    refreshHint: '重新读取工作台状态',
    retry: '重试',
    loading: '加载中…',
    loadingData: '正在读取工作台…',
    noData: '尚未载入数据集',
    errorTitle: '出错了',
    crashTitle: '这个面板渲染失败',
    crashHint: '其它面板仍可使用，可以点击重试恢复。',
    close: '关闭',
    cancel: '取消',
    copy: '复制',
    copied: '已复制',
    copyFailed: '复制失败（浏览器不允许写入剪贴板）',
    search: '搜索',
    all: '全部',
    none: '无',
    yes: '是',
    no: '否',
    unknown: '未知',
    dataset: '数据集',
    samples: '样本',
    genes: '基因',
    groups: '分组',
    batches: '批次',
    engine: '引擎',
    modules: '模块',
    sourceFormat: '来源格式',
    species: '物种',
    warnings: '警告',
    notes: '说明',
    current: '当前',
    completed: '已完成',
    pending: '待运行',
    running: '运行中',
    failed: '失败',
    plot: '出图',
    plotting: '出图中…',
    noFigure: '暂无图',
    figureFail: '图片加载失败（文件可能已被清理）',
    exportBtn: '导出',
    exported: '已导出',
    exportFail: '导出失败',
    copyPath: '复制路径',
    pathLabel: '文件路径',
    tableRows: '共 {v} 行',
    tableTruncated: '仅展示前 {v} 行',
    emptyTable: '空表',
    total: '总计',
    selectAll: '全选',
    clearSelection: '清空',
    of: '/',
    selectDataset: '切换数据集',
    closeDataset: '关闭数据集',
    workdir: '工作目录',
    engineReady: '引擎就绪',
    engineBooting: '引擎启动中',
    engineBad: '引擎不可用',
    engineGuideTitle: 'Python 引擎尚未就绪',
    engineGuideBody: '工作台依赖本机 Python 环境（numpy / pandas / scipy / sklearn / matplotlib / anndata / lifelines）。整页仍可浏览，但计算与出图需要引擎先就绪。',
    engineGuideSteps: '排查顺序：① 查看工作目录下 engine.log 的报错；② 确认插件设置里的 python3 路径可执行；③ 重启插件或宿主后点「重试」。',
    pythonLabel: 'Python',
    rLabel: 'R',
    rMissing: 'R 未安装',
    rVersion: 'R 版本',
    engineInfo: '引擎信息',
    statusOk: '正常',
    statusWarn: '注意',
    statusFail: '异常',
    pipelineTitle: '分析流程',
    pipelineHint: '按顺序运行，完成后可重跑（重跑会清空下游结果）',
    datasetsTitle: '数据集',
    stepStatusDone: '已完成',
    stepStatusRunning: '运行中',
    stepStatusReady: '待运行',
    stepStatusBlocked: '待前置',
    stepStatusFailed: '失败',
    stepRun: '运行',
    stepRerun: '重跑',
    stepRunning: '运行中…',
    stepParams: '参数',
    stepNeedData: '先导入或打开一个数据集',
    stepNeedPrereq: '需先完成：',
    stepSummaryEmpty: '暂无结果摘要',
    stepResult: '结果',
    datasetNoDataset: '还没有数据集',
    datasetOpenHint: '点「导入数据」打开表达矩阵，或用「GEO 检索」下载公开数据。',
    datasetChips: '数据集',
    datasetStepsDone: '已完成 {v} 步',
    datasetCurrentStep: '当前：{v}',
    warningsLabel: '警告',
    tabOverview: '概览',
    tabExpr: '表达矩阵',
    tabDe: '差异',
    tabEnrich: '富集',
    tabImmune: '免疫',
    tabSurvival: '生存',
    tabPrognosis: '预后',
    tabSc: '单细胞',
    tabNcrna: '非编码',
    tabRecords: '记录',
    tabLog: '日志',
    overviewTitle: '数据集概览',
    metricSamples: '样本数',
    metricGenes: '基因数',
    metricGroups: '分组',
    metricBatches: '批次',
    metricEngine: '引擎',
    metricModules: '非编码模块',
    pcaPlot: 'PCA 样本散点',
    corrPlot: '样本相关性',
    geneQuickTitle: '快速基因箱线图',
    geneQuickHint: '输入基因名（逗号分隔），点「出图」查看表达分布。',
    geneInputPlaceholder: '如：CDH1, VIM, MKI67',
    nextSteps: '建议下一步',
    noSuggestions: '暂无建议，继续按左栏流程运行即可。',
    suggestionDetail: '建议详情',
    analysisCatalog: '分析目录',
    geoSearch: 'GEO 检索',
    importData: '导入数据',
    demoLoad: '载入内置示例',
    exportExprCsv: '表达矩阵 CSV',
    exportDeCsv: '差异结果 CSV',
    exportEnrichCsv: '富集结果 CSV',
    exportImmuneCsv: '免疫浸润 CSV',
    exportSampleCsv: '样本信息 CSV',
    exportSummary: '摘要 JSON',
    exportProtocol: '分析协议',
    exportMethods: '方法描述',
    exportReport: '分析报告',
    exportScript: '复现脚本',
    exprTitle: '表达矩阵',
    exprHint: '搜索基因并添加到已选清单；箱线图/热图走引擎绘图，完整矩阵请导出 CSV。',
    exprGeneSearch: '搜索基因（回车添加）',
    exprAddGene: '添加',
    exprSelected: '已选基因',
    exprExport: '导出表达矩阵',
    exprExportHint: '导出的 CSV 在宿主本机，路径可复制到分析脚本使用。',
    exprNoGene: '请至少选择一个基因',
    exprPlotBox: '基因箱线图',
    exprPlotHeatmap: '基因热图',
    exprGeneNotInData: '该基因不在当前数据集基因列表里',
    exprEmpty: '先在顶部搜索并添加基因。',
    exprTableFallback: '表达矩阵较大时不回传整表：完整数据请用「导出表达矩阵」，表内结果来自分析接口。',
    exprExportedPath: '已导出表达矩阵',
    deTitle: '差异表达',
    deHint: '左栏运行「差异表达」后，这里显示差异表、火山图、MA 图与热图。',
    deEmpty: '还没有差异表达结果。',
    deVolcano: '火山图',
    deMa: 'MA 图',
    deHeatmap: '差异基因热图',
    deSortHint: '点表头可排序',
    deRegulation: '调控方向',
    deUp: '上调',
    deDown: '下调',
    deNs: '不显著',
    deExport: '导出差异表',
    deTableTitle: '差异基因表',
    enrichTitle: '功能富集',
    enrichHint: '先完成差异表达；选择数据库与方法后运行富集。',
    enrichDb: '数据库',
    enrichMethod: '方法',
    enrichP: 'P 值阈值',
    enrichQ: 'q 值阈值',
    enrichTopN: '展示条数',
    enrichRun: '运行富集',
    enrichBar: '富集条形图',
    enrichDot: '富集气泡图',
    enrichEmpty: '还没有富集结果。',
    enrichNoDb: '请至少选择一个数据库',
    immuneTitle: '免疫浸润',
    immuneHint: '方法越多越慢；优先 R IOBR，缺失时自动回退并在结果里标注。',
    immuneMethods: '方法',
    immuneGroupCol: '分组列',
    immuneRun: '运行免疫浸润',
    immuneStack: '细胞比例堆叠图',
    immuneCorr: '方法相关性热图',
    immuneTable: '细胞比例表',
    immuneDiff: '组间差异',
    immuneEmpty: '还没有免疫浸润结果。',
    survivalTitle: '临床生存',
    survivalHint: '需要样本信息中包含时间列与事件列（可在导入时提供临床表）。',
    survivalTimeCol: '时间列',
    survivalEventCol: '事件列',
    survivalEndpoint: '终点',
    survivalGene: '目标基因（可选）',
    survivalCutoff: 'cutoff 选择',
    survivalRun: '运行生存分析',
    survivalKm: 'KM 曲线',
    survivalForest: 'Cox 森林图',
    survivalEmpty: '还没有生存分析结果。',
    survivalNeedColumns: '未检测到时间/事件列，请先导入临床表或检查列名。',
    prognosisTitle: '高级预后',
    prognosisHint: '需先完成临床生存分析。',
    prognosisModel: '模型',
    prognosisGenes: '基因清单（留空取 DE top）',
    prognosisNfold: '交叉验证折数',
    prognosisRun: '运行预后模型',
    prognosisRisk: '风险评分分布图',
    prognosisRoc: '时间依赖 ROC',
    prognosisCindex: 'C-index',
    prognosisEmpty: '还没有预后结果。',
    scTitle: '单细胞映射',
    scHint: '选择 h5ad 参考：本机路径 / 本机已有 h5ad / 从单细胞工作台拉取。',
    scRefPath: 'h5ad 路径',
    scLocalRefs: '本机参考',
    scPull: '从单细胞工作台拉取',
    scAttach: '挂载参考',
    scGeneSet: '签名来源',
    scGenes: '自定义基因',
    scCellTypeCol: '细胞类型列',
    scScoreMethod: '打分方法',
    scRun: '运行单细胞映射',
    scScoreTable: '细胞类型打分表',
    scExprFig: '亚群表达图',
    scEmpty: '还没有单细胞映射结果。',
    scNoRef: '请先选择并挂载 h5ad 参考',
    scAttachOk: '参考已挂载',
    scAttachFail: '挂载参考失败',
    scFromWorkspace: '单细胞工作台数据集',
    scLocal: '本机 h5ad',
    scAvailable: '单细胞工作台在线',
    scUnavailable: '单细胞工作台未运行（仅可用本机 h5ad）',
    scHintNoRefs: '没有探测到可用的 h5ad，请直接输入本机绝对路径。',
    ncrnaTitle: '非编码 / 融合 / 剪接',
    ncrnaHint: '先挂载对应模态的数据文件，再运行该模态分析。',
    ncrnaKind: '模态',
    ncrnaPath: '模块文件路径',
    ncrnaAttach: '挂载模块',
    ncrnaAttached: '已挂载',
    ncrnaRun: '运行非编码分析',
    ncrnaEmpty: '还没有非编码分析结果。',
    ncrnaNotAttached: '该模态尚未挂载数据',
    ncrnaAttachOk: '模块已挂载',
    ncrnaAttachFail: '挂载模块失败',
    moduleCirc: 'circRNA',
    moduleMirna: 'miRNA',
    moduleLncrna: 'lncRNA',
    moduleFusion: '融合基因',
    moduleSplicing: '可变剪接',
    recordsTitle: '记录与自检',
    recordsHistory: '运行历史',
    recordsAudit: '确定性自检',
    recordsAuditScore: '自检得分',
    recordsRunAudit: '运行自检',
    recordsExport: '导出',
    recordsNoHistory: '暂无运行记录',
    recordsAnalysis: '高级分析目录',
    recordsExportOk: '已导出',
    recordsExportFail: '导出失败',
    recordsAuditEmpty: '尚未运行自检',
    recordsExporting: '导出中…',
    auditOk: '通过',
    auditWarn: '注意',
    auditFail: '异常',
    auditHint: '自检来自确定性规则（不调用模型），得分 = 100 - 12×注意 - 30×异常。',
    logTitle: '活动日志',
    logEmpty: '暂无活动记录',
    logRefresh: '刷新日志',
    agentTitle: 'Agent 面板',
    agentStatus: '工作台状态',
    agentSuggestions: '建议下一步',
    agentActivity: '最近活动',
    agentNoActivity: '暂无活动',
    agentNoDataset: '载入数据集后显示状态与建议',
    agentCopy: '复制指令',
    agentCopied: '指令已复制',
    agentCopyFailed: '复制失败',
    geoTitle: 'GEO 数据检索',
    geoSubtitle: '搜索 NCBI GEO Series，中文化字段；选中系列可指定本地目录下载并一键导入。',
    geoTermPh: '胃癌 顺铂 耐药 表达谱',
    geoSpecies: '物种',
    geoAll: '全部',
    geoHuman: '人',
    geoMouse: '小鼠',
    geoStudyType: '研究类型',
    geoStudyAll: '全部',
    geoStudyArray: '芯片表达谱',
    geoStudySeq: '高通量测序',
    geoMinSamples: '最小样本数',
    geoHasClinical: '有临床信息',
    geoSearchBtn: '搜索',
    geoSearching: '检索中…',
    geoResults: '检索结果',
    geoTotal: '共 {v} 条',
    geoNoResults: '没有匹配的 GEO 系列，换个关键词试试。',
    geoDetailBtn: '查看详情',
    geoDetailTitle: '系列详情',
    geoLoadingDetail: '读取系列详情…',
    geoFiles: '可下载文件',
    geoDownloadDir: '下载目录',
    geoDownloadDirHint: '默认 ~/.dsh/geobulk/GEO；可改为任意绝对路径',
    geoDownloadBtn: '下载所选',
    geoDownloading: '下载中…',
    geoImportBtn: '导入并清洗',
    geoImporting: '导入中…',
    geoSelectAll: '全选',
    geoClear: '取消全选',
    geoJobProgress: '下载进度',
    geoJobDone: '下载完成',
    geoJobFailed: '任务失败',
    geoNeedFile: '请至少选择一个文件',
    geoSummary: '摘要',
    geoDesign: '实验设计',
    geoPlatform: '平台',
    geoRelease: '发布日期',
    geoPubmed: 'PubMed',
    geoMatrix: '系列矩阵',
    geoClinicalTag: '临床标签',
    geoImportDone: '已导入，正在启动质控…',
    geoHint: '网络不可用时检索会失败；可改用「导入数据」加载本地文件。',
    geoErrorSearch: 'GEO 检索失败',
    geoErrorDetail: '读取系列详情失败',
    geoYear: '年份',
    geoOpenAfter: '下载后自动导入',
    geoSeriesMatrix: '系列矩阵',
    geoSupp: '补充文件',
    geoRaw: '原始数据',
    geoOther: '其它',
    geoNoFiles: '该系列没有解析到可下载文件',
    importTitle: '导入数据',
    importSubtitle: '支持表达矩阵、GEO series matrix、xlsx、zip/tar、RDS/RData 等',
    importSourceLocal: '本机路径',
    importSourceUrl: 'URL',
    importSourceUpload: '上传文件',
    importPathPh: '/absolute/path/expression.csv 或目录',
    importUrlPh: 'https://example.com/GSE12345_series_matrix.txt.gz',
    importChooseFile: '选择文件',
    importFormat: '数据格式',
    importAuto: '自动识别',
    importAdvanced: '高级选项',
    importSheet: 'Sheet 名',
    importGeneCol: '基因列',
    importTranspose: '转置矩阵',
    importMember: '压缩包内文件',
    importObject: 'RData 对象',
    importSpecies: '物种',
    importSampleSheet: '样本信息表路径',
    importClinical: '临床表路径',
    importStart: '导入',
    importBusy: '导入中…',
    importInspect: '仅探测',
    importInspectBusy: '探测中…',
    importCandidates: '候选文件，点击选择后重试',
    importDone: '导入成功',
    importError: '导入失败',
    importHint: '同一路径重复导入会复用已打开的数据集；样本名会在样本表中自动匹配。',
    importUploadHint: '上传文件先保存到工作台目录，再按所选格式打开。',
    importInspectTitle: '探测结果',
    importNoPath: '请先填写本机路径或 URL',
    importUploading: '上传中…',
    analysisTitle: '分析目录',
    analysisHint: '高级分析按分组展示；灰色项表示前置未满足。',
    analysisAvailable: '可用',
    analysisBlocked: '不可用',
    analysisParams: '参数',
    analysisRun: '运行分析',
    analysisRunning: '运行中…',
    analysisEmpty: '暂无可用的高级分析',
    analysisResult: '分析结果',
    analysisTables: '表格',
    analysisFigures: '图',
    analysisWarnings: '提示',
    analysisClose: '收起',
    analysisNoParams: '该分析无需参数',
    analysisNeedDataset: '先载入数据集',
    analysisNeeds: '前置',
    resultTitle: '结果',
    resultMetrics: '指标',
    resultTables: '表格',
    resultFigures: '图',
    resultWarnings: '警告',
    resultNotes: '说明',
    resultFiles: '数据文件',
    resultEmpty: '没有可展示的结果',
    resultRows: '行',
    jobRunning: '任务运行中',
    jobFailed: '任务失败',
    jobTimeout: '任务超时，请稍后刷新查看',
    jobProgress: '进度',
    activityUser: '用户',
    activityAgent: 'Agent',
    statusDot: '状态',
    engineMode: '运行模式',
    savedTo: '已保存',
    openInFolder: '在文件管理器中打开',
    noSelection: '未选择',
    defaultOption: '默认',
    runningShort: '运行',
    errorInline: '请求失败',
    confirmClose: '确定关闭这个数据集？',
    notFound: '未找到',
    stepLog: '步骤记录',
    analysisLog: '分析记录',
    at: '时间',
    actor: '操作者',
    action: '动作',
    duration: '耗时',
    resultColumn: '结果',
  }

  var DICT_EN = {
    nav: 'Bulk Workbench',
    appTitle: 'GEO · Bulk Transcriptomics Workbench',
    appSubtitle: 'GEO search · expression QC · DE / enrichment / immunity / survival · single-cell mapping',
    refresh: 'Refresh',
    refreshHint: 'Reload workbench state',
    retry: 'Retry',
    loading: 'Loading…',
    loadingData: 'Reading workbench…',
    noData: 'No dataset loaded',
    errorTitle: 'Something went wrong',
    crashTitle: 'This panel failed to render',
    crashHint: 'Other panels still work. Click retry to recover.',
    close: 'Close',
    cancel: 'Cancel',
    copy: 'Copy',
    copied: 'Copied',
    copyFailed: 'Copy failed (clipboard unavailable)',
    search: 'Search',
    all: 'All',
    none: 'None',
    yes: 'Yes',
    no: 'No',
    unknown: 'Unknown',
    dataset: 'Dataset',
    samples: 'Samples',
    genes: 'Genes',
    groups: 'Groups',
    batches: 'Batches',
    engine: 'Engine',
    modules: 'Modules',
    sourceFormat: 'Source format',
    species: 'Species',
    warnings: 'Warnings',
    notes: 'Notes',
    current: 'Current',
    completed: 'Completed',
    pending: 'Pending',
    running: 'Running',
    failed: 'Failed',
    plot: 'Plot',
    plotting: 'Plotting…',
    noFigure: 'No figure',
    figureFail: 'Failed to load image (file may be cleaned)',
    exportBtn: 'Export',
    exported: 'Exported',
    exportFail: 'Export failed',
    copyPath: 'Copy path',
    pathLabel: 'File path',
    tableRows: '{v} rows',
    tableTruncated: 'showing first {v} rows',
    emptyTable: 'Empty table',
    total: 'Total',
    selectAll: 'Select all',
    clearSelection: 'Clear',
    of: '/',
    selectDataset: 'Switch dataset',
    closeDataset: 'Close dataset',
    workdir: 'Workdir',
    engineReady: 'Engine ready',
    engineBooting: 'Engine starting',
    engineBad: 'Engine unavailable',
    engineGuideTitle: 'Python engine is not ready',
    engineGuideBody: 'The workbench needs local Python (numpy / pandas / scipy / sklearn / matplotlib / anndata / lifelines). The page still opens, but computing and plotting need the engine.',
    engineGuideSteps: 'Check: (1) errors in engine.log under workdir; (2) python3 path in plugin settings; (3) restart the plugin or host, then retry.',
    pythonLabel: 'Python',
    rLabel: 'R',
    rMissing: 'R not installed',
    rVersion: 'R version',
    engineInfo: 'Engine info',
    statusOk: 'OK',
    statusWarn: 'Warning',
    statusFail: 'Error',
    pipelineTitle: 'Pipeline',
    pipelineHint: 'Run in order; rerun when needed (downstream results are cleared)',
    datasetsTitle: 'Datasets',
    stepStatusDone: 'Done',
    stepStatusRunning: 'Running',
    stepStatusReady: 'Ready',
    stepStatusBlocked: 'Blocked',
    stepStatusFailed: 'Failed',
    stepRun: 'Run',
    stepRerun: 'Rerun',
    stepRunning: 'Running…',
    stepParams: 'Parameters',
    stepNeedData: 'Import or open a dataset first',
    stepNeedPrereq: 'Requires: ',
    stepSummaryEmpty: 'No result summary yet',
    stepResult: 'Result',
    datasetNoDataset: 'No dataset yet',
    datasetOpenHint: 'Open an expression matrix via Import, or download public data via GEO search.',
    datasetChips: 'Datasets',
    datasetStepsDone: '{v} steps done',
    datasetCurrentStep: 'Current: {v}',
    warningsLabel: 'Warnings',
    tabOverview: 'Overview',
    tabExpr: 'Expression',
    tabDe: 'DE',
    tabEnrich: 'Enrichment',
    tabImmune: 'Immunity',
    tabSurvival: 'Survival',
    tabPrognosis: 'Prognosis',
    tabSc: 'Single-cell',
    tabNcrna: 'Noncoding',
    tabRecords: 'Records',
    tabLog: 'Log',
    overviewTitle: 'Dataset overview',
    metricSamples: 'Samples',
    metricGenes: 'Genes',
    metricGroups: 'Groups',
    metricBatches: 'Batches',
    metricEngine: 'Engine',
    metricModules: 'Modules',
    pcaPlot: 'PCA samples',
    corrPlot: 'Sample correlation',
    geneQuickTitle: 'Quick gene boxplot',
    geneQuickHint: 'Enter comma-separated genes and plot expression distribution.',
    geneInputPlaceholder: 'e.g. CDH1, VIM, MKI67',
    nextSteps: 'Suggested next steps',
    noSuggestions: 'No suggestions; continue with the pipeline on the left.',
    suggestionDetail: 'Suggestion',
    analysisCatalog: 'Analysis catalog',
    geoSearch: 'GEO search',
    importData: 'Import data',
    demoLoad: 'Load demo data',
    exportExprCsv: 'Expression CSV',
    exportDeCsv: 'DE results CSV',
    exportEnrichCsv: 'Enrichment CSV',
    exportImmuneCsv: 'Immune infiltration CSV',
    exportSampleCsv: 'Sample info CSV',
    exportSummary: 'Summary JSON',
    exportProtocol: 'Analysis protocol',
    exportMethods: 'Methods text',
    exportReport: 'Analysis report',
    exportScript: 'Reproducible script',
    exprTitle: 'Expression matrix',
    exprHint: 'Search genes and add them; box/heatmap plots come from the engine. Export CSV for the full matrix.',
    exprGeneSearch: 'Search gene (Enter to add)',
    exprAddGene: 'Add',
    exprSelected: 'Selected genes',
    exprExport: 'Export expression CSV',
    exprExportHint: 'The CSV lives on the host; copy the path for downstream scripts.',
    exprNoGene: 'Select at least one gene',
    exprPlotBox: 'Gene boxplot',
    exprPlotHeatmap: 'Gene heatmap',
    exprGeneNotInData: 'Gene is not in this dataset',
    exprEmpty: 'Search and add genes above.',
    exprTableFallback: 'The full matrix is not returned inline: export the expression CSV; table results come from analysis endpoints.',
    exprExportedPath: 'Expression CSV exported',
    deTitle: 'Differential expression',
    deHint: 'Run the DE step on the left to see the table, volcano, MA plot and heatmap.',
    deEmpty: 'No DE result yet.',
    deVolcano: 'Volcano',
    deMa: 'MA plot',
    deHeatmap: 'DE heatmap',
    deSortHint: 'Click a header to sort',
    deRegulation: 'Regulation',
    deUp: 'Up',
    deDown: 'Down',
    deNs: 'NS',
    deExport: 'Export DE table',
    deTableTitle: 'DE genes',
    enrichTitle: 'Functional enrichment',
    enrichHint: 'Complete DE first; choose libraries and method, then run.',
    enrichDb: 'Libraries',
    enrichMethod: 'Method',
    enrichP: 'P threshold',
    enrichQ: 'q threshold',
    enrichTopN: 'Top N',
    enrichRun: 'Run enrichment',
    enrichBar: 'Enrichment barplot',
    enrichDot: 'Enrichment dotplot',
    enrichEmpty: 'No enrichment result yet.',
    enrichNoDb: 'Select at least one library',
    immuneTitle: 'Immune infiltration',
    immuneHint: 'More methods take longer; R IOBR is preferred with automatic Python fallback.',
    immuneMethods: 'Methods',
    immuneGroupCol: 'Group column',
    immuneRun: 'Run immune infiltration',
    immuneStack: 'Cell fraction stackplot',
    immuneCorr: 'Method correlation',
    immuneTable: 'Cell fraction table',
    immuneDiff: 'Group difference',
    immuneEmpty: 'No immune result yet.',
    survivalTitle: 'Clinical survival',
    survivalHint: 'Sample info needs time and event columns (provide a clinical table on import).',
    survivalTimeCol: 'Time column',
    survivalEventCol: 'Event column',
    survivalEndpoint: 'Endpoint',
    survivalGene: 'Gene (optional)',
    survivalCutoff: 'Cutoff',
    survivalRun: 'Run survival analysis',
    survivalKm: 'KM curves',
    survivalForest: 'Cox forest',
    survivalEmpty: 'No survival result yet.',
    survivalNeedColumns: 'No time/event columns detected; import a clinical table or check column names.',
    prognosisTitle: 'Advanced prognosis',
    prognosisHint: 'Run survival analysis first.',
    prognosisModel: 'Model',
    prognosisGenes: 'Genes (empty = DE top)',
    prognosisNfold: 'CV folds',
    prognosisRun: 'Run prognosis model',
    prognosisRisk: 'Risk score plot',
    prognosisRoc: 'Time-dependent ROC',
    prognosisCindex: 'C-index',
    prognosisEmpty: 'No prognosis result yet.',
    scTitle: 'Single-cell mapping',
    scHint: 'Choose an h5ad reference: local path / detected local h5ad / pull from single-cell workbench.',
    scRefPath: 'h5ad path',
    scLocalRefs: 'Local references',
    scPull: 'Pull from single-cell workbench',
    scAttach: 'Attach reference',
    scGeneSet: 'Signature source',
    scGenes: 'Custom genes',
    scCellTypeCol: 'Cell type column',
    scScoreMethod: 'Scoring method',
    scRun: 'Run single-cell mapping',
    scScoreTable: 'Cell type scores',
    scExprFig: 'Subcluster expression',
    scEmpty: 'No single-cell mapping result yet.',
    scNoRef: 'Choose and attach an h5ad reference first',
    scAttachOk: 'Reference attached',
    scAttachFail: 'Failed to attach reference',
    scFromWorkspace: 'Single-cell workbench datasets',
    scLocal: 'Local h5ad',
    scAvailable: 'Single-cell workbench online',
    scUnavailable: 'Single-cell workbench not running (local h5ad only)',
    scHintNoRefs: 'No h5ad detected; enter an absolute local path.',
    ncrnaTitle: 'Noncoding / fusion / splicing',
    ncrnaHint: 'Attach module data first, then run that module analysis.',
    ncrnaKind: 'Module',
    ncrnaPath: 'Module file path',
    ncrnaAttach: 'Attach module',
    ncrnaAttached: 'Attached',
    ncrnaRun: 'Run noncoding analysis',
    ncrnaEmpty: 'No noncoding result yet.',
    ncrnaNotAttached: 'No data attached for this module',
    ncrnaAttachOk: 'Module attached',
    ncrnaAttachFail: 'Failed to attach module',
    moduleCirc: 'circRNA',
    moduleMirna: 'miRNA',
    moduleLncrna: 'lncRNA',
    moduleFusion: 'Fusion',
    moduleSplicing: 'Splicing',
    recordsTitle: 'Records & audit',
    recordsHistory: 'Run history',
    recordsAudit: 'Deterministic audit',
    recordsAuditScore: 'Audit score',
    recordsRunAudit: 'Run audit',
    recordsExport: 'Export',
    recordsNoHistory: 'No runs yet',
    recordsAnalysis: 'Advanced analysis catalog',
    recordsExportOk: 'Exported',
    recordsExportFail: 'Export failed',
    recordsAuditEmpty: 'Audit has not run yet',
    recordsExporting: 'Exporting…',
    auditOk: 'Pass',
    auditWarn: 'Warn',
    auditFail: 'Fail',
    auditHint: 'Audit uses deterministic rules (no model): score = 100 - 12×warn - 30×fail.',
    logTitle: 'Activity log',
    logEmpty: 'No activity yet',
    logRefresh: 'Refresh log',
    agentTitle: 'Agent panel',
    agentStatus: 'Workbench status',
    agentSuggestions: 'Suggested next steps',
    agentActivity: 'Recent activity',
    agentNoActivity: 'No activity',
    agentNoDataset: 'Load a dataset to see status and suggestions',
    agentCopy: 'Copy instruction',
    agentCopied: 'Instruction copied',
    agentCopyFailed: 'Copy failed',
    geoTitle: 'GEO data search',
    geoSubtitle: 'Search NCBI GEO Series with localized fields; download selected series to a local directory and import in one click.',
    geoTermPh: 'gastric cancer cisplatin resistance expression profiling',
    geoSpecies: 'Species',
    geoAll: 'All',
    geoHuman: 'Human',
    geoMouse: 'Mouse',
    geoStudyType: 'Study type',
    geoStudyAll: 'All',
    geoStudyArray: 'Array profiling',
    geoStudySeq: 'High-throughput sequencing',
    geoMinSamples: 'Min samples',
    geoHasClinical: 'Has clinical data',
    geoSearchBtn: 'Search',
    geoSearching: 'Searching…',
    geoResults: 'Results',
    geoTotal: '{v} results',
    geoNoResults: 'No matching GEO series; try another term.',
    geoDetailBtn: 'Details',
    geoDetailTitle: 'Series details',
    geoLoadingDetail: 'Loading series details…',
    geoFiles: 'Downloadable files',
    geoDownloadDir: 'Download directory',
    geoDownloadDirHint: 'Default ~/.dsh/geobulk/GEO; any absolute path works',
    geoDownloadBtn: 'Download selected',
    geoDownloading: 'Downloading…',
    geoImportBtn: 'Import & QC',
    geoImporting: 'Importing…',
    geoSelectAll: 'Select all',
    geoClear: 'Clear',
    geoJobProgress: 'Download progress',
    geoJobDone: 'Download complete',
    geoJobFailed: 'Job failed',
    geoNeedFile: 'Select at least one file',
    geoSummary: 'Summary',
    geoDesign: 'Design',
    geoPlatform: 'Platform',
    geoRelease: 'Released',
    geoPubmed: 'PubMed',
    geoMatrix: 'Series matrix',
    geoClinicalTag: 'Clinical tags',
    geoImportDone: 'Imported; starting QC…',
    geoHint: 'Search fails offline; use Import to load local files instead.',
    geoErrorSearch: 'GEO search failed',
    geoErrorDetail: 'Failed to load series details',
    geoYear: 'Year',
    geoOpenAfter: 'Import after download',
    geoSeriesMatrix: 'Series matrix',
    geoSupp: 'Supplementary',
    geoRaw: 'Raw data',
    geoOther: 'Other',
    geoNoFiles: 'No downloadable files parsed for this series',
    importTitle: 'Import data',
    importSubtitle: 'Expression matrix, GEO series matrix, xlsx, zip/tar, RDS/RData and more',
    importSourceLocal: 'Local path',
    importSourceUrl: 'URL',
    importSourceUpload: 'Upload file',
    importPathPh: '/absolute/path/expression.csv or a directory',
    importUrlPh: 'https://example.com/GSE12345_series_matrix.txt.gz',
    importChooseFile: 'Choose file',
    importFormat: 'Data format',
    importAuto: 'Auto detect',
    importAdvanced: 'Advanced options',
    importSheet: 'Sheet',
    importGeneCol: 'Gene column',
    importTranspose: 'Transpose matrix',
    importMember: 'Archive member',
    importObject: 'RData object',
    importSpecies: 'Species',
    importSampleSheet: 'Sample sheet path',
    importClinical: 'Clinical table path',
    importStart: 'Import',
    importBusy: 'Importing…',
    importInspect: 'Inspect only',
    importInspectBusy: 'Inspecting…',
    importCandidates: 'Candidates; pick one and retry',
    importDone: 'Imported',
    importError: 'Import failed',
    importHint: 'Importing the same path reuses the open dataset; sample names are matched automatically.',
    importUploadHint: 'Uploads are saved to the workbench directory before opening.',
    importInspectTitle: 'Inspect result',
    importNoPath: 'Provide a local path or URL first',
    importUploading: 'Uploading…',
    analysisTitle: 'Analysis catalog',
    analysisHint: 'Advanced analyses grouped by category; grey items have unmet prerequisites.',
    analysisAvailable: 'Available',
    analysisBlocked: 'Unavailable',
    analysisParams: 'Parameters',
    analysisRun: 'Run analysis',
    analysisRunning: 'Running…',
    analysisEmpty: 'No advanced analysis available',
    analysisResult: 'Analysis result',
    analysisTables: 'Tables',
    analysisFigures: 'Figures',
    analysisWarnings: 'Notes',
    analysisClose: 'Collapse',
    analysisNoParams: 'No parameters',
    analysisNeedDataset: 'Load a dataset first',
    analysisNeeds: 'Requires',
    resultTitle: 'Result',
    resultMetrics: 'Metrics',
    resultTables: 'Tables',
    resultFigures: 'Figures',
    resultWarnings: 'Warnings',
    resultNotes: 'Notes',
    resultFiles: 'Data files',
    resultEmpty: 'Nothing to show',
    resultRows: 'rows',
    jobRunning: 'Job running',
    jobFailed: 'Job failed',
    jobTimeout: 'Job timed out; refresh later',
    jobProgress: 'Progress',
    activityUser: 'user',
    activityAgent: 'agent',
    statusDot: 'Status',
    engineMode: 'Mode',
    savedTo: 'Saved',
    openInFolder: 'Open in file manager',
    noSelection: 'Not selected',
    defaultOption: 'default',
    runningShort: 'Run',
    errorInline: 'Request failed',
    confirmClose: 'Close this dataset?',
    notFound: 'Not found',
    stepLog: 'Step history',
    analysisLog: 'Analysis history',
    at: 'Time',
    actor: 'Actor',
    action: 'Action',
    duration: 'Duration',
    resultColumn: 'Result',
  }

  // ==========================================================================
  // 纯函数层：格式化 / 取色 / 坐标 / 搜索 payload / 任务状态机 / 响应归一化
  // ==========================================================================

  function isNil(value) {
    return value === null || value === undefined
  }

  function clamp(value, min, max) {
    var number = Number(value)
    if (!isFinite(number)) number = 0
    if (isFinite(Number(min)) && number < Number(min)) number = Number(min)
    if (isFinite(Number(max)) && number > Number(max)) number = Number(max)
    return number
  }

  function asNumber(value) {
    if (isNil(value) || value === '') return null
    var number = Number(value)
    return isFinite(number) ? number : null
  }

  /** 去掉空白、把 undefined 变空串；可选截断。 */
  function plainText(value, max) {
    var text = isNil(value) ? '' : String(value).replace(/\s+/g, ' ').trim()
    if (isFinite(Number(max)) && Number(max) > 0 && text.length > Number(max)) {
      return text.slice(0, Number(max) - 1) + '…'
    }
    return text
  }

  /** 千分位/紧凑数字：12345 → 12,345；1234567 → 1.23M；null → —。 */
  function fmtNumber(value, digits) {
    if (isNil(value) || value === '') return '—'
    var number = Number(value)
    if (!isFinite(number)) return '—'
    var fixed = isFinite(Number(digits)) ? clamp(digits, 0, 6) : 2
    var abs = Math.abs(number)
    if (abs >= 1e9) return trimFixed(number / 1e9, fixed) + 'B'
    if (abs >= 1e6) return trimFixed(number / 1e6, fixed) + 'M'
    if (abs >= 1e4) return trimFixed(number / 1e3, fixed) + 'K'
    var parts = trimFixed(number, fixed).split('.')
    parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ',')
    return parts.join('.')
  }

  function trimFixed(number, digits) {
    var text = Number(number).toFixed(clamp(digits, 0, 6))
    if (text.indexOf('.') >= 0) text = text.replace(/0+$/, '').replace(/\.$/, '')
    if (text === '-0') text = '0'
    return text
  }

  /** 比例→百分比：0.1234 → 12.3%；>1 视为已经是百分数。 */
  function fmtPct(value, digits) {
    if (isNil(value) || value === '') return '—'
    var number = Number(value)
    if (!isFinite(number)) return '—'
    var fixed = isFinite(Number(digits)) ? clamp(digits, 0, 4) : 1
    var percent = Math.abs(number) <= 1 ? number * 100 : number
    return trimFixed(percent, fixed) + '%'
  }

  /** P 值：<0.001 用科学计数；其余保留 3-4 位。 */
  function fmtP(value) {
    if (isNil(value) || value === '') return '—'
    var number = Number(value)
    if (!isFinite(number)) return '—'
    if (number <= 0) return '0'
    if (number < 0.001) return number.toExponential(2).replace('e-', 'e-')
    if (number < 0.01) return number.toFixed(4)
    return number.toFixed(3)
  }

  /** 字节：1536 → 1.5 KB。 */
  function formatBytes(bytes) {
    if (isNil(bytes) || bytes === '') return '—'
    var number = Number(bytes)
    if (!isFinite(number) || number < 0) return '—'
    if (number === 0) return '0 B'
    var units = ['B', 'KB', 'MB', 'GB', 'TB']
    var exponent = Math.min(Math.floor(Math.log(number) / Math.log(1024)), units.length - 1)
    var value = number / Math.pow(1024, exponent)
    return value.toFixed(exponent === 0 ? 0 : 1) + ' ' + units[exponent]
  }

  /** 毫秒：812 → 812 ms；65000 → 1m 5s。 */
  function formatDuration(ms) {
    if (isNil(ms) || ms === '') return '—'
    var number = Number(ms)
    if (!isFinite(number) || number < 0) return '—'
    if (number < 1000) return Math.round(number) + ' ms'
    if (number < 60000) return trimFixed(number / 1000, 1) + ' s'
    var totalSec = Math.round(number / 1000)
    var minutes = Math.floor(totalSec / 60)
    if (minutes < 60) return minutes + 'm ' + (totalSec % 60) + 's'
    var hours = Math.floor(minutes / 60)
    return hours + 'h ' + (minutes % 60) + 'm'
  }

  /** ISO 时间 → 06-12 14:03（本地时区，不依赖 toLocaleString 输出格式）。 */
  function formatClock(value) {
    if (isNil(value) || value === '') return '—'
    var date = value instanceof Date ? value : new Date(value)
    if (!date || isNaN(date.getTime())) {
      var text = String(value)
      return text.length > 19 ? text.slice(5, 16).replace('T', ' ') : text
    }
    function pad(number) { return number < 10 ? '0' + number : String(number) }
    return pad(date.getMonth() + 1) + '-' + pad(date.getDate()) + ' ' + pad(date.getHours()) + ':' + pad(date.getMinutes())
  }

  function formatError(error) {
    if (!error) return ''
    if (typeof error === 'string') return error
    if (typeof error.message === 'string' && error.message) return error.message
    try { return JSON.stringify(error) } catch (ignored) { return String(error) }
  }

  function hexToRgb(hex) {
    var text = String(hex || '').replace('#', '')
    if (text.length === 3) text = text.charAt(0) + text.charAt(0) + text.charAt(1) + text.charAt(1) + text.charAt(2) + text.charAt(2)
    if (!/^[0-9a-fA-F]{6}$/.test(text)) return { r: 0, g: 0, b: 0 }
    return {
      r: parseInt(text.slice(0, 2), 16),
      g: parseInt(text.slice(2, 4), 16),
      b: parseInt(text.slice(4, 6), 16),
    }
  }

  /** 分组稳定取色：同一名称永远同一色；数字视为调色板下标。 */
  function colorForGroup(value) {
    if (typeof value === 'number' && isFinite(value)) {
      return PALETTE[Math.abs(Math.round(value)) % PALETTE.length]
    }
    var text = isNil(value) ? '' : String(value)
    var hash = 0
    for (var index = 0; index < text.length; index++) {
      hash = (hash * 31 + text.charCodeAt(index)) >>> 0
    }
    return PALETTE[hash % PALETTE.length]
  }

  /** 按 key（字符串字段名或函数）分组，返回普通对象，保持首次出现顺序。 */
  function groupBy(list, keyFn) {
    var source = Array.isArray(list) ? list : []
    var pick = typeof keyFn === 'function'
      ? keyFn
      : function (item) { return isNil(item) ? '' : item[keyFn] }
    var out = {}
    for (var index = 0; index < source.length; index++) {
      var key = pick(source[index], index)
      key = isNil(key) ? '' : String(key)
      if (!Object.prototype.hasOwnProperty.call(out, key)) out[key] = []
      out[key].push(source[index])
    }
    return out
  }

  /**
   * 火山图点：rows 为 DE 行（gene/logFC/P.Value/adj.P.Val）。
   * 过滤非数值行；y=-log10(p)；sig = |logFC|>=lfcThreshold 且 p<=padjThreshold。
   */
  function volcanoPoints(rows, options) {
    var opts = options && typeof options === 'object' ? options : {}
    var lfcThreshold = isFinite(Number(opts.lfcThreshold)) ? Number(opts.lfcThreshold) : 1
    var padjThreshold = isFinite(Number(opts.padjThreshold)) ? Number(opts.padjThreshold) : 0.05
    var useAdjusted = opts.useAdjusted !== false
    var source = Array.isArray(rows) ? rows : []
    var points = []
    for (var index = 0; index < source.length; index++) {
      var row = source[index] && typeof source[index] === 'object' ? source[index] : {}
      var gene = plainText(row.gene || row.ID || row.symbol || '')
      var logFC = asNumber(row.logFC)
      if (logFC === null) logFC = asNumber(row.log2FoldChange)
      var rawP = asNumber(row['P.Value'])
      if (rawP === null) rawP = asNumber(row.pvalue)
      if (rawP === null) rawP = asNumber(row.p_value)
      var adjusted = asNumber(row['adj.P.Val'])
      if (adjusted === null) adjusted = asNumber(row.padj)
      if (adjusted === null) adjusted = asNumber(row.p_adjust)
      var judgeP = useAdjusted && adjusted !== null ? adjusted : rawP
      if (logFC === null || judgeP === null) continue
      if (rawP !== null && rawP < 0) continue
      if (judgeP < 0) continue
      var safeP = judgeP <= 0 ? 1e-300 : judgeP
      var y = Math.min(-Math.log(safeP) / Math.LN10, 300)
      var regulation = row.regulation === 'up' || row.regulation === 'down' || row.regulation === 'ns'
        ? row.regulation
        : (Math.abs(logFC) >= lfcThreshold ? (logFC > 0 ? 'up' : 'down') : 'ns')
      var sig = Math.abs(logFC) >= lfcThreshold && judgeP <= padjThreshold
      points.push({
        gene: gene || ('#' + (index + 1)),
        x: logFC,
        y: y,
        p: rawP === null ? judgeP : rawP,
        padj: adjusted,
        regulation: regulation,
        sig: sig,
        color: VOLCANO_COLORS[regulation] || VOLCANO_COLORS.ns,
      })
    }
    return points
  }

  /**
   * PCA 坐标变换（世界坐标 → 屏幕像素）：
   *   自动按点集包围盒缩放并留 padding，y 轴翻转（PC2 向上）。
   *   view = { scale, tx, ty }：scale 在 fit 基础上再缩放，tx/ty 是屏幕像素平移。
   */
  function pcaScreen(points, view, width, height, padding) {
    var boxWidth = asNumber(width)
    var boxHeight = asNumber(height)
    if (boxWidth === null || boxWidth <= 0) boxWidth = 600
    if (boxHeight === null || boxHeight <= 0) boxHeight = 400
    var pad = asNumber(padding)
    if (pad === null || pad < 0) pad = 24
    pad = Math.min(pad, Math.min(boxWidth, boxHeight) / 4)
    var list = []
    var source = Array.isArray(points) ? points : []
    for (var index = 0; index < source.length; index++) {
      var point = source[index]
      var px = null
      var py = null
      if (Array.isArray(point)) {
        px = asNumber(point[0])
        py = asNumber(point[1])
      } else if (point && typeof point === 'object') {
        px = asNumber(point.x)
        py = asNumber(point.y)
      }
      if (px === null || py === null) continue
      list.push({ x: px, y: py, label: point && typeof point === 'object' && point.label ? String(point.label) : '', source: point })
    }
    var minX = 0
    var maxX = 0
    var minY = 0
    var maxY = 0
    for (var cursor = 0; cursor < list.length; cursor++) {
      if (cursor === 0) {
        minX = maxX = list[cursor].x
        minY = maxY = list[cursor].y
      } else {
        if (list[cursor].x < minX) minX = list[cursor].x
        if (list[cursor].x > maxX) maxX = list[cursor].x
        if (list[cursor].y < minY) minY = list[cursor].y
        if (list[cursor].y > maxY) maxY = list[cursor].y
      }
    }
    var spanX = maxX - minX
    var spanY = maxY - minY
    if (!(spanX > 1e-6)) spanX = 1
    if (!(spanY > 1e-6)) spanY = 1
    var centerX = (minX + maxX) / 2
    var centerY = (minY + maxY) / 2
    var fitScale = Math.min((boxWidth - pad * 2) / spanX, (boxHeight - pad * 2) / spanY)
    if (!(fitScale > 0) || !isFinite(fitScale)) fitScale = 1
    var viewScale = view && isFinite(Number(view.scale)) && Number(view.scale) > 0 ? Number(view.scale) : 1
    var translateX = view && isFinite(Number(view.tx)) ? Number(view.tx) : 0
    var translateY = view && isFinite(Number(view.ty)) ? Number(view.ty) : 0
    var scale = fitScale * viewScale
    var out = []
    for (var cursor2 = 0; cursor2 < list.length; cursor2++) {
      out.push({
        x: boxWidth / 2 + (list[cursor2].x - centerX) * scale + translateX,
        y: boxHeight / 2 - (list[cursor2].y - centerY) * scale + translateY,
        label: list[cursor2].label,
        raw: list[cursor2],
      })
    }
    return out
  }

  /** GEO 搜索 → POST /geo/search 的 body（与 CONTRACT §4 字段一致）。 */
  function buildSearchPayload(filters) {
    var raw = filters && typeof filters === 'object' ? filters : {}
    var term = plainText(raw.term, 500)
    var payload = {
      term: term,
      species: raw.species === 'human' || raw.species === 'mouse' ? raw.species : 'all',
      study_type: raw.study_type === 'expression profiling by array' || raw.study_type === 'high throughput sequencing'
        ? raw.study_type
        : 'all',
      limit: clamp(isFinite(Number(raw.limit)) && Number(raw.limit) > 0 ? Number(raw.limit) : 20, 1, 100),
      page: clamp(isFinite(Number(raw.page)) ? Number(raw.page) : 0, 0, 10000),
      sort: raw.sort === 'date' || raw.sort === 'samples' ? raw.sort : 'relevance',
      actor: 'user',
    }
    var minSamples = asNumber(raw.min_samples)
    if (minSamples !== null && minSamples > 0) payload.min_samples = Math.round(minSamples)
    if (raw.has_clinical === true) payload.has_clinical = true
    return payload
  }

  /** 建议 → 可直接粘给 Agent 的中文指令（含步骤名/step id）。 */
  function suggestionInstruction(suggestion) {
    var raw = suggestion && typeof suggestion === 'object' ? suggestion : {}
    var title = plainText(raw.title, 120)
    var detail = plainText(raw.detail, 400)
    var step = plainText(raw.step, 40)
    var stepMeta = step && STEP_BY_ID[step] ? STEP_BY_ID[step] : null
    if (!title) title = stepMeta ? stepMeta.title : '继续分析'
    var text = '请在工作台执行「' + title + '」'
    if (stepMeta) text += '（步骤：' + stepMeta.title + '，step=' + step + '）'
    else if (step) text += '（step=' + step + '）'
    if (detail && detail !== title) text += '：' + detail
    return text + '。'
  }

  /** 步骤默认参数（数组会复制，避免多个表单共享同一引用）。 */
  function defaultParams(step) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    var out = {}
    if (!meta || !Array.isArray(meta.params)) return out
    for (var index = 0; index < meta.params.length; index++) {
      var descriptor = meta.params[index]
      if (!descriptor || !descriptor.key) continue
      var value = descriptor.default
      if (Array.isArray(value)) value = value.slice()
      out[descriptor.key] = value === undefined ? null : value
    }
    return out
  }

  /** 表单值 → 请求 params（number 转 NaN→null；columns/genes 支持逗号字符串）。 */
  function collectStepParams(step, values) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    var merged = defaultParams(meta)
    var source = values && typeof values === 'object' ? values : {}
    var keys = Object.keys(source)
    for (var index = 0; index < keys.length; index++) merged[keys[index]] = source[keys[index]]
    if (!meta || !Array.isArray(meta.params)) return merged
    var out = {}
    for (var cursor = 0; cursor < meta.params.length; cursor++) {
      var descriptor = meta.params[cursor]
      if (!descriptor || !descriptor.key) continue
      var value = merged[descriptor.key]
      var type = descriptor.type || 'text'
      if (type === 'number') {
        var number = asNumber(value)
        out[descriptor.key] = number === null ? null : number
      } else if (type === 'boolean') {
        out[descriptor.key] = value === true
      } else if (type === 'multiselect' || type === 'columns' || type === 'genes' || type === 'array') {
        out[descriptor.key] = parseTokenList(value)
      } else if (type === 'column') {
        out[descriptor.key] = plainText(value) || null
      } else {
        out[descriptor.key] = isNil(value) ? '' : String(value)
      }
    }
    return out
  }

  function parseTokenList(value) {
    if (Array.isArray(value)) return value.slice()
    var text = plainText(value)
    if (!text) return []
    return text.split(/[,，;；\s]+/).filter(function (item) { return !!item })
  }

  function formatParamValue(descriptor, value) {
    if (isNil(value) || value === '') return '—'
    if (descriptor && descriptor.type === 'multiselect') {
      if (!Array.isArray(value)) return String(value)
      var labels = {}
      var options = Array.isArray(descriptor.options) ? descriptor.options : []
      for (var index = 0; index < options.length; index++) labels[options[index].value] = options[index].label || options[index].value
      return value.map(function (item) { return labels[item] || item }).join('、')
    }
    if (Array.isArray(value)) return value.length ? value.join(', ') : '—'
    if (typeof value === 'boolean') return value ? '是' : '否'
    return String(value)
  }

  /** 已完成步骤集合。 */
  function doneStepSet(dataset) {
    var set = {}
    var steps = dataset && Array.isArray(dataset.steps) ? dataset.steps : []
    for (var index = 0; index < steps.length; index++) set[String(steps[index])] = true
    return set
  }

  /** 前置是否全部完成（pca 的 requires 是 normalize，不是 batch）。 */
  function isStepReady(step, dataset) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    if (!meta) return false
    if (!dataset) return false
    var done = doneStepSet(dataset)
    var requires = Array.isArray(meta.requires) ? meta.requires : []
    for (var index = 0; index < requires.length; index++) {
      if (!done[requires[index]]) return false
    }
    return true
  }

  /**
   * 左栏状态：running > done > ready > blocked。
   * busyStepId 是正在运行的步骤；failedStepId 最近失败的步骤。
   */
  function stepStatusOf(step, dataset, busyStepId, failedStepId) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    if (!meta) return 'pending'
    if (busyStepId === meta.id) return 'running'
    if (failedStepId === meta.id) return 'failed'
    if (dataset && doneStepSet(dataset)[meta.id]) return 'done'
    if (dataset && isStepReady(meta, dataset)) return 'ready'
    return 'blocked'
  }

  function stepStatusKey(status) {
    if (status === 'done') return 'stepStatusDone'
    if (status === 'running') return 'stepStatusRunning'
    if (status === 'ready') return 'stepStatusReady'
    if (status === 'failed') return 'stepStatusFailed'
    return 'stepStatusBlocked'
  }

  function figureSrc(raw) {
    if (typeof raw === 'string' && raw) {
      if (raw.indexOf('data:') === 0 || raw.charAt(0) === '/') return raw
      var nameOnly = raw.split('/').pop()
      if (/^[A-Za-z0-9._-]+\.png$/.test(nameOnly)) return API + '/figure?id=' + encodeURIComponent(nameOnly)
      return ''
    }
    if (!raw || typeof raw !== 'object') return ''
    if (typeof raw.url === 'string' && raw.url) return raw.url
    if (typeof raw.file === 'string' && raw.file) {
      var name = raw.file.split('/').pop()
      if (/^[A-Za-z0-9._-]+\.png$/.test(name)) return API + '/figure?id=' + encodeURIComponent(name)
    }
    return ''
  }

  function normalizeFigure(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.figure && typeof payload.figure === 'object' ? payload.figure : payload)
      : {}
    var file = typeof raw.file === 'string' ? raw.file : ''
    var url = figureSrc({ url: typeof raw.url === 'string' ? raw.url : '', file: file })
    if (!url) return null
    return {
      id: isNil(raw.id) ? '' : String(raw.id),
      title: plainText(raw.title || raw.kind || '', 160),
      url: url,
      file: file,
      width: asNumber(raw.width) || 0,
      height: asNumber(raw.height) || 0,
      kind: typeof raw.kind === 'string' ? raw.kind : '',
      bytes: asNumber(raw.bytes),
    }
  }

  function normalizeTable(table, index) {
    var raw = table && typeof table === 'object' ? table : {}
    var columns = []
    if (Array.isArray(raw.columns)) {
      for (var columnIndex = 0; columnIndex < raw.columns.length; columnIndex++) {
        var column = raw.columns[columnIndex] && typeof raw.columns[columnIndex] === 'object' ? raw.columns[columnIndex] : {}
        var key = plainText(column.key || column.name || ('' + columnIndex))
        if (!key) continue
        columns.push({
          key: key,
          label: plainText(column.label || key),
          type: column.type === 'number' || column.numeric === true ? 'number' : 'string',
        })
      }
    }
    var rows = Array.isArray(raw.rows) ? raw.rows : []
    if (!columns.length && rows.length && rows[0] && typeof rows[0] === 'object') {
      var firstKeys = Object.keys(rows[0])
      for (var keyIndex = 0; keyIndex < firstKeys.length; keyIndex++) {
        columns.push({
          key: firstKeys[keyIndex],
          label: firstKeys[keyIndex],
          type: typeof rows[0][firstKeys[keyIndex]] === 'number' ? 'number' : 'string',
        })
      }
    }
    return {
      id: plainText(raw.id || ('table-' + index)),
      title: plainText(raw.title || ''),
      columns: columns,
      rows: rows,
      note: plainText(raw.note || '', 300),
      truncated: raw.truncated === true,
    }
  }

  /** §2.5 result 信封归一化（步骤/分析/导出共用）。 */
  function normalizeResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    if (raw.result && typeof raw.result === 'object' && !raw.summary) {
      raw = Object.assign({}, raw.result, {
        analysis_id: raw.analysis_id,
        ms: raw.ms || raw.result.ms,
        dataset: raw.dataset,
      })
    }
    var summary = raw.summary && typeof raw.summary === 'object' ? raw.summary : {}
    var bullets = []
    if (Array.isArray(summary.bullets)) {
      for (var index = 0; index < summary.bullets.length; index++) {
        var bullet = plainText(summary.bullets[index], 500)
        if (bullet) bullets.push(bullet)
      }
    }
    var metrics = {}
    if (summary.metrics && typeof summary.metrics === 'object') {
      var metricKeys = Object.keys(summary.metrics)
      for (var keyIndex = 0; keyIndex < metricKeys.length; keyIndex++) {
        metrics[metricKeys[keyIndex]] = summary.metrics[metricKeys[keyIndex]]
      }
    }
    var tables = []
    if (Array.isArray(raw.tables)) {
      for (var tableIndex = 0; tableIndex < raw.tables.length; tableIndex++) {
        tables.push(normalizeTable(raw.tables[tableIndex], tableIndex))
      }
    }
    var figures = []
    if (Array.isArray(raw.figures)) {
      for (var figureIndex = 0; figureIndex < raw.figures.length; figureIndex++) {
        var figure = normalizeFigure(raw.figures[figureIndex])
        if (figure) figures.push(figure)
      }
    }
    var dataOut = []
    if (Array.isArray(raw.data_out)) {
      for (var outIndex = 0; outIndex < raw.data_out.length; outIndex++) {
        var item = raw.data_out[outIndex] && typeof raw.data_out[outIndex] === 'object' ? raw.data_out[outIndex] : {}
        if (!item.path) continue
        dataOut.push({ kind: plainText(item.kind || ''), path: String(item.path), label: plainText(item.label || item.kind || '') })
      }
    }
    return {
      title: plainText(summary.title || raw.title || '', 200),
      summary: { title: plainText(summary.title || '', 200), bullets: bullets, metrics: metrics },
      bullets: bullets,
      metrics: metrics,
      tables: tables,
      figures: figures,
      data_out: dataOut,
      warnings: normalizeStringList(raw.warnings),
      notes: normalizeStringList(raw.notes),
      analysis_id: plainText(raw.analysis_id || ''),
      ms: asNumber(raw.ms),
      step: plainText(raw.step || ''),
      raw: raw,
    }
  }

  function normalizeStringList(value) {
    var list = Array.isArray(value) ? value : (isNil(value) || value === '' ? [] : [value])
    var out = []
    for (var index = 0; index < list.length; index++) {
      var text = plainText(list[index], 500)
      if (text) out.push(text)
    }
    return out
  }

  function normalizeDataset(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.dataset && typeof payload.dataset === 'object' ? payload.dataset : payload)
      : {}
    var id = plainText(raw.id || raw.dataset || raw.dataset_id || '')
    var steps = Array.isArray(raw.steps) ? raw.steps.map(function (item) { return plainText(item) }).filter(Boolean) : []
    var modules = []
    if (Array.isArray(raw.modules)) {
      for (var index = 0; index < raw.modules.length; index++) {
        var module = raw.modules[index] && typeof raw.modules[index] === 'object' ? raw.modules[index] : {}
        modules.push({
          kind: plainText(module.kind || module.module || ''),
          n: asNumber(module.n) || 0,
          path: typeof module.path === 'string' ? module.path : '',
        })
      }
    }
    var sampleColumns = []
    if (Array.isArray(raw.sample_columns)) {
      for (var columnIndex = 0; columnIndex < raw.sample_columns.length; columnIndex++) {
        var column = raw.sample_columns[columnIndex] && typeof raw.sample_columns[columnIndex] === 'object' ? raw.sample_columns[columnIndex] : {}
        var key = plainText(column.key || column.name || '')
        if (!key) continue
        sampleColumns.push({ key: key, kind: column.kind === 'numeric' || column.kind === 'time' || column.kind === 'event' ? column.kind : 'categorical', n_unique: asNumber(column.n_unique) })
      }
    }
    var history = []
    if (Array.isArray(raw.history)) {
      var entries = raw.history.slice(-MAX_HISTORY_ROWS)
      for (var historyIndex = 0; historyIndex < entries.length; historyIndex++) {
        var entry = entries[historyIndex] && typeof entries[historyIndex] === 'object' ? entries[historyIndex] : {}
        history.push({
          at: plainText(entry.at || entry.time || ''),
          step: plainText(entry.step || entry.action || ''),
          ok: entry.ok !== false && entry.status !== 'failed',
          ms: asNumber(entry.ms),
          summary: plainText(entry.summary || entry.message || '', 240),
          actor: plainText(entry.actor || ''),
        })
      }
    }
    var runs = []
    if (Array.isArray(raw.analyses_run)) {
      for (var runIndex = 0; runIndex < raw.analyses_run.length; runIndex++) {
        var run = raw.analyses_run[runIndex] && typeof raw.analyses_run[runIndex] === 'object' ? raw.analyses_run[runIndex] : {}
        runs.push({ at: plainText(run.at || ''), id: plainText(run.id || ''), ok: run.ok !== false, ms: asNumber(run.ms) })
      }
    }
    return {
      id: id,
      name: plainText(raw.name || id || ''),
      path: typeof raw.path === 'string' ? raw.path : '',
      source_format: plainText(raw.source_format || raw.format || ''),
      n_samples: asNumber(raw.n_samples),
      n_genes: asNumber(raw.n_genes),
      species: plainText(raw.species || ''),
      sample_columns: sampleColumns,
      steps: steps,
      current: plainText(raw.current || ''),
      modules: modules,
      warnings: normalizeStringList(raw.warnings),
      notes: normalizeStringList(raw.notes),
      history: history,
      analyses_run: runs,
      sc_ref: plainText(raw.sc_ref || ''),
      group_key: plainText(raw.group_key || ''),
      sc_ref: plainText(raw.sc_ref || ''),
    }
  }

  function normalizeHealth(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var engine = raw.engine && typeof raw.engine === 'object' ? raw.engine : raw
    var r = engine.r && typeof engine.r === 'object' ? engine.r : (raw.r && typeof raw.r === 'object' ? raw.r : {})
    var packages = engine.packages && typeof engine.packages === 'object' ? engine.packages : {}
    var rPackages = r.packages && typeof r.packages === 'object' ? r.packages : {}
    var ready = engine.ready === true
    if (isNil(engine.ready) && (typeof engine.python === 'string' || typeof engine.mode === 'string')) ready = true
    return {
      ready: ready,
      mode: plainText(engine.mode || ''),
      python: plainText(engine.python || ''),
      platform: plainText(engine.platform || ''),
      workdir: plainText(raw.workdir || engine.workdir || ''),
      packages: packages,
      geo: raw.geo !== false,
      r: {
        available: r.available !== false && !!r.version,
        version: plainText(r.version || ''),
        packages: rPackages,
      },
      error: plainText(raw.error || engine.error || ''),
    }
  }

  function engineBadge(engine) {
    var health = engine && (engine.r || engine.ready !== undefined) ? engine : normalizeHealth(engine)
    var py = health.python || ''
    var rVersion = health.r && health.r.available ? health.r.version : ''
    var detail = []
    if (py) detail.push('Python ' + py)
    detail.push(rVersion ? 'R ' + rVersion : 'R —')
    if (!health.ready) {
      return { ready: false, tone: 'bad', labelKey: health.error ? 'engineBad' : 'engineBooting', detail: detail.join(' · ') }
    }
    return { ready: true, tone: 'ok', labelKey: 'engineReady', detail: detail.join(' · '), python: py, r: rVersion }
  }

  function normalizeSuggestion(item, index) {
    var raw = item && typeof item === 'object' ? item : { title: String(item || '') }
    return {
      id: plainText(raw.id || ('suggestion-' + index)),
      title: plainText(raw.title || raw.detail || ''),
      detail: plainText(raw.detail || ''),
      step: plainText(raw.step || ''),
      params: raw.params && typeof raw.params === 'object' ? raw.params : null,
    }
  }

  function normalizeSuggestions(list) {
    if (!Array.isArray(list)) return []
    var out = []
    for (var index = 0; index < list.length; index++) out.push(normalizeSuggestion(list[index], index))
    return out
  }

  function normalizeActivityEntry(entry, index) {
    var raw = entry && typeof entry === 'object' ? entry : {}
    var ok = raw.ok !== false && raw.status !== 'failed' && raw.status !== 'error'
    var action = plainText(raw.action || raw.cmd || raw.step || raw.kind || raw.op || '')
    var text = plainText(raw.summary || raw.message || raw.detail || raw.error || (typeof entry === 'string' ? entry : ''), 300)
    return {
      id: plainText(raw.id || ('activity-' + index)),
      at: plainText(raw.at || raw.time || raw.ts || ''),
      actor: plainText(raw.actor || 'user'),
      action: action,
      text: text,
      ok: ok,
      ms: asNumber(raw.ms),
      dataset: plainText(raw.dataset || ''),
      step: plainText(raw.step || ''),
    }
  }

  function normalizeState(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var datasets = []
    if (Array.isArray(raw.datasets)) {
      for (var index = 0; index < raw.datasets.length; index++) datasets.push(normalizeDataset(raw.datasets[index]))
    }
    var activity = []
    if (Array.isArray(raw.activity)) {
      for (var activityIndex = 0; activityIndex < raw.activity.length; activityIndex++) {
        activity.push(normalizeActivityEntry(raw.activity[activityIndex], activityIndex))
      }
    }
    var active = plainText(raw.active || raw.active_id || '')
    if (!active && datasets.length) active = datasets[0].id
    var sc = raw.sc && typeof raw.sc === 'object' ? raw.sc : {}
    return {
      active: active,
      datasets: datasets,
      engine: normalizeHealth(raw.engine ? { engine: raw.engine } : raw),
      suggestions: normalizeSuggestions(raw.suggestions),
      recipes: Array.isArray(raw.recipes) ? raw.recipes : [],
      activity: activity,
      sc: {
        available: sc.available === true,
        datasets: Array.isArray(sc.datasets) ? sc.datasets : [],
        local: Array.isArray(sc.local) ? sc.local : [],
      },
    }
  }

  function normalizeAnalyses(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var list = Array.isArray(raw.analyses) ? raw.analyses : []
    var out = []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      var params = []
      if (Array.isArray(item.params)) {
        for (var paramIndex = 0; paramIndex < item.params.length; paramIndex++) {
          var descriptor = item.params[paramIndex] && typeof item.params[paramIndex] === 'object' ? item.params[paramIndex] : {}
          if (!descriptor.key) continue
          params.push({
            key: plainText(descriptor.key),
            label: plainText(descriptor.label || descriptor.key),
            type: plainText(descriptor.type || 'text'),
            default: descriptor.default,
            min: asNumber(descriptor.min),
            max: asNumber(descriptor.max),
            step: asNumber(descriptor.step),
            options: Array.isArray(descriptor.options)
              ? descriptor.options.map(function (option) {
                  if (option && typeof option === 'object') return { value: plainText(option.value), label: plainText(option.label || option.value) }
                  return { value: plainText(option), label: plainText(option) }
                })
              : null,
          })
        }
      }
      out.push({
        id: plainText(item.id || ''),
        title: plainText(item.title || item.id || ''),
        group: plainText(item.group || '其它'),
        method: plainText(item.method || ''),
        needs: Array.isArray(item.needs) ? item.needs.map(function (need) { return plainText(need) }) : [],
        params: params,
        available: item.available !== false,
        reason: plainText(item.reason || ''),
      })
    }
    var groups = Array.isArray(raw.groups) ? raw.groups.map(function (item) { return plainText(item) }).filter(Boolean) : []
    if (!groups.length) {
      var seen = {}
      for (var groupIndex = 0; groupIndex < out.length; groupIndex++) {
        if (!seen[out[groupIndex].group]) { seen[out[groupIndex].group] = true; groups.push(out[groupIndex].group) }
      }
    }
    return { analyses: out, groups: groups }
  }

  function normalizePreview(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    return {
      genes: Array.isArray(raw.genes) ? raw.genes.map(function (item) { return plainText(item) }).filter(Boolean) : [],
      sample_columns: Array.isArray(raw.sample_columns)
        ? raw.sample_columns.map(function (item) {
            var column = item && typeof item === 'object' ? item : {}
            return { key: plainText(column.key || column.name || ''), kind: plainText(column.kind || 'categorical'), n_unique: asNumber(column.n_unique) }
          }).filter(function (item) { return !!item.key })
        : [],
      n_samples: asNumber(raw.n_samples),
      n_genes: asNumber(raw.n_genes),
      steps: Array.isArray(raw.steps) ? raw.steps.map(function (item) { return plainText(item) }) : [],
      current: plainText(raw.current || ''),
      modules: Array.isArray(raw.modules) ? raw.modules : [],
    }
  }

  function normalizeAudit(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var checks = []
    var list = Array.isArray(raw.checks) ? raw.checks : []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      checks.push({
        id: plainText(item.id || ('check-' + index)),
        level: item.level === 'fail' ? 'fail' : (item.level === 'warn' ? 'warn' : 'ok'),
        title: plainText(item.title || item.id || ''),
        detail: plainText(item.detail || ''),
        hint: plainText(item.hint || ''),
      })
    }
    return {
      score: isFinite(Number(raw.score)) ? Math.round(Number(raw.score)) : 0,
      level: raw.level === 'fail' ? 'fail' : (raw.level === 'warn' ? 'warn' : 'ok'),
      checks: checks,
    }
  }

  function normalizeJob(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var job = raw.job && typeof raw.job === 'object' ? raw.job : raw
    return {
      id: plainText(job.id || ''),
      kind: plainText(job.kind || ''),
      status: plainText(job.status || ''),
      pct: clamp(asNumber(job.pct) === null ? 0 : Number(job.pct), 0, 100),
      message: plainText(job.message || job.detail || ''),
      result: job.result && typeof job.result === 'object' ? job.result : null,
      error: plainText(job.error || ''),
    }
  }

  /**
   * /jobs 响应状态机：
   *   running | done | failed | timeout（timeout 由调用方传 elapsedMs/timeoutMs 判定）
   */
  function jobPolling(response, options) {
    var job = normalizeJob(response)
    var opts = options && typeof options === 'object' ? options : {}
    var elapsedMs = asNumber(opts.elapsedMs)
    if (elapsedMs === null || elapsedMs < 0) elapsedMs = 0
    var timeoutMs = asNumber(opts.timeoutMs)
    if (timeoutMs === null || timeoutMs <= 0) timeoutMs = POLL_TIMEOUT_MS
    var done = job.status === 'done' || job.status === 'success' || job.status === 'completed'
    var failed = job.status === 'failed' || job.status === 'error'
    var timedOut = !done && !failed && elapsedMs >= timeoutMs
    var status = done ? 'done' : (failed ? 'failed' : (timedOut ? 'timeout' : 'running'))
    return {
      id: job.id,
      kind: job.kind,
      status: status,
      pct: job.pct,
      message: job.message,
      result: job.result,
      error: job.error,
      done: done,
      failed: failed,
      timedOut: timedOut,
      running: status === 'running',
    }
  }

  function jobDone(response) {
    return jobPolling(response).done
  }

  // ==========================================================================
  // API 客户端：全部走同源相对路径；fetch 可注入（单测用假 fetch）。
  // ==========================================================================

  function createApi(fetchImpl, baseUrl) {
    var base = typeof baseUrl === 'string' && baseUrl ? baseUrl.replace(/\/$/, '') : API

    function resolveFetch() {
      if (typeof fetchImpl === 'function') return fetchImpl
      if (typeof globalThis.fetch === 'function') return globalThis.fetch.bind(globalThis)
      return null
    }

    function apiError(message, data, status) {
      var error = new Error(message)
      error.__bgApiError = true
      error.code = ''
      if (data && typeof data === 'object') {
        error.data = data
        if (typeof data.code === 'string' && data.code) error.code = data.code
        if (Array.isArray(data.candidates)) error.candidates = data.candidates
      }
      if (isFinite(Number(status)) && Number(status) > 0) error.status = Number(status)
      return error
    }

    function pickError(data) {
      if (!data || typeof data !== 'object') return ''
      if (typeof data.error === 'string' && data.error) return data.error
      if (typeof data.message === 'string' && data.message) return data.message
      return ''
    }

    function parseResponse(response) {
      if (!response) return Promise.reject(apiError('Host 返回了空响应'))
      var ok = response.ok !== false
      var status = Number(response.status) || 0
      var textPromise = typeof response.text === 'function'
        ? Promise.resolve().then(function () { return response.text() }).catch(function () { return '' })
        : Promise.resolve('')
      return textPromise.then(function (text) {
        var data = null
        if (text) {
          try { data = JSON.parse(text) } catch (ignored) { data = null }
        }
        if (!ok) throw apiError(pickError(data) || ('HTTP ' + status), data, status)
        if (data && data.ok === false) throw apiError(pickError(data) || '请求失败', data, status)
        return data && typeof data === 'object' ? data : {}
      })
    }

    function request(path, settings, extra) {
      var doFetch = resolveFetch()
      if (!doFetch) return Promise.reject(apiError('当前环境不支持 fetch'))
      var options = Object.assign({}, settings || {})
      var extraOptions = extra && typeof extra === 'object' ? extra : {}
      var timeoutMs = isFinite(Number(extraOptions.timeoutMs)) && Number(extraOptions.timeoutMs) > 0
        ? Number(extraOptions.timeoutMs)
        : DEFAULT_API_TIMEOUT_MS
      var controller = null
      var timer = null
      if (extraOptions.signal) {
        options.signal = extraOptions.signal
      } else if (typeof AbortController === 'function' && timeoutMs > 0) {
        controller = new AbortController()
        options.signal = controller.signal
        timer = setTimeout(function () {
          try { controller.abort() } catch (ignored) { /* ignore */ }
        }, timeoutMs)
      }
      function clear() {
        if (timer) { clearTimeout(timer); timer = null }
      }
      var promise = Promise.resolve()
        .then(function () { return doFetch(base + path, options) })
        .then(parseResponse)
      return promise.then(function (value) {
        clear()
        return value
      }, function (error) {
        clear()
        if (error && error.__bgApiError) throw error
        var message = formatError(error)
        if (error && (error.name === 'AbortError' || /aborted/i.test(message))) {
          if (extraOptions.signal) throw error
          throw apiError('请求超时（' + Math.round(timeoutMs / 1000) + ' 秒）：引擎可能还在计算，请稍后刷新重试。')
        }
        throw apiError('网络请求失败：' + (message || '无法连接 Host'))
      })
    }

    function jsonPost(path, body, extra) {
      return request(path, {
        method: 'POST',
        headers: { 'content-type': 'application/json' },
        body: JSON.stringify(body || {}),
      }, extra)
    }

    function withActor(body) {
      return Object.assign({ actor: 'user' }, body && typeof body === 'object' ? body : {})
    }

    return {
      base: base,
      request: request,
      getHealth: function (extra) { return request('/health', { method: 'GET' }, extra) },
      getState: function (extra) { return request('/state', { method: 'GET' }, extra) },
      getActivity: function (extra) { return request('/activity', { method: 'GET' }, extra) },
      getFormats: function (extra) { return request('/formats', { method: 'GET' }, extra) },
      geoSearch: function (filters, extra) {
        return jsonPost('/geo/search', buildSearchPayload(filters), Object.assign({ timeoutMs: 90000 }, extra || {}))
      },
      geoDetail: function (gse, extra) {
        return jsonPost('/geo/detail', withActor({ gse: plainText(gse).toUpperCase() }), Object.assign({ timeoutMs: 90000 }, extra || {}))
      },
      geoDownload: function (options, extra) {
        var body = {
          gse: plainText(options && options.gse).toUpperCase(),
          files: Array.isArray(options && options.files) ? options.files.slice() : [],
          dest_dir: plainText(options && options.dest_dir) || DEFAULT_GEO_DIR,
          overwrite: options && options.overwrite === true,
        }
        if (options && options.open_after === true) body.open_after = true
        return jsonPost('/geo/download', withActor(body), Object.assign({ timeoutMs: 60000 }, extra || {}))
      },
      getJob: function (id, extra) {
        return request('/jobs?id=' + encodeURIComponent(plainText(id)), { method: 'GET' }, extra)
      },
      openDataset: function (options, extra) {
        return jsonPost('/dataset/open', withActor(buildImportPayload(options)),
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      demo: function (extra) { return jsonPost('/dataset/demo', withActor({}), extra) },
      selectDataset: function (dataset, extra) {
        return jsonPost('/dataset/select', withActor({ dataset: plainText(dataset) }), extra)
      },
      closeDataset: function (dataset, extra) {
        return jsonPost('/dataset/close', withActor({ dataset: plainText(dataset) }), extra)
      },
      uploadDataset: function (filename, data, extra) {
        return request('/dataset/upload', {
          method: 'POST',
          headers: {
            'content-type': 'application/octet-stream',
            'x-dsh-filename': encodeURIComponent(plainText(filename) || 'upload.bin'),
          },
          body: data,
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      inspectSource: function (options, extra) {
        return jsonPost('/inspect', withActor(buildImportPayload(options)),
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      preview: function (dataset, extra) {
        return jsonPost('/preview', withActor({ dataset: plainText(dataset) }), extra)
      },
      samples: function (payload, extra) {
        var body = withActor(payload && typeof payload === 'object' ? payload : {})
        return jsonPost('/samples', body, extra)
      },
      runStep: function (dataset, step, params, extra) {
        return jsonPost('/run', withActor({
          dataset: plainText(dataset),
          step: plainText(step),
          params: params && typeof params === 'object' ? params : {},
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      plot: function (dataset, kind, params, extra) {
        return jsonPost('/plot', withActor({
          dataset: plainText(dataset),
          kind: plainText(kind),
          params: params && typeof params === 'object' ? params : {},
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      analyses: function (dataset, extra) {
        var query = dataset ? '?dataset=' + encodeURIComponent(plainText(dataset)) : ''
        return request('/analyses' + query, { method: 'GET' }, extra)
      },
      analysis: function (dataset, id, params, extra) {
        return jsonPost('/analysis', withActor({
          dataset: plainText(dataset),
          id: plainText(id),
          params: params && typeof params === 'object' ? params : {},
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS * 2 }, extra || {}))
      },
      attach: function (dataset, module, path, extra) {
        return jsonPost('/attach', withActor({
          dataset: plainText(dataset),
          module: plainText(module),
          path: plainText(path),
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      exportData: function (dataset, what, path, extra) {
        var body = { dataset: plainText(dataset), what: plainText(what) }
        if (path) body.path = plainText(path)
        return jsonPost('/export', withActor(body), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      audit: function (dataset, extra) {
        return jsonPost('/audit', withActor({ dataset: plainText(dataset) }), extra)
      },
      suggest: function (dataset, extra) {
        return jsonPost('/suggest', withActor({ dataset: plainText(dataset) }), extra)
      },
      recipe: function (dataset, id, extra) {
        var body = {}
        if (dataset) body.dataset = plainText(dataset)
        if (id) body.id = plainText(id)
        return jsonPost('/recipe', withActor(body), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      scRefs: function (extra) { return request('/sc/refs', { method: 'GET' }, extra) },
      scAttach: function (payload, extra) {
        return jsonPost('/sc/attach', withActor(payload && typeof payload === 'object' ? payload : {}),
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      reset: function (extra) { return jsonPost('/reset', withActor({}), extra) },
    }
  }

  function buildImportPayload(options) {
    var raw = options && typeof options === 'object' ? options : {}
    var out = {}
    var path = plainText(raw.path)
    var url = plainText(raw.url)
    if (path) out.path = path
    if (url) out.url = url
    var name = plainText(raw.name)
    if (name) out.name = name
    var format = plainText(raw.data_format)
    if (format && format !== 'auto') out.data_format = format
    var member = plainText(raw.member)
    if (member) out.member = member
    var object = plainText(raw.object)
    if (object) out.object = object
    var sheet = plainText(raw.sheet)
    if (sheet) out.sheet = sheet
    var geneCol = plainText(raw.gene_col)
    if (geneCol) out.gene_col = geneCol
    var sampleSheet = plainText(raw.sample_sheet)
    if (sampleSheet) out.sample_sheet = sampleSheet
    var clinical = plainText(raw.clinical)
    if (clinical) out.clinical = clinical
    var species = plainText(raw.species)
    if (species && species !== 'auto') out.species = species
    if (raw.transpose === true) out.transpose = true
    return out
  }

  // ==========================================================================
  // 推荐/表格小工具
  // ==========================================================================

  /** 基因搜索：前缀优先、包含次之，最多 limit 个。 */
  function matchGenes(allGenes, query, limit) {
    var source = Array.isArray(allGenes) ? allGenes : []
    var needle = plainText(query).toUpperCase()
    var max = clamp(isFinite(Number(limit)) ? Number(limit) : 50, 1, 500)
    var prefix = []
    var contains = []
    for (var index = 0; index < source.length; index++) {
      var gene = String(source[index] || '')
      var upper = gene.toUpperCase()
      if (!needle || upper === needle) { prefix.push(gene); continue }
      if (upper.indexOf(needle) === 0) prefix.push(gene)
      else if (upper.indexOf(needle) > 0) contains.push(gene)
      if (prefix.length >= max) break
    }
    return prefix.concat(contains).slice(0, max)
  }

  /** 表格行排序（稳定；数字列按数值）。 */
  function sortRows(rows, key, direction, columns) {
    var list = Array.isArray(rows) ? rows.slice() : []
    var numeric = false
    if (Array.isArray(columns)) {
      for (var index = 0; index < columns.length; index++) {
        if (columns[index] && columns[index].key === key && columns[index].type === 'number') numeric = true
      }
    }
    var factor = direction === 'desc' ? -1 : 1
    return list.sort(function (a, b) {
      var left = a && typeof a === 'object' ? a[key] : null
      var right = b && typeof b === 'object' ? b[key] : null
      if (isNil(left) && isNil(right)) return 0
      if (isNil(left)) return 1
      if (isNil(right)) return -1
      if (numeric && typeof left === 'number' && typeof right === 'number') return (left - right) * factor
      return String(left).localeCompare(String(right)) * factor
    })
  }

  /** 表头列 → 展示列（rows 里可能缺列）。 */
  function tableColumns(table) {
    return table && Array.isArray(table.columns) ? table.columns : []
  }

  function formatCell(value, numeric) {
    if (isNil(value) || value === '') return '—'
    if (typeof value === 'number') {
      if (!numeric) return String(value)
      var rounded = Math.abs(value) >= 1000 ? value : Number(value.toPrecision(4))
      return String(rounded)
    }
    if (typeof value === 'boolean') return value ? '是' : '否'
    if (Array.isArray(value)) return value.join(', ')
    return String(value)
  }

  /** 结果表 → 渲染行（上限 MAX_TABLE_RENDER_ROWS）。 */
  function formatTableRows(table, limit) {
    var columns = tableColumns(table)
    var rows = table && Array.isArray(table.rows) ? table.rows : []
    var max = clamp(isFinite(Number(limit)) ? Number(limit) : MAX_TABLE_RENDER_ROWS, 1, 5000)
    var shown = rows.slice(0, max)
    var out = []
    for (var index = 0; index < shown.length; index++) {
      var row = shown[index] && typeof shown[index] === 'object' ? shown[index] : {}
      var cells = []
      for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) {
        var column = columns[columnIndex]
        cells.push({
          key: column.key,
          text: formatCell(row[column.key], column.type === 'number'),
          numeric: column.type === 'number',
        })
      }
      out.push({ key: 'row-' + index, cells: cells })
    }
    return {
      columns: columns,
      rows: out,
      total: rows.length,
      shown: out.length,
      truncated: rows.length > out.length || (table && table.truncated === true),
    }
  }

  function activityTone(entry) {
    return entry && entry.ok === false ? 'bad' : 'ok'
  }

  function historyTone(entry) {
    return entry && entry.ok === false ? 'bad' : 'ok'
  }

  /**
   * PlotBlock 的缓存键：kind + 稳定序列化的 params。
   * 键排序、跳过空值，保证同一组参数永远得到同一个 key（避免重复出图）。
   */
  function plotKey(kind, params) {
    var name = plainText(kind)
    var source = params && typeof params === 'object' && !Array.isArray(params) ? params : {}
    var keys = Object.keys(source).sort()
    var parts = []
    for (var index = 0; index < keys.length; index++) {
      var key = keys[index]
      var value = source[key]
      if (isNil(value) || value === '') continue
      if (Array.isArray(value)) {
        var joined = value.map(function (item) { return plainText(item) }).filter(Boolean).join(',')
        if (joined) parts.push(key + '=' + joined)
        continue
      }
      if (typeof value === 'object') {
        try { parts.push(key + '=' + JSON.stringify(value)) } catch (ignored) { /* 不可序列化就跳过 */ }
        continue
      }
      parts.push(key + '=' + String(value))
    }
    return name + (parts.length ? '|' + parts.join('&') : '')
  }

  /**
   * 宿主 locale 服务 → 取词函数。宿主不可用时静默回落 DICT_ZH，
   * 保证任何环境（含单测）下 T() 都不会返回裸 key 或抛异常。
   */
  function createTranslator(locale) {
    var bound = null
    try {
      if (locale && typeof locale.bind === 'function') bound = locale.bind(NS)
    } catch (ignored) {
      bound = null
    }
    return function translate(key) {
      if (typeof bound === 'function') {
        try {
          var value = bound(key)
          if (typeof value === 'string' && value && value !== key) return value
        } catch (ignored) { /* 宿主 locale 不可用时静默回落 */ }
      }
      var fallback = DICT_ZH[key]
      return typeof fallback === 'string' ? fallback : String(key)
    }
  }

  // ==========================================================================
  // 纯函数测试出口。node 环境（typeof window === 'undefined'）挂到
  // globalThis.__BG_TEST__；浏览器里也会挂进 factory 返回的 __test__，便于联调。
  // ==========================================================================
  var BG_TEST = {
    NS: NS,
    API: API,
    POLL_MS: POLL_MS,
    POLL_TIMEOUT_MS: POLL_TIMEOUT_MS,
    DEFAULT_GEO_DIR: DEFAULT_GEO_DIR,
    MAX_HISTORY_ROWS: MAX_HISTORY_ROWS,
    STEPS: STEPS,
    STEP_BY_ID: STEP_BY_ID,
    TAB_KEYS: TAB_KEYS,
    TAB_META: TAB_META,
    PALETTE: PALETTE,
    CONTINUOUS_STOPS: CONTINUOUS_STOPS,
    VOLCANO_COLORS: VOLCANO_COLORS,
    AUDIT_LEVEL_STYLE: AUDIT_LEVEL_STYLE,
    DE_COLUMNS: DE_COLUMNS,
    ENRICH_DBS: ENRICH_DBS,
    IMMUNE_METHODS: IMMUNE_METHODS,
    NCRNA_MODULES: NCRNA_MODULES,
    EXPORT_KINDS: EXPORT_KINDS,
    QUICK_EXPORTS: QUICK_EXPORTS,
    DICT_ZH: DICT_ZH,
    DICT_EN: DICT_EN,
    clamp: clamp,
    asNumber: asNumber,
    plainText: plainText,
    fmtNumber: fmtNumber,
    fmtPct: fmtPct,
    fmtP: fmtP,
    fmtBytes: formatBytes,
    plotKey: plotKey,
    createTranslator: createTranslator,
    formatBytes: formatBytes,
    formatDuration: formatDuration,
    formatClock: formatClock,
    formatError: formatError,
    hexToRgb: hexToRgb,
    colorForGroup: colorForGroup,
    groupBy: groupBy,
    volcanoPoints: volcanoPoints,
    pcaScreen: pcaScreen,
    buildSearchPayload: buildSearchPayload,
    buildImportPayload: buildImportPayload,
    suggestionInstruction: suggestionInstruction,
    defaultParams: defaultParams,
    collectStepParams: collectStepParams,
    parseTokenList: parseTokenList,
    formatParamValue: formatParamValue,
    doneStepSet: doneStepSet,
    isStepReady: isStepReady,
    stepStatusOf: stepStatusOf,
    stepStatusKey: stepStatusKey,
    stepById: stepById,
    figureSrc: figureSrc,
    normalizeFigure: normalizeFigure,
    normalizeTable: normalizeTable,
    normalizeResult: normalizeResult,
    normalizeDataset: normalizeDataset,
    normalizeHealth: normalizeHealth,
    normalizeState: normalizeState,
    normalizeSuggestions: normalizeSuggestions,
    normalizeActivityEntry: normalizeActivityEntry,
    normalizeAnalyses: normalizeAnalyses,
    normalizePreview: normalizePreview,
    normalizeAudit: normalizeAudit,
    normalizeJob: normalizeJob,
    normalizeStringList: normalizeStringList,
    engineBadge: engineBadge,
    jobPolling: jobPolling,
    jobDone: jobDone,
    apiClient: createApi,
    createApi: createApi,
    matchGenes: matchGenes,
    sortRows: sortRows,
    formatTableRows: formatTableRows,
    tableColumns: tableColumns,
    formatCell: formatCell,
    activityTone: activityTone,
    historyTone: historyTone,
  }

  // ==========================================================================
  // 浏览器层：React 组件 + 样式 + 槽注册。只 require('react')。
  // ==========================================================================

  function factory(require) {
    var React = require('react')
    var useState = React.useState
    var useEffect = React.useEffect
    var useRef = React.useRef
    var useMemo = React.useMemo
    var useCallback = React.useCallback

    /** React.createElement 简写（不写 JSX，避免宿主构建链依赖）。 */
    function h(type, props) {
      var children = Array.prototype.slice.call(arguments, 2)
      return React.createElement.apply(React, [type, props].concat(children))
    }

    /** t()：apply() 里用宿主 locale 服务重新绑定；默认回落中文词典。 */
    var T = function (key) {
      var fallback = DICT_ZH[key]
      return typeof fallback === 'string' ? fallback : String(key)
    }

    function warnPlugin(message, error) {
      try { console.warn('[dsh-geobulk] ' + message, error || '') } catch (ignored) { /* ignore */ }
    }

    /** 模板取词：字典里用 `{v}` 占位，避免依赖宿主 t() 的插值实现。 */
    function tf(translate, key, value) {
      var raw = typeof translate === 'function' ? translate(key) : ''
      var template = typeof raw === 'string' && raw && raw !== key ? raw : (DICT_ZH[key] || '')
      return String(template).replace(/\{v\}/g, isNil(value) ? '' : String(value))
    }

    // ==========================================================================
    // 样式：一个 CSS 字符串，由 <style> 渲染在组件内部（生命周期跟随组件）。
    // 颜色只用 --dsw-alias-* 主题变量，派生色统一走 color-mix，深浅色自适应。
    // 类名一律 bg- 前缀。
    // ==========================================================================
    var CSS = [
      '.bg-root{display:flex;flex-direction:column;height:100%;min-height:72vh;overflow:hidden;position:relative;isolation:isolate;color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-base);border-radius:12px;font-size:13px;line-height:1.5;font-variant-numeric:tabular-nums}',
      '.bg-root *,.bg-root *::before,.bg-root *::after{box-sizing:border-box}',
      '.bg-root button,.bg-root input,.bg-root select,.bg-root textarea{font-family:inherit;font-variant-numeric:tabular-nums}',
      '.bg-root ::selection{background:color-mix(in srgb,var(--dsw-alias-brand-primary) 24%,transparent)}',

      /* 顶栏 + 1px 渐变线 */
      '.bg-top{position:relative;display:flex;align-items:center;gap:10px;flex-wrap:wrap;flex:none;min-height:52px;padding:9px 14px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.bg-top::after{content:"";position:absolute;left:0;right:0;bottom:-1px;height:1px;background:linear-gradient(90deg,transparent,color-mix(in srgb,var(--dsw-alias-brand-primary) 55%,transparent) 26%,color-mix(in srgb,var(--dsw-alias-state-success-primary) 45%,transparent) 58%,transparent);pointer-events:none}',
      '.bg-brand{display:flex;align-items:center;gap:9px;min-width:0}',
      '.bg-brand-mark{display:inline-flex;align-items:center;justify-content:center;flex:none;width:28px;height:28px;border-radius:9px;border:1px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 34%,var(--dsw-alias-border-l1));background:linear-gradient(135deg,color-mix(in srgb,var(--dsw-alias-brand-primary) 16%,var(--dsw-alias-bg-layer-1)),var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.bg-title-block{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.bg-title{font-size:14px;font-weight:600;letter-spacing:.01em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-subtitle{font-size:11px;color:var(--dsw-alias-label-secondary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-badge{display:inline-flex;align-items:center;gap:6px;height:22px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:11px;line-height:1;cursor:default;white-space:nowrap}',
      '.bg-badge-dot{width:6px;height:6px;border-radius:50%;flex:none;background:currentColor;box-shadow:0 0 0 2px color-mix(in srgb,currentColor 18%,transparent)}',
      '.bg-badge-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-success-primary)}',
      '.bg-badge-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 42%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-warn-primary)}',
      '.bg-badge-bad{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 42%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-error-primary)}',
      '.bg-chips{display:flex;align-items:center;gap:5px;flex-wrap:wrap;min-width:0;flex:1 1 220px}',
      '.bg-chip{display:inline-flex;align-items:center;gap:5px;max-width:200px;height:23px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font:inherit;font-size:11.5px;cursor:pointer;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;transition:border-color .14s ease,background .14s ease,color .14s ease}',
      '.bg-chip:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-primary)}',
      '.bg-chip-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.bg-chip-x{display:inline-flex;align-items:center;justify-content:center;width:14px;height:14px;margin-right:-2px;border-radius:50%;font-size:10px;color:var(--dsw-alias-label-tertiary)}',
      '.bg-chip-x:hover{background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 18%,transparent);color:var(--dsw-alias-state-error-primary)}',
      '.bg-top-actions{display:flex;align-items:center;gap:6px;flex-wrap:wrap;margin-left:auto;justify-content:flex-end}',
      '.bg-top-sep{width:1px;height:18px;background:var(--dsw-alias-border-l1);flex:none}',

      /* 按钮 / 输入 */
      '.bg-btn{display:inline-flex;align-items:center;justify-content:center;gap:5px;height:28px;padding:0 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;line-height:1;cursor:pointer;white-space:nowrap;transition:border-color .14s ease,background .14s ease,opacity .14s ease,box-shadow .14s ease}',
      '.bg-btn:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-overlay)}',
      '.bg-btn:active:not(:disabled){transform:translateY(.5px)}',
      '.bg-btn:disabled{opacity:.45;cursor:default}',
      '.bg-btn-primary{background:var(--dsw-alias-brand-primary);border-color:var(--dsw-alias-brand-primary);color:var(--dsw-alias-bg-base);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-brand-primary) 32%,transparent)}',
      '.bg-btn-primary:hover:not(:disabled){opacity:.9;background:var(--dsw-alias-brand-primary);border-color:var(--dsw-alias-brand-primary)}',
      '.bg-btn-ghost{background:transparent;border-color:transparent;color:var(--dsw-alias-label-secondary)}',
      '.bg-btn-ghost:hover:not(:disabled){background:var(--dsw-alias-bg-layer-2);border-color:var(--dsw-alias-border-l1);color:var(--dsw-alias-label-primary)}',
      '.bg-btn-sm{height:24px;padding:0 8px;font-size:11.5px;border-radius:6px}',
      '.bg-btn-danger{color:var(--dsw-alias-state-error-primary);border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 36%,var(--dsw-alias-border-l1))}',
      '.bg-input{height:28px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;outline:none;min-width:0;transition:border-color .14s ease,box-shadow .14s ease}',
      '.bg-input:focus{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 48%,var(--dsw-alias-border-l1));box-shadow:0 0 0 2px color-mix(in srgb,var(--dsw-alias-brand-primary) 12%,transparent)}',
      '.bg-input:disabled{opacity:.55}',
      '.bg-input::placeholder{color:var(--dsw-alias-label-tertiary)}',
      '.bg-select{cursor:pointer;padding-right:5px}',
      '.bg-textarea{min-height:54px;padding:6px 9px;resize:vertical;line-height:1.45}',
      '.bg-check{display:inline-flex;align-items:center;gap:6px;font-size:12px;color:var(--dsw-alias-label-secondary);cursor:pointer;user-select:none}',
      '.bg-check input{accent-color:var(--dsw-alias-brand-primary);width:14px;height:14px;cursor:pointer}',
      '.bg-form{display:flex;flex-direction:column;gap:7px}',
      '.bg-field{display:flex;align-items:center;justify-content:space-between;gap:8px;min-width:0}',
      '.bg-field-col{flex-direction:column;align-items:stretch;gap:4px}',
      '.bg-field-label{flex:0 1 auto;min-width:0;font-size:11.5px;color:var(--dsw-alias-label-secondary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-field .bg-input{flex:1 1 auto;max-width:150px;height:25px;font-size:11.5px}',
      '.bg-field-col .bg-input{max-width:none;width:100%}',
      '.bg-multi{display:flex;flex-wrap:wrap;gap:4px}',
      '.bg-multi-opt{display:inline-flex;align-items:center;gap:4px;height:22px;padding:0 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:6px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font:inherit;font-size:11px;cursor:pointer;transition:border-color .14s ease,background .14s ease,color .14s ease}',
      '.bg-multi-opt:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-primary)}',
      '.bg-multi-opt-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 46%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',

      /* 三栏骨架 */
      '.bg-body{display:flex;flex:1 1 auto;min-height:0;overflow:hidden}',
      '.bg-left{display:flex;flex-direction:column;flex:0 0 268px;width:268px;min-height:0;border-right:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.bg-center{display:flex;flex-direction:column;flex:1 1 auto;min-width:0;min-height:0;background:var(--dsw-alias-bg-base)}',
      '.bg-right{display:flex;flex-direction:column;flex:0 0 300px;width:300px;min-height:0;overflow-y:auto;border-left:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.bg-panel-head{display:flex;align-items:center;justify-content:space-between;gap:8px;flex:none;padding:9px 12px;border-bottom:1px solid var(--dsw-alias-border-l1);font-size:11px;font-weight:600;letter-spacing:.04em;color:var(--dsw-alias-label-secondary)}',
      '.bg-panel-note{font-weight:400;font-size:10.5px;letter-spacing:0;color:var(--dsw-alias-label-tertiary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-banner{display:flex;align-items:center;justify-content:space-between;gap:10px;flex:none;padding:6px 14px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-state-warn-primary) 32%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 7%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-label-primary);font-size:11.5px;animation:bg-fade .16s ease}',
      '.bg-banner-bad{border-bottom-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 34%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 6%,var(--dsw-alias-bg-layer-1))}',
      '.bg-banner-text{flex:1 1 auto;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.bg-panel-hint{padding:7px 12px 0;font-size:11px;color:var(--dsw-alias-label-tertiary);line-height:1.45}',
      '.bg-panel-body{display:flex;flex-direction:column;gap:9px;padding:10px 12px 16px;overflow-y:auto;min-height:0}',

      /* 左栏：数据集卡 + 流程卡 */
      '.bg-dataset-card{display:flex;flex-direction:column;gap:5px;padding:10px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:linear-gradient(160deg,color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,var(--dsw-alias-bg-base)),var(--dsw-alias-bg-base) 58%);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-label-primary) 5%,transparent)}',
      '.bg-dataset-name{font-size:12.5px;font-weight:600;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-dataset-meta{display:flex;flex-wrap:wrap;gap:4px 8px;font-size:11px;color:var(--dsw-alias-label-tertiary)}',
      '.bg-dataset-meta b{font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.bg-steps{display:flex;flex-direction:column;gap:8px;padding:10px 12px 16px;overflow-y:auto;min-height:0}',
      '.bg-step{display:flex;flex-direction:column;gap:5px;padding:8px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:var(--dsw-alias-bg-base);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-label-primary) 4%,transparent);transition:border-color .16s ease,background .16s ease,opacity .16s ease}',
      '.bg-step-done{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 28%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 4%,var(--dsw-alias-bg-base))}',
      '.bg-step-running{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 48%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,var(--dsw-alias-bg-base))}',
      '.bg-step-blocked{opacity:.8}',
      '.bg-step-head{display:flex;align-items:center;gap:7px;min-width:0}',
      '.bg-step-dot{position:relative;display:inline-flex;align-items:center;justify-content:center;flex:none;width:20px;height:20px;border-radius:50%;border:1px solid var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:10.5px;font-weight:600}',
      '.bg-step-dot::after{content:"";position:absolute;right:-1px;bottom:-1px;width:6px;height:6px;border-radius:50%;background:currentColor;opacity:0;box-shadow:0 0 0 2px color-mix(in srgb,currentColor 20%,transparent)}',
      '.bg-step-dot-done{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 55%,transparent);background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 14%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-success-primary)}',
      '.bg-step-dot-done::after{opacity:1}',
      '.bg-step-dot-running{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 58%,transparent);background:color-mix(in srgb,var(--dsw-alias-brand-primary) 14%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary);animation:bg-pulse 1.6s ease-in-out infinite}',
      '.bg-step-dot-pending::after{opacity:0}',
      '.bg-step-dot-blocked{border-color:var(--dsw-alias-border-l1);color:var(--dsw-alias-label-tertiary)}',
      '.bg-step-title{flex:1 1 auto;min-width:0;font-size:12.5px;font-weight:500;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-step-status{flex:none;font-size:10.5px;color:var(--dsw-alias-label-tertiary);white-space:nowrap}',
      '.bg-status-done{color:var(--dsw-alias-state-success-primary)}',
      '.bg-status-running{color:var(--dsw-alias-brand-primary)}',
      '.bg-status-failed{color:var(--dsw-alias-state-error-primary)}',
      '.bg-step-caret{flex:none;display:inline-flex;align-items:center;justify-content:center;width:18px;height:18px;padding:0;border:none;border-radius:6px;background:transparent;color:var(--dsw-alias-label-tertiary);cursor:pointer;transition:background .14s ease,color .14s ease}',
      '.bg-step-caret:hover{background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary)}',
      '.bg-step-params{display:flex;flex-direction:column;gap:6px;padding:7px 1px 2px;border-top:1px dashed var(--dsw-alias-border-l1);margin-top:2px}',
      '.bg-step-foot{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.bg-step-summary{flex:1 1 auto;min-width:0;font-size:11px;color:var(--dsw-alias-label-tertiary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-step-summary-full{white-space:normal;line-height:1.4;color:var(--dsw-alias-label-secondary)}',
      '.bg-step .bg-inline-err{margin-top:1px}',

      /* 中栏 tabs */
      '.bg-tabs{display:flex;align-items:center;gap:2px;flex:none;padding:6px 10px 0;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1);overflow-x:auto;scrollbar-width:thin}',
      '.bg-tab{display:inline-flex;align-items:center;justify-content:center;gap:5px;height:30px;padding:0 11px;border:none;border-bottom:2px solid transparent;border-radius:7px 7px 0 0;background:transparent;color:var(--dsw-alias-label-secondary);font:inherit;font-size:12px;cursor:pointer;white-space:nowrap;transition:color .14s ease,background .14s ease,border-color .14s ease}',
      '.bg-tab:hover{color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-layer-2)}',
      '.bg-tab-on{color:var(--dsw-alias-brand-primary);border-bottom-color:var(--dsw-alias-brand-primary);font-weight:600}',
      '.bg-pane{display:flex;flex:1 1 auto;flex-direction:column;gap:11px;min-height:0;padding:12px 14px 18px;overflow-y:auto}',
      '.bg-toolbar{display:flex;align-items:center;gap:8px;flex-wrap:wrap;flex:none;min-height:28px}',
      '.bg-toolbar-label{font-size:11px;color:var(--dsw-alias-label-tertiary);white-space:nowrap}',

      /* 卡片 / 指标 / 区块 */
      '.bg-card{display:flex;flex-direction:column;gap:8px;padding:11px 13px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:var(--dsw-alias-bg-layer-1);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-label-primary) 4%,transparent)}',
      '.bg-card-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-secondary);letter-spacing:.02em}',
      '.bg-card-title small{font-weight:400;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.bg-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(108px,1fr));gap:8px}',
      '.bg-metric{display:flex;flex-direction:column;gap:1px;padding:8px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.bg-metric-k{font-size:10.5px;color:var(--dsw-alias-label-tertiary);letter-spacing:.03em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-metric-v{font-size:20px;font-weight:300;line-height:1.2;color:var(--dsw-alias-label-primary);font-variant-numeric:tabular-nums}',
      '.bg-metric-v small{font-size:11px;font-weight:500;color:var(--dsw-alias-label-secondary);margin-left:3px}',
      '.bg-metric-unit{font-size:11px;color:var(--dsw-alias-label-tertiary);font-weight:400}',
      '.bg-grid-2{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:11px;align-items:start}',
      '.bg-section-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11px;font-weight:600;letter-spacing:.04em;color:var(--dsw-alias-label-secondary);margin-top:2px}',
      '.bg-muted{font-size:11.5px;color:var(--dsw-alias-label-secondary);line-height:1.5}',
      '.bg-inline-row{display:flex;align-items:center;gap:6px;flex-wrap:wrap;min-width:0}',
      '.bg-inline-err{padding:5px 8px;border-radius:6px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,transparent);background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 8%,transparent);color:var(--dsw-alias-state-error-primary);font-size:11.5px;word-break:break-word}',
      '.bg-warn-text{color:var(--dsw-alias-state-warn-primary);font-size:11.5px}',
      '.bg-ok-text{color:var(--dsw-alias-state-success-primary);font-size:11.5px}',
      '.bg-divider{flex:none;height:1px;margin:2px 0;background:var(--dsw-alias-border-l1)}',
      '.bg-path{font-size:10.5px;color:var(--dsw-alias-label-tertiary);word-break:break-all;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}',

      /* 图 */
      '.bg-figure{display:flex;flex-direction:column;gap:6px;min-width:0}',
      '.bg-figure-title{font-size:11px;color:var(--dsw-alias-label-secondary);font-weight:500}',
      '.bg-figure img{display:block;max-width:100%;height:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-2)}',
      '.bg-figure-meta{font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.bg-figure-empty{display:flex;align-items:center;justify-content:center;min-height:120px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-tertiary);font-size:11.5px;padding:12px;text-align:center}',
      '.bg-figure img.bg-figure-broken{display:none}',

      /* 表格 */
      '.bg-table-wrap{max-height:420px;overflow:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.bg-table{width:100%;border-collapse:separate;border-spacing:0;font-size:11.5px}',
      '.bg-table th{position:sticky;top:0;z-index:1;padding:6px 10px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-weight:600;text-align:left;white-space:nowrap}',
      '.bg-table th.bg-th-sort{cursor:pointer;user-select:none}',
      '.bg-table th.bg-th-sort:hover{color:var(--dsw-alias-label-primary)}',
      '.bg-table td{padding:4px 10px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);color:var(--dsw-alias-label-primary);white-space:nowrap;max-width:340px;overflow:hidden;text-overflow:ellipsis}',
      '.bg-table tbody tr:hover td{background:color-mix(in srgb,var(--dsw-alias-brand-primary) 4%,transparent)}',
      '.bg-td-num{text-align:right;font-variant-numeric:tabular-nums}',
      '.bg-table-note{padding:5px 2px 0;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.bg-empty-table{padding:18px 12px;text-align:center;color:var(--dsw-alias-label-tertiary);font-size:11.5px}',

      /* 结果视图 */
      '.bg-result{display:flex;flex-direction:column;gap:10px}',
      '.bg-bullets{margin:0;padding-left:16px;display:flex;flex-direction:column;gap:3px;font-size:12px;color:var(--dsw-alias-label-primary)}',
      '.bg-bullets li::marker{color:var(--dsw-alias-brand-primary)}',
      '.bg-result-block{display:flex;flex-direction:column;gap:6px}',
      '.bg-result-block-title{font-size:11px;font-weight:600;color:var(--dsw-alias-label-secondary);letter-spacing:.03em}',
      '.bg-warn-list{display:flex;flex-direction:column;gap:3px;font-size:11.5px;color:var(--dsw-alias-state-warn-primary)}',
      '.bg-note-list{display:flex;flex-direction:column;gap:3px;font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.bg-file-chip{display:inline-flex;align-items:center;gap:6px;max-width:100%;padding:4px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-2);font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.bg-file-chip .bg-path{flex:1 1 auto;min-width:0}',

      /* 引导 / 空态 / 崩溃 */
      '.bg-guide{display:flex;flex-direction:column;align-items:flex-start;gap:9px;justify-content:center;padding:22px 20px;border:1.5px dashed var(--dsw-alias-border-l2);border-radius:12px;background:var(--dsw-alias-bg-layer-1);text-align:left}',
      '.bg-guide-title{font-size:13.5px;font-weight:600}',
      '.bg-guide-body{font-size:12px;color:var(--dsw-alias-label-secondary);max-width:640px;line-height:1.55}',
      '.bg-guide-hint{font-size:11.5px;color:var(--dsw-alias-label-tertiary);max-width:660px;line-height:1.5}',
      '.bg-notice{display:flex;align-items:center;gap:8px;padding:9px 12px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-size:11.5px}',
      '.bg-crash{display:flex;flex-direction:column;align-items:flex-start;gap:8px;margin:14px;padding:16px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,var(--dsw-alias-border-l1));border-radius:12px;background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 5%,var(--dsw-alias-bg-layer-1))}',
      '.bg-crash-title{font-size:13px;font-weight:600;color:var(--dsw-alias-state-error-primary)}',
      '.bg-crash-hint{font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.bg-crash-detail{font-size:10.5px;color:var(--dsw-alias-label-tertiary);font-family:ui-monospace,SFMono-Regular,Menlo,monospace;word-break:break-all}',
      '.bg-crash-mini{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;border-radius:50%;background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 16%,transparent);color:var(--dsw-alias-state-error-primary);font-size:11px;font-weight:700}',
      '.bg-loading{display:inline-flex;align-items:center;gap:7px;font-size:12px;color:var(--dsw-alias-label-secondary)}',
      '.bg-spinner{display:inline-block;width:13px;height:13px;border-radius:50%;border:1.6px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 30%,transparent);border-top-color:var(--dsw-alias-brand-primary);animation:bg-spin .8s linear infinite;flex:none}',
      '@keyframes bg-spin{to{transform:rotate(360deg)}}',
      '@keyframes bg-pulse{0%,100%{box-shadow:0 0 0 0 color-mix(in srgb,var(--dsw-alias-brand-primary) 30%,transparent)}50%{box-shadow:0 0 0 4px color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,transparent)}}',
      '@keyframes bg-fade{from{opacity:0;transform:translateY(2px)}to{opacity:1;transform:none}}',

      /* 弹层 */
      '.bg-modal-mask{position:fixed;inset:0;z-index:60;display:flex;align-items:center;justify-content:center;padding:24px;background:color-mix(in srgb,var(--dsw-alias-label-primary) 34%,transparent);backdrop-filter:blur(3px);animation:bg-fade .14s ease}',
      '.bg-modal{display:flex;flex-direction:column;width:min(920px,94vw);max-height:88vh;border:1px solid var(--dsw-alias-border-l1);border-radius:14px;background:var(--dsw-alias-bg-base);box-shadow:0 24px 80px color-mix(in srgb,var(--dsw-alias-label-primary) 24%,transparent);overflow:hidden}',
      '.bg-modal-wide{width:min(1080px,95vw)}',
      '.bg-modal-head{display:flex;align-items:flex-start;justify-content:space-between;gap:10px;flex:none;padding:13px 16px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.bg-modal-title{font-size:14px;font-weight:600}',
      '.bg-modal-subtitle{font-size:11.5px;color:var(--dsw-alias-label-secondary);margin-top:1px}',
      '.bg-modal-body{display:flex;flex-direction:column;gap:11px;padding:13px 16px 18px;overflow-y:auto;min-height:0}',
      '.bg-modal-foot{display:flex;align-items:center;justify-content:flex-end;gap:8px;flex:none;padding:10px 16px;border-top:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.bg-modal-input{width:100%;max-width:none}',
      '.bg-advanced{margin-top:1px;padding:9px 12px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.bg-advanced summary{cursor:pointer;font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.bg-advanced[open] summary{margin-bottom:7px;color:var(--dsw-alias-label-primary)}',
      '.bg-form-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:8px 10px}',
      '.bg-form-grid .bg-field{flex-direction:column;align-items:stretch;gap:3px}',
      '.bg-form-grid .bg-input{max-width:none;width:100%}',

      /* GEO 检索 */
      '.bg-geo-search{display:flex;align-items:center;gap:7px;flex-wrap:wrap}',
      '.bg-geo-term{flex:1 1 300px;min-width:220px;height:32px;font-size:12.5px}',
      '.bg-geo-filters{display:flex;align-items:center;gap:10px;flex-wrap:wrap;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.bg-geo-layout{display:grid;grid-template-columns:minmax(260px,340px) 1fr;gap:12px;min-height:0;flex:1 1 auto}',
      '.bg-geo-list{display:flex;flex-direction:column;gap:7px;overflow-y:auto;max-height:52vh;padding-right:2px}',
      '.bg-geo-card{display:flex;flex-direction:column;gap:4px;padding:9px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1);cursor:pointer;text-align:left;font:inherit;transition:border-color .14s ease,background .14s ease,box-shadow .14s ease}',
      '.bg-geo-card:hover{border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2)}',
      '.bg-geo-card-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 50%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,var(--dsw-alias-bg-layer-1));box-shadow:0 0 0 1px color-mix(in srgb,var(--dsw-alias-brand-primary) 24%,transparent)}',
      '.bg-geo-gse{font-size:12px;font-weight:600;color:var(--dsw-alias-brand-primary)}',
      '.bg-geo-title{font-size:12px;color:var(--dsw-alias-label-primary);line-height:1.4}',
      '.bg-geo-title-en{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.4;display:-webkit-box;-webkit-line-clamp:2;-webkit-box-orient:vertical;overflow:hidden}',
      '.bg-geo-meta{display:flex;flex-wrap:wrap;align-items:center;gap:4px 7px;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.bg-geo-tag{padding:0 6px;height:17px;display:inline-flex;align-items:center;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);font-size:10px;color:var(--dsw-alias-label-secondary)}',
      '.bg-geo-tag-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-success-primary)}',
      '.bg-geo-detail{display:flex;flex-direction:column;gap:9px;min-height:0;overflow-y:auto;padding:1px 2px 4px;max-height:56vh}',
      '.bg-geo-detail-title{font-size:13px;font-weight:600;line-height:1.4}',
      '.bg-geo-kv{display:grid;grid-template-columns:84px 1fr;gap:3px 8px;font-size:11.5px}',
      '.bg-geo-kv dt{color:var(--dsw-alias-label-tertiary)}',
      '.bg-geo-kv dd{margin:0;color:var(--dsw-alias-label-primary);word-break:break-word}',
      '.bg-geo-file{display:flex;align-items:center;gap:7px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);font-size:11.5px}',
      '.bg-geo-file-name{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.bg-progress{display:flex;flex-direction:column;gap:4px;min-width:180px;flex:1 1 auto}',
      '.bg-progress-bar{position:relative;height:7px;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);overflow:hidden}',
      '.bg-progress-fill{position:absolute;left:0;top:0;bottom:0;border-radius:999px;background:linear-gradient(90deg,color-mix(in srgb,var(--dsw-alias-brand-primary) 70%,transparent),var(--dsw-alias-brand-primary));transition:width .2s ease}',
      '.bg-progress-text{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',

      /* 分析抽屉 */
      '.bg-drawer-mask{position:fixed;inset:0;z-index:55;display:flex;justify-content:flex-end;background:color-mix(in srgb,var(--dsw-alias-label-primary) 22%,transparent);animation:bg-fade .14s ease}',
      '.bg-drawer{display:flex;flex-direction:column;width:min(560px,96vw);height:100%;border-left:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-base);box-shadow:-18px 0 60px color-mix(in srgb,var(--dsw-alias-label-primary) 18%,transparent)}',
      '.bg-drawer-body{display:flex;flex-direction:column;gap:10px;padding:12px 14px 18px;overflow-y:auto;min-height:0}',
      '.bg-analysis-group{display:flex;flex-direction:column;gap:5px}',
      '.bg-analysis-group-title{font-size:10.5px;font-weight:600;letter-spacing:.05em;color:var(--dsw-alias-label-tertiary);text-transform:none}',
      '.bg-analysis-item{display:flex;align-items:center;justify-content:space-between;gap:8px;padding:7px 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-layer-1);cursor:pointer;text-align:left;font:inherit;transition:border-color .14s ease,background .14s ease}',
      '.bg-analysis-item:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2)}',
      '.bg-analysis-item:disabled{opacity:.55;cursor:not-allowed}',
      '.bg-analysis-main{display:flex;flex-direction:column;gap:1px;min-width:0;flex:1 1 auto}',
      '.bg-analysis-item-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 48%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,var(--dsw-alias-bg-layer-1))}',
      '.bg-analysis-name{font-size:12px;font-weight:500;color:var(--dsw-alias-label-primary)}',
      '.bg-analysis-method{font-size:10.5px;color:var(--dsw-alias-label-tertiary);margin-top:1px;line-height:1.4}',
      '.bg-analysis-forms{display:flex;flex-direction:column;gap:8px;padding:10px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:11px;background:var(--dsw-alias-bg-layer-1)}',

      /* 右栏 Agent 面板 */
      '.bg-agent-section{display:flex;flex-direction:column;gap:6px;padding:10px 12px;border-bottom:1px solid var(--dsw-alias-border-l1)}',
      '.bg-agent-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:10.5px;font-weight:600;letter-spacing:.05em;color:var(--dsw-alias-label-tertiary)}',
      '.bg-suggestion{display:flex;flex-direction:column;gap:4px;padding:8px 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-layer-2);animation:bg-fade .16s ease}',
      '.bg-suggestion-title{font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-primary);line-height:1.4}',
      '.bg-suggestion-detail{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.45}',
      '.bg-activity{display:flex;flex-direction:column;gap:6px}',
      '.bg-activity-item{display:flex;gap:7px;font-size:11px;line-height:1.4}',
      '.bg-activity-dot{position:relative;flex:none;width:6px;height:6px;margin-top:4px;border-radius:50%;background:var(--dsw-alias-state-success-primary);box-shadow:0 0 0 2px color-mix(in srgb,var(--dsw-alias-state-success-primary) 16%,transparent)}',
      '.bg-activity-dot-bad{background:var(--dsw-alias-state-error-primary);box-shadow:0 0 0 2px color-mix(in srgb,var(--dsw-alias-state-error-primary) 16%,transparent)}',
      '.bg-activity-main{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.bg-activity-head{display:flex;align-items:center;gap:5px;flex-wrap:wrap;color:var(--dsw-alias-label-secondary)}',
      '.bg-activity-actor{font-size:9.5px;padding:0 5px;height:15px;display:inline-flex;align-items:center;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);color:var(--dsw-alias-label-tertiary)}',
      '.bg-activity-text{color:var(--dsw-alias-label-tertiary);font-size:10.5px;word-break:break-word}',
      '.bg-stat-row{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.bg-stat-row b{font-weight:600;color:var(--dsw-alias-label-primary);font-variant-numeric:tabular-nums}',
      '.bg-history{display:flex;flex-direction:column;gap:4px}',
      '.bg-history-item{display:flex;align-items:center;gap:7px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1);font-size:11px}',
      '.bg-history-at{flex:none;font-size:10.5px;color:var(--dsw-alias-label-tertiary);font-variant-numeric:tabular-nums}',
      '.bg-history-step{flex:none;font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.bg-history-summary{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:var(--dsw-alias-label-tertiary)}',
      '.bg-history-ms{flex:none;color:var(--dsw-alias-label-tertiary);font-size:10.5px}',
      '.bg-audit-item{display:flex;align-items:flex-start;gap:7px;padding:6px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1);font-size:11px}',
      '.bg-audit-item-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 36%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 5%,var(--dsw-alias-bg-layer-1))}',
      '.bg-audit-item-fail{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 36%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 5%,var(--dsw-alias-bg-layer-1))}',
      '.bg-audit-mark{flex:none;width:16px;height:16px;display:inline-flex;align-items:center;justify-content:center;border-radius:50%;font-size:9.5px;font-weight:700;background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 16%,transparent);color:var(--dsw-alias-state-success-primary)}',
      '.bg-audit-mark-warn{background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 16%,transparent);color:var(--dsw-alias-state-warn-primary)}',
      '.bg-audit-mark-fail{background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 16%,transparent);color:var(--dsw-alias-state-error-primary)}',
      '.bg-audit-main{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.bg-audit-title{font-weight:600;color:var(--dsw-alias-label-primary)}',
      '.bg-audit-detail{color:var(--dsw-alias-label-tertiary);font-size:10.5px;line-height:1.45}',
      '.bg-score{display:flex;align-items:baseline;gap:6px}',
      '.bg-score-v{font-size:30px;font-weight:250;line-height:1;color:var(--dsw-alias-label-primary)}',
      '.bg-score-label{font-size:11px;color:var(--dsw-alias-label-tertiary)}',

      /* 响应式 */
      '@media (max-width:1180px){.bg-left{flex-basis:240px;width:240px}.bg-right{flex-basis:262px;width:262px}}',
      '@media (max-width:980px){',
      '  .bg-left{flex-basis:216px;width:216px}.bg-right{display:none}',
      '  .bg-geo-layout{grid-template-columns:1fr}',
      '  .bg-field .bg-input{max-width:120px}',
      '}',
      '@media (prefers-reduced-motion:reduce){.bg-spinner{animation:none}.bg-step-dot-running{animation:none}}',
    ].join('')

    // ==========================================================================
    // 通用小组件
    // ==========================================================================

    function Spinner() {
      return h('span', { className: 'bg-spinner', 'aria-hidden': 'true' })
    }

    function LoadingLine(props) {
      return h('div', { className: 'bg-loading' }, h(Spinner), h('span', null, (props && props.text) || T('loading')))
    }

    function ErrorBanner(props) {
      return h('div', { className: 'bg-crash', role: 'alert' },
        h('div', { className: 'bg-crash-title' }, (props && props.title) || T('errorTitle')),
        props && props.message ? h('div', { className: 'bg-crash-detail' }, String(props.message)) : null,
        props && props.onRetry
          ? h('button', { type: 'button', className: 'bg-btn bg-btn-sm', onClick: props.onRetry }, (props && props.retryLabel) || T('retry'))
          : null)
    }

    function InlineError(props) {
      if (!props || !props.message) return null
      return h('div', { className: 'bg-inline-err', role: 'alert' }, String(props.message))
    }

    function Notice(props) {
      return h('div', { className: 'bg-notice' }, (props && props.text) || '')
    }

    /** 内联 SVG 图标（不用 emoji）。 */
    function Icon(props) {
      var name = (props && props.name) || 'dot'
      var size = (props && props.size) || 14
      var paths = {
        refresh: 'M20 11a8 8 0 1 0-2.3 5.7M20 4v7h-7',
        search: 'M10.5 17a6.5 6.5 0 1 0 0-13 6.5 6.5 0 0 0 0 13zM15.2 15.2 21 21',
        import: 'M12 3v12m0 0 4-4m-4 4-4-4M4 17v3h16v-3',
        download: 'M12 3v12m0 0 4-4m-4 4-4-4M4 20h16',
        close: 'M6 6l12 12M18 6 6 18',
        play: 'M7 5l12 7-12 7z',
        chart: 'M4 20V10M10 20V4M16 20v-7M22 20H2',
        gear: 'M12 15.2a3.2 3.2 0 1 0 0-6.4 3.2 3.2 0 0 0 0 6.4zM19.4 15a7.6 7.6 0 0 0 .1-1.2 7.6 7.6 0 0 0-.1-1.2l2-1.6-2-3.4-2.3.9a7.7 7.7 0 0 0-2-1.2L14.6 4h-4l-.5 2.3a7.7 7.7 0 0 0-2 1.2l-2.3-.9-2 3.4 2 1.6a7.6 7.6 0 0 0-.1 1.2c0 .4 0 .8.1 1.2l-2 1.6 2 3.4 2.3-.9c.6.5 1.3.9 2 1.2l.5 2.3h4l.5-2.3c.7-.3 1.4-.7 2-1.2l2.3.9 2-3.4z',
        chevron: 'M9 6l6 6-6 6',
        chevronDown: 'M6 9l6 6 6-6',
        copy: 'M9 9h11v11H9zM5 15V4h11',
        check: 'M4 12.5l5 5L20 6.5',
        external: 'M14 4h6v6M20 4l-9 9M18 14v6H4V6h6',
        folder: 'M3 6h6l2 2h10v11H3z',
        plus: 'M12 5v14M5 12h14',
        flask: 'M9 3h6M10 3v6l-5 9a2 2 0 0 0 1.7 3h10.6a2 2 0 0 0 1.7-3l-5-9V3',
        shield: 'M12 3l8 3v6c0 5-3.5 8-8 9-4.5-1-8-4-8-9V6z',
        clock: 'M12 21a9 9 0 1 0 0-18 9 9 0 0 0 0 18zM12 7v5l3.2 2',
        database: 'M12 3c4.4 0 8 1.3 8 3s-3.6 3-8 3-8-1.3-8-3 3.6-3 8-3zM4 6v6c0 1.7 3.6 3 8 3s8-1.3 8-3V6M4 12v6c0 1.7 3.6 3 8 3s8-1.3 8-3v-6',
        dot: 'M12 12h.01',
      }
      return h('svg', {
        width: size,
        height: size,
        viewBox: '0 0 24 24',
        fill: name === 'play' ? 'currentColor' : 'none',
        stroke: 'currentColor',
        strokeWidth: props && props.strokeWidth ? props.strokeWidth : 1.7,
        strokeLinecap: 'round',
        strokeLinejoin: 'round',
        'aria-hidden': 'true',
        focusable: 'false',
        style: { flex: 'none', display: 'inline-block' },
      }, h('path', { d: paths[name] || paths.dot }))
    }

    function Btn(props) {
      var children = Array.prototype.slice.call(arguments, 1)
      var extra = props || {}
      // 两种调用方式都支持：Btn(props, ...children) 直呼；h(Btn, props, ...children) 经 React。
      if (!children.length && extra.children !== undefined && extra.children !== null && extra.children !== false) {
        children = Array.isArray(extra.children) ? extra.children : [extra.children]
      }
      var className = 'bg-btn' +
        (extra.variant === 'primary' ? ' bg-btn-primary' : '') +
        (extra.variant === 'ghost' ? ' bg-btn-ghost' : '') +
        (extra.variant === 'danger' ? ' bg-btn-danger' : '') +
        (extra.small ? ' bg-btn-sm' : '') +
        (extra.className ? ' ' + extra.className : '')
      return h('button', {
        type: 'button',
        className: className,
        title: extra.title || undefined,
        disabled: extra.disabled === true,
        onClick: extra.onClick,
      }, extra.icon ? h(Icon, { name: extra.icon, size: extra.small ? 12 : 14 }) : null, children.filter(Boolean))
    }

    /** 渲染崩溃兜底：每个主要区块都包一层，避免整个页面白屏。 */
    class PluginErrorBoundary extends React.Component {
      constructor(props) {
        super(props)
        this.state = { error: null }
        this.handleReset = this.handleReset.bind(this)
      }

      static getDerivedStateFromError(error) {
        return { error: error || new Error('unknown render error') }
      }

      componentDidCatch(error, info) {
        warnPlugin('UI crashed: ' + formatError(error), info)
      }

      handleReset() {
        this.setState({ error: null })
      }

      render() {
        if (this.state.error) {
          if (this.props && this.props.mini) {
            return h('span', { className: 'bg-crash-mini', title: T('crashTitle') }, '!')
          }
          return h('div', { className: 'bg-crash', role: 'alert' },
            h('div', { className: 'bg-crash-title' }, T('crashTitle')),
            h('div', { className: 'bg-crash-hint' }, T('crashHint')),
            h('div', { className: 'bg-crash-detail' }, formatError(this.state.error)),
            h('button', { type: 'button', className: 'bg-btn bg-btn-primary bg-btn-sm', onClick: this.handleReset }, T('retry')))
        }
        return this.props ? this.props.children : null
      }
    }

    // ==========================================================================
    // 侧边栏图标（纯 SVG）
    // ==========================================================================

    function SidebarGlyph(props) {
      var size = (props && props.size) || 18
      return h('svg', {
        width: size,
        height: size,
        viewBox: '0 0 24 24',
        fill: 'none',
        stroke: 'currentColor',
        strokeWidth: 1.6,
        strokeLinecap: 'round',
        strokeLinejoin: 'round',
        'aria-hidden': 'true',
        focusable: 'false',
      },
        h('path', { d: 'M4 19.5V10M9.3 19.5V5M14.7 19.5v-6.5M20 19.5V8' }),
        h('circle', { cx: 9.3, cy: 5, r: 1.4, fill: 'currentColor', stroke: 'none' }),
        h('circle', { cx: 14.7, cy: 13, r: 1.4, fill: 'currentColor', stroke: 'none' }),
        h('path', { d: 'M3 21.5h18' }))
    }

    function SidebarIcon(props) {
      return h(PluginErrorBoundary, { mini: true }, h(SidebarGlyph, { size: props && props.size }))
    }

    // ==========================================================================
    // 通用展示组件：指标卡 / 图 / 表 / 结果信封 / 参数表单 / 弹层
    // ==========================================================================

    function MetricCard(props) {
      return h('div', { className: 'bg-metric', title: props.hint || undefined },
        h('div', { className: 'bg-metric-k' }, props.label),
        h('div', { className: 'bg-metric-v' }, props.value, props.unit ? h('small', null, props.unit) : null))
    }

    function SectionCard(props) {
      return h('div', { className: 'bg-card' },
        props.title || props.extra
          ? h('div', { className: 'bg-card-title' },
              h('span', null, props.title || ''),
              props.extra ? h('span', null, props.extra) : null)
          : null,
        props.children)
    }

    function FigureBlock(props) {
      var failed = useState(false)
      var broken = failed[0]
      var setBroken = failed[1]
      var figure = props.figure || null
      useEffect(function () {
        setBroken(false)
      }, [figure && figure.url ? figure.url : ''])
      if (props.loading) {
        return h('div', { className: 'bg-figure' },
          props.title ? h('div', { className: 'bg-figure-title' }, props.title) : null,
          h('div', { className: 'bg-figure-empty' }, h(LoadingLine, { text: props.loadingText || T('plotting') })))
      }
      if (props.error) {
        return h('div', { className: 'bg-figure' },
          props.title ? h('div', { className: 'bg-figure-title' }, props.title) : null,
          h('div', { className: 'bg-inline-err' }, props.error),
          props.onReload ? h('div', null, Btn({ small: true, icon: 'refresh', onClick: props.onReload }, T('retry'))) : null)
      }
      if (!figure || !figure.url) {
        return h('div', { className: 'bg-figure' },
          props.title ? h('div', { className: 'bg-figure-title' }, props.title) : null,
          h('div', { className: 'bg-figure-empty' }, props.emptyText || T('noFigure')))
      }
      return h('div', { className: 'bg-figure' },
        props.title ? h('div', { className: 'bg-figure-title' }, props.title) : null,
        broken
          ? h('div', { className: 'bg-figure-empty' }, T('figureFail'))
          : h('img', {
              src: figure.url,
              alt: figure.title || props.title || 'figure',
              loading: 'lazy',
              onError: function () { setBroken(true) },
            }),
        figure.file ? h('div', { className: 'bg-figure-meta' }, figure.file) : null,
        props.onReload && !broken ? h('div', null, Btn({ small: true, variant: 'ghost', icon: 'refresh', onClick: props.onReload }, T('retry'))) : null)
    }

    /** 自动加载 + 缓存展示一张服务端图（pca/volcano/...）。 */
    function PlotBlock(props) {
      var key = plotKey(props.kind, props.params)
      var entry = props.plots[key] || null
      var loadPlot = props.loadPlot
      var enabled = props.enabled !== false
      useEffect(function () {
        if (!enabled) return
        if (entry && (entry.loading || entry.figure || entry.error)) return
        loadPlot(props.kind, props.params)
      }, [enabled, key, entry ? (entry.loading ? 2 : (entry.figure ? 1 : 3)) : 0])
      var reload = function () { loadPlot(props.kind, props.params, true) }
      return h(FigureBlock, {
        title: props.title,
        figure: entry && entry.figure ? entry.figure : null,
        loading: !!(entry && entry.loading),
        error: entry && entry.error ? entry.error : '',
        emptyText: props.emptyText || T('noFigure'),
        onReload: reload,
      })
    }

    function DataTable(props) {
      var table = props.table || { columns: [], rows: [] }
      var sortState = useState({ key: '', dir: 'asc' })
      var sort = sortState[0]
      var setSort = sortState[1]
      var limit = isFinite(Number(props.maxRows)) ? Number(props.maxRows) : MAX_TABLE_RENDER_ROWS
      var rows = table.rows || []
      var view
      if (props.sortable !== false) {
        var ordered = sort.key ? sortRows(rows, sort.key, sort.dir, table.columns) : rows
        view = formatTableRows({ columns: table.columns, rows: ordered, truncated: table.truncated }, limit)
      } else {
        view = formatTableRows(table, limit)
      }
      function clickHeader(column) {
        if (!props.sortable) return
        setSort(function (previous) {
          if (previous.key === column.key) return { key: column.key, dir: previous.dir === 'asc' ? 'desc' : 'asc' }
          return { key: column.key, dir: 'asc' }
        })
      }
      if (!view.columns.length) {
        return h('div', { className: 'bg-empty-table' }, props.emptyText || T('emptyTable'))
      }
      return h('div', null,
        h('div', { className: 'bg-table-wrap', style: props.maxHeight ? { maxHeight: props.maxHeight } : undefined },
          h('table', { className: 'bg-table' },
            h('thead', null,
              h('tr', null, view.columns.map(function (column) {
                var active = sort.key === column.key
                return h('th', {
                  key: column.key,
                  className: (props.sortable !== false ? 'bg-th-sort' : '') + (column.type === 'number' ? ' bg-td-num' : ''),
                  title: props.sortable !== false ? T('deSortHint') : undefined,
                  onClick: function () { clickHeader(column) },
                }, column.label + (active ? (sort.dir === 'asc' ? ' ↑' : ' ↓') : ''))
              }))),
            h('tbody', null, view.rows.map(function (row) {
              return h('tr', { key: row.key }, row.cells.map(function (cell) {
                return h('td', { key: cell.key, className: cell.numeric ? 'bg-td-num' : '', title: cell.text }, cell.text)
              }))
            })))),
        h('div', { className: 'bg-table-note' },
          tf(T, 'tableRows', view.total) +
          (view.truncated ? ' · ' + tf(T, 'tableTruncated', view.shown) : '') +
          (table.note ? ' · ' + table.note : '') +
          (props.sortable !== false && sort.key ? ' · ' + T('deSortHint') : '')))
    }

    /** 通用结果信封渲染（§2.5）：bullets / metrics / tables / figures / warnings / data_out。 */
    function ResultView(props) {
      var result = props.result || null
      if (!result) {
        return props.emptyText ? h('div', { className: 'bg-notice' }, props.emptyText) : null
      }
      var blocks = []
      if (result.bullets && result.bullets.length) {
        blocks.push(h('div', { key: 'bullets', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, result.title || T('resultTitle')),
          h('ul', { className: 'bg-bullets' }, result.bullets.map(function (bullet, index) {
            return h('li', { key: 'b-' + index }, bullet)
          }))))
      }
      var metricKeys = Object.keys(result.metrics || {})
      if (metricKeys.length) {
        blocks.push(h('div', { key: 'metrics', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, T('resultMetrics')),
          h('div', { className: 'bg-metrics' }, metricKeys.map(function (key) {
            var value = result.metrics[key]
            return h(MetricCard, {
              key: key,
              label: key,
              value: typeof value === 'number' ? fmtNumber(value) : formatCell(value, false),
            })
          }))))
      }
      if (result.tables && result.tables.length) {
        blocks.push(h('div', { key: 'tables', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, T('resultTables') + ' · ' + result.tables.length),
          result.tables.map(function (table, index) {
            return h('div', { key: 'table-' + (table.id || index), className: 'bg-result-block' },
              table.title ? h('div', { className: 'bg-figure-title' }, table.title) : null,
              h(DataTable, { table: table, maxRows: props.maxTableRows || MAX_TABLE_RENDER_ROWS }))
          })))
      }
      if (result.figures && result.figures.length) {
        blocks.push(h('div', { key: 'figures', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, T('resultFigures') + ' · ' + result.figures.length),
          h('div', { className: 'bg-grid-2' }, result.figures.map(function (figure, index) {
            return h('div', { key: 'figure-' + (figure.id || index), className: 'bg-result-block' },
              h(FigureBlock, { title: figure.title, figure: figure }))
          }))))
      }
      if (result.warnings && result.warnings.length) {
        blocks.push(h('div', { key: 'warnings', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, T('resultWarnings')),
          h('div', { className: 'bg-warn-list' }, result.warnings.map(function (warning, index) {
            return h('div', { key: 'w-' + index }, '· ' + warning)
          }))))
      }
      if (result.notes && result.notes.length) {
        blocks.push(h('div', { key: 'notes', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, T('resultNotes')),
          h('div', { className: 'bg-note-list' }, result.notes.map(function (note, index) {
            return h('div', { key: 'n-' + index }, '· ' + note)
          }))))
      }
      if (result.data_out && result.data_out.length) {
        blocks.push(h('div', { key: 'files', className: 'bg-result-block' },
          h('div', { className: 'bg-result-block-title' }, T('resultFiles')),
          h('div', { className: 'bg-inline-row' }, result.data_out.map(function (item, index) {
            return h('span', { key: 'f-' + index, className: 'bg-file-chip', title: item.path },
              h(Icon, { name: 'download', size: 11 }),
              h('span', null, item.label || item.kind),
              h('span', { className: 'bg-path' }, item.path))
          }))))
      }
      if (!blocks.length) return h('div', { className: 'bg-notice' }, T('resultEmpty'))
      return h('div', { className: 'bg-result' }, blocks)
    }

    // ==========================================================================
    // 参数表单：按 params schema 通用渲染
    //   number / select / boolean / column / columns / text / genes / file / multiselect
    // ==========================================================================

    function ParamControl(props) {
      var descriptor = props.descriptor || {}
      var value = props.value
      var disabled = props.disabled === true
      var onChange = props.onChange || function () {}
      var type = descriptor.type || 'text'
      if (type === 'number') {
        return h('input', {
          className: 'bg-input',
          type: 'number',
          min: isNil(descriptor.min) ? undefined : descriptor.min,
          max: isNil(descriptor.max) ? undefined : descriptor.max,
          step: isNil(descriptor.step) ? 'any' : descriptor.step,
          value: isNil(value) ? '' : String(value),
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value === '' ? null : Number(event.target.value)) },
        })
      }
      if (type === 'boolean') {
        return h('label', { className: 'bg-check' },
          h('input', {
            type: 'checkbox',
            checked: value === true,
            disabled: disabled,
            onChange: function (event) { onChange(event.target.checked) },
          }),
          h('span', null, descriptor.checkLabel || (value === true ? T('yes') : T('no'))))
      }
      if (type === 'select') {
        var options = Array.isArray(descriptor.options) ? descriptor.options : []
        var current = isNil(value) ? '' : String(value)
        var known = false
        for (var optionIndex = 0; optionIndex < options.length; optionIndex++) {
          if (String(options[optionIndex].value) === current) known = true
        }
        var list = options.slice()
        if (current && !known) list.unshift({ value: current, label: current })
        return h('select', {
          className: 'bg-input bg-select',
          value: current,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value) },
        }, list.map(function (option) {
          return h('option', { key: String(option.value), value: String(option.value) }, option.label)
        }))
      }
      if (type === 'column') {
        var columns = Array.isArray(props.columns) ? props.columns : []
        var columnOptions = [{ value: '', label: T('noSelection') }]
        for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) {
          var column = columns[columnIndex]
          columnOptions.push({
            value: column.key,
            label: column.key + (column.kind && column.kind !== 'categorical' ? ' · ' + column.kind : ''),
          })
        }
        var columnValue = isNil(value) ? '' : String(value)
        var hasColumn = false
        for (var checkIndex = 0; checkIndex < columnOptions.length; checkIndex++) {
          if (columnOptions[checkIndex].value === columnValue) hasColumn = true
        }
        if (columnValue && !hasColumn) columnOptions.push({ value: columnValue, label: columnValue })
        return h('select', {
          className: 'bg-input bg-select',
          value: columnValue,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value || null) },
        }, columnOptions.map(function (option) {
          return h('option', { key: option.value || 'none', value: option.value }, option.label)
        }))
      }
      if (type === 'multiselect') {
        var optionList = Array.isArray(descriptor.options) ? descriptor.options : []
        var selected = Array.isArray(value) ? value : parseTokenList(value)
        return h('div', { className: 'bg-multi' }, optionList.map(function (option) {
          var on = selected.indexOf(option.value) >= 0
          return h('button', {
            key: String(option.value),
            type: 'button',
            disabled: disabled,
            className: 'bg-multi-opt' + (on ? ' bg-multi-opt-on' : ''),
            onClick: function () {
              var next = selected.slice()
              var index = next.indexOf(option.value)
              if (index >= 0) next.splice(index, 1)
              else next.push(option.value)
              onChange(next)
            },
          }, (on ? '✓ ' : '') + (option.label || option.value))
        }))
      }
      var textValue = isNil(value) ? '' : (Array.isArray(value) ? value.join(', ') : String(value))
      var isList = type === 'genes' || type === 'columns' || type === 'array'
      if (isList || type === 'file' && props.multiline) {
        return h('textarea', {
          className: 'bg-input bg-textarea',
          value: textValue,
          disabled: disabled,
          spellCheck: false,
          placeholder: descriptor.placeholder || (type === 'genes' ? 'CDH1, VIM, MKI67' : ''),
          onChange: function (event) { onChange(parseTokenList(event.target.value)) },
        })
      }
      return h('input', {
        className: 'bg-input',
        type: 'text',
        value: textValue,
        disabled: disabled,
        spellCheck: false,
        placeholder: descriptor.placeholder || (type === 'file' || type === 'text' && descriptor.path ? '/absolute/path' : ''),
        onChange: function (event) { onChange(event.target.value) },
      })
    }

    function ParamForm(props) {
      var params = Array.isArray(props.params) ? props.params : []
      var values = props.values && typeof props.values === 'object' ? props.values : {}
      var onChange = props.onChange || function () {}
      var columns = Array.isArray(props.columns) ? props.columns : []
      if (!params.length) return h('div', { className: 'bg-muted' }, T('analysisNoParams'))
      return h('div', { className: 'bg-form' }, params.map(function (descriptor) {
        if (!descriptor || !descriptor.key) return null
        var label = descriptor.label || descriptor.labelKey || descriptor.key
        var control = h(ParamControl, {
          key: 'control-' + descriptor.key,
          descriptor: descriptor,
          columns: columns,
          value: values[descriptor.key],
          disabled: props.disabled,
          onChange: function (next) {
            var patch = {}
            patch[descriptor.key] = next
            onChange(patch)
          },
        })
        var isBoolean = descriptor.type === 'boolean'
        var isFull = descriptor.type === 'multiselect' || descriptor.type === 'genes' || descriptor.type === 'columns' || descriptor.type === 'array' || (descriptor.type === 'text' && descriptor.wide)
        return h('div', { key: descriptor.key, className: 'bg-field' + (isFull || isBoolean ? ' bg-field-col' : '') },
          h('span', { className: 'bg-field-label', title: label }, label),
          control)
      }))
    }

    function Modal(props) {
      function onMask(event) {
        if (event.target === event.currentTarget && !props.busy && props.onClose) props.onClose()
      }
      return h('div', { className: 'bg-modal-mask', onMouseDown: onMask },
        h('div', { className: 'bg-modal' + (props.wide ? ' bg-modal-wide' : '') },
          h('div', { className: 'bg-modal-head' },
            h('div', null,
              h('div', { className: 'bg-modal-title' }, props.title || ''),
              props.subtitle ? h('div', { className: 'bg-modal-subtitle' }, props.subtitle) : null),
            h('button', {
              type: 'button',
              className: 'bg-btn bg-btn-sm bg-btn-ghost',
              disabled: props.busy === true,
              title: T('close'),
              onClick: props.onClose,
            }, h(Icon, { name: 'close', size: 13 }))),
          props.children,
          props.footer ? h('div', { className: 'bg-modal-foot' }, props.footer) : null))
    }

    function GuideBlock(props) {
      return h('div', { className: 'bg-guide' },
        h('div', { className: 'bg-guide-title' }, props.title || ''),
        props.body ? h('div', { className: 'bg-guide-body' }, props.body) : null,
        props.hint ? h('div', { className: 'bg-guide-hint' }, props.hint) : null,
        props.actions ? h('div', { className: 'bg-inline-row' }, props.actions) : null)
    }

    /** 弹层里的字段行。 */
    function FormField(props) {
      return h('label', { className: 'bg-field' + (props.full ? ' bg-field-col' : '') },
        h('span', { className: 'bg-field-label' }, props.label),
        props.children)
    }

    // ==========================================================================
    // 顶栏 / 数据集 chips
    // ==========================================================================

    function TopBar(props) {
      var badge = engineBadge(props.health)
      var datasets = Array.isArray(props.datasets) ? props.datasets : []
      return h('header', { className: 'bg-top' },
        h('div', { className: 'bg-brand' },
          h('span', { className: 'bg-brand-mark' }, h(Icon, { name: 'chart', size: 15, strokeWidth: 1.8 })),
          h('div', { className: 'bg-title-block' },
            h('div', { className: 'bg-title', title: T('appTitle') }, T('appTitle')),
            h('div', { className: 'bg-subtitle', title: T('appSubtitle') }, T('appSubtitle')))),
        h('span', { className: 'bg-badge bg-badge-' + (props.health && props.health.ready ? 'ok' : (props.health && props.health.error ? 'bad' : 'warn')), title: badge.detail },
          h('span', { className: 'bg-badge-dot' }),
          badge.ready ? T(badge.labelKey) : T(badge.labelKey) + (badge.detail ? ' · ' + badge.detail : '')),
        props.health && props.health.python
          ? h('span', { className: 'bg-badge', title: T('engineInfo') }, 'Python ' + props.health.python)
          : null,
        props.health
          ? h('span', { className: 'bg-badge', title: T('rVersion') },
              props.health.r && props.health.r.available ? 'R ' + props.health.r.version : T('rMissing'))
          : null,
        h('div', { className: 'bg-chips' }, datasets.map(function (dataset) {
          var on = dataset.id === props.activeId
          return h('span', {
            key: dataset.id,
            className: 'bg-chip' + (on ? ' bg-chip-on' : ''),
            title: (dataset.path || dataset.name) + (dataset.n_samples ? ' · ' + dataset.n_samples + ' ' + T('samples') : ''),
          },
            h('button', {
              type: 'button',
              className: 'bg-btn bg-btn-ghost bg-btn-sm',
              style: { height: 'auto', padding: '0 2px', border: 'none', color: 'inherit', fontSize: '11.5px' },
              disabled: props.busy,
              onClick: function () { if (!on) props.onSelect(dataset.id) },
            }, dataset.name || dataset.id),
            on
              ? h('button', {
                  type: 'button',
                  className: 'bg-chip-x',
                  title: T('closeDataset'),
                  disabled: props.busy,
                  onClick: function (event) { event.stopPropagation(); props.onCloseDataset(dataset.id) },
                }, '×')
              : null)
        })),
        h('div', { className: 'bg-top-actions' },
          Btn({ small: true, icon: 'flask', onClick: props.onAnalysis, disabled: !props.activeId, title: T('analysisHint') }, T('analysisCatalog')),
          Btn({ small: true, icon: 'search', onClick: props.onGeo }, T('geoSearch')),
          Btn({ small: true, icon: 'import', onClick: props.onImport }, T('importData')),
          Btn({ small: true, icon: 'refresh', onClick: props.onRefresh, disabled: props.refreshing, title: T('refreshHint') }, T('refresh'))))
    }

    // ==========================================================================
    // 左栏：数据集摘要 + 12 张流程卡
    // ==========================================================================

    function columnSummary(dataset, keys) {
      var columns = dataset && Array.isArray(dataset.sample_columns) ? dataset.sample_columns : []
      for (var index = 0; index < columns.length; index++) {
        for (var keyIndex = 0; keyIndex < keys.length; keyIndex++) {
          if (columns[index].key === keys[keyIndex]) return columns[index]
        }
      }
      return null
    }

    function DatasetCard(props) {
      var dataset = props.dataset
      if (!dataset) {
        return h('div', { className: 'bg-dataset-card' },
          h('div', { className: 'bg-dataset-name' }, T('datasetNoDataset')),
          h('div', { className: 'bg-dataset-meta' }, T('datasetOpenHint')))
      }
      var groupColumn = columnSummary(dataset, [dataset.group_key || 'group'])
      var batchColumn = columnSummary(dataset, ['batch', 'batch_key'])
      var meta = []
      if (!isNil(dataset.n_samples)) meta.push(h('span', { key: 's' }, T('samples') + ' ', h('b', null, fmtNumber(dataset.n_samples))))
      if (!isNil(dataset.n_genes)) meta.push(h('span', { key: 'g' }, T('genes') + ' ', h('b', null, fmtNumber(dataset.n_genes))))
      if (groupColumn) meta.push(h('span', { key: 'grp' }, T('groups') + ' ', h('b', null, groupColumn.n_unique || '—')))
      if (batchColumn) meta.push(h('span', { key: 'bat' }, T('batches') + ' ', h('b', null, batchColumn.n_unique || '—')))
      if (dataset.source_format) meta.push(h('span', { key: 'fmt' }, h('b', null, dataset.source_format)))
      return h('div', { className: 'bg-dataset-card' },
        h('div', { className: 'bg-dataset-name', title: dataset.path || dataset.name }, dataset.name || dataset.id),
        h('div', { className: 'bg-dataset-meta' }, meta),
        h('div', { className: 'bg-dataset-meta' },
          h('span', null, tf(T, 'datasetStepsDone', dataset.steps.length)),
          dataset.current ? h('span', null, tf(T, 'datasetCurrentStep', dataset.current)) : null),
        dataset.warnings.length
          ? h('div', { className: 'bg-warn-text', title: dataset.warnings.join('\n') }, T('warnings') + ' · ' + dataset.warnings.length)
          : null)
    }

    function StepCard(props) {
      var step = props.step
      var status = props.status
      var expandedState = useState(status === 'running' || status === 'failed')
      var expanded = expandedState[0]
      var setExpanded = expandedState[1]
      var statusKey = stepStatusKey(status)
      var dotClass = status === 'done' ? ' bg-step-dot-done' : (status === 'running' ? ' bg-step-dot-running' : (status === 'blocked' ? ' bg-step-dot-blocked' : ' bg-step-dot-pending'))
      var cardClass = 'bg-step' + (status === 'done' ? ' bg-step-done' : (status === 'running' ? ' bg-step-running' : (status === 'blocked' ? ' bg-step-blocked' : '')))
      var blockedRequires = []
      var done = doneStepSet(props.dataset)
      var requires = Array.isArray(step.requires) ? step.requires : []
      for (var index = 0; index < requires.length; index++) {
        if (!done[requires[index]]) blockedRequires.push(stepById(requires[index]) ? stepById(requires[index]).title : requires[index])
      }
      var disabled = props.busy || !props.dataset || status === 'blocked'
      var summary = props.summary || ''
      return h('div', { className: cardClass, 'data-step': step.id },
        h('div', { className: 'bg-step-head' },
          h('span', { className: 'bg-step-dot' + dotClass, title: T(statusKey) }, String(props.index + 1)),
          h('span', { className: 'bg-step-title', title: step.title }, step.title),
          h('span', { className: 'bg-step-status bg-status-' + status }, T(statusKey)),
          h('button', {
            type: 'button',
            className: 'bg-step-caret',
            title: expanded ? T('analysisClose') : T('stepParams'),
            'aria-expanded': expanded ? 'true' : 'false',
            onClick: function () { setExpanded(!expanded) },
          }, h(Icon, { name: expanded ? 'chevronDown' : 'chevron', size: 12 }))),
        expanded
          ? h('div', { className: 'bg-step-params' },
              h(ParamForm, {
                params: step.params,
                values: props.values || {},
                columns: props.columns || [],
                disabled: props.busy,
                onChange: props.onParamChange,
              }))
          : null,
        h('div', { className: 'bg-step-foot' },
          Btn({
            small: true,
            variant: status === 'done' ? undefined : 'primary',
            icon: status === 'done' ? 'refresh' : 'play',
            disabled: disabled,
            title: status === 'blocked' ? T('stepNeedPrereq') + blockedRequires.join('、') : T('stepRun'),
            onClick: function () { props.onRun(step.id) },
          }, props.busy ? T('stepRunning') : (status === 'done' ? T('stepRerun') : T('stepRun'))),
          h('span', { className: 'bg-step-summary', title: summary || '' }, summary || (props.busy ? T('stepRunning') : (status === 'blocked' ? T('stepNeedPrereq') + blockedRequires.join('、') : T('stepSummaryEmpty'))))),
        props.error ? h(InlineError, { message: props.error }) : null)
    }

    function PipelinePanel(props) {
      var steps = STEPS
      return h('aside', { className: 'bg-left' },
        h('div', { className: 'bg-panel-head' },
          h('span', null, T('datasetsTitle')),
          h('span', { className: 'bg-panel-note' }, T('pipelineHint'))),
        h('div', { className: 'bg-panel-body', style: { paddingBottom: '8px' } },
          h(DatasetCard, { dataset: props.dataset })),
        h('div', { className: 'bg-steps' }, steps.map(function (step, index) {
          return h(StepCard, {
            key: step.id,
            step: step,
            index: index,
            status: props.statuses[step.id],
            dataset: props.dataset,
            values: props.values[step.id],
            columns: props.columns,
            summary: props.summaries[step.id],
            error: props.errors[step.id],
            busy: props.busyStep === step.id,
            onParamChange: function (patch) { props.onParamChange(step.id, patch) },
            onRun: props.onRun,
          })
        })))
    }

    // ==========================================================================
    // 右栏：Agent 面板（状态摘要 / 建议 + 复制指令 / 活动时间线）
    // ==========================================================================

    function AgentPanel(props) {
      var dataset = props.dataset
      var suggestions = Array.isArray(props.suggestions) ? props.suggestions : []
      var activity = Array.isArray(props.activity) ? props.activity : []
      var sections = []
      if (!dataset) {
        sections.push(h('div', { key: 'empty', className: 'bg-agent-section' },
          h('div', { className: 'bg-agent-title' }, T('agentStatus')),
          h('div', { className: 'bg-muted' }, T('agentNoDataset'))))
      } else {
        var currentStep = stepById(dataset.current)
        sections.push(h('div', { key: 'status', className: 'bg-agent-section' },
          h('div', { className: 'bg-agent-title' }, T('agentStatus')),
          h('div', { className: 'bg-stat-row' }, h('span', null, dataset.name || dataset.id), h('b', null, (dataset.n_samples || '—') + ' × ' + (dataset.n_genes || '—'))),
          h('div', { className: 'bg-stat-row' }, h('span', null, T('completed')), h('b', null, dataset.steps.length + '/' + STEPS.length)),
          h('div', { className: 'bg-stat-row' }, h('span', null, T('current')), h('b', null, currentStep ? currentStep.title : '—')),
          h('div', { className: 'bg-stat-row' }, h('span', null, T('engine')), h('b', null, props.health && props.health.ready ? T('engineReady') : T('engineBooting')))))
      }
      sections.push(h('div', { key: 'suggestions', className: 'bg-agent-section' },
        h('div', { className: 'bg-agent-title' }, T('agentSuggestions')),
        suggestions.length
          ? suggestions.map(function (suggestion) {
              var stepMeta = suggestion.step ? stepById(suggestion.step) : null
              var copied = props.copiedKey === ('suggestion-' + suggestion.id)
              return h('div', { key: suggestion.id, className: 'bg-suggestion' },
                h('div', { className: 'bg-suggestion-title' }, suggestion.title || T('nextSteps')),
                suggestion.detail ? h('div', { className: 'bg-suggestion-detail' }, suggestion.detail) : null,
                h('div', { className: 'bg-inline-row' },
                  stepMeta ? h('span', { className: 'bg-geo-tag' }, stepMeta.title) : null,
                  Btn({
                    small: true,
                    variant: 'ghost',
                    icon: copied ? 'check' : 'copy',
                    onClick: function () { props.onCopy(suggestionInstruction(suggestion), 'suggestion-' + suggestion.id) },
                  }, copied ? T('copied') : T('agentCopy'))))
            })
          : h('div', { className: 'bg-muted' }, T('noSuggestions'))))
      sections.push(h('div', { key: 'activity', className: 'bg-agent-section', style: { borderBottom: 'none' } },
        h('div', { className: 'bg-agent-title' }, T('agentActivity')),
        activity.length
          ? h('div', { className: 'bg-activity' }, activity.slice(0, 10).map(function (entry) {
              return h('div', { key: entry.id, className: 'bg-activity-item' },
                h('span', { className: 'bg-activity-dot' + (entry.ok ? '' : ' bg-activity-dot-bad') }),
                h('div', { className: 'bg-activity-main' },
                  h('div', { className: 'bg-activity-head' },
                    h('span', { className: 'bg-activity-actor' }, entry.actor === 'agent' ? T('activityAgent') : T('activityUser')),
                    h('span', null, entry.action || T('statusDot')),
                    h('span', { className: 'bg-history-at' }, formatClock(entry.at))),
                  entry.text ? h('div', { className: 'bg-activity-text' }, entry.text) : null))
            }))
          : h('div', { className: 'bg-muted' }, T('agentNoActivity'))))
      return h('aside', { className: 'bg-right' },
        h('div', { className: 'bg-panel-head' },
          h('span', null, T('agentTitle')),
          h('span', { className: 'bg-panel-note' }, 'Agent')),
        sections)
    }

    // ==========================================================================
    // 中栏：Tab 框架 + 概览 / 表达矩阵 / 差异 / 富集
    // ==========================================================================

    function TabsBar(props) {
      return h('nav', { className: 'bg-tabs', role: 'tablist' }, TAB_KEYS.map(function (key) {
        return h('button', {
          key: key,
          type: 'button',
          role: 'tab',
          'aria-selected': props.active === key ? 'true' : 'false',
          className: 'bg-tab' + (props.active === key ? ' bg-tab-on' : ''),
          onClick: function () { props.onChange(key) },
        }, T(TAB_META[key].titleKey))
      }))
    }

    function TabGuide(props) {
      return h(GuideBlock, {
        title: props.title,
        body: props.body,
        hint: props.hint,
        actions: props.action ? Btn({ small: true, variant: 'primary', icon: 'play', onClick: props.action.onClick, disabled: props.action.disabled }, props.action.label) : null,
      })
    }

    function overviewStat(dataset, key) {
      var column = columnSummary(dataset, [key, key + '_key'])
      return column ? (column.n_unique || 0) : 0
    }

    function OverviewTab(props) {
      var dataset = props.dataset
      if (!dataset) {
        return h('div', { className: 'bg-pane' },
          h(TabGuide, {
            title: T('datasetNoDataset'),
            body: T('datasetOpenHint'),
            action: { label: T('importData'), onClick: props.onOpenImport },
          }))
      }
      var preview = props.preview || normalizePreview({})
      var geneState = useState('')
      var geneInput = geneState[0]
      var setGeneInput = geneState[1]
      var defaultGene = (preview.genes && preview.genes.length ? preview.genes[0] : '')
      var effectiveGenes = parseTokenList(geneInput || defaultGene).slice(0, 4)
      var engineInfo = props.health && props.health.ready
        ? (props.health.python ? 'Python ' + props.health.python : T('engineReady')) + (props.health.r && props.health.r.available ? ' · R ' + props.health.r.version : '')
        : T('engineBooting')
      var metrics = [
        { label: T('metricSamples'), value: fmtNumber(dataset.n_samples) },
        { label: T('metricGenes'), value: fmtNumber(dataset.n_genes) },
        { label: T('metricGroups'), value: fmtNumber(overviewStat(dataset, 'group')) },
        { label: T('metricBatches'), value: fmtNumber(overviewStat(dataset, 'batch')) },
        { label: T('metricModules'), value: fmtNumber(dataset.modules.length) },
        { label: T('metricEngine'), value: engineInfo },
      ]
      var suggestions = Array.isArray(props.suggestions) ? props.suggestions : []
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-toolbar' },
          h('span', { className: 'bg-section-title' }, T('overviewTitle')),
          props.busyStep ? h('span', { className: 'bg-loading' }, h(Spinner), h('span', null, props.busyStep)) : null),
        h('div', { className: 'bg-metrics' }, metrics.map(function (metric) {
          return h(MetricCard, { key: metric.label, label: metric.label, value: metric.value })
        })),
        h('div', { className: 'bg-grid-2' },
          h(SectionCard, { title: T('pcaPlot') },
            h(PlotBlock, { plots: props.plots, loadPlot: props.loadPlot, kind: 'pca_scatter', params: {}, enabled: !!dataset })),
          h(SectionCard, { title: T('corrPlot') },
            h(PlotBlock, { plots: props.plots, loadPlot: props.loadPlot, kind: 'sample_corr', params: {}, enabled: !!dataset }))),
        h('div', { className: 'bg-grid-2' },
          h(SectionCard, {
            title: T('geneQuickTitle'),
            extra: h(Btn, { small: true, variant: 'ghost', icon: 'chart', onClick: function () { setGeneInput('') } }, T('clearSelection')),
          },
            h('div', { className: 'bg-muted' }, T('geneQuickHint')),
            h('div', { className: 'bg-inline-row' },
              h('input', {
                className: 'bg-input',
                style: { flex: '1 1 auto', minWidth: '180px' },
                value: geneInput,
                placeholder: T('geneInputPlaceholder'),
                onChange: function (event) { setGeneInput(event.target.value) },
              })),
            h(PlotBlock, {
              plots: props.plots,
              loadPlot: props.loadPlot,
              kind: 'gene_box',
              params: { genes: effectiveGenes },
              enabled: effectiveGenes.length > 0,
              emptyText: T('exprEmpty'),
            })),
          h(SectionCard, { title: T('nextSteps'), extra: T('agentSuggestions') },
            suggestions.length
              ? h('div', { style: { display: 'flex', flexDirection: 'column', gap: '6px' } }, suggestions.slice(0, 5).map(function (suggestion) {
                  return h('div', { key: suggestion.id, className: 'bg-suggestion' },
                    h('div', { className: 'bg-suggestion-title' }, suggestion.title),
                    suggestion.detail ? h('div', { className: 'bg-suggestion-detail' }, suggestion.detail) : null,
                    h('div', { className: 'bg-inline-row' },
                      suggestion.step && stepById(suggestion.step)
                        ? Btn({ small: true, variant: 'ghost', icon: 'play', onClick: function () { props.onRunStep(suggestion.step) } }, stepById(suggestion.step).title)
                        : null))
                }))
              : h('div', { className: 'bg-muted' }, T('noSuggestions')))),
        dataset.warnings.length
          ? h(SectionCard, { title: T('warnings'), extra: dataset.warnings.length },
              h('div', { className: 'bg-warn-list' }, dataset.warnings.map(function (warning, index) {
                return h('div', { key: 'w-' + index }, '· ' + warning)
              })))
          : null)
    }

    function GenePicker(props) {
      var queryState = useState('')
      var query = queryState[0]
      var setQuery = queryState[1]
      var allGenes = Array.isArray(props.allGenes) ? props.allGenes : []
      var selected = Array.isArray(props.value) ? props.value : []
      var matches = useMemo(function () { return matchGenes(allGenes, query, 8) }, [allGenes.length, query])
      function addGene(gene) {
        var name = plainText(gene)
        if (!name || selected.indexOf(name) >= 0) return
        props.onChange(selected.concat([name]))
        setQuery('')
      }
      return h('div', { style: { display: 'flex', flexDirection: 'column', gap: '6px' } },
        h('div', { className: 'bg-inline-row' },
          h('input', {
            className: 'bg-input',
            style: { flex: '1 1 auto', minWidth: '160px' },
            value: query,
            placeholder: props.placeholder || T('exprGeneSearch'),
            onChange: function (event) { setQuery(event.target.value) },
            onKeyDown: function (event) {
              if (event.key === 'Enter') {
                event.preventDefault()
                addGene(matches.length ? matches[0] : query)
              }
            },
          }),
          Btn({ small: true, icon: 'plus', onClick: function () { addGene(matches.length ? matches[0] : query) } }, T('exprAddGene'))),
        query
          ? h('div', { className: 'bg-multi' }, matches.length
              ? matches.map(function (gene) {
                  return h('button', { key: gene, type: 'button', className: 'bg-multi-opt', onClick: function () { addGene(gene) } }, gene)
                })
              : h('span', { className: 'bg-muted' }, T('exprGeneNotInData')))
          : null,
        selected.length
          ? h('div', { className: 'bg-multi' }, selected.map(function (gene) {
              return h('button', {
                key: gene,
                type: 'button',
                className: 'bg-multi-opt bg-multi-opt-on',
                title: T('clearSelection'),
                onClick: function () { props.onChange(selected.filter(function (item) { return item !== gene })) },
              }, gene + ' ×')
            }))
          : h('div', { className: 'bg-muted' }, T('exprEmpty')))
    }

    function ExprTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var preview = props.preview || normalizePreview({})
      var genesState = useState([])
      var genes = genesState[0]
      var setGenes = genesState[1]
      var analysisState = useState(null)
      var exprAnalysis = analysisState[0]
      var setExprAnalysis = analysisState[1]
      var analysisBusyState = useState(false)
      var analysisBusy = analysisBusyState[0]
      var setAnalysisBusy = analysisBusyState[1]
      var analysisErrorState = useState('')
      var analysisError = analysisErrorState[0]
      var setAnalysisError = analysisErrorState[1]
      var exported = props.exports.expr_csv || {}
      function runHeatmap() {
        if (!genes.length) return
        setAnalysisBusy(true)
        setAnalysisError('')
        props.callAnalysis('gene_heatmap', { genes: genes })
          .then(function (result) { setExprAnalysis(result) })
          .catch(function (error) { setAnalysisError(formatError(error)) })
          .then(function () { setAnalysisBusy(false) })
      }
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-section-title' }, T('exprTitle')),
        h(SectionCard, { title: T('exprSelected'), extra: h('span', { className: 'bg-muted' }, T('exprHint')) },
          h(GenePicker, { allGenes: preview.genes, value: genes, onChange: setGenes })),
        h('div', { className: 'bg-toolbar' },
          Btn({ small: true, variant: 'primary', icon: 'chart', disabled: !genes.length, onClick: runHeatmap }, analysisBusy ? T('running') : T('exprPlotHeatmap')),
          Btn({ small: true, icon: 'download', disabled: exported.busy === true, onClick: function () { props.runExport('expr_csv') } }, exported.busy ? T('recordsExporting') : T('exprExport')),
          exported.file ? h('span', { className: 'bg-file-chip' }, h(Icon, { name: 'check', size: 11 }), h('span', { className: 'bg-path' }, exported.file)) : null,
          exported.error ? h(InlineError, { message: exported.error }) : null),
        h(SectionCard, { title: T('exprPlotBox') },
          h(PlotBlock, {
            plots: props.plots,
            loadPlot: props.loadPlot,
            kind: 'gene_box',
            params: { genes: genes },
            enabled: genes.length > 0,
            emptyText: T('exprEmpty'),
          })),
        analysisError ? h(InlineError, { message: analysisError }) : null,
        exprAnalysis ? h(SectionCard, { title: T('exprPlotHeatmap') }, h(ResultView, { result: exprAnalysis })) : null,
        h(Notice, { text: T('exprTableFallback') }),
        h('div', { className: 'bg-muted' }, T('exprExportHint')))
    }

    function DeTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var result = props.stepResults.de || null
      var figureState = useState('volcano')
      var figureKind = figureState[0]
      var setFigureKind = figureState[1]
      var exported = props.exports.de_csv || {}
      var figureKinds = [
        { key: 'volcano', labelKey: 'deVolcano' },
        { key: 'ma_plot', labelKey: 'deMa' },
        { key: 'de_heatmap', labelKey: 'deHeatmap' },
      ]
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-toolbar' },
          h('span', { className: 'bg-section-title' }, T('deTitle')),
          h('span', { className: 'bg-muted' }, T('deHint')),
          h('span', { style: { flex: '1 1 auto' } }),
          Btn({ small: true, icon: 'download', disabled: exported.busy === true || !result, onClick: function () { props.runExport('de_csv') } }, exported.busy ? T('recordsExporting') : T('deExport'))),
        !result
          ? h(TabGuide, {
              title: T('deEmpty'),
              body: T('deHint'),
              action: { label: T('stepRun') + ' · ' + stepById('de').title, onClick: function () { props.onRunStep('de') }, disabled: props.stepBusy },
            })
          : h('div', { className: 'bg-pane', style: { padding: 0, overflow: 'visible' } },
              exported.file ? h('div', { className: 'bg-file-chip' }, h(Icon, { name: 'check', size: 11 }), h('span', { className: 'bg-path' }, exported.file)) : null,
              exported.error ? h(InlineError, { message: exported.error }) : null,
              h(ResultView, { result: result, maxTableRows: 200 }),
              h('div', { className: 'bg-toolbar' },
                figureKinds.map(function (item) {
                  return Btn({
                    key: item.key,
                    small: true,
                    variant: figureKind === item.key ? 'primary' : undefined,
                    icon: 'chart',
                    onClick: function () { setFigureKind(item.key) },
                  }, T(item.labelKey))
                })),
              h(SectionCard, null,
                h(PlotBlock, {
                  plots: props.plots,
                  loadPlot: props.loadPlot,
                  kind: figureKind,
                  params: {},
                  enabled: !!result,
                  title: T(figureKinds[0].key === figureKind ? 'deVolcano' : (figureKind === 'ma_plot' ? 'deMa' : 'deHeatmap')),
                }))))
    }

    function EnrichTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var step = STEP_BY_ID.enrich
      var valuesState = useState(function () { return defaultParams(step) })
      var values = valuesState[0]
      var setValues = valuesState[1]
      var result = props.stepResults.enrich || null
      var figureState = useState('enrich_bar')
      var figureKind = figureState[0]
      var setFigureKind = figureState[1]
      function patch(next) { setValues(function (previous) { return Object.assign({}, previous, next) }) }
      function run() {
        var params = collectStepParams(step, values)
        if (!params.db || !params.db.length) {
          setValues(function (previous) { return previous })
          return
        }
        props.runStep('enrich', params)
      }
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-section-title' }, T('enrichTitle')),
        h('div', { className: 'bg-muted' }, T('enrichHint')),
        h(SectionCard, { title: T('stepParams') },
          h(ParamForm, {
            params: step.params,
            values: values,
            columns: (props.preview || {}).sample_columns,
            disabled: props.stepBusy === 'enrich',
            onChange: patch,
          }),
          h('div', { className: 'bg-inline-row' },
            Btn({
              small: true,
              variant: 'primary',
              icon: 'play',
              disabled: props.stepBusy === 'enrich',
              onClick: run,
            }, props.stepBusy === 'enrich' ? T('stepRunning') : T('enrichRun')))),
        result
          ? h('div', { style: { display: 'flex', flexDirection: 'column', gap: '10px' } },
              h('div', { className: 'bg-toolbar' },
                Btn({ small: true, variant: figureKind === 'enrich_bar' ? 'primary' : undefined, icon: 'chart', onClick: function () { setFigureKind('enrich_bar') } }, T('enrichBar')),
                Btn({ small: true, variant: figureKind === 'enrich_dot' ? 'primary' : undefined, icon: 'chart', onClick: function () { setFigureKind('enrich_dot') } }, T('enrichDot'))),
              h(SectionCard, null,
                h(PlotBlock, {
                  plots: props.plots,
                  loadPlot: props.loadPlot,
                  kind: figureKind,
                  params: {},
                  enabled: !!result,
                  title: figureKind === 'enrich_bar' ? T('enrichBar') : T('enrichDot'),
                })),
              h(ResultView, { result: result }))
          : h(Notice, { text: T('enrichEmpty') }))
    }

    // ==========================================================================
    // 中栏：免疫 / 生存 / 预后 / 单细胞 / 非编码 / 记录 / 日志
    // ==========================================================================

    function StepFormTab(props) {
      var step = STEP_BY_ID[props.stepId]
      var valuesState = useState(function () { return defaultParams(step) })
      var values = valuesState[0]
      var setValues = valuesState[1]
      function patch(next) { setValues(function (previous) { return Object.assign({}, previous, next) }) }
      var busy = props.stepBusy === props.stepId
      var result = props.stepResults[props.stepId] || null
      var columns = (props.preview || {}).sample_columns || []
      var extraToolbar = props.guideAction && !result
        ? null
        : null
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-section-title' }, props.title),
        h('div', { className: 'bg-muted' }, props.hint),
        props.notice ? h(Notice, { text: props.notice }) : null,
        h(SectionCard, { title: T('stepParams') },
          h(ParamForm, {
            params: step.params,
            values: values,
            columns: columns,
            disabled: busy,
            onChange: patch,
          }),
          h('div', { className: 'bg-inline-row' },
            Btn({
              small: true,
              variant: 'primary',
              icon: 'play',
              disabled: busy || !props.dataset,
              onClick: function () { props.runStep(props.stepId, collectStepParams(step, values)) },
            }, busy ? T('stepRunning') : props.runLabel),
            extraToolbar)),
        result
          ? h(React.Fragment, null,
              h('div', { className: 'bg-toolbar' }, (props.figureKinds || []).map(function (item) {
                return Btn({
                  key: item.kind,
                  small: true,
                  variant: props.figureKind === item.kind ? 'primary' : undefined,
                  icon: 'chart',
                  onClick: function () { props.setFigureKind(item.kind) },
                }, T(item.labelKey))
              })),
              props.figureKind
                ? h(SectionCard, null,
                    h(PlotBlock, {
                      plots: props.plots,
                      loadPlot: props.loadPlot,
                      kind: props.figureKind,
                      params: {},
                      enabled: !!result,
                      title: T((props.figureKinds || []).filter(function (item) { return item.kind === props.figureKind })[0].labelKey),
                    }))
                : null,
              h(ResultView, { result: result }))
          : h(Notice, { text: props.emptyText }))
    }

    function ImmuneTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var figureState = useState('immune_stack')
      return h(StepFormTab, {
        stepId: 'immune',
        title: T('immuneTitle'),
        hint: T('immuneHint'),
        emptyText: T('immuneEmpty'),
        runLabel: T('immuneRun'),
        figureKinds: [
          { kind: 'immune_stack', labelKey: 'immuneStack' },
          { kind: 'immune_corr', labelKey: 'immuneCorr' },
        ],
        figureKind: figureState[0],
        setFigureKind: figureState[1],
        dataset: props.dataset,
        preview: props.preview,
        stepResults: props.stepResults,
        stepBusy: props.stepBusy,
        runStep: props.runStep,
        plots: props.plots,
        loadPlot: props.loadPlot,
      })
    }

    function SurvivalTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var figureState = useState('km_curve')
      var columns = (props.preview || {}).sample_columns || []
      var timeColumn = columns.filter(function (column) { return column.kind === 'time' })[0]
      var eventColumn = columns.filter(function (column) { return column.kind === 'event' })[0]
      if (!timeColumn && !eventColumn) {
        return h('div', { className: 'bg-pane' },
          h(TabGuide, {
            title: T('survivalTitle'),
            body: T('survivalHint'),
            hint: T('survivalNeedColumns'),
            action: { label: T('importData'), onClick: props.onOpenImport },
          }))
      }
      return h(StepFormTab, {
        stepId: 'survival',
        title: T('survivalTitle'),
        hint: T('survivalHint'),
        emptyText: T('survivalEmpty'),
        runLabel: T('survivalRun'),
        notice: T('survivalNeedColumns'),
        figureKinds: [
          { kind: 'km_curve', labelKey: 'survivalKm' },
          { kind: 'forest', labelKey: 'survivalForest' },
        ],
        figureKind: figureState[0],
        setFigureKind: figureState[1],
        dataset: props.dataset,
        preview: props.preview,
        stepResults: props.stepResults,
        stepBusy: props.stepBusy,
        runStep: props.runStep,
        plots: props.plots,
        loadPlot: props.loadPlot,
      })
    }

    function PrognosisTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var figureState = useState('risk_plot')
      return h(StepFormTab, {
        stepId: 'prognosis',
        title: T('prognosisTitle'),
        hint: T('prognosisHint'),
        emptyText: T('prognosisEmpty'),
        runLabel: T('prognosisRun'),
        figureKinds: [
          { kind: 'risk_plot', labelKey: 'prognosisRisk' },
          { kind: 'roc_time', labelKey: 'prognosisRoc' },
        ],
        figureKind: figureState[0],
        setFigureKind: figureState[1],
        dataset: props.dataset,
        preview: props.preview,
        stepResults: props.stepResults,
        stepBusy: props.stepBusy,
        runStep: props.runStep,
        plots: props.plots,
        loadPlot: props.loadPlot,
      })
    }

    function ScTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var step = STEP_BY_ID.scmap
      var refsState = useState({ loading: true, error: '', available: false, datasets: [], local: [] })
      var refs = refsState[0]
      var setRefs = refsState[1]
      var valuesState = useState(function () { return defaultParams(step) })
      var values = valuesState[0]
      var setValues = valuesState[1]
      var attachState = useState({ busy: false, error: '', ok: '' })
      var attach = attachState[0]
      var setAttach = attachState[1]
      var attachRef = useRef(false)
      useEffect(function () {
        if (attachRef.current) return
        attachRef.current = true
        setRefs(function (previous) { return Object.assign({}, previous, { loading: true }) })
        props.api.scRefs()
          .then(function (payload) {
            var datasets = Array.isArray(payload.datasets) ? payload.datasets : []
            var local = Array.isArray(payload.local) ? payload.local : []
            setRefs({ loading: false, error: '', available: payload.available === true, datasets: datasets, local: local })
          })
          .catch(function (error) {
            setRefs({ loading: false, error: formatError(error), available: false, datasets: [], local: [] })
          })
      }, [])
      function patch(next) { setValues(function (previous) { return Object.assign({}, previous, next) }) }
      function attachRefPath(path, scDataset) {
        setAttach({ busy: true, error: '', ok: '' })
        var payload = { dataset: props.dataset.id }
        if (scDataset) payload.sc_dataset = scDataset
        else payload.ref_path = path
        props.api.scAttach(payload)
          .then(function (payload2) {
            var resolved = payload2.ref_path || payload2.path || path || ''
            patch({ ref_path: resolved })
            setAttach({ busy: false, error: '', ok: T('scAttachOk') })
            props.refreshState()
          })
          .catch(function (error) { setAttach({ busy: false, error: formatError(error), ok: '' }) })
      }
      var busy = props.stepBusy === 'scmap'
      var result = props.stepResults.scmap || null
      var workspaceOptions = refs.datasets || []
      var localOptions = refs.local || []
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-section-title' }, T('scTitle')),
        h('div', { className: 'bg-muted' }, T('scHint')),
        h(SectionCard, {
          title: T('scLocalRefs') + ' / ' + T('scPull'),
          extra: refs.loading ? h('span', { className: 'bg-loading' }, h(Spinner)) : (refs.available ? h('span', { className: 'bg-ok-text' }, T('scAvailable')) : h('span', { className: 'bg-muted' }, T('scUnavailable'))),
        },
          refs.error ? h(InlineError, { message: refs.error }) : null,
          !refs.loading && !refs.error && workspaceOptions.length
            ? h('div', { className: 'bg-field bg-field-col' },
                h('span', { className: 'bg-field-label' }, T('scFromWorkspace')),
                h('div', { className: 'bg-inline-row' }, h('select', {
                  className: 'bg-input bg-select',
                  style: { flex: '1 1 200px' },
                  value: '',
                  onChange: function (event) {
                    var value = event.target.value
                    if (!value) return
                    var matched = workspaceOptions.filter(function (item) { return (item.id || item.name) === value })[0]
                    attachRefPath(matched && matched.path ? matched.path : '', matched ? (matched.id || value) : value)
                  },
                }, h('option', { value: '' }, T('noSelection')),
                  workspaceOptions.map(function (item) {
                    return h('option', { key: item.id || item.name, value: item.id || item.name },
                      (item.name || item.id || '') + (item.n_obs ? ' · ' + item.n_obs + ' cells' : ''))
                  }))))
            : null,
          !refs.loading && localOptions.length
            ? h('div', { className: 'bg-field bg-field-col' },
                h('span', { className: 'bg-field-label' }, T('scLocal')),
                h('div', { className: 'bg-inline-row' }, h('select', {
                  className: 'bg-input bg-select',
                  style: { flex: '1 1 200px' },
                  value: '',
                  onChange: function (event) {
                    var value = event.target.value
                    if (value) patch({ ref_path: value })
                  },
                }, h('option', { value: '' }, T('noSelection')),
                  localOptions.map(function (item) {
                    return h('option', { key: item.path, value: item.path }, (item.name || item.path) + (item.bytes ? ' · ' + formatBytes(item.bytes) : ''))
                  }))))
            : (!refs.loading && !refs.error ? h('div', { className: 'bg-muted' }, T('scHintNoRefs')) : null),
          attach.ok ? h('div', { className: 'bg-ok-text' }, attach.ok) : null,
          attach.error ? h(InlineError, { message: attach.error }) : null),
        h(SectionCard, { title: T('stepParams') },
          h(ParamForm, {
            params: [
              { key: 'ref_path', label: T('scRefPath'), type: 'text', default: '' },
              { key: 'gene_set', label: T('scGeneSet'), type: 'select', default: 'de_up', options: [
                { value: 'de_up', label: 'DE up' }, { value: 'de_down', label: 'DE down' }, { value: 'custom', label: 'Custom' },
              ] },
              { key: 'genes', label: T('scGenes'), type: 'genes', default: [] },
              { key: 'cell_type_col', label: T('scCellTypeCol'), type: 'text', default: 'cell_type' },
              { key: 'score_method', label: T('scScoreMethod'), type: 'select', default: 'mean', options: [
                { value: 'mean', label: 'mean' }, { value: 'aucell', label: 'AUCell' }, { value: 'addmodule', label: 'AddModule' },
              ] },
            ],
            values: values,
            columns: [],
            disabled: busy,
            onChange: patch,
          }),
          h('div', { className: 'bg-inline-row' },
            Btn({ small: true, variant: 'primary', icon: 'play', disabled: busy, onClick: function () { props.runStep('scmap', collectStepParams(step, values)) } }, busy ? T('stepRunning') : T('scRun')),
            Btn({
              small: true,
              icon: 'database',
              disabled: attach.busy || !values.ref_path,
              onClick: function () { attachRefPath(values.ref_path, '') },
            }, attach.busy ? T('running') : T('scAttach')))),
        result
          ? h(React.Fragment, null,
              h(SectionCard, null,
                h(PlotBlock, {
                  plots: props.plots,
                  loadPlot: props.loadPlot,
                  kind: 'gene_box',
                  params: { genes: (values.genes || []).slice(0, 4) },
                  enabled: !!result && (values.genes || []).length > 0,
                  title: T('scExprFig'),
                })),
              h(ResultView, { result: result }))
          : h(Notice, { text: T('scEmpty') }))
    }

    function NcrnaTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var valuesState = useState(function () { return defaultParams(STEP_BY_ID.ncrna) })
      var values = valuesState[0]
      var setValues = valuesState[1]
      var attachState = useState({ busy: false, error: '', ok: '' })
      var attach = attachState[0]
      var setAttach = attachState[1]
      function patch(next) { setValues(function (previous) { return Object.assign({}, previous, next) }) }
      var modules = props.dataset.modules || []
      var attached = modules.filter(function (module) { return module.kind === values.kind })[0]
      var busy = props.stepBusy === 'ncrna'
      var result = props.stepResults.ncrna || null
      function attachModule() {
        if (!values.module_path) {
          setAttach({ busy: false, error: T('importNoPath'), ok: '' })
          return
        }
        setAttach({ busy: true, error: '', ok: '' })
        props.api.attach(props.dataset.id, values.kind, values.module_path)
          .then(function () {
            setAttach({ busy: false, error: '', ok: T('ncrnaAttachOk') })
            props.refreshState()
          })
          .catch(function (error) { setAttach({ busy: false, error: formatError(error), ok: '' }) })
      }
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-section-title' }, T('ncrnaTitle')),
        h('div', { className: 'bg-muted' }, T('ncrnaHint')),
        h('div', { className: 'bg-multi' }, modules.map(function (module, index) {
          return h('span', { key: 'module-' + index, className: 'bg-geo-tag bg-geo-tag-ok' },
            module.kind + (module.n ? ' · ' + module.n : ''))
        })),
        h(SectionCard, { title: T('stepParams') },
          h(ParamForm, {
            params: [
              { key: 'kind', label: T('ncrnaKind'), type: 'select', default: 'circ', options: NCRNA_MODULES.map(function (item) {
                return { value: item.value, label: T(item.labelKey) }
              }) },
              { key: 'module_path', label: T('ncrnaPath'), type: 'text', default: '' },
            ],
            values: values,
            columns: [],
            disabled: attach.busy,
            onChange: patch,
          }),
          attached
            ? h('div', { className: 'bg-ok-text' }, T('ncrnaAttached') + ' · ' + values.kind)
            : h('div', { className: 'bg-muted' }, T('ncrnaNotAttached')),
          attach.ok ? h('div', { className: 'bg-ok-text' }, attach.ok) : null,
          attach.error ? h(InlineError, { message: attach.error }) : null,
          h('div', { className: 'bg-inline-row' },
            Btn({ small: true, icon: 'import', disabled: attach.busy, onClick: attachModule }, attach.busy ? T('running') : T('ncrnaAttach')),
            Btn({
              small: true,
              variant: 'primary',
              icon: 'play',
              disabled: busy,
              onClick: function () { props.runStep('ncrna', collectStepParams(STEP_BY_ID.ncrna, values)) },
            }, busy ? T('stepRunning') : T('ncrnaRun')))),
        result ? h(ResultView, { result: result }) : h(Notice, { text: T('ncrnaEmpty') }))
    }

    function RecordsTab(props) {
      if (!props.dataset) return h('div', { className: 'bg-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      var auditState = useState({ busy: false, error: '', audit: null })
      var audit = auditState[0]
      var setAudit = auditState[1]
      function runAudit() {
        setAudit({ busy: true, error: '', audit: null })
        props.api.audit(props.dataset.id)
          .then(function (payload) { setAudit({ busy: false, error: '', audit: normalizeAudit(payload) }) })
          .catch(function (error) { setAudit({ busy: false, error: formatError(error), audit: null }) })
      }
      var history = (props.dataset.history || []).slice().reverse()
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-toolbar' },
          h('span', { className: 'bg-section-title' }, T('recordsTitle')),
          h('span', { style: { flex: '1 1 auto' } }),
          Btn({ small: true, icon: 'flask', onClick: props.onOpenAnalysis }, T('recordsAnalysis'))),
        h('div', { className: 'bg-grid-2' },
          h(SectionCard, { title: T('recordsHistory'), extra: tf(T, 'tableRows', history.length) },
            history.length
              ? h('div', { className: 'bg-history' }, history.map(function (entry, index) {
                  return h('div', { key: 'h-' + index, className: 'bg-history-item' },
                    h('span', { className: 'bg-history-at' }, formatClock(entry.at)),
                    h('span', { className: 'bg-history-step' }, (stepById(entry.step) || {}).title || entry.step || '—'),
                    h('span', { className: 'bg-history-summary', title: entry.summary || '' }, entry.summary || ''),
                    h('span', { className: 'bg-history-ms' }, entry.ms ? formatDuration(entry.ms) : ''),
                    h('span', { className: 'bg-history-ms', className2: '' }, entry.ok ? '' : T('failed')))
                }))
              : h('div', { className: 'bg-muted' }, T('recordsNoHistory'))),
          h(SectionCard, {
            title: T('recordsAudit'),
            extra: h(Btn, { small: true, icon: 'shield', disabled: audit.busy, onClick: runAudit }, audit.busy ? T('running') : T('recordsRunAudit')),
          },
            audit.error ? h(InlineError, { message: audit.error }) : null,
            audit.audit
              ? h('div', { style: { display: 'flex', flexDirection: 'column', gap: '7px' } },
                  h('div', { className: 'bg-score' },
                    h('span', { className: 'bg-score-v' }, audit.audit.score),
                    h('span', { className: 'bg-score-label' }, T('recordsAuditScore') + ' · ' + T(AUDIT_LEVEL_STYLE[audit.audit.level] ? AUDIT_LEVEL_STYLE[audit.audit.level].labelKey : 'auditOk'))),
                  h('div', { className: 'bg-muted' }, T('auditHint')),
                  audit.audit.checks.map(function (check) {
                    var style = AUDIT_LEVEL_STYLE[check.level] || AUDIT_LEVEL_STYLE.ok
                    return h('div', { key: check.id, className: 'bg-audit-item' + (check.level === 'warn' ? ' bg-audit-item-warn' : (check.level === 'fail' ? ' bg-audit-item-fail' : '')) },
                      h('span', { className: 'bg-audit-mark' + (check.level === 'warn' ? ' bg-audit-mark-warn' : (check.level === 'fail' ? ' bg-audit-mark-fail' : '')) }, check.level === 'ok' ? '✓' : (check.level === 'warn' ? '!' : '×')),
                      h('div', { className: 'bg-audit-main' },
                        h('span', { className: 'bg-audit-title' }, check.title || check.id),
                        check.detail ? h('span', { className: 'bg-audit-detail' }, check.detail) : null,
                        check.hint ? h('span', { className: 'bg-audit-detail' }, check.hint) : null))
                  }))
              : h('div', { className: 'bg-muted' }, T('recordsAuditEmpty')))),
        h(SectionCard, { title: T('recordsExport'), extra: h('span', { className: 'bg-muted' }, T('exportBtn')) },
          h('div', { className: 'bg-inline-row' }, EXPORT_KINDS.map(function (item) {
            var state = props.exports[item.id] || {}
            return h('span', { key: item.id, style: { display: 'inline-flex', alignItems: 'center', gap: '5px' } },
              Btn({
                small: true,
                icon: 'download',
                disabled: state.busy === true,
                title: item.hint,
                onClick: function () { props.runExport(item.id) },
              }, T(item.labelKey)),
              state.file ? h('span', { className: 'bg-path', title: state.file }, '✓ ' + state.file) : null,
              state.error ? h(InlineError, { message: state.error }) : null)
          }))))
    }

    function LogTab(props) {
      var state = useState({ loading: true, error: '', entries: [] })
      var log = state[0]
      var setLog = state[1]
      var mountedRef = useRef(false)
      function load() {
        setLog(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        props.api.getActivity()
          .then(function (payload) {
            var list = Array.isArray(payload.activity) ? payload.activity : (Array.isArray(payload.entries) ? payload.entries : [])
            setLog({ loading: false, error: '', entries: list.map(function (entry, index) { return normalizeActivityEntry(entry, index) }) })
          })
          .catch(function (error) { setLog({ loading: false, error: formatError(error), entries: [] }) })
      }
      useEffect(function () {
        if (mountedRef.current) return
        mountedRef.current = true
        load()
      }, [])
      return h('div', { className: 'bg-pane' },
        h('div', { className: 'bg-toolbar' },
          h('span', { className: 'bg-section-title' }, T('logTitle')),
          h('span', { style: { flex: '1 1 auto' } }),
          Btn({ small: true, icon: 'refresh', onClick: load, disabled: log.loading }, T('logRefresh'))),
        log.error ? h(InlineError, { message: log.error }) : null,
        log.loading
          ? h(LoadingLine, { text: T('loading') })
          : (log.entries.length
              ? h('div', { className: 'bg-table-wrap' },
                  h('table', { className: 'bg-table' },
                    h('thead', null, h('tr', null,
                      h('th', null, T('at')), h('th', null, T('actor')), h('th', null, T('action')),
                      h('th', null, T('duration')), h('th', null, T('resultColumn')))),
                    h('tbody', null, log.entries.map(function (entry) {
                      return h('tr', { key: entry.id },
                        h('td', null, formatClock(entry.at)),
                        h('td', null, entry.actor === 'agent' ? T('activityAgent') : T('activityUser')),
                        h('td', null, entry.action || '—'),
                        h('td', { className: 'bg-td-num' }, entry.ms ? formatDuration(entry.ms) : '—'),
                        h('td', null, h('span', { className: entry.ok ? 'bg-ok-text' : 'bg-warn-text' }, entry.ok ? T('completed') : T('failed'))))
                    }))))
              : h(Notice, { text: T('logEmpty') })))
    }

    // ==========================================================================
    // GEO 检索弹层：搜索 → 详情 → 指定目录下载（轮询）→ 导入并质控
    // ==========================================================================

    function normalizeGeoItem(item) {
      var raw = item && typeof item === 'object' ? item : {}
      return {
        gse: plainText(raw.gse || raw.id || ''),
        title: plainText(raw.title || ''),
        titleZh: plainText(raw.title_zh || raw.titleZh || ''),
        summary: plainText(raw.summary || '', 600),
        species: plainText(raw.species || ''),
        nSamples: asNumber(raw.n_samples),
        platform: plainText(raw.platform || (Array.isArray(raw.platforms) ? raw.platforms.join(', ') : '')),
        studyType: plainText(raw.study_type || ''),
        release: plainText(raw.release || ''),
        hasMatrix: raw.has_matrix === true,
        hasClinical: raw.has_clinical === true,
        disease: plainText(raw.disease || ''),
        tags: Array.isArray(raw.tags) ? raw.tags.map(function (tag) { return plainText(tag) }).filter(Boolean) : [],
      }
    }

    function normalizeGeoDetail(payload) {
      var raw = payload && typeof payload === 'object' ? payload : {}
      var files = []
      var source = Array.isArray(raw.files) ? raw.files : []
      for (var index = 0; index < source.length; index++) {
        var file = source[index] && typeof source[index] === 'object' ? source[index] : {}
        var name = plainText(file.name || '')
        if (!name) continue
        files.push({
          name: name,
          size: asNumber(file.size),
          url: plainText(file.url || ''),
          kind: plainText(file.kind || 'other'),
        })
      }
      var clinical = raw.clinical && typeof raw.clinical === 'object' ? raw.clinical : {}
      var fields = []
      if (Array.isArray(clinical.fields)) {
        for (var fieldIndex = 0; fieldIndex < clinical.fields.length; fieldIndex++) {
          var field = clinical.fields[fieldIndex] && typeof clinical.fields[fieldIndex] === 'object' ? clinical.fields[fieldIndex] : {}
          fields.push({ key: plainText(field.key || ''), n: asNumber(field.n), example: plainText(field.example || '') })
        }
      }
      var suggested = plainText(raw.suggested || '')
      if (!suggested) {
        var matrix = files.filter(function (file) { return file.kind === 'series_matrix' || /series_matrix/i.test(file.name) })[0]
        if (matrix) suggested = matrix.name
      }
      return {
        gse: plainText(raw.gse || ''),
        title: plainText(raw.title || ''),
        summary: plainText(raw.summary || '', 1200),
        design: plainText(raw.design || '', 600),
        species: plainText(raw.species || ''),
        nSamples: asNumber(raw.n_samples),
        platforms: Array.isArray(raw.platforms) ? raw.platforms.map(function (item) { return plainText(item) }) : [],
        pubmed: plainText(raw.pubmed || ''),
        release: plainText(raw.release || ''),
        contacts: plainText(raw.contacts || ''),
        clinical: { fromSeriesMatrix: clinical.from_series_matrix === true, fields: fields },
        files: files,
        suggested: suggested,
      }
    }

    function GeoCard(props) {
      var item = props.item
      return h('button', {
        type: 'button',
        className: 'bg-geo-card' + (props.active ? ' bg-geo-card-on' : ''),
        onClick: function () { props.onSelect(item) },
      },
        h('div', { className: 'bg-inline-row' },
          h('span', { className: 'bg-geo-gse' }, item.gse),
          item.hasMatrix ? h('span', { className: 'bg-geo-tag bg-geo-tag-ok' }, T('geoMatrix')) : null,
          item.hasClinical ? h('span', { className: 'bg-geo-tag bg-geo-tag-ok' }, T('geoClinicalTag')) : null),
        h('div', { className: 'bg-geo-title' }, item.titleZh || item.title),
        item.titleZh ? h('div', { className: 'bg-geo-title-en' }, item.title) : null,
        h('div', { className: 'bg-geo-meta' },
          item.nSamples !== null ? h('span', null, item.nSamples + ' ' + T('samples')) : null,
          item.platform ? h('span', null, item.platform) : null,
          item.release ? h('span', null, item.release.slice(0, 4)) : null,
          item.species ? h('span', null, item.species) : null),
        item.tags.length
          ? h('div', { className: 'bg-geo-meta' }, item.tags.slice(0, 5).map(function (tag) {
              return h('span', { key: tag, className: 'bg-geo-tag' }, tag)
            }))
          : null)
    }

    function GeoDialog(props) {
      var filtersState = useState({ term: '', species: 'all', study_type: 'all', min_samples: '', has_clinical: false, limit: 20, page: 0 })
      var filters = filtersState[0]
      var setFilters = filtersState[1]
      var searchState = useState({ busy: false, error: '', total: 0, items: [] })
      var search = searchState[0]
      var setSearch = searchState[1]
      var detailState = useState({ gse: '', busy: false, error: '', detail: null })
      var detail = detailState[0]
      var setDetail = detailState[1]
      var selectedState = useState({})
      var selectedFiles = selectedState[0]
      var setSelectedFiles = selectedState[1]
      var dirState = useState(DEFAULT_GEO_DIR)
      var destDir = dirState[0]
      var setDestDir = dirState[1]
      var downloadState = useState({ busy: false, pct: 0, message: '', status: '', error: '', result: null })
      var download = downloadState[0]
      var setDownload = downloadState[1]
      var importState = useState({ busy: false, error: '' })
      var importInfo = importState[0]
      var setImportInfo = importState[1]
      var timerRef = useRef(null)
      var startRef = useRef(0)
      useEffect(function () {
        return function () {
          if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null }
        }
      }, [])
      function patchFilters(next) { setFilters(function (previous) { return Object.assign({}, previous, next) }) }
      function runSearch(event) {
        if (event && event.preventDefault) event.preventDefault()
        if (!plainText(filters.term)) return
        setSearch({ busy: true, error: '', total: 0, items: [] })
        setDetail({ gse: '', busy: false, error: '', detail: null })
        var payload = buildSearchPayload(filters)
        props.api.geoSearch(payload)
          .then(function (response) {
            var items = (Array.isArray(response.items) ? response.items : []).map(normalizeGeoItem)
            setSearch({ busy: false, error: '', total: asNumber(response.total) || items.length, items: items })
          })
          .catch(function (error) { setSearch({ busy: false, error: formatError(error), total: 0, items: [] }) })
      }
      function openDetail(item) {
        setDetail({ gse: item.gse, busy: true, error: '', detail: null })
        setSelectedFiles({})
        setDownload({ busy: false, pct: 0, message: '', status: '', error: '', result: null })
        props.api.geoDetail(item.gse)
          .then(function (payload) {
            var normalized = normalizeGeoDetail(payload)
            var defaults = {}
            if (normalized.suggested) defaults[normalized.suggested] = true
            setSelectedFiles(defaults)
            setDetail({ gse: item.gse, busy: false, error: '', detail: normalized })
          })
          .catch(function (error) { setDetail({ gse: item.gse, busy: false, error: formatError(error), detail: null }) })
      }
      function pollJob(jobId) {
        startRef.current = Date.now()
        if (timerRef.current) clearInterval(timerRef.current)
        timerRef.current = setInterval(function () {
          props.api.getJob(jobId)
            .then(function (payload) {
              var state = jobPolling(payload, { elapsedMs: Date.now() - startRef.current })
              if (state.running) {
                setDownload(function (previous) { return Object.assign({}, previous, { pct: state.pct, message: state.message }) })
                return
              }
              if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null }
              if (state.done) {
                setDownload({ busy: false, pct: 100, message: state.message || T('geoJobDone'), status: 'done', error: '', result: state.result })
                // 一键「下载并导入」：后端在 job.result.opened 里带回已打开的数据集（已跑 qc），直接切过去。
                if (state.result && state.result.opened && typeof props.onImported === 'function') props.onImported(state.result.opened)
              } else {
                setDownload({
                  busy: false,
                  pct: state.pct,
                  message: state.message,
                  status: state.failed ? 'failed' : 'timeout',
                  error: state.error || (state.failed ? T('geoJobFailed') : T('jobTimeout')),
                  result: null,
                })
              }
            })
            .catch(function (error) {
              if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null }
              setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: formatError(error), result: null })
            })
        }, POLL_MS)
      }
      function startDownload(openAfter) {
        if (!detail.detail) return
        var files = Object.keys(selectedFiles).filter(function (name) { return selectedFiles[name] })
        if (!files.length) {
          setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: T('geoNeedFile'), result: null })
          return
        }
        setDownload({ busy: true, pct: 0, message: '', status: 'running', error: '', result: null })
        props.api.geoDownload({ gse: detail.gse, files: files, dest_dir: destDir, open_after: openAfter === true })
          .then(function (payload) {
            var job = payload.job && typeof payload.job === 'object' ? payload.job : payload
            if (job && job.id) pollJob(job.id)
            else if (job && job.status === 'done') {
              setDownload({ busy: false, pct: 100, message: T('geoJobDone'), status: 'done', error: '', result: job.result || null })
              if (job.result && job.result.opened && typeof props.onImported === 'function') props.onImported(job.result.opened)
            } else {
              setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: T('geoJobFailed'), result: null })
            }
          })
          .catch(function (error) { setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: formatError(error), result: null }) })
      }
      function importedFileFrom(result) {
        var files = result && Array.isArray(result.files) ? result.files : []
        var candidates = files.filter(function (file) { return file && file.path })
        var matrix = candidates.filter(function (file) { return /series_matrix/i.test(String(file.name || file.path)) })[0]
        return matrix || candidates[0] || null
      }
      function importAndQc() {
        var target = importedFileFrom(download.result)
        if (!target) {
          setImportInfo({ busy: false, error: T('geoNeedFile') })
          return
        }
        setImportInfo({ busy: true, error: '' })
        var options = { path: target.path, name: detail.gse || '' }
        if (/series_matrix/i.test(String(target.name || target.path))) options.data_format = 'geo_series'
        props.api.openDataset(options)
          .then(function (payload) {
            var dataset = payload.dataset && typeof payload.dataset === 'object' ? payload.dataset : payload
            var datasetId = plainText(dataset.id || payload.active || '')
            var qcParams = defaultParams(STEP_BY_ID.qc)
            if (!datasetId) throw new Error(T('importError'))
            return props.api.runStep(datasetId, 'qc', qcParams).then(function () {
              setImportInfo({ busy: false, error: '' })
              setDownload(function (previous) { return Object.assign({}, previous, { message: T('geoImportDone') }) })
              props.onImported(datasetId)
            })
          })
          .catch(function (error) { setImportInfo({ busy: false, error: formatError(error) }) })
      }
      var items = search.items
      var activeDetail = detail.detail
      var selectedCount = Object.keys(selectedFiles).filter(function (name) { return selectedFiles[name] }).length
      var searchBody = h('div', { className: 'bg-modal-body' },
        h('form', { className: 'bg-geo-search', onSubmit: runSearch },
          h('input', {
            className: 'bg-input bg-geo-term',
            value: filters.term,
            placeholder: T('geoTermPh'),
            spellCheck: false,
            onChange: function (event) { patchFilters({ term: event.target.value }) },
          }),
          Btn({ small: false, variant: 'primary', icon: 'search', disabled: search.busy, onClick: runSearch },
            search.busy ? T('geoSearching') : T('geoSearchBtn'))),
        h('div', { className: 'bg-geo-filters' },
          h('label', { className: 'bg-field' },
            h('span', { className: 'bg-field-label' }, T('geoSpecies')),
            h('select', {
              className: 'bg-input bg-select',
              value: filters.species,
              onChange: function (event) { patchFilters({ species: event.target.value }) },
            },
              h('option', { value: 'all' }, T('geoAll')),
              h('option', { value: 'human' }, T('geoHuman')),
              h('option', { value: 'mouse' }, T('geoMouse')))),
          h('label', { className: 'bg-field' },
            h('span', { className: 'bg-field-label' }, T('geoStudyType')),
            h('select', {
              className: 'bg-input bg-select',
              value: filters.study_type,
              onChange: function (event) { patchFilters({ study_type: event.target.value }) },
            },
              h('option', { value: 'all' }, T('geoStudyAll')),
              h('option', { value: 'expression profiling by array' }, T('geoStudyArray')),
              h('option', { value: 'high throughput sequencing' }, T('geoStudySeq')))),
          h('label', { className: 'bg-field' },
            h('span', { className: 'bg-field-label' }, T('geoMinSamples')),
            h('input', {
              className: 'bg-input',
              type: 'number',
              min: 0,
              style: { width: '76px' },
              value: filters.min_samples,
              onChange: function (event) { patchFilters({ min_samples: event.target.value }) },
            })),
          h('label', { className: 'bg-check' },
            h('input', { type: 'checkbox', checked: filters.has_clinical, onChange: function (event) { patchFilters({ has_clinical: event.target.checked }) } }),
            h('span', null, T('geoHasClinical')))),
        h('div', { className: 'bg-muted' }, T('geoHint')),
        search.error ? h(InlineError, { message: T('geoErrorSearch') + '：' + search.error }) : null,
        search.busy
          ? h(LoadingLine, { text: T('geoSearching') })
          : h('div', { className: 'bg-geo-layout' },
              h('div', null,
                h('div', { className: 'bg-section-title' },
                  h('span', null, T('geoResults')),
                  h('span', { className: 'bg-muted' }, items.length ? tf(T, 'geoTotal', search.total) : '')),
                items.length
                  ? h('div', { className: 'bg-geo-list' }, items.map(function (item) {
                      return h(GeoCard, {
                        key: item.gse,
                        item: item,
                        active: detail.gse === item.gse,
                        onSelect: openDetail,
                      })
                    }))
                  : h('div', { className: 'bg-muted' }, T('geoNoResults'))),
              h('div', { className: 'bg-geo-detail' },
                detail.busy ? h(LoadingLine, { text: T('geoLoadingDetail') }) : null,
                detail.error ? h(InlineError, { message: T('geoErrorDetail') + '：' + detail.error }) : null,
                activeDetail
                  ? h(React.Fragment, null,
                      h('div', { className: 'bg-geo-detail-title' }, activeDetail.title || detail.gse),
                      activeDetail.title && activeDetail.titleZh ? h('div', { className: 'bg-geo-title-en' }, activeDetail.titleZh) : null,
                      h('dl', { className: 'bg-geo-kv' },
                        h('dt', null, T('dataset')), h('dd', null, activeDetail.gse),
                        h('dt', null, T('geoPlatform')), h('dd', null, (activeDetail.platforms || []).join(', ') || '—'),
                        h('dt', null, T('samples')), h('dd', null, activeDetail.nSamples === null ? '—' : String(activeDetail.nSamples)),
                        h('dt', null, T('geoRelease')), h('dd', null, activeDetail.release || '—'),
                        h('dt', null, T('geoPubmed')), h('dd', null, activeDetail.pubmed || '—')),
                      activeDetail.summary ? h('div', null, h('div', { className: 'bg-field-label' }, T('geoSummary')), h('div', { className: 'bg-muted' }, activeDetail.summary)) : null,
                      activeDetail.design ? h('div', null, h('div', { className: 'bg-field-label' }, T('geoDesign')), h('div', { className: 'bg-muted' }, activeDetail.design)) : null,
                      activeDetail.clinical.fields.length
                        ? h('div', null,
                            h('div', { className: 'bg-field-label' }, T('geoClinicalTag') + ' · ' + activeDetail.clinical.fields.length),
                            h('div', { className: 'bg-multi' }, activeDetail.clinical.fields.slice(0, 12).map(function (field) {
                              return h('span', { key: field.key, className: 'bg-geo-tag' }, field.key + (field.n ? ' · ' + field.n : ''))
                            })))
                        : null,
                      h('div', { className: 'bg-field-label' }, T('geoFiles') + ' · ' + activeDetail.files.length),
                      activeDetail.files.length
                        ? h('div', { style: { display: 'flex', flexDirection: 'column', gap: '4px' } }, activeDetail.files.map(function (file) {
                            var on = selectedFiles[file.name] === true
                            return h('label', { key: file.name, className: 'bg-geo-file' },
                              h('input', {
                                type: 'checkbox',
                                checked: on,
                                onChange: function (event) {
                                  var next = Object.assign({}, selectedFiles)
                                  next[file.name] = event.target.checked
                                  setSelectedFiles(next)
                                },
                              }),
                              h('span', { className: 'bg-geo-file-name', title: file.name }, file.name),
                              h('span', { className: 'bg-muted' }, file.kind || ''),
                              h('span', { className: 'bg-muted' }, isNil(file.size) ? '' : formatBytes(file.size)))
                          }))
                        : h('div', { className: 'bg-muted' }, T('geoNoFiles')),
                      h('div', { className: 'bg-inline-row' },
                        h('span', { className: 'bg-field-label' }, T('geoSelectAll')),
                        Btn({ small: true, icon: 'check', onClick: function () {
                          var all = {}
                          activeDetail.files.forEach(function (file) { all[file.name] = true })
                          setSelectedFiles(all)
                        } }, T('geoSelectAll')),
                        Btn({ small: true, icon: 'close', onClick: function () { setSelectedFiles({}) } }, T('geoClear')),
                        h('span', { className: 'bg-muted' }, selectedCount + ' / ' + activeDetail.files.length)))
                  : h('div', { className: 'bg-muted' }, T('geoDetailTitle') + ' · ' + T('geoDetailBtn')))),
        h('div', { className: 'bg-divider' }),
        h('div', { className: 'bg-inline-row' },
          h('span', { className: 'bg-field-label', style: { whiteSpace: 'nowrap' } }, T('geoDownloadDir')),
          h('input', {
            className: 'bg-input',
            style: { flex: '1 1 260px' },
            value: destDir,
            spellCheck: false,
            placeholder: DEFAULT_GEO_DIR,
            onChange: function (event) { setDestDir(event.target.value) },
          }),
          h('span', { className: 'bg-muted' }, T('geoDownloadDirHint'))),
        download.busy || download.status
          ? h('div', { className: 'bg-progress' },
              h('div', { className: 'bg-progress-bar' },
                h('div', { className: 'bg-progress-fill', style: { width: clamp(download.pct, 0, 100) + '%' } })),
              h('div', { className: 'bg-progress-text' },
                h('span', null, download.message || (download.busy ? T('geoDownloading') : (download.status === 'done' ? T('geoJobDone') : ''))),
                h('span', null, Math.round(clamp(download.pct, 0, 100)) + '%')))
          : null,
        download.error ? h(InlineError, { message: download.error }) : null,
        importInfo.error ? h(InlineError, { message: importInfo.error }) : null,
        h('div', { className: 'bg-inline-row' },
          Btn({
            small: false,
            variant: 'primary',
            icon: 'download',
            disabled: !activeDetail || download.busy || selectedCount === 0 || !plainText(destDir),
            onClick: function () { startDownload(false) },
          }, download.busy ? T('geoDownloading') : T('geoDownloadBtn')),
          Btn({
            small: false,
            icon: 'import',
            title: T('geoOpenAfter'),
            disabled: !activeDetail || download.busy || selectedCount === 0 || !plainText(destDir),
            onClick: function () { startDownload(true) },
          }, T('geoOpenAfter')),
          download.result
            ? Btn({ small: false, icon: 'import', disabled: importInfo.busy, onClick: importAndQc },
                importInfo.busy ? T('geoImporting') : T('geoImportBtn'))
            : null,
          Btn({ small: false, icon: 'close', onClick: props.onClose }, T('close'))))
      return h(Modal, {
        title: T('geoTitle'),
        subtitle: T('geoSubtitle'),
        wide: true,
        busy: download.busy || importInfo.busy,
        onClose: props.onClose,
      }, searchBody)
    }

    // ==========================================================================
    // 导入弹层：本机路径 / URL / 上传 + 格式 + 高级项 + 候选点选重试
    // ==========================================================================

    function formatOptions(formatsPayload) {
      var raw = formatsPayload && typeof formatsPayload === 'object' ? formatsPayload : {}
      var list = Array.isArray(raw) ? raw : (Array.isArray(raw.formats) ? raw.formats : [])
      return list.map(function (item) {
        if (typeof item === 'string') return { value: item, label: item }
        var format = item && typeof item === 'object' ? item : {}
        var value = plainText(format.id || format.value || format.name || '')
        var label = plainText(format.label || format.title || format.id || format.name || value)
        return { value: value, label: label }
      }).filter(function (item) { return !!item.value })
    }

    function ImportDialog(props) {
      var sourceState = useState('local')
      var source = sourceState[0]
      var setSource = sourceState[1]
      var formState = useState({
        path: '', url: '', data_format: 'auto', sheet: '', gene_col: '', transpose: false,
        member: '', object: '', species: 'auto', sample_sheet: '', clinical: '',
      })
      var form = formState[0]
      var setForm = formState[1]
      var fileState = useState(null)
      var file = fileState[0]
      var setFile = fileState[1]
      var busyState = useState({ import: false, inspect: false, error: '', inspectError: '', result: null, candidates: [] })
      var busy = busyState[0]
      var setBusy = busyState[1]
      function patchForm(next) { setForm(function (previous) { return Object.assign({}, previous, next) }) }
      var formats = formatOptions(props.formats)
      var options = function () {
        var payload = Object.assign({}, form)
        if (source === 'url') payload.path = ''
        if (source === 'local') payload.url = ''
        return payload
      }
      function handleAmbiguous(error) {
        var candidates = Array.isArray(error.candidates) ? error.candidates : []
        var names = candidates.map(function (candidate) {
          if (typeof candidate === 'string') return candidate
          if (candidate && typeof candidate === 'object') return plainText(candidate.name || candidate.path || candidate.sheet || candidate.id || '')
          return ''
        }).filter(Boolean)
        setBusy(function (previous) { return Object.assign({}, previous, { error: formatError(error), candidates: names }) })
      }
      function submit() {
        if (source === 'upload') {
          if (!file) {
            setBusy(function (previous) { return Object.assign({}, previous, { error: T('importNoPath'), candidates: [] }) })
            return
          }
          setBusy(function (previous) { return Object.assign({}, previous, { import: true, error: '', candidates: [] }) })
          Promise.resolve(typeof file.arrayBuffer === 'function' ? file.arrayBuffer() : file)
            .then(function (buffer) { return props.api.uploadDataset(file.name || 'upload.bin', buffer) })
            .then(function (payload) {
              setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: '', result: payload }) })
              props.onImported(payload)
            })
            .catch(function (error) {
              setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: formatError(error), candidates: [] }) })
            })
          return
        }
        var payload = options()
        if (!plainText(payload.path) && !plainText(payload.url)) {
          setBusy(function (previous) { return Object.assign({}, previous, { error: T('importNoPath'), candidates: [] }) })
          return
        }
        setBusy(function (previous) { return Object.assign({}, previous, { import: true, error: '', candidates: [] }) })
        var call = source === 'url'
          ? props.api.openDataset({ url: plainText(payload.url), name: payload.name, data_format: payload.data_format === 'auto' ? '' : payload.data_format, species: payload.species })
          : props.api.openDataset(payload)
        call.then(function (response) {
          setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: '', result: response }) })
          props.onImported(response)
        }).catch(function (error) {
          if (error && error.code === 'ambiguous') { handleAmbiguous(error); return }
          setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: formatError(error), candidates: [] }) })
        })
      }
      function inspect() {
        var payload = options()
        if (!plainText(payload.path) && !plainText(payload.url)) {
          setBusy(function (previous) { return Object.assign({}, previous, { inspectError: T('importNoPath') }) })
          return
        }
        setBusy(function (previous) { return Object.assign({}, previous, { inspect: true, inspectError: '', result: null }) })
        props.api.inspectSource(payload)
          .then(function (response) {
            setBusy(function (previous) { return Object.assign({}, previous, { inspect: false, inspectError: '', result: response }) })
          })
          .catch(function (error) {
            setBusy(function (previous) { return Object.assign({}, previous, { inspect: false, inspectError: formatError(error) }) })
          })
      }
      function pickCandidate(candidate) {
        if (source === 'local' || source === 'url') {
          if (/\.xlsx?$/i.test(candidate)) patchForm({ sheet: candidate })
          else patchForm({ member: candidate })
        }
        setBusy(function (previous) { return Object.assign({}, previous, { candidates: [] }) })
      }
      var sourceTabs = [
        { id: 'local', labelKey: 'importSourceLocal' },
        { id: 'url', labelKey: 'importSourceUrl' },
        { id: 'upload', labelKey: 'importSourceUpload' },
      ]
      return h(Modal, {
        title: T('importTitle'),
        subtitle: T('importSubtitle'),
        busy: busy.import,
        onClose: props.onClose,
      },
        h('div', { className: 'bg-modal-body' },
          h('div', { className: 'bg-inline-row' }, sourceTabs.map(function (item) {
            return h('button', {
              key: item.id,
              type: 'button',
              className: 'bg-multi-opt' + (source === item.id ? ' bg-multi-opt-on' : ''),
              onClick: function () { setSource(item.id) },
            }, T(item.labelKey))
          })),
          source === 'upload'
            ? h('div', { className: 'bg-field bg-field-col' },
                h('span', { className: 'bg-field-label' }, T('importChooseFile')),
                h('input', {
                  type: 'file',
                  className: 'bg-input',
                  style: { height: 'auto', padding: '5px 8px' },
                  onChange: function (event) { setFile(event.target.files && event.target.files[0] ? event.target.files[0] : null) },
                }),
                h('div', { className: 'bg-muted' }, T('importUploadHint')))
            : h('div', { className: 'bg-field bg-field-col' },
                h('span', { className: 'bg-field-label' }, source === 'url' ? T('importSourceUrl') : T('importSourceLocal')),
                h('input', {
                  className: 'bg-input',
                  value: source === 'url' ? form.url : form.path,
                  spellCheck: false,
                  placeholder: source === 'url' ? T('importUrlPh') : T('importPathPh'),
                  onChange: function (event) { patchForm(source === 'url' ? { url: event.target.value } : { path: event.target.value }) },
                })),
          h('div', { className: 'bg-form-grid' },
            h(FormField, { label: T('importFormat') },
              h('select', {
                className: 'bg-input bg-select',
                value: form.data_format,
                onChange: function (event) { patchForm({ data_format: event.target.value }) },
              }, h('option', { value: 'auto' }, T('importAuto')), formats.map(function (item) {
                return h('option', { key: item.value, value: item.value }, item.label)
              }))),
            h(FormField, { label: T('importSpecies') },
              h('select', {
                className: 'bg-input bg-select',
                value: form.species,
                onChange: function (event) { patchForm({ species: event.target.value }) },
              },
                h('option', { value: 'auto' }, T('importAuto')),
                h('option', { value: 'human' }, T('geoHuman')),
                h('option', { value: 'mouse' }, T('geoMouse'))))),
          h('details', { className: 'bg-advanced' },
            h('summary', null, T('importAdvanced')),
            h('div', { className: 'bg-form-grid', style: { marginTop: '8px' } },
              h(FormField, { label: T('importSheet') },
                h('input', { className: 'bg-input', value: form.sheet, spellCheck: false, onChange: function (event) { patchForm({ sheet: event.target.value }) } })),
              h(FormField, { label: T('importGeneCol') },
                h('input', { className: 'bg-input', value: form.gene_col, spellCheck: false, onChange: function (event) { patchForm({ gene_col: event.target.value }) } })),
              h(FormField, { label: T('importMember') },
                h('input', { className: 'bg-input', value: form.member, spellCheck: false, onChange: function (event) { patchForm({ member: event.target.value }) } })),
              h(FormField, { label: T('importObject') },
                h('input', { className: 'bg-input', value: form.object, spellCheck: false, onChange: function (event) { patchForm({ object: event.target.value }) } })),
              h(FormField, { label: T('importSampleSheet') },
                h('input', { className: 'bg-input', value: form.sample_sheet, spellCheck: false, placeholder: '/absolute/path/samples.csv', onChange: function (event) { patchForm({ sample_sheet: event.target.value }) } })),
              h(FormField, { label: T('importClinical') },
                h('input', { className: 'bg-input', value: form.clinical, spellCheck: false, placeholder: '/absolute/path/clinical.csv', onChange: function (event) { patchForm({ clinical: event.target.value }) } })),
              h('label', { className: 'bg-check' },
                h('input', { type: 'checkbox', checked: form.transpose, onChange: function (event) { patchForm({ transpose: event.target.checked }) } }),
                h('span', null, T('importTranspose'))))),
          h('div', { className: 'bg-muted' }, T('importHint')),
          busy.error ? h(InlineError, { message: T('importError') + '：' + busy.error }) : null,
          busy.candidates.length
            ? h('div', { className: 'bg-card' },
                h('div', { className: 'bg-card-title' }, T('importCandidates')),
                h('div', { className: 'bg-multi' }, busy.candidates.map(function (candidate) {
                  return h('button', {
                    key: candidate,
                    type: 'button',
                    className: 'bg-multi-opt',
                    onClick: function () { pickCandidate(candidate) },
                  }, candidate)
                })))
            : null,
          busy.inspectError ? h(InlineError, { message: busy.inspectError }) : null,
          busy.result && busy.result.format
            ? h('div', { className: 'bg-card' },
                h('div', { className: 'bg-card-title' }, T('importInspectTitle')),
                h('div', { className: 'bg-inline-row' },
                  h('span', { className: 'bg-geo-tag' }, String(busy.result.format)),
                  !isNil(busy.result.n_samples) ? h('span', { className: 'bg-geo-tag' }, busy.result.n_samples + ' ' + T('samples')) : null,
                  !isNil(busy.result.n_genes) ? h('span', { className: 'bg-geo-tag' }, busy.result.n_genes + ' ' + T('genes')) : null,
                  !isNil(busy.result.size) ? h('span', { className: 'bg-geo-tag' }, formatBytes(busy.result.size)) : null))
            : null),
        h('div', { className: 'bg-modal-foot' },
          Btn({
            small: true,
            icon: 'flask',
            disabled: busy.import || busy.inspect,
            title: T('demoLoad'),
            onClick: function () {
              setBusy(function (previous) { return Object.assign({}, previous, { import: true, error: '', candidates: [] }) })
              props.api.demo()
                .then(function (payload) {
                  setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: '', result: payload }) })
                  props.onImported(payload)
                })
                .catch(function (error) {
                  setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: formatError(error), candidates: [] }) })
                })
            },
          }, T('demoLoad')),
          Btn({ small: true, icon: 'search', disabled: busy.import || busy.inspect, onClick: inspect }, busy.inspect ? T('importInspectBusy') : T('importInspect')),
          Btn({ variant: 'primary', icon: 'import', disabled: busy.import || (source === 'upload' ? !file : (!plainText(form.path) && !plainText(form.url))), onClick: submit }, busy.import ? T('importBusy') : T('importStart')),
          Btn({ icon: 'close', onClick: props.onClose, disabled: busy.import }, T('close'))))
    }

    // ==========================================================================
    // 分析目录抽屉：GET /geobulk/analyses 分组 + 通用参数表单 + POST /geobulk/analysis
    // ==========================================================================

    /** 分析参数默认值（数组复制，避免多个分析共享同一引用）。 */
    function analysisDefaults(descriptors) {
      var out = {}
      var list = Array.isArray(descriptors) ? descriptors : []
      for (var index = 0; index < list.length; index++) {
        var descriptor = list[index]
        if (!descriptor || !descriptor.key) continue
        var value = descriptor.default
        if (Array.isArray(value)) value = value.slice()
        out[descriptor.key] = value === undefined ? null : value
      }
      return out
    }

    function AnalysisDrawer(props) {
      var api = props.api
      var dataset = props.dataset
      var datasetId = dataset ? dataset.id : ''
      var catalogState = useState({ loading: true, error: '', groups: [], analyses: [] })
      var catalog = catalogState[0]
      var setCatalog = catalogState[1]
      var selectState = useState('')
      var selectedId = selectState[0]
      var setSelectedId = selectState[1]
      var valuesState = useState({})
      var values = valuesState[0]
      var setValues = valuesState[1]
      var runState = useState({ busy: false, error: '', result: null, id: '' })
      var run = runState[0]
      var setRun = runState[1]

      useEffect(function () {
        if (!datasetId) {
          setCatalog({ loading: false, error: '', groups: [], analyses: [] })
          return undefined
        }
        var cancelled = false
        setCatalog(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        api.analyses(datasetId)
          .then(function (payload) {
            if (cancelled) return
            var normalized = normalizeAnalyses(payload)
            setCatalog({ loading: false, error: '', groups: normalized.groups, analyses: normalized.analyses })
          })
          .catch(function (error) {
            if (cancelled) return
            setCatalog({ loading: false, error: formatError(error), groups: [], analyses: [] })
          })
        return function () { cancelled = true }
      }, [datasetId])

      var metaList = catalog.analyses.filter(function (item) { return item.id === selectedId })
      var meta = metaList.length ? metaList[0] : null

      function selectAnalysis(item) {
        if (!item || item.available === false) return
        setSelectedId(item.id)
        setValues(analysisDefaults(item.params))
        setRun({ busy: false, error: '', result: null, id: '' })
      }

      function patchValues(patch) {
        setValues(function (previous) { return Object.assign({}, previous, patch || {}) })
      }

      function runAnalysis() {
        if (!meta || !datasetId || run.busy) return
        setRun({ busy: true, error: '', result: null, id: meta.id })
        api.analysis(datasetId, meta.id, values || {})
          .then(function (payload) {
            setRun({ busy: false, error: '', result: normalizeResult(payload), id: meta.id })
            if (typeof props.onRan === 'function') props.onRan(meta.id)
          })
          .catch(function (error) {
            setRun({ busy: false, error: formatError(error), result: null, id: meta.id })
          })
      }

      var grouped = groupBy(catalog.analyses, 'group')
      var groupNames = Array.isArray(catalog.groups) ? catalog.groups.slice() : []
      var extraGroups = Object.keys(grouped)
      for (var index = 0; index < extraGroups.length; index++) {
        if (groupNames.indexOf(extraGroups[index]) < 0) groupNames.push(extraGroups[index])
      }
      var groupNodes = groupNames.map(function (group) {
        var items = grouped[group] || []
        if (!items.length) return null
        return h('div', { key: 'group-' + group, className: 'bg-analysis-group' },
          h('div', { className: 'bg-analysis-group-title' }, group),
          items.map(function (item) {
            var disabled = item.available === false
            return h('button', {
              key: item.id,
              type: 'button',
              className: 'bg-analysis-item' + (selectedId === item.id ? ' bg-analysis-item-on' : ''),
              disabled: disabled,
              title: disabled ? (item.reason || T('analysisBlocked')) : (item.method || item.title),
              onClick: function () { selectAnalysis(item) },
            },
              h('div', { className: 'bg-analysis-main' },
                h('div', { className: 'bg-analysis-name' }, item.title || item.id),
                item.method ? h('div', { className: 'bg-analysis-method' }, item.method) : null),
              h('span', { className: 'bg-geo-tag' + (disabled ? '' : ' bg-geo-tag-ok') },
                disabled ? T('analysisBlocked') : T('analysisAvailable')))
          }))
      })

      var detailNode = null
      if (meta) {
        detailNode = h('div', { className: 'bg-analysis-forms' },
          h('div', { className: 'bg-card-title' },
            h('span', null, meta.title || meta.id),
            meta.method ? h('small', null, meta.method) : null),
          meta.needs && meta.needs.length
            ? h('div', { className: 'bg-inline-row' },
                h('span', { className: 'bg-field-label' }, T('analysisNeeds')),
                meta.needs.map(function (need) {
                  var step = stepById(need)
                  return h('span', { key: 'need-' + need, className: 'bg-geo-tag' }, step ? step.title : need)
                }))
            : null,
          meta.available === false && meta.reason ? h('div', { className: 'bg-warn-text' }, meta.reason) : null,
          h(ParamForm, {
            params: meta.params,
            values: values,
            columns: props.columns || [],
            disabled: run.busy,
            onChange: patchValues,
          }),
          h('div', { className: 'bg-inline-row' },
            Btn({
              small: true,
              variant: 'primary',
              icon: 'play',
              disabled: run.busy || !datasetId,
              onClick: runAnalysis,
            }, run.busy ? T('analysisRunning') : T('analysisRun'))),
          run.error ? h(InlineError, { message: run.error }) : null)
      }

      return h('div', {
        className: 'bg-drawer-mask',
        onMouseDown: function (event) {
          if (event.target === event.currentTarget && !run.busy && props.onClose) props.onClose()
        },
      },
        h('aside', { className: 'bg-drawer' },
          h('div', { className: 'bg-modal-head' },
            h('div', null,
              h('div', { className: 'bg-modal-title' }, T('analysisTitle')),
              h('div', { className: 'bg-modal-subtitle' }, T('analysisHint'))),
            Btn({ small: true, variant: 'ghost', icon: 'close', disabled: run.busy, onClick: props.onClose }, T('close'))),
          h('div', { className: 'bg-drawer-body' },
            !datasetId ? h(Notice, { text: T('analysisNeedDataset') }) : null,
            catalog.loading ? h(LoadingLine, { text: T('loading') }) : null,
            catalog.error ? h(InlineError, { message: catalog.error }) : null,
            !catalog.loading && !catalog.error && !catalog.analyses.length ? h(Notice, { text: T('analysisEmpty') }) : null,
            groupNodes,
            detailNode,
            run.result ? h(SectionCard, { title: T('analysisResult') + (run.id ? ' · ' + run.id : '') }, h(ResultView, { result: run.result })) : null)))
    }

    // ==========================================================================
    // 整页工作台：状态容器 + 三栏布局 + 两个弹层 + 分析抽屉
    //   · 左栏状态由 STEPS + dataset.steps 派生（stepStatusOf）
    //   · 结果信封只保存在本会话（切数据集会清空，后端 steps 仍保留）
    //   · 所有故障单点包 PluginErrorBoundary，不黑整页
    // ==========================================================================

    /** 步骤 → 中栏 tab（跑完自动跳转，方便立刻看到图表）。 */
    var TAB_FOR_STEP = {
      de: 'de',
      enrich: 'enrich',
      immune: 'immune',
      survival: 'survival',
      prognosis: 'prognosis',
      scmap: 'sc',
      ncrna: 'ncrna',
    }

    function pickDataset(datasets, id) {
      var list = Array.isArray(datasets) ? datasets : []
      for (var index = 0; index < list.length; index++) {
        if (list[index].id === id) return list[index]
      }
      return null
    }

    function datasetIdOf(payload) {
      if (!payload) return ''
      if (typeof payload === 'string') return plainText(payload)
      if (typeof payload !== 'object') return ''
      if (payload.dataset && typeof payload.dataset === 'object') return plainText(payload.dataset.id || '')
      return plainText(payload.id || payload.active || '')
    }

    /** 步骤卡上的结果摘要：优先 bullets，其次 metrics，最后图/表计数。 */
    function resultSummaryText(result) {
      if (!result) return ''
      var bullets = Array.isArray(result.bullets) ? result.bullets : []
      if (bullets.length) return bullets.slice(0, 2).join('；')
      var metrics = result.metrics && typeof result.metrics === 'object' ? result.metrics : {}
      var keys = Object.keys(metrics)
      if (keys.length) {
        var parts = []
        for (var index = 0; index < keys.length && parts.length < 4; index++) {
          var value = metrics[keys[index]]
          parts.push(keys[index] + ' ' + (typeof value === 'number' ? fmtNumber(value) : formatCell(value, false)))
        }
        return parts.join(' · ')
      }
      var counts = []
      if (result.tables && result.tables.length) counts.push(result.tables.length + ' ' + T('analysisTables'))
      if (result.figures && result.figures.length) counts.push(result.figures.length + ' ' + T('analysisFigures'))
      return counts.join(' · ')
    }

    function WorkbenchBody() {
      var api = useMemo(function () { return createApi() }, [])
      var bootState = useState({ loading: true, error: '', health: normalizeHealth({}), formats: null })
      var boot = bootState[0]
      var setBoot = bootState[1]
      var bootNonceState = useState(0)
      var bootNonce = bootNonceState[0]
      var setBootNonce = bootNonceState[1]
      var stateState = useState(function () { return normalizeState({}) })
      var snapshot = stateState[0]
      var setSnapshot = stateState[1]
      var previewState = useState(function () { return normalizePreview({}) })
      var preview = previewState[0]
      var setPreview = previewState[1]
      var valuesState = useState({})
      var stepValues = valuesState[0]
      var setStepValues = valuesState[1]
      var resultsState = useState({})
      var stepResults = resultsState[0]
      var setStepResults = resultsState[1]
      var summariesState = useState({})
      var stepSummaries = summariesState[0]
      var setStepSummaries = summariesState[1]
      var errorsState = useState({})
      var stepErrors = errorsState[0]
      var setStepErrors = errorsState[1]
      var busyStepState = useState('')
      var stepBusy = busyStepState[0]
      var setStepBusy = busyStepState[1]
      var failedState = useState('')
      var failedStep = failedState[0]
      var setFailedStep = failedState[1]
      var plotsState = useState({})
      var plots = plotsState[0]
      var setPlots = plotsState[1]
      var exportsState = useState({})
      var exportsMap = exportsState[0]
      var setExports = exportsState[1]
      var tabState = useState('overview')
      var tab = tabState[0]
      var setTab = tabState[1]
      var dialogState = useState({ geo: false, import: false, analysis: false })
      var dialogs = dialogState[0]
      var setDialogs = dialogState[1]
      var copiedState = useState('')
      var copiedKey = copiedState[0]
      var setCopiedKey = copiedState[1]
      var refreshingState = useState(false)
      var refreshing = refreshingState[0]
      var setRefreshing = refreshingState[1]
      var busyState = useState(false)
      var busy = busyState[0]
      var setBusy = busyState[1]
      var bannerState = useState({ kind: '', text: '' })
      var banner = bannerState[0]
      var setBanner = bannerState[1]
      var activeRef = useRef('')
      var copyTimerRef = useRef(null)

      // 首屏：health + state + formats 并发拉取；bootNonce 用于「重试」。
      useEffect(function () {
        var cancelled = false
        setBoot(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        Promise.all([
          api.getHealth().catch(function () { return null }),
          api.getState(),
          api.getFormats().catch(function () { return null }),
        ]).then(function (list) {
          if (cancelled) return
          var next = normalizeState(list[1])
          var health = normalizeHealth(list[0] || {})
          if (!health.ready && next.engine.ready) health = next.engine
          setSnapshot(next)
          setBoot({ loading: false, error: '', health: health, formats: list[2] })
        }).catch(function (error) {
          if (cancelled) return
          setBoot({ loading: false, error: formatError(error), health: normalizeHealth({}), formats: null })
        })
        return function () { cancelled = true }
      }, [bootNonce, api])

      // 切换数据集：清掉本会话的图/结果，重新拉 preview（列清单）。
      useEffect(function () {
        if (activeRef.current === snapshot.active) return undefined
        activeRef.current = snapshot.active
        setStepResults({})
        setStepSummaries({})
        setStepErrors({})
        setStepBusy('')
        setFailedStep('')
        setExports({})
        setPlots({})
        setStepValues({})
        setCopiedKey('')
        if (!snapshot.active) {
          setPreview(normalizePreview({}))
          return undefined
        }
        var cancelled = false
        api.preview(snapshot.active)
          .then(function (payload) {
            if (!cancelled) setPreview(normalizePreview(payload))
          })
          .catch(function () {
            if (!cancelled) setPreview(normalizePreview({}))
          })
        return function () { cancelled = true }
      }, [snapshot.active, api])

      useEffect(function () {
        return function () {
          if (copyTimerRef.current) clearTimeout(copyTimerRef.current)
        }
      }, [])

      var dataset = pickDataset(snapshot.datasets, snapshot.active)
      var health = snapshot.engine && snapshot.engine.ready ? snapshot.engine : boot.health
      var columns = preview.sample_columns && preview.sample_columns.length
        ? preview.sample_columns
        : (dataset && dataset.sample_columns ? dataset.sample_columns : [])

      var statuses = {}
      var paramValues = {}
      var stepsLoop = 0
      for (stepsLoop = 0; stepsLoop < STEPS.length; stepsLoop++) {
        var stepMeta = STEPS[stepsLoop]
        statuses[stepMeta.id] = stepStatusOf(stepMeta, dataset, stepBusy, failedStep)
        paramValues[stepMeta.id] = Object.assign({}, defaultParams(stepMeta), stepValues[stepMeta.id] || {})
      }

      function refreshState(options) {
        var opts = options && typeof options === 'object' ? options : {}
        if (!opts.silent) setRefreshing(true)
        return api.getState().then(function (payload) {
          var next = normalizeState(payload)
          setSnapshot(next)
          setBoot(function (previous) {
            var merged = previous.health
            if (next.engine.ready && !merged.ready) merged = next.engine
            return Object.assign({}, previous, { health: merged })
          })
          if (!opts.silent) setRefreshing(false)
          return next
        }).catch(function (error) {
          if (!opts.silent) {
            setRefreshing(false)
            setBanner({ kind: 'bad', text: formatError(error) })
          }
          return null
        })
      }

      function refreshPreview(targetId) {
        var id = plainText(targetId) || snapshot.active
        if (!id) {
          setPreview(normalizePreview({}))
          return Promise.resolve(null)
        }
        return api.preview(id).then(function (payload) {
          var normalized = normalizePreview(payload)
          setPreview(normalized)
          return normalized
        }).catch(function () { return null })
      }

      function loadPlot(kind, params, force) {
        var id = snapshot.active
        if (!id) return
        var safeParams = params && typeof params === 'object' ? params : {}
        var key = plotKey(kind, safeParams)
        var entry = plots[key]
        if (!force && entry && (entry.loading || entry.figure)) return
        setPlots(function (previous) {
          var next = Object.assign({}, previous)
          next[key] = { loading: true, error: '', figure: null }
          return next
        })
        api.plot(id, kind, safeParams).then(function (payload) {
          var raw = payload && payload.figure ? payload.figure : payload
          var figure = normalizeFigure(raw)
          setPlots(function (previous) {
            var next = Object.assign({}, previous)
            next[key] = {
              loading: false,
              error: figure ? '' : T('noFigure'),
              figure: figure,
            }
            return next
          })
        }).catch(function (error) {
          setPlots(function (previous) {
            var next = Object.assign({}, previous)
            next[key] = { loading: false, error: formatError(error), figure: null }
            return next
          })
        })
      }

      function callAnalysis(analysisId, params) {
        var id = snapshot.active
        if (!id) return Promise.reject(new Error(T('analysisNeedDataset')))
        return api.analysis(id, analysisId, params || {}).then(function (payload) {
          refreshState({ silent: true })
          return normalizeResult(payload)
        })
      }

      function runExport(what) {
        var id = snapshot.active
        if (!id) return
        setExports(function (previous) {
          var next = Object.assign({}, previous)
          next[what] = Object.assign({}, previous[what] || {}, { busy: true, error: '' })
          return next
        })
        api.exportData(id, what).then(function (payload) {
          var file = plainText(payload && (payload.file || payload.path)) ||
            plainText(payload && payload.result && payload.result.file)
          setExports(function (previous) {
            var next = Object.assign({}, previous)
            next[what] = { busy: false, file: file, error: '' }
            return next
          })
          refreshState({ silent: true })
        }).catch(function (error) {
          setExports(function (previous) {
            var next = Object.assign({}, previous)
            next[what] = { busy: false, file: '', error: formatError(error) }
            return next
          })
        })
      }

      function handleRunStep(stepId, params) {
        var step = stepById(stepId)
        var id = snapshot.active
        if (!step || !id) return Promise.resolve(null)
        var requestParams = params && typeof params === 'object' && !Array.isArray(params)
          ? params
          : collectStepParams(step, stepValues[stepId] || {})
        setStepBusy(stepId)
        setFailedStep('')
        setBanner({ kind: '', text: '' })
        setStepErrors(function (previous) {
          var next = Object.assign({}, previous)
          delete next[stepId]
          return next
        })
        return api.runStep(id, stepId, requestParams).then(function (payload) {
          var result = normalizeResult(payload)
          setStepResults(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = result
            return next
          })
          setStepSummaries(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = resultSummaryText(result)
            return next
          })
          setStepValues(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = Object.assign({}, requestParams)
            return next
          })
          setStepBusy('')
          if (TAB_FOR_STEP[stepId]) setTab(TAB_FOR_STEP[stepId])
          refreshState({ silent: true })
          if (stepId === 'qc' || stepId === 'filter' || stepId === 'normalize' || stepId === 'batch') refreshPreview()
          return result
        }).catch(function (error) {
          setStepBusy('')
          setFailedStep(stepId)
          setStepErrors(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = formatError(error)
            return next
          })
          return null
        })
      }

      function handleParamChange(stepId, patch) {
        var step = stepById(stepId)
        setStepValues(function (previous) {
          var next = Object.assign({}, previous)
          next[stepId] = Object.assign({}, step ? defaultParams(step) : {}, previous[stepId] || {}, patch || {})
          return next
        })
      }

      function handleSelectDataset(id) {
        if (!id || id === snapshot.active || busy) return
        setBusy(true)
        setBanner({ kind: '', text: '' })
        api.selectDataset(id).then(function () {
          return refreshState({ silent: true })
        }).then(function () {
          setBusy(false)
        }).catch(function (error) {
          setBusy(false)
          setBanner({ kind: 'bad', text: formatError(error) })
        })
      }

      function handleCloseDataset(id) {
        if (!id || busy) return
        setBusy(true)
        api.closeDataset(id).then(function () {
          return refreshState({ silent: true })
        }).then(function () {
          setBusy(false)
        }).catch(function (error) {
          setBusy(false)
          setBanner({ kind: 'bad', text: formatError(error) })
        })
      }

      function handleImported(payload) {
        var targetId = datasetIdOf(payload)
        setDialogs({ geo: false, import: false, analysis: false })
        setTab('overview')
        setBusy(true)
        refreshState({ silent: true }).then(function (next) {
          var active = next ? next.active : ''
          if (targetId && active && targetId !== active) {
            return api.selectDataset(targetId).then(function () { return refreshState({ silent: true }) })
          }
          return null
        }).then(function () {
          setBusy(false)
          refreshPreview(targetId)
        }).catch(function (error) {
          setBusy(false)
          setBanner({ kind: 'bad', text: formatError(error) })
        })
      }

      function copyText(text, key) {
        var value = plainText(text, 4000)
        if (!value) return
        var markCopied = function () {
          setCopiedKey(key || 'copy')
          if (copyTimerRef.current) clearTimeout(copyTimerRef.current)
          copyTimerRef.current = setTimeout(function () { setCopiedKey('') }, 1800)
        }
        var fallback = function () {
          try {
            if (typeof document === 'undefined' || !document.body) return false
            var area = document.createElement('textarea')
            area.value = value
            area.setAttribute('readonly', 'readonly')
            area.style.position = 'fixed'
            area.style.top = '-1000px'
            area.style.opacity = '0'
            document.body.appendChild(area)
            area.select()
            var ok = typeof document.execCommand === 'function' ? document.execCommand('copy') : false
            document.body.removeChild(area)
            return ok === true
          } catch (ignored) {
            return false
          }
        }
        var clipboard = null
        try { clipboard = typeof navigator === 'undefined' ? null : navigator.clipboard } catch (ignored) { clipboard = null }
        if (clipboard && typeof clipboard.writeText === 'function') {
          Promise.resolve(clipboard.writeText(value)).then(markCopied).catch(function () {
            if (fallback()) markCopied()
            else setBanner({ kind: 'bad', text: T('copyFailed') })
          })
          return
        }
        if (fallback()) markCopied()
        else setBanner({ kind: 'bad', text: T('copyFailed') })
      }

      function renderTab() {
        if (tab === 'expr') {
          return h(ExprTab, {
            dataset: dataset,
            preview: preview,
            plots: plots,
            loadPlot: loadPlot,
            exports: exportsMap,
            callAnalysis: callAnalysis,
            runExport: runExport,
          })
        }
        if (tab === 'de') {
          return h(DeTab, {
            dataset: dataset,
            stepResults: stepResults,
            plots: plots,
            loadPlot: loadPlot,
            stepBusy: stepBusy,
            onRunStep: handleRunStep,
            exports: exportsMap,
            runExport: runExport,
          })
        }
        if (tab === 'enrich') {
          return h(EnrichTab, {
            dataset: dataset,
            stepResults: stepResults,
            plots: plots,
            loadPlot: loadPlot,
            stepBusy: stepBusy,
            preview: preview,
            runStep: handleRunStep,
          })
        }
        if (tab === 'immune') {
          return h(ImmuneTab, {
            dataset: dataset,
            preview: preview,
            stepResults: stepResults,
            stepBusy: stepBusy,
            runStep: handleRunStep,
            plots: plots,
            loadPlot: loadPlot,
          })
        }
        if (tab === 'survival') {
          return h(SurvivalTab, {
            dataset: dataset,
            preview: preview,
            stepResults: stepResults,
            stepBusy: stepBusy,
            runStep: handleRunStep,
            plots: plots,
            loadPlot: loadPlot,
            onOpenImport: function () { setDialogs(function (previous) { return Object.assign({}, previous, { import: true }) }) },
          })
        }
        if (tab === 'prognosis') {
          return h(PrognosisTab, {
            dataset: dataset,
            preview: preview,
            stepResults: stepResults,
            stepBusy: stepBusy,
            runStep: handleRunStep,
            plots: plots,
            loadPlot: loadPlot,
          })
        }
        if (tab === 'sc') {
          return h(ScTab, {
            dataset: dataset,
            api: api,
            plots: plots,
            loadPlot: loadPlot,
            stepBusy: stepBusy,
            stepResults: stepResults,
            runStep: handleRunStep,
            refreshState: function () { refreshState({ silent: true }) },
          })
        }
        if (tab === 'ncrna') {
          return h(NcrnaTab, {
            dataset: dataset,
            api: api,
            stepBusy: stepBusy,
            stepResults: stepResults,
            runStep: handleRunStep,
            refreshState: function () { refreshState({ silent: true }) },
          })
        }
        if (tab === 'records') {
          return h(RecordsTab, {
            dataset: dataset,
            api: api,
            exports: exportsMap,
            runExport: runExport,
            onOpenAnalysis: function () { setDialogs(function (previous) { return Object.assign({}, previous, { analysis: true }) }) },
          })
        }
        if (tab === 'log') {
          return h(LogTab, { api: api })
        }
        return h(OverviewTab, {
          dataset: dataset,
          preview: preview,
          health: health,
          suggestions: snapshot.suggestions,
          busyStep: stepBusy ? (stepById(stepBusy) || { title: stepBusy }).title : '',
          plots: plots,
          loadPlot: loadPlot,
          onRunStep: handleRunStep,
          onOpenImport: function () { setDialogs(function (previous) { return Object.assign({}, previous, { import: true }) }) },
        })
      }

      var bannerNode = banner.text
        ? h('div', { className: 'bg-banner' + (banner.kind === 'bad' ? ' bg-banner-bad' : '') },
            h('span', { className: 'bg-banner-text', title: banner.text }, banner.text),
            Btn({
              small: true,
              variant: 'ghost',
              icon: 'close',
              onClick: function () { setBanner({ kind: '', text: '' }) },
            }, T('close')))
        : null

      var centerNode = boot.loading
        ? h('div', { className: 'bg-pane' }, h(LoadingLine, { text: T('loadingData') }))
        : (boot.error
            ? h('div', { className: 'bg-pane' },
                h(ErrorBanner, {
                  title: T('errorTitle'),
                  message: boot.error,
                  retryLabel: T('retry'),
                  onRetry: function () { setBootNonce(function (previous) { return previous + 1 }) },
                }))
            : h(PluginErrorBoundary, { key: 'pane-' + tab }, renderTab()))

      return h('div', { className: 'bg-root' },
        h('style', { 'data-plugin': NS }, CSS),
        h(PluginErrorBoundary, null,
          h(TopBar, {
            health: health,
            datasets: snapshot.datasets,
            activeId: snapshot.active,
            busy: busy || refreshing,
            refreshing: refreshing,
            onSelect: handleSelectDataset,
            onCloseDataset: handleCloseDataset,
            onGeo: function () { setDialogs(function (previous) { return Object.assign({}, previous, { geo: true }) }) },
            onImport: function () { setDialogs(function (previous) { return Object.assign({}, previous, { import: true }) }) },
            onRefresh: function () {
              refreshState({ silent: false }).then(function (next) { if (next) refreshPreview() })
            },
            onAnalysis: function () { setDialogs(function (previous) { return Object.assign({}, previous, { analysis: true }) }) },
          })),
        bannerNode,
        h('div', { className: 'bg-body' },
          h(PluginErrorBoundary, null,
            h(PipelinePanel, {
              dataset: dataset,
              statuses: statuses,
              values: paramValues,
              columns: columns,
              summaries: stepSummaries,
              errors: stepErrors,
              busyStep: stepBusy,
              onParamChange: handleParamChange,
              onRun: handleRunStep,
            })),
          h('div', { className: 'bg-center' },
            h(PluginErrorBoundary, null, h(TabsBar, { active: tab, onChange: setTab })),
            centerNode),
          h(PluginErrorBoundary, null,
            h(AgentPanel, {
              dataset: dataset,
              health: health,
              suggestions: snapshot.suggestions,
              activity: snapshot.activity,
              copiedKey: copiedKey,
              onCopy: copyText,
            }))),
        dialogs.geo
          ? h(PluginErrorBoundary, { key: 'dialog-geo' },
              h(GeoDialog, {
                api: api,
                onImported: handleImported,
                onClose: function () { setDialogs(function (previous) { return Object.assign({}, previous, { geo: false }) }) },
              }))
          : null,
        dialogs.import
          ? h(PluginErrorBoundary, { key: 'dialog-import' },
              h(ImportDialog, {
                api: api,
                formats: boot.formats,
                onImported: handleImported,
                onClose: function () { setDialogs(function (previous) { return Object.assign({}, previous, { import: false }) }) },
              }))
          : null,
        dialogs.analysis
          ? h(PluginErrorBoundary, { key: 'dialog-analysis' },
              h(AnalysisDrawer, {
                api: api,
                dataset: dataset,
                columns: columns,
                onRan: function () { refreshState({ silent: true }) },
                onClose: function () { setDialogs(function (previous) { return Object.assign({}, previous, { analysis: false }) }) },
              }))
          : null)
    }

    /** 整页入口：整体再包一层错误边界（单点异常不黑整页）。 */
    function WorkbenchPage() {
      return h(PluginErrorBoundary, null, h(WorkbenchBody))
    }

    // ==========================================================================
    // 插件入口：locale 字典 + 两个槽（整页 main + 侧边栏入口 sidebar.panellist）
    // ==========================================================================

    function apply(ctx) {
      var locale = ctx.locale
      if (locale && typeof locale.register === 'function') {
        var registerDictionaries = function () {
          try {
            return locale.register(NS, { zh: DICT_ZH, en: DICT_EN })
          } catch (error) {
            warnPlugin('failed to register locale dictionaries', error)
            return undefined
          }
        }
        try {
          if (typeof ctx.effect === 'function') ctx.effect(registerDictionaries, NS + ': locale')
          else registerDictionaries()
        } catch (error) {
          warnPlugin('failed to install locale effect', error)
        }
        T = createTranslator(locale)
      }

      var slots = ctx.slots
      if (!slots && typeof ctx.get === 'function') {
        try { slots = ctx.get('slots') } catch (error) { slots = null }
      }
      if (!slots || typeof slots.register !== 'function' || typeof slots.inject !== 'function') {
        warnPlugin('slots service unavailable; UI not registered')
        return
      }

      var disposers = []
      function keep(handle) {
        if (typeof handle === 'function') disposers.push(handle)
        return handle
      }

      try {
        if (typeof ctx.effect === 'function') {
          ctx.effect(function () {
            return function () {
              while (disposers.length) {
                var dispose = disposers.pop()
                try { dispose() } catch (ignored) { /* 单个 disposer 失败不影响其它 */ }
              }
            }
          }, NS + ': slot disposers')
        }
      } catch (error) {
        warnPlugin('failed to install slot cleanup', error)
      }

      try {
        keep(slots.inject('main', function () {
          return keep(slots.register({ name: 'main', key: 'geobulk', locale: NS }, WorkbenchPage))
        }))
      } catch (error) {
        warnPlugin('failed to register main slot', error)
      }

      try {
        keep(slots.inject('sidebar.panellist', function () {
          return keep(slots.register({
            name: 'sidebar.panellist',
            id: 'geobulk',
            order: 4,
            label: function () { return T('nav') },
            locale: NS,
          }, SidebarIcon))
        }))
      } catch (error) {
        warnPlugin('failed to register sidebar.panellist slot', error)
      }
    }

    // __test__ 只给自测/联调用（宿主不消费）：纯函数 + 组件引用。
    return {
      name: NS,
      inject: ['slots', 'locale'],
      apply: apply,
      __test__: Object.assign({}, BG_TEST, {
        components: {
          WorkbenchPage: WorkbenchPage,
          WorkbenchBody: WorkbenchBody,
          AnalysisDrawer: AnalysisDrawer,
          PipelinePanel: PipelinePanel,
          StepCard: StepCard,
          DatasetCard: DatasetCard,
          TopBar: TopBar,
          AgentPanel: AgentPanel,
          TabsBar: TabsBar,
          OverviewTab: OverviewTab,
          ExprTab: ExprTab,
          DeTab: DeTab,
          EnrichTab: EnrichTab,
          ImmuneTab: ImmuneTab,
          SurvivalTab: SurvivalTab,
          PrognosisTab: PrognosisTab,
          ScTab: ScTab,
          NcrnaTab: NcrnaTab,
          RecordsTab: RecordsTab,
          LogTab: LogTab,
          GeoDialog: GeoDialog,
          ImportDialog: ImportDialog,
          ResultView: ResultView,
          ParamForm: ParamForm,
          PluginErrorBoundary: PluginErrorBoundary,
          SidebarIcon: SidebarIcon,
        },
      }),
    }
  }

  // ==========================================================================
  // 顶层守卫：node 里只暴露纯函数给单测；浏览器里注册 Client 模块。
  // ==========================================================================

  if (typeof window === 'undefined') {
    globalThis.__BG_TEST__ = BG_TEST
  } else if (window.__ModuleLoader__ && typeof window.__ModuleLoader__.load === 'function') {
    window.__ModuleLoader__.load({ id: NS, factory: factory })
  }
})()
