/**
 * dsh-geobulk 浏览器端纯函数单测（不依赖浏览器、不装 npm 依赖）。
 *
 * 运行：node test/client-logic.test.mjs
 *
 * client.js 是浏览器模块（window.__ModuleLoader__.load），但在 node 里
 * （typeof window === 'undefined'）只会把纯函数挂到 globalThis.__BG_TEST__，
 * 所以这里直接 import 它、断言纯函数行为 + 静态合规约束（§6.1/§6.3/§6.4）。
 */
import assert from 'node:assert/strict'
import { readFile } from 'node:fs/promises'
import { test } from 'node:test'

const CLIENT_URL = new URL('../client.js', import.meta.url)
await import(CLIENT_URL.href)

const BG = globalThis.__BG_TEST__
assert.ok(BG, 'node 环境应通过 globalThis.__BG_TEST__ 暴露纯函数')

const CLIENT_SOURCE = await readFile(CLIENT_URL, 'utf8')

// ---------------------------------------------------------------------------
// 0. 加载与导出形状
// ---------------------------------------------------------------------------

test('import：无 window 时不炸，且暴露 NS/API/纯函数层', () => {
  assert.equal(typeof globalThis.window, 'undefined')
  assert.equal(BG.NS, 'dsh-geobulk')
  assert.equal(BG.API, '/geobulk')
  assert.equal(BG.DEFAULT_GEO_DIR, '~/.dsh/geobulk/GEO')
  assert.ok(BG.POLL_MS > 0)
  assert.ok(BG.DICT_ZH.appTitle && BG.DICT_EN.appTitle)
  assert.ok(BG.DICT_ZH.nav && BG.DICT_EN.nav)
  assert.ok(BG.DICT_ZH.analysisTitle && BG.DICT_EN.analysisTitle)
  assert.equal(typeof BG.apiClient, 'function')
  assert.equal(typeof BG.jobPolling, 'function')
  assert.equal(typeof BG.jobDone, 'function')
  assert.equal(typeof BG.pcaScreen, 'function')
  assert.equal(typeof BG.volcanoPoints, 'function')
  assert.equal(typeof BG.buildSearchPayload, 'function')
  assert.equal(typeof BG.suggestionInstruction, 'function')
  assert.equal(typeof BG.plotKey, 'function')
  assert.equal(typeof BG.createTranslator, 'function')
})

test('字典：中英文 key 集完全一致，且所有 T()/tf()/labelKey 都有词条', () => {
  const zhStart = CLIENT_SOURCE.indexOf('var DICT_ZH = {')
  const enStart = CLIENT_SOURCE.indexOf('var DICT_EN = {')
  const enEnd = CLIENT_SOURCE.indexOf('\n  }\n', enStart)
  assert.ok(zhStart > 0 && enStart > zhStart && enEnd > enStart, '两个字典块应存在')
  const zhBlock = CLIENT_SOURCE.slice(zhStart, enStart)
  const enBlock = CLIENT_SOURCE.slice(enStart, enEnd)
  const zhKeys = new Set()
  for (const match of zhBlock.matchAll(/^ {4}([A-Za-z0-9_]+):/gm)) zhKeys.add(match[1])
  const enKeys = new Set()
  for (const match of enBlock.matchAll(/^ {4}([A-Za-z0-9_]+):/gm)) enKeys.add(match[1])
  assert.ok(zhKeys.size > 300, `字典词条太少：${zhKeys.size}`)
  assert.deepEqual([...zhKeys].filter((key) => !enKeys.has(key)), [], 'EN 缺中文侧的 key')
  assert.deepEqual([...enKeys].filter((key) => !zhKeys.has(key)), [], 'ZH 缺英文侧的 key')

  const used = new Set()
  for (const match of CLIENT_SOURCE.matchAll(/\bT\('([A-Za-z0-9_]+)'\)/g)) used.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/tf\(T, '([A-Za-z0-9_]+)'/g)) used.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/labelKey: '([A-Za-z0-9_]+)'/g)) used.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/titleKey: '([A-Za-z0-9_]+)'/g)) used.add(match[1])
  // 动态 key：STEPS 状态 / 引擎徽章 / EXPORT_KINDS / NCRNA_MODULES / 审计等级 / TAB_META
  for (const key of [
    'stepStatusDone', 'stepStatusRunning', 'stepStatusReady', 'stepStatusBlocked', 'stepStatusFailed',
    'engineReady', 'engineBooting', 'engineBad',
    'exportExprCsv', 'exportDeCsv', 'exportEnrichCsv', 'exportImmuneCsv', 'exportSampleCsv',
    'exportSummary', 'exportProtocol', 'exportMethods', 'exportReport', 'exportScript',
    'moduleCirc', 'moduleMirna', 'moduleLncrna', 'moduleFusion', 'moduleSplicing',
    'auditOk', 'auditWarn', 'auditFail',
    'tabOverview', 'tabExpr', 'tabDe', 'tabEnrich', 'tabImmune', 'tabSurvival',
    'tabPrognosis', 'tabSc', 'tabNcrna', 'tabRecords', 'tabLog',
  ]) used.add(key)
  const missing = [...used].filter((key) => !zhKeys.has(key)).sort()
  assert.deepEqual(missing, [], `字典缺 key：${missing.join(', ')}`)
  const missingEn = [...used].filter((key) => !enKeys.has(key)).sort()
  assert.deepEqual(missingEn, [], `EN 字典缺 key：${missingEn.join(', ')}`)
})

// ---------------------------------------------------------------------------
// 1. STEPS / TAB（§2.6 / §6.2 顺序）
// ---------------------------------------------------------------------------

test('STEPS：12 个步骤顺序与 CONTRACT §2.6 一致', () => {
  assert.deepEqual(BG.STEPS.map((step) => step.id), [
    'qc', 'filter', 'normalize', 'batch', 'pca', 'de', 'enrich',
    'immune', 'survival', 'prognosis', 'scmap', 'ncrna',
  ])
  const requires = {}
  for (const step of BG.STEPS) requires[step.id] = step.requires
  assert.deepEqual(requires, {
    qc: [],
    filter: ['qc'],
    normalize: ['filter'],
    batch: ['normalize'],
    pca: ['normalize'],
    de: ['normalize'],
    enrich: ['de'],
    immune: ['normalize'],
    survival: ['normalize'],
    prognosis: ['survival'],
    scmap: ['de'],
    ncrna: ['normalize'],
  })
  for (const step of BG.STEPS) {
    assert.ok(step.title, `${step.id} 应有中文标题`)
    assert.ok(Array.isArray(step.params) && step.params.length > 0, `${step.id} 应有参数`)
    for (const descriptor of step.params) {
      assert.ok(descriptor.key, `${step.id} 参数缺 key`)
      assert.ok([
        'number', 'select', 'boolean', 'column', 'columns', 'text', 'genes', 'multiselect', 'array',
      ].includes(descriptor.type), `${step.id}.${descriptor.key} 用了不支持的控件类型：${descriptor.type}`)
    }
  }
})

test('TAB_KEYS：11 个 tab 与 CONTRACT §6.2 一致，且都有 TAB_META 标题', () => {
  assert.deepEqual(BG.TAB_KEYS, [
    'overview', 'expr', 'de', 'enrich', 'immune', 'survival',
    'prognosis', 'sc', 'ncrna', 'records', 'log',
  ])
  assert.equal(BG.TAB_KEYS.length, 11)
  for (const key of BG.TAB_KEYS) {
    assert.ok(BG.TAB_META[key] && BG.TAB_META[key].titleKey, `${key} 缺 TAB_META`)
  }
  assert.equal(BG.STEPS.length, 12)
  assert.equal(BG.ENRICH_DBS.length, 6)
  assert.equal(BG.IMMUNE_METHODS.length, 8)
  assert.equal(BG.NCRNA_MODULES.length, 5)
  assert.equal(BG.EXPORT_KINDS.length, 10)
})

// ---------------------------------------------------------------------------
// 2. 数值 / P / 字节 / 时间格式化
// ---------------------------------------------------------------------------

test('fmtNumber：千分位与紧凑单位', () => {
  assert.equal(BG.fmtNumber(0), '0')
  assert.equal(BG.fmtNumber(999), '999')
  assert.equal(BG.fmtNumber(1234), '1,234')
  assert.equal(BG.fmtNumber(25000), '25K')
  assert.equal(BG.fmtNumber(1500000), '1.5M')
  assert.equal(BG.fmtNumber(2000000000), '2B')
  assert.equal(BG.fmtNumber(null), '—')
  assert.equal(BG.fmtNumber(undefined), '—')
  assert.equal(BG.fmtNumber(''), '—')
  assert.equal(BG.fmtNumber(Number.NaN), '—')
  assert.equal(BG.fmtBytes, BG.formatBytes)
})

test('fmtPct / fmtP：比例、百分数与 P 值', () => {
  assert.equal(BG.fmtPct(0.5), '50%')
  assert.equal(BG.fmtPct(1), '100%')
  assert.equal(BG.fmtPct(0), '0%')
  assert.equal(BG.fmtPct(37.5), '37.5%')
  assert.equal(BG.fmtPct(null), '—')

  assert.equal(BG.fmtP(0), '0')
  assert.equal(BG.fmtP(0.05), '0.050')
  assert.equal(BG.fmtP(0.005), '0.0050')
  assert.equal(BG.fmtP(0.0005), '5.00e-4')
  assert.match(BG.fmtP(1.2e-30), /^1\.20e-30$/)
  assert.equal(BG.fmtP(null), '—')
})

test('formatBytes / formatDuration / formatClock', () => {
  assert.equal(BG.formatBytes(0), '0 B')
  assert.equal(BG.formatBytes(512), '512 B')
  assert.equal(BG.formatBytes(1536), '1.5 KB')
  assert.equal(BG.formatBytes(1048576), '1.0 MB')
  assert.equal(BG.formatBytes(-1), '—')
  assert.equal(BG.formatBytes(null), '—')

  assert.equal(BG.formatDuration(812), '812 ms')
  assert.equal(BG.formatDuration(1500), '1.5 s')
  assert.equal(BG.formatDuration(65000), '1m 5s')
  assert.equal(BG.formatDuration(null), '—')

  const date = new Date(2024, 0, 2, 3, 4, 5)
  assert.equal(BG.formatClock(date.toISOString()), '01-02 03:04')
  assert.equal(BG.formatClock(null), '—')
  assert.equal(BG.formatError(new Error('引擎未就绪')), '引擎未就绪')
  assert.equal(BG.formatError('失败了'), '失败了')
  assert.equal(BG.formatError(null), '')
})

test('colorForGroup：同名同色、数字当调色板下标', () => {
  assert.equal(BG.colorForGroup('case'), BG.colorForGroup('case'))
  assert.ok(BG.PALETTE.includes(BG.colorForGroup('case')))
  assert.ok(BG.PALETTE.includes(BG.colorForGroup('control')))
  assert.equal(BG.colorForGroup(0), BG.PALETTE[0])
  assert.equal(BG.colorForGroup(2), BG.PALETTE[2])
  assert.equal(BG.colorForGroup(21), BG.PALETTE[1])
  assert.ok(BG.PALETTE.includes(BG.colorForGroup(null)))
})

// ---------------------------------------------------------------------------
// 3. 坐标变换 / 火山图 / 图 URL
// ---------------------------------------------------------------------------

test('pcaScreen：包围盒缩放 + y 轴翻转 + view 平移/缩放', () => {
  const points = [{ x: -1, y: -1 }, { x: 1, y: 1 }]
  const fitted = BG.pcaScreen(points, {}, 100, 100, 10)
  assert.equal(fitted.length, 2)
  assert.equal(fitted[0].x, 10)
  assert.equal(fitted[0].y, 90)
  assert.equal(fitted[1].x, 90)
  assert.equal(fitted[1].y, 10)

  const zoomed = BG.pcaScreen(points, { scale: 2 }, 100, 100, 10)
  assert.equal(zoomed[0].x, -30)
  assert.equal(zoomed[1].x, 130)

  const panned = BG.pcaScreen(points, { tx: 5, ty: -6 }, 100, 100, 10)
  assert.equal(panned[0].x, 15)
  assert.equal(panned[0].y, 84)

  const invalid = BG.pcaScreen([{ x: 1 }, { x: 'a', y: 2 }, [3, 4]], {}, 100, 100, 10)
  assert.equal(invalid.length, 1)
  assert.equal(invalid[0].x, 50)
  assert.equal(invalid[0].y, 50)

  assert.equal(BG.pcaScreen([], {}, 100, 100, 10).length, 0)
  assert.equal(BG.pcaScreen(null).length, 0)
})

test('volcanoPoints：logFC/P 过滤、上下调着色、-log10 纵坐标', () => {
  const rows = [
    { gene: 'A', logFC: 2, 'P.Value': 0.001, 'adj.P.Val': 0.001 },
    { gene: 'B', logFC: -2.5, 'P.Value': 0.01, 'adj.P.Val': 0.01 },
    { gene: 'C', logFC: 0.2, 'P.Value': 0.5, 'adj.P.Val': 0.5 },
    { gene: 'D', logFC: 'bad', 'P.Value': 0.01 },
    { gene: 'E', logFC: 3 },
    { gene: 'F', logFC: 1.5, padj: 0.2 },
  ]
  const points = BG.volcanoPoints(rows, { lfcThreshold: 1, padjThreshold: 0.05 })
  assert.deepEqual(points.map((point) => point.gene), ['A', 'B', 'C', 'F'])
  assert.equal(points[0].regulation, 'up')
  assert.equal(points[0].sig, true)
  assert.equal(points[0].color, BG.VOLCANO_COLORS.up)
  assert.equal(points[1].regulation, 'down')
  assert.equal(points[2].regulation, 'ns')
  assert.equal(points[2].sig, false)
  assert.equal(points[3].sig, false)
  assert.ok(Math.abs(points[0].y - 3) < 1e-9)
  assert.equal(BG.volcanoPoints(null, {}).length, 0)
})

test('figureSrc / normalizeFigure：裸 PNG 名走 /geobulk/figure，绝对路径原样', () => {
  assert.equal(BG.figureSrc('/geobulk/figure?id=a.png'), '/geobulk/figure?id=a.png')
  assert.equal(BG.figureSrc('data:image/png;base64,xx'), 'data:image/png;base64,xx')
  assert.equal(BG.figureSrc('/tmp/analyses/volcano.png'), '/tmp/analyses/volcano.png')
  assert.equal(BG.figureSrc('volcano.png'), '/geobulk/figure?id=volcano.png')
  assert.equal(BG.figureSrc('bad.svg'), '')
  assert.equal(BG.figureSrc(''), '')
  const figure = BG.normalizeFigure({ file: '/tmp/x/pca.png', title: 'PCA' })
  assert.equal(figure.url, '/geobulk/figure?id=pca.png')
  assert.equal(figure.title, 'PCA')
  assert.equal(BG.normalizeFigure({ file: 'bad.svg' }), null)
})

// ---------------------------------------------------------------------------
// 4. GEO 搜索 payload / 建议指令 / 参数收集
// ---------------------------------------------------------------------------

test('buildSearchPayload：字段归一化、limit/page 夹取、actor=user', () => {
  assert.deepEqual(BG.buildSearchPayload({}), {
    term: '', species: 'all', study_type: 'all', limit: 20, page: 0, sort: 'relevance', actor: 'user',
  })
  const payload = BG.buildSearchPayload({
    term: '  胃癌 顺铂 耐药  ',
    species: 'human',
    study_type: 'high throughput sequencing',
    min_samples: '8',
    has_clinical: true,
    limit: 1000,
    page: -3,
    sort: 'samples',
  })
  assert.equal(payload.term, '胃癌 顺铂 耐药')
  assert.equal(payload.species, 'human')
  assert.equal(payload.study_type, 'high throughput sequencing')
  assert.equal(payload.min_samples, 8)
  assert.equal(payload.has_clinical, true)
  assert.equal(payload.limit, 100)
  assert.equal(payload.page, 0)
  assert.equal(payload.sort, 'samples')

  const invalid = BG.buildSearchPayload({ species: 'pig', study_type: 'x', min_samples: 0 })
  assert.equal(invalid.species, 'all')
  assert.equal(invalid.study_type, 'all')
  assert.equal(invalid.min_samples, undefined)
  assert.equal(invalid.has_clinical, undefined)
})

test('suggestionInstruction：生成可直接粘给 Agent 的中文指令', () => {
  const text = BG.suggestionInstruction({ title: '运行差异表达', detail: '按分组做 limma', step: 'de' })
  assert.match(text, /^请在工作台执行「运行差异表达」/)
  assert.match(text, /（步骤：差异表达，step=de）/)
  assert.match(text, /：按分组做 limma/)
  assert.match(text, /。$/)

  const noStep = BG.suggestionInstruction({ title: '先导入数据' })
  assert.match(noStep, /请在工作台执行「先导入数据」/)

  const empty = BG.suggestionInstruction(null)
  assert.match(empty, /^请在工作台执行「继续分析」/)

  const unknownStep = BG.suggestionInstruction({ step: 'mystery' })
  assert.match(unknownStep, /（step=mystery）/)
})

test('defaultParams / collectStepParams：类型收敛与数组隔离', () => {
  const deDefaults = BG.defaultParams('de')
  assert.equal(deDefaults.group_col, 'group')
  assert.equal(deDefaults.method, 'limma')
  assert.equal(deDefaults.paired, false)
  deDefaults.covars.push('age')
  assert.deepEqual(BG.defaultParams('de').covars, [], '默认数组不能共享引用')

  const enrich = BG.defaultParams('enrich')
  assert.deepEqual(enrich.db, ['GO_BP', 'KEGG'])
  assert.ok(Array.isArray(BG.defaultParams('enrich').db))

  const collected = BG.collectStepParams('de', {
    group_col: 'group', case: 'tumor', control: 'normal', method: 'limma',
    covars: 'age, sex', lfc: '1', padj: '0.01', paired: true,
  })
  assert.equal(collected.case, 'tumor')
  assert.deepEqual(collected.covars, ['age', 'sex'])
  assert.equal(collected.lfc, 1)
  assert.equal(collected.padj, 0.01)
  assert.equal(collected.paired, true)
  assert.equal(collected.group_col, 'group')

  const blank = BG.collectStepParams('de', { group_col: '', lfc: '', padj: '' })
  assert.equal(blank.group_col, null)
  assert.equal(blank.lfc, null)
  assert.equal(blank.padj, null)

  assert.deepEqual(BG.parseTokenList('a, b；c d'), ['a', 'b', 'c', 'd'])
  assert.deepEqual(BG.parseTokenList(['x', 'y']), ['x', 'y'])
  assert.deepEqual(BG.parseTokenList(''), [])
})

test('doneStepSet / isStepReady / stepStatusOf：依赖判定与 running 优先', () => {
  const dataset = { id: 'ds1', steps: ['qc', 'filter', 'normalize', 'de'] }
  assert.deepEqual(BG.doneStepSet(dataset), { qc: true, filter: true, normalize: true, de: true })
  assert.equal(BG.isStepReady('batch', dataset), true)
  assert.equal(BG.isStepReady('pca', dataset), true)
  assert.equal(BG.isStepReady('enrich', dataset), true, 'enrich 只依赖 de')
  assert.equal(BG.isStepReady('survival', dataset), true)
  assert.equal(BG.isStepReady('prognosis', dataset), false, 'prognosis 依赖 survival')
  assert.equal(BG.isStepReady('qc', null), false)

  assert.equal(BG.stepStatusOf('qc', dataset, '', ''), 'done')
  assert.equal(BG.stepStatusOf('qc', dataset, 'qc', ''), 'running')
  assert.equal(BG.stepStatusOf('qc', dataset, '', 'qc'), 'failed')
  assert.equal(BG.stepStatusOf('survival', dataset, '', ''), 'ready')
  assert.equal(BG.stepStatusOf('prognosis', dataset, '', ''), 'blocked')
  assert.equal(BG.stepStatusKey('blocked'), 'stepStatusBlocked')
  assert.equal(BG.stepStatusKey('ready'), 'stepStatusReady')
  assert.equal(BG.stepStatusKey('weird'), 'stepStatusBlocked')
})

// ---------------------------------------------------------------------------
// 5. API 客户端：method / url / body / 错误信封
// ---------------------------------------------------------------------------

function makeFetch(payload, options) {
  const opts = options || {}
  const calls = []
  const impl = function (url, init) {
    calls.push({ url: String(url), init: init || {} })
    if (opts.reject) return Promise.resolve(impl).then(() => { throw opts.reject })
    const body = typeof payload === 'function' ? payload(url, init) : payload
    return Promise.resolve({
      ok: opts.httpOk === false ? false : true,
      status: opts.status || 200,
      text: function () { return Promise.resolve(typeof body === 'string' ? body : JSON.stringify(body)) },
    })
  }
  impl.calls = calls
  return impl
}

test('apiClient：GET/POST、同源前缀、JSON body 与 actor', async () => {
  const fake = makeFetch({ ok: true, engine: { ready: true } })
  const api = BG.apiClient(fake)

  await api.getHealth()
  assert.equal(fake.calls.length, 1)
  assert.equal(fake.calls[0].url, '/geobulk/health')
  assert.equal(fake.calls[0].init.method, 'GET')

  await api.runStep('ds1', 'de', { group_col: 'group' })
  const runCall = fake.calls[1]
  assert.equal(runCall.url, '/geobulk/run')
  assert.equal(runCall.init.method, 'POST')
  assert.equal(runCall.init.headers['content-type'], 'application/json')
  assert.deepEqual(JSON.parse(runCall.init.body), {
    dataset: 'ds1', step: 'de', params: { group_col: 'group' }, actor: 'user',
  })

  await api.geoSearch({ term: '胃癌 顺铂', limit: 5 })
  const searchCall = fake.calls[2]
  assert.equal(searchCall.url, '/geobulk/geo/search')
  const searchBody = JSON.parse(searchCall.init.body)
  assert.equal(searchBody.term, '胃癌 顺铂')
  assert.equal(searchBody.limit, 5)
  assert.equal(searchBody.actor, 'user')

  await api.geoDownload({ gse: 'gse13861', files: ['a.txt'], dest_dir: '~/.dsh/geobulk/GEO' })
  const downloadBody = JSON.parse(fake.calls[3].init.body)
  assert.equal(downloadBody.gse, 'GSE13861')
  assert.deepEqual(downloadBody.files, ['a.txt'])
  assert.equal(downloadBody.dest_dir, '~/.dsh/geobulk/GEO')
  assert.equal(fake.calls[3].url, '/geobulk/geo/download')

  await api.getJob('job 1')
  assert.equal(fake.calls[4].url, '/geobulk/jobs?id=job%201')
  assert.equal(fake.calls[4].init.method, 'GET')

  await api.analyses('ds1')
  assert.equal(fake.calls[5].url, '/geobulk/analyses?dataset=ds1')

  await api.selectDataset('ds 2')
  assert.deepEqual(JSON.parse(fake.calls[6].init.body), { dataset: 'ds 2', actor: 'user' })

  await api.exportData('ds1', 'de_csv')
  assert.deepEqual(JSON.parse(fake.calls[7].init.body), { dataset: 'ds1', what: 'de_csv', actor: 'user' })

  await api.scAttach({ dataset: 'ds1', ref_path: '/tmp/ref.h5ad' })
  assert.deepEqual(JSON.parse(fake.calls[8].init.body), { dataset: 'ds1', ref_path: '/tmp/ref.h5ad', actor: 'user' })
})

test('apiClient：ok:false 信封与 HTTP 失败都抛带 code 的错误', async () => {
  const soft = makeFetch({ ok: false, code: 'bad_request', error: '参数错误' })
  const apiOdd = BG.apiClient(soft)
  await assert.rejects(apiOdd.getState(), (error) => {
    assert.equal(error.message, '参数错误')
    assert.equal(error.code, 'bad_request')
    return true
  })

  const hard = makeFetch({ ok: false, code: 'internal', error: '引擎崩了' }, { httpOk: false, status: 500 })
  const apiBad = BG.apiClient(hard)
  await assert.rejects(apiBad.getHealth(), (error) => {
    assert.equal(error.message, '引擎崩了')
    assert.equal(error.status, 500)
    return true
  })
})

test('apiClient：inspect/samples/attach/plot/export 的路径与 body 合法', async () => {
  const fake = makeFetch({ ok: true })
  const api = BG.apiClient(fake)
  await api.inspectSource({ path: '/tmp/m.csv', data_format: 'matrix_text', transpose: true })
  assert.deepEqual(JSON.parse(fake.calls[0].init.body), {
    path: '/tmp/m.csv', data_format: 'matrix_text', transpose: true, actor: 'user',
  })
  await api.plot('ds1', 'volcano', {})
  assert.equal(fake.calls[1].url, '/geobulk/plot')
  assert.deepEqual(JSON.parse(fake.calls[1].init.body), { dataset: 'ds1', kind: 'volcano', params: {}, actor: 'user' })
  await api.analysis('ds1', 'immune', { methods: ['epic'] })
  assert.deepEqual(JSON.parse(fake.calls[2].init.body), {
    dataset: 'ds1', id: 'immune', params: { methods: ['epic'] }, actor: 'user',
  })
  await api.attach('ds1', 'circ', '/tmp/circ.txt')
  assert.deepEqual(JSON.parse(fake.calls[3].init.body), {
    dataset: 'ds1', module: 'circ', path: '/tmp/circ.txt', actor: 'user',
  })
  await api.recipe('ds1', 'standard')
  assert.deepEqual(JSON.parse(fake.calls[4].init.body), { dataset: 'ds1', id: 'standard', actor: 'user' })
})

test('plotKey：参数顺序无关、空值忽略、不同参数不同 key', () => {
  assert.equal(BG.plotKey('volcano', { lfc: 1, padj: 0.05 }), BG.plotKey('volcano', { padj: 0.05, lfc: 1 }))
  assert.equal(BG.plotKey('gene_box', { genes: ['A', 'B'] }), BG.plotKey('gene_box', { genes: ['A', 'B'] }))
  assert.notEqual(BG.plotKey('gene_box', { genes: ['A'] }), BG.plotKey('gene_box', { genes: ['A', 'B'] }))
  assert.equal(BG.plotKey('pca_scatter', { group: null }), 'pca_scatter')
  assert.equal(BG.plotKey('pca_scatter', {}), 'pca_scatter')
  assert.equal(BG.plotKey('', {}), '')
})

// ---------------------------------------------------------------------------
// 6. 异步任务轮询状态机
// ---------------------------------------------------------------------------

test('jobPolling：running / done / failed / timeout 四态', () => {
  const running = BG.jobPolling({ job: { id: 'j1', status: 'running', pct: 42, message: '下载中' } })
  assert.equal(running.status, 'running')
  assert.equal(running.running, true)
  assert.equal(running.pct, 42)
  assert.equal(running.done, false)

  const done = BG.jobPolling({ job: { id: 'j1', status: 'done', pct: 100, result: { dir: '/tmp/geo' } } })
  assert.equal(done.status, 'done')
  assert.equal(done.done, true)
  assert.deepEqual(done.result, { dir: '/tmp/geo' })
  assert.equal(BG.jobDone({ job: { status: 'done' } }), true)
  assert.equal(BG.jobDone({ job: { status: 'running' } }), false)

  const failed = BG.jobPolling({ ok: false, job: { id: 'j2', status: 'failed', error: '网络错误' } })
  assert.equal(failed.status, 'failed')
  assert.equal(failed.failed, true)
  assert.equal(failed.error, '网络错误')

  const timeout = BG.jobPolling({ job: { status: 'running', pct: 10 } }, { elapsedMs: 999999, timeoutMs: 1000 })
  assert.equal(timeout.status, 'timeout')
  assert.equal(timeout.timedOut, true)
  assert.equal(timeout.running, false)

  const doneWins = BG.jobPolling({ job: { status: 'done' } }, { elapsedMs: 999999, timeoutMs: 1 })
  assert.equal(doneWins.status, 'done')
})

// ---------------------------------------------------------------------------
// 7. 响应归一化（§2.5 信封 / state / analyses）
// ---------------------------------------------------------------------------

test('normalizeResult：summary/tables/figures/warnings 信封', () => {
  const result = BG.normalizeResult({
    ok: true,
    summary: { title: '差异表达', bullets: ['检出 120 个上调基因', '  '], metrics: { up: 120, down: 80 } },
    tables: [{
      id: 'de',
      title: 'DE 表',
      columns: [{ key: 'gene', label: 'gene' }, { key: 'logFC', type: 'number' }],
      rows: [{ gene: 'CDH1', logFC: 2.1 }],
      total: 1,
    }],
    figures: [{ id: 'f1', title: '火山图', file: '/tmp/analyses/ds1/de/volcano.png' }],
    warnings: ['样本数偏少'],
    notes: ['测试'],
    data_out: [{ kind: 'de_csv', path: '/tmp/de.csv', label: 'DE CSV' }],
  })
  assert.equal(result.title, '差异表达')
  assert.deepEqual(result.bullets, ['检出 120 个上调基因'])
  assert.deepEqual(result.metrics, { up: 120, down: 80 })
  assert.equal(result.tables.length, 1)
  assert.equal(result.tables[0].rows.length, 1)
  assert.equal(result.tables[0].columns[1].type, 'number')
  assert.equal(result.figures.length, 1)
  assert.equal(result.figures[0].url, '/geobulk/figure?id=volcano.png')
  assert.deepEqual(result.warnings, ['样本数偏少'])
  assert.deepEqual(result.notes, ['测试'])
  assert.equal(result.data_out[0].path, '/tmp/de.csv')

  const empty = BG.normalizeResult(null)
  assert.deepEqual(empty.bullets, [])
  assert.deepEqual(empty.tables, [])
  assert.deepEqual(empty.figures, [])
})

test('normalizeResult：接受嵌套 result 包装', () => {
  const result = BG.normalizeResult({
    ok: true,
    result: { summary: { title: '富集', bullets: ['KEGG 命中 12 条'] }, tables: [], figures: [] },
    ms: 812,
  })
  assert.equal(result.title, '富集')
  assert.equal(result.ms, 812)
})

test('normalizeState / normalizeDataset / normalizeJob', () => {
  const state = BG.normalizeState({
    active: '',
    datasets: [{
      id: 'ds1', name: 'demo', n_samples: 20, n_genes: 2000,
      steps: ['qc', 'de'], sample_columns: [{ key: 'group', kind: 'categorical', n_unique: 2 }],
      modules: [{ kind: 'circ', n: 12 }], history: [{ at: '2026-01-01T00:00:00Z', step: 'qc', ok: true, ms: 100 }],
      warnings: ['x'], analyses_run: [],
    }],
    engine: { ready: true, mode: 'builtin', python: '3.11', r: { version: '4.3' } },
    suggestions: [{ id: 's1', title: '下一步', step: 'de' }],
    activity: [{ at: '2026-01-01T00:00:00Z', actor: 'agent', action: 'run:qc', detail: '完成', ok: true }],
    sc: { available: false, datasets: [] },
  })
  assert.equal(state.active, 'ds1', '未显式给 active 时回落到第一个数据集')
  assert.equal(state.datasets.length, 1)
  assert.equal(state.datasets[0].n_samples, 20)
  assert.equal(state.datasets[0].sample_columns[0].n_unique, 2)
  assert.equal(state.datasets[0].modules[0].kind, 'circ')
  assert.equal(state.datasets[0].history.length, 1)
  assert.equal(state.engine.ready, true)
  assert.equal(state.engine.r.available, true)
  assert.equal(state.suggestions[0].title, '下一步')
  assert.equal(state.activity[0].action, 'run:qc')
  assert.equal(state.activity[0].ok, true)
  assert.equal(state.sc.available, false)

  const dataset = BG.normalizeDataset({ dataset: { id: 'ds2', steps: ['qc'], warnings: 'a' } })
  assert.equal(dataset.id, 'ds2')
  assert.deepEqual(dataset.warnings, ['a'])

  const job = BG.normalizeJob({ job: { id: 'j', status: 'running', pct: '30', message: '下载' } })
  assert.equal(job.id, 'j')
  assert.equal(job.pct, 30)
  assert.equal(job.message, '下载')
})

test('normalizeAnalyses：分组、参数 schema、可用性', () => {
  const catalog = BG.normalizeAnalyses({
    analyses: [
      { id: 'immune', title: '免疫浸润', group: '临床', method: 'CIBERSORT', needs: ['normalize'], params: [{ key: 'methods', type: 'multiselect', default: ['epic'] }] },
      { id: 'cnv', title: 'CNV', group: '基因组', available: false, reason: '缺数据' },
    ],
  })
  assert.equal(catalog.analyses.length, 2)
  assert.deepEqual(catalog.groups, ['临床', '基因组'])
  assert.equal(catalog.analyses[0].needs[0], 'normalize')
  assert.equal(catalog.analyses[0].params[0].type, 'multiselect')
  assert.equal(catalog.analyses[1].available, false)
  assert.equal(catalog.analyses[1].reason, '缺数据')
})

test('normalizePreview / normalizeAudit / normalizeHealth / engineBadge', () => {
  const preview = BG.normalizePreview({ genes: ['A', '', 'B'], sample_columns: [{ name: 'group', kind: 'categorical' }] })
  assert.deepEqual(preview.genes, ['A', 'B'])
  assert.equal(preview.sample_columns[0].key, 'group')

  const audit = BG.normalizeAudit({ score: 87.6, level: 'warn', checks: [{ id: 'mt', level: 'warn', title: '线粒体', detail: '偏高' }] })
  assert.equal(audit.score, 88)
  assert.equal(audit.checks[0].level, 'warn')

  const health = BG.normalizeHealth({ engine: { ready: true, python: '3.11.9', r: { version: '4.3.1' } } })
  assert.equal(health.ready, true)
  assert.equal(health.python, '3.11.9')
  assert.equal(health.r.available, true)
  assert.equal(BG.engineBadge(health).labelKey, 'engineReady')
  assert.equal(BG.engineBadge({ ready: false }).labelKey, 'engineBooting')
  assert.equal(BG.engineBadge({ ready: false, error: '没有 python' }).labelKey, 'engineBad')
})

// ---------------------------------------------------------------------------
// 8. 静态合规（§6.1 骨架 / §6.3 样式 / §6.4 可测性）
// ---------------------------------------------------------------------------

test('模块骨架：IIFE + ModuleLoader + inject 两槽 + PluginErrorBoundary', () => {
  assert.match(CLIENT_SOURCE, /^\/\*\*/)
  assert.match(CLIENT_SOURCE, /;\(function \(\) \{/)
  assert.ok(CLIENT_SOURCE.trimEnd().endsWith('})()'), '整文件应以 })() 收尾')
  assert.match(CLIENT_SOURCE, /window\.__ModuleLoader__\.load\(\{ id: NS, factory: factory \}\)/)
  assert.match(CLIENT_SOURCE, /globalThis\.__BG_TEST__ = BG_TEST/)
  assert.match(CLIENT_SOURCE, /inject: \['slots', 'locale'\]/)

  // main 槽（key=geobulk）
  assert.match(CLIENT_SOURCE, /slots\.inject\('main',/)
  assert.match(CLIENT_SOURCE, /name: 'main', key: 'geobulk', locale: NS/)
  // 侧边栏入口（id=geobulk, order=4）
  assert.match(CLIENT_SOURCE, /slots\.inject\('sidebar\.panellist',/)
  assert.match(CLIENT_SOURCE, /name: 'sidebar\.panellist',\s*\n\s*id: 'geobulk',\s*\n\s*order: 4,/)
  // locale 注册 zh/en
  assert.match(CLIENT_SOURCE, /locale\.register\(NS, \{ zh: DICT_ZH, en: DICT_EN \}\)/)
  assert.match(CLIENT_SOURCE, /class PluginErrorBoundary extends React\.Component/)
  assert.match(CLIENT_SOURCE, /function WorkbenchPage/)
  assert.match(CLIENT_SOURCE, /function apply\(ctx\)/)

  // 只 require('react')，不碰 @deepseek-ai Client 包
  const requires = [...CLIENT_SOURCE.matchAll(/\brequire\((['"])([^'"]+)\1\)/g)].map((match) => match[2])
  assert.deepEqual([...new Set(requires)], ['react'], '唯一允许的 require 目标是 react')
  assert.ok(!CLIENT_SOURCE.includes('@deepseek-ai/dsh-client'), '不得直接 import 宿主 Client 包')
})

test('样式：类名 bg- 前缀、颜色只走 --dsw-alias-*、圆角/数值/过渡符合 §6.3', () => {
  const cssStart = CLIENT_SOURCE.indexOf('var CSS = [')
  const cssEnd = CLIENT_SOURCE.indexOf('].join(\'\')', cssStart)
  assert.ok(cssStart > 0 && cssEnd > cssStart, 'CSS 段落应存在')
  const css = CLIENT_SOURCE.slice(cssStart, cssEnd)
  assert.ok(css.length > 12000, `CSS 太短：${css.length}`)

  // 选择器一律 .bg- 前缀
  const selectors = [...css.matchAll(/\.([A-Za-z][A-Za-z0-9_-]*)\s*[{,:]/g)].map((match) => match[1])
  const badSelectors = [...new Set(selectors)].filter((name) => !name.startsWith('bg-'))
  assert.deepEqual(badSelectors, [], `CSS 选择器不是 bg- 前缀：${badSelectors.join(', ')}`)

  // 不用固定颜色（数据色在 JS 常量里，不进 CSS）
  const hexColors = [...css.matchAll(/#[0-9a-fA-F]{3,8}\b/g)].map((match) => match[0])
  assert.deepEqual(hexColors, [], `CSS 不应出现硬编码颜色：${hexColors.join(', ')}`)
  assert.ok(!/\brgba?\(/.test(css), 'CSS 不应出现 rgb()/rgba()，请用 color-mix + 主题变量')
  const aliasUses = [...css.matchAll(/--dsw-alias-[a-z0-9-]+/g)].map((match) => match[0])
  assert.ok(aliasUses.length > 120, `主题变量用得不够：${aliasUses.length}`)

  assert.match(css, /border-radius:12px/)
  assert.match(css, /font-variant-numeric:tabular-nums/)
  assert.match(css, /@media \(max-width:980px\)/)
  assert.match(css, /transition:[^;}]*\.1[2-8]s/)
  // 字号 12/13px 为主：9-11px 只允许出现在辅助/标签上
  assert.match(css, /font-size:13px/)
  assert.match(css, /font-size:12px/)
})

test('组件里使用的 bg- 类名都有样式定义（动态拼接前缀除外）', () => {
  const cssStart = CLIENT_SOURCE.indexOf('var CSS = [')
  const cssEnd = CLIENT_SOURCE.indexOf('].join(\'\')', cssStart)
  const css = CLIENT_SOURCE.slice(cssStart, cssEnd)
  const defined = new Set([...css.matchAll(/\.(bg-[A-Za-z0-9_-]+)/g)].map((match) => match[1]))
  const used = new Set()
  for (const match of CLIENT_SOURCE.matchAll(/className: '([^']*)'/g)) {
    for (const token of match[1].split(/\s+/)) {
      if (token.startsWith('bg-') && !token.endsWith('-')) used.add(token)
    }
  }
  const missing = [...used].filter((name) => !defined.has(name)).sort()
  assert.deepEqual(missing, [], `用了未定义的样式类：${missing.join(', ')}`)
})

test('不用 emoji 当图标，图标走内联 SVG / 文字', () => {
  const emoji = CLIENT_SOURCE.match(/[\u{1F000}-\u{1FAFF}\u{1F1E6}-\u{1F1FF}]/gu)
  assert.equal(emoji, null, `出现 emoji：${emoji ? emoji.join('') : ''}`)
  assert.match(CLIENT_SOURCE, /function Icon\(props\)/)
  assert.match(CLIENT_SOURCE, /function SidebarGlyph\(props\)/)
})

test('数据面：客户端只打同源 /geobulk 路径', () => {
  const apiBlock = CLIENT_SOURCE.slice(CLIENT_SOURCE.indexOf('function createApi(fetchImpl, baseUrl)'), CLIENT_SOURCE.indexOf('function buildImportPayload'))
  assert.ok(apiBlock.includes('base + path'), 'apiClient 应统一拼 base 前缀')
  const hardcoded = [...apiBlock.matchAll(/fetch\(\s*['"`]([^'"`]+)['"`]/g)]
  assert.deepEqual(hardcoded, [], '不要硬编码绝对 URL')
  assert.match(apiBlock, /var base = typeof baseUrl === 'string' && baseUrl \? baseUrl\.replace\(\/\\\/\$\/, ''\) : API/)
})

test('纯函数层不碰 DOM / 不 require react', () => {
  const pure = CLIENT_SOURCE.slice(0, CLIENT_SOURCE.indexOf('function factory(require)'))
  assert.ok(pure.length > 50000, `纯函数层太小：${pure.length}`)
  assert.ok(pure.length > CLIENT_SOURCE.length * 0.25, '纯函数层应占可观的比重')
  // 去掉注释再查：头注释里会提到 window.__ModuleLoader__ / factory(require)
  const pureCode = pure.replace(/\/\*[\s\S]*?\*\//g, '').replace(/^\s*\/\/.*$/gm, '')
  assert.ok(!/\bdocument\./.test(pureCode), '纯函数层不应碰 document')
  assert.ok(!/\bwindow\./.test(pureCode), '纯函数层不应碰 window')
  assert.ok(!/\brequire\(/.test(pureCode), '纯函数层不应 require')
  // 浏览器层只在 factory 里
  const browser = CLIENT_SOURCE.slice(CLIENT_SOURCE.indexOf('function factory(require)'))
  assert.match(browser, /require\('react'\)/)
})
