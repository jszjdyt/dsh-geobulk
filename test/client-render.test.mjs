/**
 * 临时烟测（不属于交付物，验证完删除）：用极简 React 运行时把整棵组件树跑一遍，
 * 覆盖 factory 执行、apply 槽注册、首屏 bootstrap、各 tab 渲染。
 */
import assert from 'node:assert/strict'

// ---------- 捕获 console.warn ----------
const warnings = []
const originalWarn = console.warn
console.warn = function (...args) { warnings.push(args.map(String).join(' ')) }

// ---------- 极简 React ----------
let currentHooks = null
let cursor = 0
let currentStore = null
let currentId = ''
let effectQueue = []
const hookStore = new Map()
const pathStack = []

function depsEqual(a, b) {
  if (!a || !b) return a === b
  if (a.length !== b.length) return false
  for (let i = 0; i < a.length; i += 1) if (a[i] !== b[i]) return false
  return true
}

function createElement(type, props, ...children) {
  const merged = Object.assign({}, props || {})
  if (children.length === 1) merged.children = children[0]
  else if (children.length > 1) merged.children = children
  return { type, props: merged }
}

class Component {
  constructor(props) { this.props = props || {}; this.state = {} }
  setState(update) {
    const previous = this.state || {}
    this.state = Object.assign({}, previous, typeof update === 'function' ? update(previous) : update)
  }
}

const React = {
  createElement,
  Component,
  Fragment: '__Fragment__',
  useState(initial) {
    const index = cursor++
    if (!(index in currentHooks)) {
      const cell = { value: typeof initial === 'function' ? initial() : initial }
      cell.set = function (update) { cell.value = typeof update === 'function' ? update(cell.value) : update }
      currentHooks[index] = cell
    }
    return [currentHooks[index].value, currentHooks[index].set]
  },
  useEffect(fn, deps) {
    const index = cursor++
    const key = currentId + '|e' + index
    const previous = currentStore.effects.get(key)
    const changed = !previous || !depsEqual(previous.deps, deps || null)
    currentStore.effects.set(key, { deps: deps || null, fn })
    if (changed) effectQueue.push({ key, fn })
  },
  useRef(initial) {
    const index = cursor++
    if (!(index in currentHooks)) currentHooks[index] = { current: initial }
    return currentHooks[index]
  },
  useMemo(fn) {
    const index = cursor++
    if (!(index in currentHooks)) currentHooks[index] = { value: fn() }
    return currentHooks[index].value
  },
  useCallback(fn) { return fn },
}

// ---------- 载入浏览器模块 ----------
const loaded = []
globalThis.window = { __ModuleLoader__: { load(spec) { loaded.push(spec) } } }
await import(new URL('../client.js', import.meta.url).href)
assert.equal(loaded.length, 1, '__ModuleLoader__.load 应被调用一次')
assert.equal(loaded[0].id, 'dsh-geobulk')
const mod = loaded[0].factory(function (name) {
  assert.equal(name, 'react', 'factory 只允许 require react')
  return React
})
assert.equal(mod.name, 'dsh-geobulk')
assert.deepEqual(mod.inject, ['slots', 'locale'])
assert.equal(typeof mod.apply, 'function')
const components = mod.__test__.components
assert.ok(components.WorkbenchPage && components.AnalysisDrawer)

// ---------- apply()：locale + 两个槽 ----------
const registered = []
const ctx = {
  effect(fn) { const dispose = fn(); return typeof dispose === 'function' ? dispose : function () {} },
  locale: {
    register(id, dicts) {
      assert.equal(id, 'dsh-geobulk')
      assert.ok(dicts.zh && dicts.en)
      this._dicts = dicts
      return function () {}
    },
    bind() { const dicts = this._dicts; return function (key) { return dicts.zh[key] } },
  },
  slots: {
    inject(owner, fn) { fn(); return function () {} },
    register(spec, component) {
      registered.push({ spec, component })
      return function () {}
    },
  },
}
mod.apply(ctx)
assert.equal(registered.length, 2)
assert.equal(registered[0].spec.name, 'main')
assert.equal(registered[0].spec.key, 'geobulk')
assert.equal(registered[1].spec.name, 'sidebar.panellist')
assert.equal(registered[1].spec.id, 'geobulk')
assert.equal(registered[1].spec.order, 4)
assert.equal(registered[1].spec.label(), 'Bulk 工作台')

// ---------- fetch 桩 ----------
const fetchCalls = []
globalThis.fetch = function (url, init) {
  fetchCalls.push(String(url))
  const path = String(url)
  const respond = (payload) => Promise.resolve({
    ok: true,
    status: 200,
    text: () => Promise.resolve(JSON.stringify(payload)),
  })
  if (path.endsWith('/health')) return respond({ ok: true, engine: { ready: true, python: '3.11.9', r: { version: '4.3.1' }, packages: {} }, geo: true, workdir: '/tmp' })
  if (path.endsWith('/state')) {
    return respond({
      ok: true,
      active: 'ds1',
      datasets: [{
        id: 'ds1', name: 'demo-bulk', n_samples: 20, n_genes: 2000, steps: ['qc', 'normalize', 'de'],
        sample_columns: [{ key: 'group', kind: 'categorical', n_unique: 2 }, { key: 'batch', kind: 'categorical', n_unique: 2 }],
        modules: [{ kind: 'circ', n: 12 }], history: [{ at: '2026-01-01T00:00:00Z', step: 'qc', ok: true, ms: 100, summary: 'ok' }],
        warnings: [], analyses_run: [],
      }],
      engine: { ready: true, python: '3.11.9', r: { version: '4.3.1' } },
      suggestions: [{ id: 's1', title: '运行差异表达', detail: '按分组做 limma', step: 'de' }],
      activity: [{ at: '2026-01-01T00:00:00Z', actor: 'agent', action: 'run:qc', detail: '完成', ok: true }],
      sc: { available: false, datasets: [], local: [] },
    })
  }
  if (path.endsWith('/formats')) return respond({ ok: true, formats: ['matrix_text', 'xlsx', 'geo_series'] })
  if (path.endsWith('/preview')) return respond({ ok: true, genes: ['CDH1', 'VIM', 'MKI67'], sample_columns: [{ key: 'group', kind: 'categorical', n_unique: 2 }], n_samples: 20, n_genes: 2000 })
  if (path.includes('/plot')) return respond({ ok: true, figure: { id: 'f1', url: '/geobulk/figure?id=f1.png', file: '/tmp/f1.png', title: 'plot' } })
  if (path.includes('/analyses')) return respond({ ok: true, groups: ['临床'], analyses: [{ id: 'immune', title: '免疫浸润', group: '临床', method: 'CIBERSORT', needs: ['normalize'], params: [{ key: 'methods', type: 'multiselect', default: ['epic'] }] }] })
  if (path.includes('/activity')) return respond({ entries: [{ at: '2026-01-01T00:00:00Z', actor: 'user', action: 'run:qc', detail: 'ok', ok: true }] })
  if (path.includes('/sc/refs')) return respond({ ok: true, available: false, datasets: [], local: [] })
  return respond({ ok: true })
}

// ---------- 渲染 ----------
const renderErrors = []
function renderNode(node, index) {
  if (node === null || node === undefined || node === false || node === true) return null
  if (typeof node === 'string' || typeof node === 'number') return node
  if (Array.isArray(node)) return node.map((child, childIndex) => renderNode(child, childIndex))
  if (typeof node !== 'object') return null
  const type = node.type
  const props = node.props || {}
  if (type === React.Fragment) return renderNode(props.children, index)
  if (typeof type === 'string') return { tag: type, children: renderNode(props.children, index) }
  if (typeof type !== 'function') return null
  if (type.prototype && typeof type.prototype.render === 'function') {
    const instance = new type(props)
    return renderNode(instance.render(), index)
  }
  pathStack.push((type.name || 'anon') + '#' + index)
  currentId = pathStack.join('/')
  let store = hookStore.get(currentId)
  if (!store) { store = { hooks: [], effects: new Map() }; hookStore.set(currentId, store) }
  const savedHooks = currentHooks
  const savedCursor = cursor
  const savedStore = currentStore
  currentHooks = store.hooks
  cursor = 0
  currentStore = store
  let output
  try { output = type(props) } catch (error) {
    renderErrors.push({ id: currentId, error })
    throw error
  } finally {
    currentHooks = savedHooks
    cursor = savedCursor
    currentStore = savedStore
    pathStack.pop()
  }
  return renderNode(output, index)
}

async function renderPasses(label, element, passes) {
  let tree = null
  for (let pass = 0; pass < passes; pass += 1) {
    effectQueue = []
    pathStack.length = 0
    try { tree = renderNode(createElement(element.type, element.props), 0) } catch (error) {
      renderErrors.push({ id: label, error })
      break
    }
    for (const item of effectQueue) {
      try { item.fn() } catch (error) { renderErrors.push({ id: label + ':' + item.key, error }) }
    }
    await new Promise((resolve) => setTimeout(resolve, 25))
  }
  return tree
}

const h = createElement
const ok = await renderPasses('WorkbenchPage', h(components.WorkbenchPage), 5)
assert.ok(ok, 'WorkbenchPage 应渲染出树')
assert.equal(renderErrors.length, 0, '渲染/effect 不应抛错：' + renderErrors.map((item) => item.id + ' → ' + (item.error && item.error.stack ? item.error.stack : item.error)).join(' | '))

// 主界面关键结构出现
function collectTags(node, tags) {
  if (!node) return tags
  if (typeof node === 'string') return tags
  if (Array.isArray(node)) { for (const child of node) collectTags(child, tags); return tags }
  if (node.tag) tags.add(node.tag)
  collectTags(node.children, tags)
  return tags
}
const tags = collectTags(ok, new Set())
assert.ok(tags.has('header'), '应有顶栏 header')
assert.ok(tags.has('aside'), '应有左右 aside')
assert.ok(tags.has('nav'), '应有 tab nav')
assert.ok(tags.has('div'), '应有结构化内容')

// 各 tab / 弹层 / 抽屉组件单独渲染
const statePayload = JSON.parse(await (await fetch('/geobulk/state')).text())
const dataset = mod.__test__.normalizeState(statePayload).datasets[0]
const base = {
  dataset,
  preview: mod.__test__.normalizePreview({ genes: ['CDH1'], sample_columns: [{ key: 'group', kind: 'categorical', n_unique: 2 }, { key: 'os_time', kind: 'time' }, { key: 'os_event', kind: 'event' }] }),
  health: mod.__test__.normalizeHealth({ engine: { ready: true, python: '3.11.9', r: { version: '4.3.1' } } }),
  suggestions: [], activity: [], plots: {}, exports: {}, stepResults: {},
  loadPlot() {}, runExport() {}, callAnalysis() {}, runStep() {}, onRunStep() {}, onOpenImport() {},
  stepBusy: '', api: mod.__test__.apiClient(), refreshState() {}, copiedKey: '', onCopy() {},
  onOpenAnalysis() {}, onClose() {}, onImported() {}, formats: { ok: true, formats: ['matrix_text'] },
  columns: [],
}
const deResult = mod.__test__.normalizeResult({ summary: { title: '差异表达', bullets: ['检出 120 个上调基因'], metrics: { up: 120, down: 80 } }, tables: [{ id: 'de', columns: [{ key: 'gene' }, { key: 'logFC', type: 'number' }], rows: [{ gene: 'CDH1', logFC: 2.1 }] }], figures: [{ id: 'f1', title: '火山图', file: '/tmp/f1.png' }] })
const stepResults = { de: deResult, enrich: deResult, immune: deResult, survival: deResult, prognosis: deResult, scmap: deResult, ncrna: deResult }
const cases = [
  ['OverviewTab', components.OverviewTab, base],
  ['ExprTab', components.ExprTab, base],
  ['DeTab', components.DeTab, Object.assign({}, base, { stepResults })],
  ['EnrichTab', components.EnrichTab, Object.assign({}, base, { stepResults })],
  ['ImmuneTab', components.ImmuneTab, Object.assign({}, base, { stepResults })],
  ['SurvivalTab', components.SurvivalTab, Object.assign({}, base, { stepResults })],
  ['PrognosisTab', components.PrognosisTab, Object.assign({}, base, { stepResults })],
  ['ScTab', components.ScTab, Object.assign({}, base, { stepResults })],
  ['NcrnaTab', components.NcrnaTab, Object.assign({}, base, { stepResults })],
  ['RecordsTab', components.RecordsTab, base],
  ['LogTab', components.LogTab, base],
  ['AgentPanel', components.AgentPanel, base],
  ['PipelinePanel', components.PipelinePanel, Object.assign({}, base, { statuses: { qc: 'done', normalize: 'done', de: 'running' }, values: { de: { group_col: 'group' } }, summaries: { de: '检出 120 个上调基因' }, errors: {} })],
  ['GeoDialog', components.GeoDialog, base],
  ['ImportDialog', components.ImportDialog, base],
  ['AnalysisDrawer', components.AnalysisDrawer, base],
  ['EmptyOverview', components.OverviewTab, Object.assign({}, base, { dataset: null })],
  ['EmptyDe', components.DeTab, Object.assign({}, base, { dataset: null })],
]
for (const [label, component, props] of cases) {
  await renderPasses(label, h(component, props), 2)
}
assert.equal(renderErrors.length, 0, '各 tab 渲染不应抛错：' + renderErrors.map((item) => item.id + ' → ' + (item.error && item.error.stack ? item.error.stack : item.error)).join(' | '))

// 首屏应请求过 health / state / formats / preview
for (const expected of ['/geobulk/health', '/geobulk/state', '/geobulk/formats']) {
  assert.ok(fetchCalls.some((url) => url.endsWith(expected)), '应请求 ' + expected)
}
assert.ok(fetchCalls.some((url) => url.endsWith('/preview')), '切数据集后应请求 preview')
assert.ok(fetchCalls.some((url) => url.includes('/plot?') || url.endsWith('/plot')), '概览应出图')

// 不应有插件自身警告
const pluginWarnings = warnings.filter((line) => line.includes('[dsh-geobulk]'))
assert.deepEqual(pluginWarnings, [], '不应有插件警告：' + pluginWarnings.join(' | '))

console.warn = originalWarn
console.log('SMOKE OK · fetch calls=' + fetchCalls.length + ' · renders=' + (cases.length + 1) + ' · components=' + Object.keys(components).length)
