/**
 * dsh-geobulk — Host 侧联调测试（不需要宿主进程）
 *
 * 用假 ctx / 假 req / 假 res 直接调 lib/impl.mjs 的 handle()，验证：
 *   ping/health、state、demo、preview、plot、suggest、audit、export、jobs、图路由、工具 schema、prompt 段落
 * 运行：node test/host_ops.test.mjs
 */

import { mkdtempSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { Writable } from 'node:stream'

const HERE = path.dirname(new URL(import.meta.url).pathname)
const impl = await import(new URL('../lib/impl.mjs', import.meta.url).href)

let passed = 0
let failed = 0
const failures = []

function check(name, condition, detail = '') {
  if (condition) {
    passed += 1
  } else {
    failed += 1
    failures.push(`${name}${detail ? `：${detail}` : ''}`)
    console.log(`  ✗ ${name}${detail ? `：${detail}` : ''}`)
  }
}

function fakeRes() {
  const chunks = []
  const res = new Writable({
    write(chunk, _enc, callback) {
      chunks.push(Buffer.from(chunk))
      callback()
    },
  })
  res.status = 0
  res.headers = {}
  res.headersSent = false
  res.writeHead = (status, headers) => {
    res.status = status
    res.headers = headers || {}
    res.headersSent = true
  }
  res.getBody = () => (chunks.length ? Buffer.concat(chunks) : null)
  return res
}

function fakeReq(method, url, body) {
  const payload = body === undefined ? [] : [Buffer.from(JSON.stringify(body))]
  return {
    method,
    url,
    headers: { host: '127.0.0.1:49180', 'content-type': 'application/json' },
    on(event, handler) {
      if (event === 'data') for (const chunk of payload) handler(chunk)
      if (event === 'end') handler()
      if (event === 'error') { /* noop */ }
      return this
    },
    destroy() {},
  }
}

async function call(runtime, method, route, body) {
  const res = fakeRes()
  await impl.handle(runtime, fakeReq(method, `/geobulk${route}`, body), res)
  const payload = res.getBody()
  if (!payload) return { status: res.status, json: null, raw: null }
  const text = payload.toString('utf8')
  if ((res.headers['content-type'] || '').includes('image/png')) {
    return { status: res.status, buffer: payload }
  }
  let json = null
  try { json = JSON.parse(text) } catch { json = null }
  return { status: res.status, json, raw: text }
}

const workdir = mkdtempSync(path.join(tmpdir(), 'dsh-geobulk-test-'))
const runtime = { ctx: { get: () => undefined }, config: { workdir }, activity: [], jobs: new Map() }

console.log('== host ops ==')

const health = await call(runtime, 'GET', '/health')
check('health 200', health.status === 200, `status=${health.status}`)
check('health.engine 字段', health.json && health.json.engine && typeof health.json.engine.ready === 'boolean')
check('health 有 R 信息', health.json && health.json.engine && health.json.engine.r !== undefined)

const demo = await call(runtime, 'POST', '/dataset/demo', {})
check('demo 200', demo.status === 200 && demo.json && demo.json.ok === true, JSON.stringify(demo.json && demo.json.error))
const datasetId = demo.json && demo.json.dataset ? demo.json.dataset.id : null
check('demo 数据集有样本/基因', demo.json && demo.json.dataset && demo.json.dataset.n_samples > 10 && demo.json.dataset.n_genes > 100,
  demo.json && demo.json.dataset ? `${demo.json.dataset.n_samples}x${demo.json.dataset.n_genes}` : '')
check('demo 有分组', demo.json && demo.json.dataset && Array.isArray(demo.json.dataset.groups) && demo.json.dataset.groups.length >= 2)

const state = await call(runtime, 'GET', '/state')
check('state 200', state.status === 200 && state.json && state.json.ok === true)
check('state 有 datasets', state.json && Array.isArray(state.json.datasets) && state.json.datasets.length === 1)
check('state 有 recipes', state.json && Array.isArray(state.json.recipes) && state.json.recipes.length >= 3)
check('state 有 suggestions', state.json && Array.isArray(state.json.suggestions))

const preview = await call(runtime, 'POST', '/preview', { dataset: datasetId })
check('preview 200', preview.status === 200 && preview.json && preview.json.ok === true)
check('preview 基因列表', preview.json && Array.isArray(preview.json.genes) && preview.json.genes.length > 0)
check('preview 样本列含 group', preview.json && (preview.json.sample_columns || []).some((c) => c.key === 'group'))

const plot = await call(runtime, 'POST', '/plot', { dataset: datasetId, kind: 'pca_scatter', params: {} })
check('plot 200', plot.status === 200 && plot.json && plot.json.ok === true, plot.json && plot.json.error)
check('plot 返回 PNG 绝对路径', plot.json && typeof plot.json.file === 'string' && plot.json.file.endsWith('.png'))
const figureId = plot.json && plot.json.figure ? path.basename(plot.json.file) : null
const figureRes = await call(runtime, 'GET', `/figure?id=${encodeURIComponent(figureId || '')}`)
check('figure 路由回 PNG', Boolean(figureRes.buffer && figureRes.buffer.slice(0, 8).equals(Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]))))

const badPlot = await call(runtime, 'POST', '/plot', { dataset: datasetId, kind: 'nope', params: {} })
check('未知图类型报错', badPlot.json && badPlot.json.ok === false)

const suggest = await call(runtime, 'POST', '/suggest', { dataset: datasetId })
check('suggest 200', suggest.status === 200 && suggest.json && Array.isArray(suggest.json.suggestions))

const audit = await call(runtime, 'POST', '/audit', { dataset: datasetId })
check('audit 200', audit.status === 200 && audit.json && typeof audit.json.score === 'number' && Array.isArray(audit.json.checks))

const exportRes = await call(runtime, 'POST', '/export', { dataset: datasetId, what: 'report' })
check('export report 200', exportRes.status === 200 && exportRes.json && exportRes.json.ok === true, exportRes.json && exportRes.json.error)
check('export report 文件存在', exportRes.json && typeof exportRes.json.file === 'string')

const formats = await call(runtime, 'GET', '/formats')
check('formats 200', formats.status === 200 && formats.json && Array.isArray(formats.json.formats) && formats.json.formats.length > 3)

const analyses = await call(runtime, 'GET', `/analyses?dataset=${datasetId}`)
check('analyses 200', analyses.status === 200 && analyses.json && Array.isArray(analyses.json.analyses))

const jobs = await call(runtime, 'GET', '/jobs')
check('jobs 200', jobs.status === 200 && jobs.json && Array.isArray(jobs.json.jobs))

const scRefs = await call(runtime, 'GET', '/sc/refs')
check('sc/refs 200', scRefs.status === 200 && scRefs.json && scRefs.json.ok === true)
check('sc/refs 无 singlecell 时不报错', scRefs.json && typeof scRefs.json.available === 'boolean')

const loopback = await call(runtime, 'GET', '/health')
check('回环 Host 放行', loopback.status === 200)
const foreign = fakeRes()
await impl.handle(runtime, { ...fakeReq('GET', '/geobulk/health'), headers: { host: '10.0.0.5:1234' } }, foreign)
check('非回环 Host 拒绝', foreign.status === 403)

// 工具 schema
const toolDef = await import(new URL('../lib/impl.mjs', import.meta.url).href)
check('impl 导出 STEPS 12 个', Array.isArray(toolDef.STEPS) && toolDef.STEPS.length === 12, String(toolDef.STEPS.length))
check('STEP_IDS 顺序正确', toolDef.STEP_IDS[0] === 'qc' && toolDef.STEP_IDS[11] === 'ncrna')
check('PLOT_KINDS 有 20 种', toolDef.PLOT_KINDS.length === 20, String(toolDef.PLOT_KINDS.length))

// runTool（未实现的步骤应给出可读错误）
let runToolError = null
try {
  await impl.runTool(runtime, { action: 'status' })
} catch (error) {
  runToolError = error
}
check('runTool status 可用', runToolError === null, runToolError && runToolError.message)

let badActionThrew = false
try { await impl.runTool(runtime, { action: 'nope' }) } catch { badActionThrew = true }
check('runTool 未知 action 抛错', badActionThrew)

// 提示词段落
const promptText = impl.renderPromptSection(runtime)
check('renderPromptSection 不抛错', typeof promptText === 'string')

console.log(`\n== 结果：${passed} 通过 / ${failed} 失败 ==`)
if (failed) {
  console.log(failures.map((f) => ` - ${f}`).join('\n'))
}
try { rmSync(workdir, { recursive: true, force: true }) } catch { /* ignore */ }
process.exit(failed ? 1 : 0)
