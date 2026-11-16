/**
 * dsh-geobulk — Host 半边入口（薄且稳定）
 *
 * 只做四件事：
 *   1. 把 `/geobulk/*` 的 exact 路由挂到 webServer，全部转给 lib/impl.mjs；
 *   2. 注册 `geobulk` agent 工具（和 UI 同一套操作实现）；
 *   3. 注册一段 system prompt section，让 Agent 知道工作台与当前数据集；
 *   4. 插件卸载时关掉 Python worker。
 *
 * 业务逻辑都在 lib/impl.mjs / lib/engine.mjs 里。为了改完 lib 代码不用重启宿主进程，
 * impl 用带 `?rev=<mtime>` 的 URL 动态 import；**同一个文件版本只 import 一次**，
 * 否则每个 HTTP 请求都会新建一份模块实例，内存会随轮询持续增长（HOST-04）。
 */

import { statSync } from 'node:fs'

const PREFIX = '/geobulk'
const VERSION = '0.1.0'

/** 每次调用都拿一份新的 impl（带 rev 查询参数绕开 ESM 缓存）。 */
const IMPL_URL = new URL('./lib/impl.mjs', import.meta.url)

/** 当前 impl.mjs 的版本标识：文件没改就不换模块实例。 */
function implVersion() {
  try { return String(statSync(IMPL_URL).mtimeMs) } catch { return 'static' }
}

let implCache = { version: '', promise: null }

/**
 * 取 impl 模块：**按文件版本缓存**。
 * 以前每次请求都生成 `?rev=Date.now()-随机`，等于每次请求都导入一份完整模块，
 * 长期打开会持续占用内存（HOST-04：500 次导入即增长约 95 MiB）。
 * 现在同一版本复用同一个 import Promise；开发时改了 impl.mjs，mtime 变化会自动换新模块。
 */
function loadImpl() {
  const version = implVersion()
  if (implCache.promise && implCache.version === version) return implCache.promise
  const url = new URL(`?rev=${encodeURIComponent(version)}`, IMPL_URL)
  const promise = import(url.href)
  implCache = { version, promise }
  promise.catch(() => {
    // 导入失败不要把坏 Promise 缓存住
    if (implCache.promise === promise) implCache = { version: '', promise: null }
  })
  return promise
}

export const name = 'dsh-geobulk'

export function apply(ctx, config) {
  const runtime = {
    version: VERSION,
    config: config && typeof config === 'object' ? config : {},
    ctx,
    invocation: 'apply',
  }

  // ---- HTTP 路由：页面通过同源 fetch 调这里 ----
  // 只注册一条前缀路由（'/geobulk'），具体分发在 lib/impl.mjs 的 handle() 里做。
  // 这样以后新增 /geobulk/* 接口不需要重新加载 Host 插件（index.js 只在激活时执行一次）。
  ctx.inject(['webServer'], (scoped) => {
    const webServer = scoped.webServer
    const dispose = webServer.register({
      kind: 'prefix',
      path: PREFIX,
      handler: (req, res) => dispatch(runtime, req, res),
    })
    scoped.effect(() => () => {
      try { dispose() } catch { /* ignore */ }
    }, 'dsh-geobulk: routes')
  })

  // ---- agent 工具 ----
  ctx.inject(['tools'], () => {
    loadImpl()
      .then((impl) => impl.ensureToolRegistered(runtime))
      .catch((error) => {
        try { console.warn('[dsh-geobulk] agent tool unavailable:', error && error.message ? error.message : error) } catch { /* ignore */ }
      })
  })

  // ---- system prompt section ----
  ctx.inject(['systemPrompt'], () => {
    loadImpl()
      .then((impl) => impl.ensurePromptSection(runtime))
      .catch((error) => {
        try { console.warn('[dsh-geobulk] prompt section unavailable:', error && error.message ? error.message : error) } catch { /* ignore */ }
      })
  })

  // ---- 卸载时回收 Python worker ----
  ctx.effect(() => () => {
    const engine = runtime.engine
    if (engine && typeof engine.dispose === 'function') {
      try { engine.dispose() } catch { /* ignore */ }
    }
  }, 'dsh-geobulk: engine teardown')
}

async function dispatch(runtime, req, res) {
  try {
    const impl = await loadImpl()
    await impl.handle(runtime, req, res)
  } catch (error) {
    const message = error && error.stack ? String(error.stack) : String(error)
    try {
      if (!res.headersSent) {
        const body = JSON.stringify({ ok: false, code: 'internal', error: message })
        res.writeHead(500, { 'content-type': 'application/json; charset=utf-8', 'cache-control': 'no-store' })
        res.end(body)
      } else {
        res.end()
      }
    } catch { /* ignore */ }
    try { console.error('[dsh-geobulk] dispatch failed:', error) } catch { /* ignore */ }
  }
}
