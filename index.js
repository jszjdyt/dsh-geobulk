/**
 * dsh-geobulk — Host 半边入口（薄且稳定）
 *
 * 只做四件事：
 *   1. 把 `/geobulk/*` 的 exact 路由挂到 webServer，全部转给 lib/impl.mjs；
 *   2. 注册 `geobulk` agent 工具（和 UI 同一套操作实现）；
 *   3. 注册一段 system prompt section，让 Agent 知道工作台与当前数据集；
 *   4. 插件卸载时关掉 Python worker。
 *
 * 业务逻辑都在 lib/impl.mjs / lib/engine.mjs 里，按 `?rev=` 动态 import，
 * 这样改完 lib 代码不用重启宿主进程（ESM 缓存绕开的做法见 dsh-paper-reader）。
 */

const PREFIX = '/geobulk'
const VERSION = '0.1.0'

/** 每次调用都拿一份新的 impl（带 rev 查询参数绕开 ESM 缓存）。 */
async function loadImpl() {
  const url = new URL(`./lib/impl.mjs?rev=${Date.now()}-${Math.random().toString(36).slice(2)}`, import.meta.url)
  return import(url.href)
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
