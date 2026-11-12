/**
 * dsh-geobulk — Python worker 管理器（Host 侧）
 *
 * 职责：
 *   1. 探测可用的 Python 解释器（宿主 Electron 进程的 PATH 往往很窄，不能只靠 `python3`）；
 *   2. 以 `python -u bg_engine.py --workdir <dir>` 拉起长驻 worker，JSON-lines over stdio；
 *   3. 请求排队、超时、进度缓存、stderr 日志、崩溃/超时后自动重启；
 *   4. 向 UI 暴露 status()：解释器、引擎模式（builtin/scanpy）、包版本、错误原因。
 *
 * 本文件只依赖 node 内置模块，可以被 `test/engine.test.mjs` 直接 import 测试。
 */

import { spawn } from 'node:child_process'
import { createWriteStream, mkdirSync, readdirSync, realpathSync, unlinkSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'
import readline from 'node:readline'

const ENGINE_SCRIPT = fileURLToPath(new URL('./bg_engine.py', import.meta.url))
const DEFAULT_TIMEOUT = 240_000
const PROBE_TIMEOUT = 12_000
const STDERR_TAIL = 16_000
const MAX_QUEUE = 64

/** 探测用的 import 检查：这些包是 builtin 引擎的硬需求。 */
const PROBE_CODE =
  'import sys,numpy,pandas,scipy,sklearn,matplotlib,anndata,lifelines;' +
  'print("PROBE_OK", sys.version.split()[0])'

function tail(text, n = 800) {
  const s = String(text || '')
  return s.length <= n ? s : '…' + s.slice(-n)
}

function uniq(list) {
  const seen = new Set()
  const out = []
  for (const item of list) {
    if (!item) continue
    let key = item
    try { key = realpathSync(item) } catch { /* 不存在的候选保留原值 */ }
    if (seen.has(key)) continue
    seen.add(key)
    out.push(item)
  }
  return out
}

/**
 * 候选解释器列表。顺序即优先级：
 *   显式配置 / 环境变量 > PATH 上的 python3 > macOS 框架版 > Homebrew > 系统版。
 */
export function pythonCandidates(opts = {}) {
  const list = []
  if (opts.python) list.push(opts.python)
  if (process.env.DSH_GEOBULK_PYTHON) list.push(process.env.DSH_GEOBULK_PYTHON)
  if (process.env.DSH_BG_PYTHON) list.push(process.env.DSH_BG_PYTHON)
  list.push('python3')
  if (process.platform === 'win32') {
    // Windows 上通常没有 python3 这个别名，用 python / py launcher
    list.push('python', 'py')
  }

  // macOS 官方框架版：把已装版本按版本号从高到低排进去
  const frameworkRoot = '/Library/Frameworks/Python.framework/Versions'
  try {
    const versions = readdirSync(frameworkRoot)
      .filter((v) => /^\d+\.\d+$/.test(v))
      .sort((a, b) => {
        const [a1, a2] = a.split('.').map(Number)
        const [b1, b2] = b.split('.').map(Number)
        return b1 - a1 || b2 - a2
      })
    for (const v of versions) list.push(path.join(frameworkRoot, v, 'bin', 'python3'))
  } catch { /* 非 macOS 或没有框架版 */ }

  list.push('/opt/homebrew/bin/python3', '/usr/local/bin/python3', '/usr/bin/python3')
  return uniq(list)
}

/** 跑一次 `python -c PROBE_CODE`，判断解释器是否带齐依赖。 */
export function probePython(cmd, timeout = PROBE_TIMEOUT) {
  return new Promise((resolve) => {
    let settled = false
    let out = ''
    let err = ''
    let child
    const finish = (result) => {
      if (settled) return
      settled = true
      clearTimeout(timer)
      try { child && child.kill('SIGKILL') } catch { /* ignore */ }
      resolve(result)
    }
    const timer = setTimeout(() => finish({ ok: false, reason: '探测超时' }), timeout)
    try {
      child = spawn(cmd, ['-u', '-c', PROBE_CODE], { stdio: ['ignore', 'pipe', 'pipe'] })
    } catch (error) {
      return finish({ ok: false, reason: String(error && error.message ? error.message : error) })
    }
    child.stdout.on('data', (b) => { out += b.toString('utf8') })
    child.stderr.on('data', (b) => { err += b.toString('utf8') })
    child.on('error', (error) => finish({ ok: false, reason: String(error && error.message ? error.message : error) }))
    child.on('close', (code) => {
      if (code === 0 && /PROBE_OK\s+(\S+)/.test(out)) {
        const version = out.match(/PROBE_OK\s+(\S+)/)[1]
        finish({ ok: true, python: cmd, version })
      } else {
        const detail = tail(err || out, 300).replace(/\s+/g, ' ').trim()
        finish({ ok: false, reason: detail || `退出码 ${code}` })
      }
    })
  })
}

/** 依次探测候选项，返回第一个可用的；失败时带上每个候选的原因。 */
export async function findPython(opts = {}) {
  const candidates = pythonCandidates(opts)
  const failures = []
  for (const cmd of candidates) {
    const result = await probePython(cmd, opts.probeTimeout ?? PROBE_TIMEOUT)
    if (result.ok) return { ok: true, python: cmd, version: result.version, candidates }
    failures.push(`${cmd}: ${result.reason}`)
  }
  return { ok: false, candidates, failures }
}

/** 确保工作目录骨架存在，并清掉上一轮残留的图。 */
export function prepareWorkdir(root) {
  const dirs = {
    root,
    figures: path.join(root, 'figures'),
    uploads: path.join(root, 'uploads'),
    exports: path.join(root, 'exports'),
    mpl: path.join(root, 'mpl-cache'),
    log: path.join(root, 'engine.log'),
  }
  for (const key of ['root', 'figures', 'uploads', 'exports', 'mpl']) mkdirSync(dirs[key], { recursive: true })
  // 旧图会在重启后失效（figure 表在内存里），清掉避免磁盘无限增长
  try {
    for (const name of readdirSync(dirs.figures)) {
      if (!name.endsWith('.png')) continue
      try { unlinkSync(path.join(dirs.figures, name)) } catch { /* ignore */ }
    }
  } catch { /* ignore */ }
  return dirs
}

export function defaultWorkdir() {
  return path.join(tmpdir(), 'dsh-geobulk')
}

/**
 * 创建一个引擎实例。惰性启动：第一次 request 时才探测 + spawn。
 *
 * @param {{python?:string, workdir?:string, engineScript?:string, timeout?:number, onLog?:(line:string)=>void}} opts
 */
export function createEngine(opts = {}) {
  const workdir = opts.workdir || defaultWorkdir()
  const dirs = prepareWorkdir(workdir)
  const engineScript = opts.engineScript || ENGINE_SCRIPT
  const defaultTimeout = opts.timeout || DEFAULT_TIMEOUT

  const state = {
    child: null,
    rl: null,
    starting: null,       // Promise<{python,version,mode}>
    pythonInfo: null,     // {python, version, mode, packages, layouts, warnings}
    stderrTail: '',
    pending: new Map(),   // id -> {resolve, reject, timer, cmd, progress}
    seq: 0,
    disposed: false,
    logStream: null,
    lastProgress: null,
    lastError: null,
    startedAt: null,
  }

  const logStream = () => {
    if (state.logStream) return state.logStream
    try {
      state.logStream = createWriteStream(dirs.log, { flags: 'a' })
      state.logStream.on('error', () => { state.logStream = null })
    } catch { state.logStream = null }
    return state.logStream
  }
  const log = (scope, message) => {
    const line = `[${new Date().toISOString()}] [${scope}] ${message}`
    try { logStream() && state.logStream.write(line + '\n') } catch { /* ignore */ }
    if (typeof opts.onLog === 'function') {
      try { opts.onLog(line) } catch { /* ignore */ }
    }
  }

  const rejectAll = (error) => {
    for (const [, entry] of state.pending) {
      clearTimeout(entry.timer)
      try { entry.reject(error) } catch { /* ignore */ }
    }
    state.pending.clear()
  }

  const killChild = (reason) => {
    const child = state.child
    state.child = null
    if (state.rl) { try { state.rl.close() } catch { /* ignore */ } state.rl = null }
    if (!child) return
    try { child.stdin.end() } catch { /* ignore */ }
    try { child.kill('SIGTERM') } catch { /* ignore */ }
    setTimeout(() => { try { child.kill('SIGKILL') } catch { /* ignore */ } }, 800).unref?.()
    if (reason) log('engine', `worker 关闭：${reason}`)
  }

  const spawnWorker = (python) => {
    log('engine', `spawn ${python} -u ${engineScript} --workdir ${workdir}`)
    let configJson = '{}'
    try { configJson = JSON.stringify(opts.config && typeof opts.config === 'object' ? opts.config : {}) } catch { configJson = '{}' }
    const child = spawn(python, ['-u', engineScript, '--workdir', workdir], {
      cwd: path.dirname(engineScript),
      env: {
        ...process.env,
        MPLBACKEND: 'Agg',
        MPLCONFIGDIR: dirs.mpl,
        PYTHONUNBUFFERED: '1',
        PYTHONIOENCODING: 'utf-8',
        DS_BG_WORKDIR: workdir,
        DS_BG_CONFIG: configJson,
      },
      stdio: ['pipe', 'pipe', 'pipe'],
    })
    state.child = child
    state.lastError = null

    state.rl = readline.createInterface({ input: child.stdout })
    state.rl.on('line', (line) => {
      const text = String(line || '').trim()
      if (!text) return
      let msg
      try { msg = JSON.parse(text) } catch {
        log('stdout', `非 JSON 行：${tail(text, 300)}`)
        return
      }
      if (msg && msg.event === 'progress') {
        const entry = state.pending.get(String(msg.id))
        if (entry) entry.progress = { message: msg.message || '', pct: msg.pct }
        state.lastProgress = { at: Date.now(), message: msg.message || '', pct: msg.pct }
        return
      }
      const id = String(msg && msg.id)
      const entry = state.pending.get(id)
      if (!entry) return
      state.pending.delete(id)
      clearTimeout(entry.timer)
      if (msg.ok) entry.resolve(msg.result)
      else {
        const error = new Error((msg.error && msg.error.message) || '引擎返回失败')
        error.code = (msg.error && msg.error.code) || 'internal'
        if (msg.error && msg.error.trace) error.trace = msg.error.trace
        if (msg.error && msg.error.details) error.details = msg.error.details
        entry.reject(error)
      }
    })

    child.stderr.on('data', (buf) => {
      const text = buf.toString('utf8')
      state.stderrTail = (state.stderrTail + text).slice(-STDERR_TAIL)
      log('stderr', tail(text, 2000).replace(/\s+$/, ''))
    })

    child.on('error', (error) => {
      state.lastError = `worker 启动失败：${error && error.message ? error.message : error}`
      rejectAll(Object.assign(new Error(state.lastError), { code: 'engine_exited' }))
    })

    child.on('exit', (code, signal) => {
      state.child = null
      if (state.disposed) return
      const detail = `worker 退出 code=${code} signal=${signal} :: ${tail(state.stderrTail, 500).replace(/\s+/g, ' ')}`
      log('engine', detail)
      const error = Object.assign(new Error(`单细胞引擎异常退出（${code ?? signal}）：${tail(state.stderrTail, 400).trim()}`), {
        code: 'engine_exited',
      })
      rejectAll(error)
      // 引擎崩溃意味着内存里的数据集没了，下一次请求会重新启动
      state.pythonInfo = null
      state.starting = null
    })
    return child
  }

  /** 启动（幂等）：探测 Python → spawn → ping 握手。 */
  const start = async () => {
    if (state.disposed) throw Object.assign(new Error('插件已卸载'), { code: 'disposed' })
    if (state.starting) return state.starting
    state.starting = (async () => {
      const found = await findPython({ python: opts.python })
      if (!found.ok) {
        const detail = found.failures.slice(0, 4).join('；')
        throw Object.assign(
          new Error(
            `找不到可用的 Python 解释器（需要 numpy/pandas/scipy/sklearn/anndata/matplotlib）。` +
              `已尝试：${detail}`,
          ),
          { code: 'no_python' },
        )
      }
      log('engine', `使用 Python：${found.python} (${found.version})`)
      spawnWorker(found.python)
      const info = await request('ping', {}, { timeout: 30_000, skipStart: true })
      state.pythonInfo = {
        python: found.python,
        version: found.version,
        mode: info.mode || 'builtin',
        packages: info.packages || {},
        layouts: info.layouts || ['tsne', 'pca'],
        warnings: info.warnings || [],
        platform: info.platform || process.platform,
      }
      state.startedAt = Date.now()
      log('engine', `握手成功：mode=${state.pythonInfo.mode} python=${found.version}`)
      return state.pythonInfo
    })().catch((error) => {
      state.starting = null
      state.lastError = String(error && error.message ? error.message : error)
      throw error
    })
    return state.starting
  }

  /** 发一条命令；排队 + 超时 + 崩溃时自动重试一次（仅幂等命令）。 */
  const request = (cmd, args = {}, reqOpts = {}) => {
    if (state.disposed) return Promise.reject(Object.assign(new Error('插件已卸载'), { code: 'disposed' }))
    const timeout = reqOpts.timeout || defaultTimeout
    const run = async () => {
      if (!reqOpts.skipStart) await start()
      if (!state.child) throw Object.assign(new Error('引擎未运行'), { code: 'engine_exited' })
      if (state.pending.size >= MAX_QUEUE) {
        throw Object.assign(new Error('引擎请求过多，请稍后再试'), { code: 'busy' })
      }
      return new Promise((resolve, reject) => {
        const id = String(++state.seq)
        const timer = setTimeout(() => {
          state.pending.delete(id)
          const detail = tail(state.stderrTail, 400).replace(/\s+/g, ' ').trim()
          reject(Object.assign(
            new Error(`单细胞引擎超时（${cmd}，${Math.round(timeout / 1000)}s）${detail ? '：' + detail : ''}`),
            { code: 'timeout' },
          ))
          // 超时后 worker 可能卡住，杀掉让它下次重启
          killChild(`请求超时 ${cmd}`)
          rejectAll(Object.assign(new Error('引擎因超时重启，内存中的数据集需要重新载入'), { code: 'engine_exited' }))
        }, timeout)
        state.pending.set(id, { resolve, reject, timer, cmd, progress: null })
        try {
          state.child.stdin.write(JSON.stringify({ id, cmd, args }) + '\n')
        } catch (error) {
          clearTimeout(timer)
          state.pending.delete(id)
          reject(Object.assign(new Error(`写入引擎失败：${error && error.message ? error.message : error}`), { code: 'engine_exited' }))
        }
      })
    }
    const promise = state.queue.then(run, run)
    state.queue = promise.then(() => undefined, () => undefined)
    return promise
  }
  state.queue = Promise.resolve()

  const dispose = async () => {
    state.disposed = true
    if (state.child) {
      try {
        state.child.stdin.write(JSON.stringify({ id: 'shutdown', cmd: 'shutdown', args: {} }) + '\n')
      } catch { /* ignore */ }
      await new Promise((resolve) => setTimeout(resolve, 150))
    }
    killChild('插件卸载')
    rejectAll(Object.assign(new Error('插件已卸载'), { code: 'disposed' }))
    try { state.logStream && state.logStream.end() } catch { /* ignore */ }
    state.logStream = null
  }

  const status = () => ({
    ready: Boolean(state.child),
    mode: state.pythonInfo ? state.pythonInfo.mode : state.lastError ? 'none' : 'unknown',
    python: state.pythonInfo ? state.pythonInfo.python : null,
    version: state.pythonInfo ? state.pythonInfo.version : null,
    packages: state.pythonInfo ? state.pythonInfo.packages : null,
    layouts: state.pythonInfo ? state.pythonInfo.layouts : null,
    warnings: state.pythonInfo ? state.pythonInfo.warnings : [],
    error: state.lastError,
    startedAt: state.startedAt,
    lastProgress: state.lastProgress,
    stderrTail: tail(state.stderrTail, 1200),
  })

  return {
    workdir,
    dirs,
    request,
    start,
    dispose,
    status,
    get progress() { return state.lastProgress },
    get info() { return state.pythonInfo },
    get running() { return Boolean(state.child) },
    /** 引擎脚本是否存在（测试/诊断用） */
    scriptPath: engineScript,
    _internal: state,
  }
}
