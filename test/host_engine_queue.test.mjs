/**
 * dsh-geobulk / dsh-singlecell — worker 生命周期回归测试（任务单 HOST-02）
 *
 * 复现的问题：start() 里的 ping 握手和业务请求共用同一条串行队列，worker 崩溃后
 * 排队中的请求会永久等待（超时计时器在 start 之后才创建）。
 *
 * 这里用一个只讲 JSON-lines 的假 worker 驱动真实的 engine.mjs：
 *   1. 冷启动（第一次请求就是业务请求）必须能结束；
 *   2. 先 start()，再同时提交 crash + list，list 必须在明确时限内结束；
 *   3. crash 之后的后续请求仍能正常执行；
 *   4. 不留下悬空 Promise（测试结束时没有 pending 项）。
 *
 * 运行：node test/host_engine_queue.test.mjs
 */

import { mkdtempSync, writeFileSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { createEngine } from '../lib/engine.mjs'

const HERE = path.dirname(new URL(import.meta.url).pathname)
const TMP = mkdtempSync(path.join(tmpdir(), 'dsh-engine-queue-'))

const FAKE_WORKER = `#!/usr/bin/env python3
import json, sys, os, time

def reply(msg_id, result=None, error=None):
    payload = {"id": msg_id, "ok": error is None}
    if error is None:
        payload["result"] = result if result is not None else {}
    else:
        payload["error"] = {"code": "internal", "message": error}
    sys.stdout.write(json.dumps(payload) + "\\n")
    sys.stdout.flush()

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        msg = json.loads(line)
    except Exception:
        continue
    cmd = msg.get("cmd")
    if cmd == "ping":
        reply(msg.get("id"), {"mode": "builtin", "packages": {}, "layouts": ["umap"], "platform": "test"})
    elif cmd == "list":
        reply(msg.get("id"), {"datasets": [{"id": "dsA", "name": "A", "n_obs": 1, "n_vars": 1, "steps": []}]})
    elif cmd == "crash":
        reply(msg.get("id"), {"bye": True})
        sys.stdout.flush()
        os._exit(9)
    elif cmd == "slow":
        time.sleep(float(msg.get("args", {}).get("seconds") or 1.0))
        reply(msg.get("id"), {"slept": True})
    elif cmd == "shutdown":
        reply(msg.get("id"), {"bye": True})
        break
    else:
        reply(msg.get("id"), error="unknown cmd " + str(cmd))
`

const workerPath = path.join(TMP, 'fake_worker.py')
writeFileSync(workerPath, FAKE_WORKER, { mode: 0o755 })

let passed = 0
let failed = 0
const failures = []

function check(name, condition, detail = '') {
  if (condition) {
    passed += 1
    console.log(`  PASS ${name}${detail ? `  [${detail}]` : ''}`)
  } else {
    failed += 1
    failures.push(`${name}${detail ? `：${detail}` : ''}`)
    console.log(`  FAIL ${name}${detail ? `  [${detail}]` : ''}`)
  }
}

function withTimeout(promise, ms, label) {
  return Promise.race([
    promise.then((value) => ({ settled: 'ok', value })).catch((error) => ({ settled: 'error', error })),
    new Promise((resolve) => setTimeout(() => resolve({ settled: 'hung', label }), ms)),
  ])
}

async function main() {
  const workdir = path.join(TMP, 'work')
  const engine = createEngine({ workdir, engineScript: workerPath, python: process.env.DSH_TEST_PYTHON || 'python3' })

  // 1. 冷启动：第一条就是业务请求
  const cold = await withTimeout(engine.request('list', {}, { timeout: 20_000 }), 25_000, 'cold-list')
  check('HOST-02a 冷启动的业务请求能结束', cold.settled === 'ok',
    cold.settled === 'ok' ? 'ok' : String(cold.error && cold.error.message || cold.label))

  // 2. 先 start()，再崩溃 + 排队请求
  await engine.start()
  const crashPromise = withTimeout(engine.request('crash', {}, { timeout: 20_000 }), 20_000, 'crash')
  // 让 crash 先进入队列执行
  await new Promise((resolve) => setTimeout(resolve, 120))
  const queued = await withTimeout(
    engine.request('list', {}, { timeout: 6_000 }),
    20_000,
    'queued-after-crash',
  )
  await crashPromise
  check('HOST-02b 崩溃后排队请求在时限内结束（修复前会永久挂住）', queued.settled !== 'hung',
    queued.settled === 'ok' ? 'recovered' : queued.settled === 'error' ? `明确失败：${queued.error && queued.error.code}` : 'hung')

  // 3. 崩溃之后的新请求仍能执行
  const after = await withTimeout(engine.request('list', {}, { timeout: 20_000 }), 25_000, 'after-crash')
  check('HOST-02c 崩溃后新请求可继续执行', after.settled === 'ok',
    after.settled === 'ok' ? 'ok' : String(after.error && after.error.message))

  // 4. 没有悬空 pending
  const pending = engine._internal ? engine._internal.pending.size : -1
  check('HOST-02d 没有残留的 pending 请求', pending === 0, `pending=${pending}`)

  // 5. 排队阶段的超时兜底（慢请求占住队列时，后面的请求也能超时结束）
  const slow = withTimeout(engine.request('slow', { seconds: 3.0 }, { timeout: 10_000 }), 20_000, 'slow')
  const impatient = await withTimeout(engine.request('list', {}, { timeout: 300 }), 12_000, 'impatient')
  await slow
  check('HOST-02e 排队等待也受超时保护', impatient.settled !== 'hung',
    impatient.settled === 'ok' ? 'ok' : `${impatient.settled}:${impatient.error && impatient.error.code}`)

  await engine.dispose()
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }

  console.log(`\n== 结果：${passed} 通过 / ${failed} 失败 ==`)
  if (failed) {
    for (const item of failures) console.log(`  ✗ ${item}`)
    process.exit(1)
  }
}

main().catch((error) => {
  console.error('测试本身异常：', error)
  process.exit(1)
})
