/**
 * dsh-geobulk / dsh-singlecell — 超时清理回归测试（任务单 HOST-05）
 *
 * 复现的问题：请求超时只 kill 直接子进程（Python worker），worker 自己拉起的
 * R / 外部 Python 子进程会继续跑。
 *
 * 这里用假 worker 自己再 spawn 一个长时间 sleep 的子进程，把它的 pid 写进文件；
 * 超时后断言该孙进程也退出。
 *
 * 运行：node test/host_process_tree.test.mjs
 */

import { mkdtempSync, writeFileSync, readFileSync, existsSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath } from 'node:url'
import { createEngine } from '../lib/engine.mjs'

const HERE = path.dirname(fileURLToPath(import.meta.url))
const TMP = mkdtempSync(path.join(tmpdir(), 'dsh-proc-tree-'))

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

const pidFile = path.join(TMP, 'grandchild.pid')
const FAKE_WORKER = `#!/usr/bin/env python3
import json, os, subprocess, sys, time

PIDFILE = ${JSON.stringify(pidFile)}

def reply(msg_id, result):
    sys.stdout.write(json.dumps({"id": msg_id, "ok": True, "result": result}) + "\\n")
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
        reply(msg.get("id"), {"mode": "builtin", "packages": {}, "layouts": ["umap"]})
    elif cmd == "spawnchild":
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
        with open(PIDFILE, "w") as handle:
            handle.write(str(child.pid))
        reply(msg.get("id"), {"child_pid": child.pid})
        # 停在这里不回复下一个请求，模拟卡住
        time.sleep(300)
    elif cmd == "shutdown":
        reply(msg.get("id"), {"bye": True})
        break
    else:
        reply(msg.get("id"), {})
`

const workerPath = path.join(TMP, 'fake_worker.py')
writeFileSync(workerPath, FAKE_WORKER, { mode: 0o755 })

function alive(pid) {
  try {
    process.kill(pid, 0)
    return true
  } catch {
    return false
  }
}

async function main() {
  const engine = createEngine({
    workdir: path.join(TMP, 'work'),
    engineScript: workerPath,
    python: process.env.DSH_TEST_PYTHON || 'python3',
  })

  const first = await engine.request('spawnchild', {}, { timeout: 10_000 })
  const childPid = Number(first && first.child_pid)
  check('HOST-05a 假 worker 成功拉起子进程', Number.isFinite(childPid) && childPid > 0, `pid=${childPid}`)
  check('HOST-05b 子进程初始存活', alive(childPid))

  // 这个请求会超时 → 宿主应杀掉整个进程组
  let timeoutError = null
  try {
    await engine.request('list', {}, { timeout: 800 })
  } catch (error) {
    timeoutError = error
  }
  check('HOST-05c 请求按预期超时', timeoutError !== null && (timeoutError.code === 'timeout' || timeoutError.code === 'engine_exited'),
    timeoutError ? String(timeoutError.code) : 'no error')

  // 给 SIGTERM/SIGKILL 一点时间
  let gone = false
  for (let i = 0; i < 40; i += 1) {
    if (!alive(childPid)) { gone = true; break }
    await new Promise((resolve) => setTimeout(resolve, 100))
  }
  check('HOST-05d 超时后孙进程也被清理（修复前仍存活）', gone, gone ? 'exited' : 'still alive')

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
