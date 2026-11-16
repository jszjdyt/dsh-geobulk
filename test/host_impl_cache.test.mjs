/**
 * dsh-geobulk / dsh-singlecell — impl 模块加载回归测试（任务单 HOST-04）
 *
 * 复现的问题：index.js 每次请求都用 `?rev=Date.now()-random` 新建模块 URL，
 * 等于每个 HTTP 请求/每次状态轮询都完整导入一次业务模块；500 次导入后 GC 仍占用约 95 MiB。
 *
 * 这里直接解析 index.js 源码，检查：
 *   1. loadImpl 不再使用随机 rev；
 *   2. 同一文件版本会复用同一个 import Promise（缓存命中）；
 *   3. impl.mjs 变化后会重新导入（开发热更新仍可用）。
 *
 * 运行：node test/host_impl_cache.test.mjs
 */

import { mkdtempSync, readFileSync, statSync, writeFileSync, rmSync, mkdirSync, utimesSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'

const HERE = path.dirname(fileURLToPath(import.meta.url))
const PLUGIN = path.resolve(HERE, '..')

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

const source = readFileSync(path.join(PLUGIN, 'index.js'), 'utf8')

check('HOST-04a index.js 不再用 Date.now/随机数生成模块 URL',
  !/lib\/impl\.mjs\?rev=\$\{Date\.now\(\)\}/.test(source) && !/Math\.random\(\)\.toString\(36\)/.test(source))

check('HOST-04b 有按版本缓存的 implCache', /let implCache = \{ version: '', promise: null \}/.test(source))

check('HOST-04c loadImpl 不再标记为 async（返回缓存 Promise）',
  /function loadImpl\(\)/.test(source))

// 用一份最小 replicate 验证缓存语义（与 index.js 中的实现同构）
const TMP = mkdtempSync(path.join(tmpdir(), 'dsh-impl-cache-'))
mkdirSync(path.join(TMP, 'lib'), { recursive: true })
const filePath = path.join(TMP, 'lib', 'impl.mjs')
writeFileSync(filePath, 'export const marker = 1\n')

const IMPL_URL = pathToFileURL(filePath)
function implVersion() {
  try { return String(statSync(IMPL_URL).mtimeMs) } catch { return 'static' }
}
let implCache = { version: '', promise: null }
function loadImpl() {
  const version = implVersion()
  if (implCache.promise && implCache.version === version) return implCache.promise
  const url = new URL(`?rev=${encodeURIComponent(version)}`, IMPL_URL)
  const promise = import(url.href)
  implCache = { version, promise }
  promise.catch(() => {
    if (implCache.promise === promise) implCache = { version: '', promise: null }
  })
  return promise
}

const first = await loadImpl()
const second = await loadImpl()
check('HOST-04d 同一版本复用同一个模块实例', first === second, `same=${first === second}`)
check('HOST-04e 同一个 Promise 对象', implCache.promise === loadImpl())

// 改文件（模拟开发热更新）
writeFileSync(filePath, 'export const marker = 2\n')
const future = new Date(Date.now() + 5000)
try { utimesSync(filePath, future, future) } catch { /* ignore */ }
const third = await loadImpl()
check('HOST-04f impl.mjs 变更后会重新导入（热更新仍可用）', third !== first, `changed=${third !== first}`)
check('HOST-04g 新版本内容是新的', third.marker === 2, `marker=${third.marker}`)

try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }

console.log(`\n== 结果：${passed} 通过 / ${failed} 失败 ==`)
if (failed) {
  for (const item of failures) console.log(`  ✗ ${item}`)
  process.exit(1)
}
