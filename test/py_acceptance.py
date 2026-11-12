#!/usr/bin/env python3
"""
dsh-geobulk — 端到端验收脚本（直接驱动 Python worker，不经过宿主）

流程：demo 数据集 → qc→filter→normalize→batch→pca→de→enrich→immune→survival→prognosis→scmap→ncrna
      → 挂载 circ 模块 → 跑若干高级分析 → audit → export(report/methods/protocol/script)

用法：
  cd plugins/dsh-geobulk
  python3 test/py_acceptance.py            # 跑全部（缺模块的步骤记为 SKIP 并继续）
  python3 test/py_acceptance.py --strict   # 任一步失败即非 0 退出
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ENGINE = os.path.join(ROOT, "lib", "bg_engine.py")
FIXTURES = os.path.join(HERE, "fixtures")

STRICT = "--strict" in sys.argv

results: list[dict] = []


class Worker:
    def __init__(self, workdir: str):
        self.workdir = workdir
        env = dict(os.environ)
        env.update({
            "MPLBACKEND": "Agg",
            "MPLCONFIGDIR": os.path.join(workdir, "mpl-cache"),
            "PYTHONUNBUFFERED": "1",
            "PYTHONIOENCODING": "utf-8",
        })
        self.proc = subprocess.Popen(
            [sys.executable, "-u", ENGINE, "--workdir", workdir],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            text=True,
        )
        self.seq = 0
        self.log = []

    def call(self, cmd: str, args: dict | None = None, timeout: float = 600) -> dict:
        if self.proc.poll() is not None:
            stderr = self.proc.stderr.read() if self.proc.stderr else ""
            raise RuntimeError(f"engine 已退出：{stderr[-2000:]}")
        self.seq += 1
        request_id = str(self.seq)
        payload = json.dumps({"id": request_id, "cmd": cmd, "args": args or {}}, ensure_ascii=False)
        self.proc.stdin.write(payload + "\n")
        self.proc.stdin.flush()
        started = time.time()
        while True:
            if time.time() - started > timeout:
                raise RuntimeError(f"{cmd} 超时（{timeout}s）")
            line = self.proc.stdout.readline()
            if not line:
                stderr = self.proc.stderr.read() if self.proc.stderr else ""
                raise RuntimeError(f"engine 输出结束（{cmd}）：{stderr[-2000:]}")
            line = line.strip()
            if not line:
                continue
            try:
                msg = json.loads(line)
            except json.JSONDecodeError:
                self.log.append(line)
                continue
            if msg.get("event") == "progress":
                continue
            if str(msg.get("id")) != request_id:
                continue
            if msg.get("ok"):
                return msg["result"]
            error = msg.get("error") or {}
            raise EngineError(error.get("code") or "internal", error.get("message") or "未知错误")

    def close(self):
        try:
            self.call("shutdown", {}, timeout=30)
        except Exception:  # noqa: BLE001
            pass
        try:
            self.proc.terminate()
        except Exception:  # noqa: BLE001
            pass


class EngineError(Exception):
    def __init__(self, code, message):
        super().__init__(f"[{code}] {message}")
        self.code = code
        self.message = message


def record(name: str, status: str, detail: str = "", seconds: float = 0.0):
    results.append({"name": name, "status": status, "detail": detail, "seconds": round(seconds, 1)})
    icon = {"PASS": "✓", "SKIP": "○", "FAIL": "✗"}[status]
    print(f"  {icon} {name:34s} {status:4s} {detail[:110]}" + (f"  [{seconds:.1f}s]" if seconds else ""))


def step(worker: Worker, dataset: str, step_id: str, params: dict):
    started = time.time()
    try:
        result = worker.call("run", {"dataset": dataset, "step": step_id, "params": params}, timeout=900)
    except EngineError as error:
        if error.code in ("unsupported", "prerequisite", "not_found"):
            record(f"step {step_id}", "SKIP", error.message)
            return None
        record(f"step {step_id}", "FAIL", error.message)
        if STRICT:
            raise
        return None
    elapsed = time.time() - started
    summary = (result.get("summary") or {})
    bullets = summary.get("bullets") or []
    figures = result.get("figures") or []
    tables = result.get("tables") or []
    detail = (bullets[0] if bullets else summary.get("title") or "")
    ok = bool(bullets)
    record(f"step {step_id}", "PASS" if ok else "FAIL",
           f"{detail[:80]} | {len(tables)} 表 {len(figures)} 图", elapsed)
    if not ok and STRICT:
        raise RuntimeError(f"{step_id} 没有 summary.bullets")
    return result


def analysis(worker: Worker, dataset: str, analysis_id: str, params: dict | None = None):
    started = time.time()
    try:
        result = worker.call("analysis", {"dataset": dataset, "id": analysis_id, "params": params or {}}, timeout=900)
    except EngineError as error:
        if error.code in ("unsupported", "prerequisite", "not_found"):
            record(f"analysis {analysis_id}", "SKIP", error.message)
            return None
        record(f"analysis {analysis_id}", "FAIL", error.message)
        if STRICT:
            raise
        return None
    elapsed = time.time() - started
    bullets = (result.get("summary") or {}).get("bullets") or []
    figures = result.get("figures") or []
    tables = result.get("tables") or []
    record(f"analysis {analysis_id}", "PASS" if (figures or tables) else "FAIL",
           f"{len(tables)} 表 {len(figures)} 图" + (f" | {bullets[0][:60]}" if bullets else ""), elapsed)
    return result


def main() -> int:
    workdir = tempfile.mkdtemp(prefix="bg-accept-")
    worker = Worker(workdir)
    print("== dsh-geobulk 端到端验收 ==")
    print(f"workdir: {workdir}")
    try:
        # 0. 引擎
        info = worker.call("ping", {}, timeout=180)
        r_ok = bool((info.get("r") or {}).get("available"))
        record("ping", "PASS",
               f"Python {info.get('python')} | R {'可用' if r_ok else '不可用'} | "
               f"模块 {sum(1 for v in (info.get('modules') or {}).values() if v.get('ok'))}/{len(info.get('modules') or {})}")

        # 1. demo
        demo = worker.call("open", {"demo": True}, timeout=180)
        dataset = demo["dataset"]["id"]
        record("demo 数据集", "PASS", f"{demo['dataset']['n_samples']} 样本 × {demo['dataset']['n_genes']} 基因，分组 {demo['dataset']['groups']}")

        # 2. 主流程
        step(worker, dataset, "qc", {})
        step(worker, dataset, "filter", {})
        step(worker, dataset, "normalize", {"method": "auto"})
        step(worker, dataset, "batch", {"batch_key": "batch", "method": "combat"})
        step(worker, dataset, "pca", {"n_pcs": 20})
        step(worker, dataset, "de", {"group_col": "group", "case": "case", "control": "control", "method": "limma"})
        step(worker, dataset, "enrich", {"db": ["GO_BP", "KEGG", "HALLMARK"], "method": "ora"})
        step(worker, dataset, "immune", {"methods": ["estimate", "ssgsea", "mcpcounter"]})
        step(worker, dataset, "survival", {"time_col": "os_time", "event_col": "os_event", "cutoff": "median"})
        step(worker, dataset, "prognosis", {"model": "cox"})

        sc_ref = os.path.join(FIXTURES, "demo_sc.h5ad")
        if os.path.exists(sc_ref):
            step(worker, dataset, "scmap", {"ref_path": sc_ref, "gene_set": "de_up", "score_method": "mean"})
        else:
            record("step scmap", "SKIP", "缺少 test/fixtures/demo_sc.h5ad")

        circ = os.path.join(FIXTURES, "circ_expr.tsv")
        if os.path.exists(circ):
            try:
                attached = worker.call("attach", {"dataset": dataset, "module": "circ", "path": circ}, timeout=300)
                record("attach circ", "PASS", f"{attached.get('module', {}).get('n_features')} 个 circRNA")
            except EngineError as error:
                record("attach circ", "SKIP", error.message)
        else:
            record("attach circ", "SKIP", "缺少 circ fixture")
        step(worker, dataset, "ncrna", {"kind": "circ"})

        # 3. 高级分析抽样
        catalogue = worker.call("analyses", {"dataset": dataset}, timeout=300)
        available = [item for item in catalogue.get("analyses", []) if item.get("available")]
        record("analyses 目录", "PASS" if available else "FAIL",
               f"{len(catalogue.get('analyses', []))} 项，可用 {len(available)} 项")
        sample_ids = []
        for preferred in ("ora_go", "gsea_hallmark", "immune_group_diff", "km_curve", "cox_uni",
                          "risk_score_cox", "de_summary", "sc_de_map", "circ_de", "sample_correlation",
                          "volcano_custom", "gsva_pathway"):
            if any(item["id"] == preferred for item in available):
                sample_ids.append(preferred)
        for analysis_id in sample_ids[:8]:
            analysis(worker, dataset, analysis_id)

        # 4. 自检 / 建议 / 预设
        audit = worker.call("audit", {"dataset": dataset}, timeout=120)
        record("audit", "PASS" if isinstance(audit.get("score"), int) else "FAIL",
               f"score={audit.get('score')} level={audit.get('level')} checks={len(audit.get('checks', []))}")
        suggest = worker.call("suggest", {"dataset": dataset}, timeout=120)
        record("suggest", "PASS" if suggest.get("suggestions") is not None else "FAIL",
               f"{len(suggest.get('suggestions', []))} 条建议")
        recipes = worker.call("recipe", {}, timeout=60)
        record("recipe", "PASS" if len(recipes.get("recipes", [])) >= 3 else "FAIL",
               f"{len(recipes.get('recipes', []))} 个预设")

        # 5. 导出
        for what in ("summary_json", "methods", "report", "protocol", "script", "de_csv", "expr_csv"):
            try:
                exported = worker.call("export", {"dataset": dataset, "what": what}, timeout=300)
                size = exported.get("bytes") or 0
                ok = size > 0 and os.path.exists(exported.get("file", ""))
                record(f"export {what}", "PASS" if ok else "FAIL", f"{size} bytes", 0)
            except EngineError as error:
                record(f"export {what}", "SKIP" if error.code == "prerequisite" else "FAIL", error.message)

        # 6. 图
        for kind in ("pca_scatter", "volcano", "sample_corr", "km_curve", "immune_stack", "gene_box"):
            try:
                figure = worker.call("plot", {"dataset": dataset, "kind": kind, "params": {"genes": ["CD8A", "MKI67", "EPCAM"]}}, timeout=300)
                ok = os.path.exists(figure.get("file", "")) and os.path.getsize(figure.get("file", "")) > 1000
                record(f"plot {kind}", "PASS" if ok else "FAIL", os.path.basename(figure.get("file") or ""), 0)
            except EngineError as error:
                record(f"plot {kind}", "SKIP" if error.code == "prerequisite" else "FAIL", error.message)
    finally:
        worker.close()

    total = len(results)
    passed = sum(1 for item in results if item["status"] == "PASS")
    skipped = sum(1 for item in results if item["status"] == "SKIP")
    failed = sum(1 for item in results if item["status"] == "FAIL")
    print("\n" + "=" * 78)
    print(f"总计 {total}：PASS {passed} / SKIP {skipped} / FAIL {failed}")
    if failed:
        print("失败项：")
        for item in results:
            if item["status"] == "FAIL":
                print(f"  - {item['name']}: {item['detail']}")
    print("=" * 78)
    if failed:
        return 1
    if STRICT and skipped:
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:  # noqa: BLE001
        print(f"\n未捕获错误：{error}")
        sys.exit(1)
