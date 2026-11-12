"""
dsh-geobulk — 通用 Rscript 桥

用途：把 Python 侧数据打包成 bundle，调 `Rscript --vanilla lib/<script>.R --in <bundle> --out <out> --params <json>`，
再把 R 写出的 result.json / *.csv / *.png 读回来。

约定（与 CONTRACT.md §3.2 一致）：
  * bundle 目录：expr.tsv（样本×基因；首行基因名、首列样本名）、info.csv（首列样本名）、params.json、
    extra/<name>.csv（可选附加表）
  * R 脚本必须写 <out>/result.json（§2.5 信封字段可少不可错），表写 <out>/tables/<name>.csv，
    图写 <out>/figures/*.png（英文文字、白底、dpi≥140）
  * 退出码：0=成功；3=缺前置（转 prerequisite）；其它=失败（stderr 尾部作消息）
  * R 不可用 / 脚本缺失 / 超时 → 调用方负责回退，本文件只抛 EngineError

本文件只依赖标准库 + numpy/pandas + bg_common。禁止新增第三方依赖。
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import threading
import time

import numpy as np
import pandas as pd

from bg_common import EngineError, jsonable, png_size

# Rscript 候选（顺序即优先级）
R_CANDIDATES = [
    os.environ.get("DS_BG_RSCRIPT") or "",
    "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    "/usr/local/bin/Rscript",
    "/opt/homebrew/bin/Rscript",
    "/usr/bin/Rscript",
]

# ping/探测时检查的 R 包（都装好则官方实现可用）
R_PACKAGES = [
    "limma", "edgeR", "DESeq2", "survival", "clusterProfiler", "DOSE", "fgsea",
    "msigdbr", "GSVA", "WGCNA", "preprocessCore", "ComplexHeatmap", "pheatmap",
    "timeROC", "rms", "glmnet", "survminer", "org.Hs.eg.db", "org.Mm.eg.db",
    "IOBR", "Biobase", "GSEABase", "AnnotationDbi", "matrixStats",
]

_LOCK = threading.Lock()
_STATE: dict = {
    "checked_at": 0.0,
    "available": False,
    "script": None,
    "version": None,
    "packages": {},
    "error": None,
}
CACHE_SECONDS = 300


def candidate_scripts() -> list[str]:
    out = []
    for cand in R_CANDIDATES:
        if cand and os.path.exists(cand):
            out.append(cand)
    # 最后兜底：PATH 上的 Rscript
    which = shutil.which("Rscript")
    if which:
        out.append(which)
    # 去重保序
    seen = set()
    uniq = []
    for item in out:
        if item in seen:
            continue
        seen.add(item)
        uniq.append(item)
    return uniq


def _run(cmd: list[str], timeout: int) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env={**os.environ, "R_LIBS_USER": os.environ.get("R_LIBS_USER", "")},
    )


def probe(force: bool = False) -> dict:
    """探测 Rscript 与关键包，结果缓存 5 分钟。永不抛异常。"""
    now = time.time()
    with _LOCK:
        if not force and _STATE["checked_at"] and now - _STATE["checked_at"] < CACHE_SECONDS:
            return dict(_STATE)
        state = {
            "checked_at": now,
            "available": False,
            "script": None,
            "version": None,
            "packages": {},
            "error": None,
        }
        for script in candidate_scripts():
            try:
                res = _run([script, "--version"], timeout=20)
                text = (res.stdout or "") + (res.stderr or "")
                if "Rscript" not in text and "R version" not in text:
                    continue
                state["script"] = script
                state["available"] = True
                raw = text.strip()
                import re as _re

                match = _re.search(r"version\s+([0-9]+\.[0-9]+\.[0-9]+)", raw)
                state["version"] = match.group(1) if match else (raw.splitlines()[0] if raw else "R")
                break
            except Exception as error:  # noqa: BLE001
                state["error"] = str(error)
                continue
        if state["available"]:
            code = "cat(paste(vapply(c(%s), function(p) requireNamespace(p, quietly=TRUE), logical(1)), collapse=','))" % (
                ",".join('"%s"' % p for p in R_PACKAGES)
            )
            try:
                res = _run([state["script"], "--vanilla", "-e", code], timeout=120)
                flags = (res.stdout or "").strip().split(",")
                if len(flags) == len(R_PACKAGES):
                    state["packages"] = {p: flags[i].strip() == "TRUE" for i, p in enumerate(R_PACKAGES)}
                else:
                    state["packages"] = {}
            except Exception as error:  # noqa: BLE001
                state["packages"] = {}
                state["error"] = str(error)
        _STATE.update(state)
        return dict(_STATE)


def r_available() -> bool:
    return bool(probe().get("available"))


def r_info() -> dict:
    state = probe()
    return {
        "available": bool(state.get("available")),
        "version": state.get("version"),
        "script": state.get("script"),
        "packages": state.get("packages") or {},
        "error": state.get("error"),
    }


def has_pkg(name: str) -> bool:
    return bool((probe().get("packages") or {}).get(name))


def _write_extra(df: pd.DataFrame, path: str) -> None:
    frame = df.copy()
    frame.to_csv(path, sep="\t", index=True, encoding="utf-8")


def build_bundle(payload: dict, bundle_dir: str) -> str:
    os.makedirs(bundle_dir, exist_ok=True)
    expr = payload.get("expr")
    genes = payload.get("genes")
    samples = payload.get("samples")
    info = payload.get("info")
    if expr is not None:
        arr = np.asarray(expr, dtype=np.float64)
        if arr.ndim != 2:
            raise EngineError("internal", "run_r: expr 必须是二维数组")
        if samples is None:
            samples = [f"S{i + 1}" for i in range(arr.shape[0])]
        if genes is None:
            genes = [f"G{i + 1}" for i in range(arr.shape[1])]
        frame = pd.DataFrame(arr, index=[str(s) for s in samples], columns=[str(g) for g in genes])
        # 样本×基因 → TSV（首行基因名，首列样本名；避免 R 自动改列名）
        frame.to_csv(os.path.join(bundle_dir, "expr.tsv"), sep="\t", index=True, encoding="utf-8")
    if info is not None:
        frame = info if isinstance(info, pd.DataFrame) else pd.DataFrame(info)
        frame.to_csv(os.path.join(bundle_dir, "info.csv"), sep="\t", index=True, encoding="utf-8")
    params = payload.get("params") or {}
    with open(os.path.join(bundle_dir, "params.json"), "w", encoding="utf-8") as handle:
        json.dump(jsonable(params), handle, ensure_ascii=False, indent=2)
    extra = payload.get("extra") or {}
    if extra:
        extra_dir = os.path.join(bundle_dir, "extra")
        os.makedirs(extra_dir, exist_ok=True)
        for name, value in extra.items():
            frame = value if isinstance(value, pd.DataFrame) else pd.DataFrame(value)
            frame.to_csv(os.path.join(extra_dir, f"{name}.tsv"), sep="\t", index=True, encoding="utf-8")
    return bundle_dir


def _read_result(outdir: str) -> dict:
    result_path = os.path.join(outdir, "result.json")
    payload: dict = {}
    if os.path.exists(result_path):
        try:
            with open(result_path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except Exception as error:  # noqa: BLE001
            raise EngineError("internal", f"R 输出 result.json 解析失败：{error}")
    tables: dict[str, pd.DataFrame] = {}
    tables_dir = os.path.join(outdir, "tables")
    if os.path.isdir(tables_dir):
        for name in sorted(os.listdir(tables_dir)):
            if not name.endswith((".csv", ".tsv")):
                continue
            path = os.path.join(tables_dir, name)
            try:
                if name.endswith(".tsv"):
                    tables[os.path.splitext(name)[0]] = pd.read_csv(path, sep="\t", index_col=0)
                else:
                    tables[os.path.splitext(name)[0]] = pd.read_csv(path)
            except Exception:  # noqa: BLE001
                continue
    figures = []
    fig_dir = os.path.join(outdir, "figures")
    if os.path.isdir(fig_dir):
        for name in sorted(os.listdir(fig_dir)):
            if not name.lower().endswith(".png"):
                continue
            path = os.path.join(fig_dir, name)
            width, height = png_size(path)
            figures.append({
                "id": os.path.splitext(name)[0],
                "title": os.path.splitext(name)[0],
                "file": path,
                "width": width,
                "height": height,
            })
    if isinstance(payload, dict):
        if figures and not payload.get("figures"):
            payload["figures"] = figures
        if payload.get("tables") and tables:
            # R 侧可能只给表描述，用真实 DataFrame 回填
            for item in payload["tables"]:
                if isinstance(item, dict) and item.get("id") in tables and not item.get("rows"):
                    frame = tables[item["id"]]
                    item["columns"] = item.get("columns") or [
                        {"key": str(c), "label": str(c), "type": "number" if pd.api.types.is_numeric_dtype(frame[c]) else "string"}
                        for c in frame.columns
                    ]
                    item["rows"] = jsonable(frame.reset_index().to_dict(orient="records"))[:500]
                    item.setdefault("total", int(len(frame)))
        if not payload.get("tables") and tables:
            # R 脚本只写了 tables/*.csv 没写 result.json 的 tables：自动生成表项
            generated = []
            for name, frame in tables.items():
                generated.append({
                    "id": str(name),
                    "title": str(name),
                    "columns": [
                        {"key": str(c), "label": str(c),
                         "type": "number" if pd.api.types.is_numeric_dtype(frame[c]) else "string"}
                        for c in frame.columns
                    ],
                    "rows": jsonable(frame.head(500).to_dict(orient="records")),
                    "total": int(len(frame)),
                })
            payload["tables"] = generated
        payload.setdefault("tables", [])
        payload.setdefault("figures", figures)
        payload.setdefault("warnings", [])
        payload.setdefault("summary", {"title": "", "bullets": [], "metrics": {}})
    return payload


def run_r(script_name: str, payload: dict, outdir: str, timeout: int = 900) -> dict:
    """调 R 脚本并读回结果。

    返回 {"ok":True,"result":dict,"tables":{name:DataFrame},"figures":[...],"warnings":[],"stdout":str}
    失败抛 EngineError（prerequisite / internal / timeout / no_r）。
    """
    state = probe()
    if not state.get("available"):
        raise EngineError("no_r", "本机没有可用的 Rscript（R 官方实现不可用）")
    lib_dir = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(lib_dir, script_name)
    if not os.path.exists(script):
        raise EngineError("unsupported", f"R 脚本不存在：{script_name}")
    outdir = os.path.abspath(outdir)
    bundle = os.path.join(outdir, "r_bundle")
    r_out = os.path.join(outdir, "r_out")
    for path in (bundle, r_out):
        os.makedirs(path, exist_ok=True)
    build_bundle(payload, bundle)
    params_json = json.dumps(jsonable(payload.get("params") or {}), ensure_ascii=False)
    cmd = [state["script"], "--vanilla", script, "--in", bundle, "--out", r_out, "--params", params_json]
    t0 = time.time()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise EngineError("timeout", f"R 计算超时（>{timeout}s）：{script_name}")
    except Exception as error:  # noqa: BLE001
        raise EngineError("internal", f"R 调用失败：{error}")
    elapsed = int((time.time() - t0) * 1000)
    stderr_tail = (res.stderr or "")[-2000:]
    stdout_tail = (res.stdout or "")[-2000:]
    if res.returncode == 3:
        reason = stderr_tail.strip().splitlines()[-1] if stderr_tail.strip() else "R 侧数据不满足前置条件"
        raise EngineError("prerequisite", reason)
    if res.returncode != 0:
        message = stderr_tail.strip() or stdout_tail.strip() or f"Rscript 退出码 {res.returncode}"
        raise EngineError("internal", f"R 计算失败：{message[-800:]}")
    parsed = _read_result(r_out)
    tables = {}
    tables_dir = os.path.join(r_out, "tables")
    if os.path.isdir(tables_dir):
        for name in sorted(os.listdir(tables_dir)):
            if name.endswith(".csv"):
                try:
                    tables[os.path.splitext(name)[0]] = pd.read_csv(os.path.join(tables_dir, name))
                except Exception:  # noqa: BLE001
                    continue
            elif name.endswith(".tsv"):
                try:
                    tables[os.path.splitext(name)[0]] = pd.read_csv(os.path.join(tables_dir, name), sep="\t")
                except Exception:  # noqa: BLE001
                    continue
    return {
        "ok": True,
        "result": parsed,
        "tables": tables,
        "figures": parsed.get("figures") or [],
        "warnings": parsed.get("warnings") or [],
        "ms": elapsed,
        "stdout": stdout_tail,
    }


def try_run_r(script_name: str, payload: dict, outdir: str, timeout: int = 900) -> dict | None:
    """run_r 的软版本：R 不可用/脚本失败时返回 None，由调用方回退 Python。

    返回 {"__failed__": True, "error": "..."} 表示 R 尝试过但失败（调用方可记 warning）。
    """
    state = probe()
    if not state.get("available"):
        return None
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), script_name)
    if not os.path.exists(script):
        return None
    try:
        return run_r(script_name, payload, outdir, timeout=timeout)
    except EngineError as error:
        if error.code == "prerequisite":
            raise
        return {"__failed__": True, "error": error.message, "code": error.code}
