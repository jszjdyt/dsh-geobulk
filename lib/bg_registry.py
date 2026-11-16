"""
dsh-geobulk — 模块注册表与调度层

职责：
  * 动态导入全部步骤/分析模块（缺文件只跳过，保证增量交付可用）
  * 汇总 STEPS / ANALYSES 目录（含可用性判断）
  * 构造 ctx（§2.3）、执行步骤与分析、写回样本列/矩阵替换/结果缓存
  * 样本表操作、模块挂载、预览、导出、自检（audit）、建议（suggest）、预设流程（recipe）

本文件由 Lead 维护。模块只依赖 bg_common / rbridge，不 import 本文件（避免循环）。
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import pkgutil
import re
import sys
import time
import traceback

import numpy as np
import pandas as pd

from bg_common import (
    EngineError,
    data_file,
    ensure_dir,
    export_csv,
    export_json,
    jsonable,
    result,
    table,
)

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

MODULES = ["bg_prep", "bg_de", "bg_enrich", "bg_immune", "bg_survival", "bg_prognosis", "bg_scmap", "bg_ncrna"]

# 12 个流程步骤：UI 左栏固定顺序（模块缺失时仍显示为"未就绪"）
STEP_ORDER = [
    {"id": "qc", "title": "样本质控", "requires": []},
    {"id": "filter", "title": "过滤", "requires": ["qc"]},
    {"id": "normalize", "title": "标准化", "requires": ["filter"]},
    {"id": "batch", "title": "去批次", "requires": ["normalize"]},
    {"id": "pca", "title": "降维聚类", "requires": ["normalize"]},
    {"id": "de", "title": "差异表达", "requires": ["normalize"]},
    {"id": "enrich", "title": "功能富集", "requires": ["de"]},
    {"id": "immune", "title": "免疫浸润", "requires": ["normalize"]},
    {"id": "survival", "title": "临床生存", "requires": ["normalize"]},
    {"id": "prognosis", "title": "高级预后", "requires": ["survival"]},
    {"id": "scmap", "title": "单细胞映射", "requires": ["de"]},
    {"id": "ncrna", "title": "非编码/融合/剪接", "requires": ["normalize"]},
]
STEP_INDEX = {item["id"]: i for i, item in enumerate(STEP_ORDER)}

_RUNTIME: dict = {"modules": {}, "errors": {}, "loaded_at": 0.0}


# ---------------------------------------------------------------------------
# 模块加载
# ---------------------------------------------------------------------------

def load_modules(force: bool = False) -> dict:
    """导入所有模块，收集 STEPS/ANALYSES；单模块失败不影响其它模块。"""
    if not force and _RUNTIME["modules"]:
        return _RUNTIME
    modules: dict[str, dict] = {}
    errors: dict[str, str] = {}
    for name in MODULES:
        # 不跳过上次失败的模块：文件可能刚落盘，retry 一次（失败的 import 不会被 Python 缓存）
        try:
            mod = importlib.import_module(name)
            steps = []
            for item in getattr(mod, "STEPS", []) or []:
                if not isinstance(item, dict) or "id" not in item or "run" not in item:
                    continue
                steps.append({
                    "id": str(item["id"]),
                    "title": str(item.get("title") or item["id"]),
                    "requires": [str(x) for x in (item.get("requires") or [])],
                    "params": _norm_params(item.get("params") or []),
                    "run": item["run"],
                    "module": name,
                })
            analyses = []
            for item in getattr(mod, "ANALYSES", []) or []:
                if not isinstance(item, dict) or "id" not in item or "run" not in item:
                    continue
                analyses.append({
                    "id": str(item["id"]),
                    "title": str(item.get("title") or item["id"]),
                    "group": str(item.get("group") or "其它"),
                    "method": str(item.get("method") or ""),
                    "needs": [str(x) for x in (item.get("needs") or [])],
                    "params": _norm_params(item.get("params") or []),
                    "run": item["run"],
                    "module": name,
                })
            modules[name] = {"steps": steps, "analyses": analyses, "ok": True, "error": None}
        except Exception as error:  # noqa: BLE001
            modules[name] = {"steps": [], "analyses": [], "ok": False, "error": f"{type(error).__name__}: {error}"}
            errors[name] = modules[name]["error"]
    _RUNTIME["modules"] = modules
    _RUNTIME["errors"] = errors
    _RUNTIME["loaded_at"] = time.time()
    return _RUNTIME


def _norm_params(params) -> list[dict]:
    out = []
    for item in params or []:
        if not isinstance(item, dict):
            continue
        param = dict(item)
        param.setdefault("key", "")
        param.setdefault("label", param["key"])
        param.setdefault("type", "text")
        param.setdefault("default", None)
        out.append(param)
    return out


def module_status() -> dict:
    load_modules()
    return {
        name: {"ok": bool(info.get("ok")), "error": info.get("error"),
               "steps": [s["id"] for s in info.get("steps", [])],
               "analyses": [a["id"] for a in info.get("analyses", [])]}
        for name, info in _RUNTIME["modules"].items()
    }


def _module_steps() -> dict[str, dict]:
    found: dict[str, dict] = {}
    for info in load_modules()["modules"].values():
        for step in info.get("steps", []):
            found[step["id"]] = step
    return found


def _module_analyses() -> dict[str, dict]:
    found: dict[str, dict] = {}
    for info in load_modules()["modules"].values():
        for analysis in info.get("analyses", []):
            found.setdefault(analysis["id"], analysis)
    return found


# ---------------------------------------------------------------------------
# 目录
# ---------------------------------------------------------------------------

def steps_catalog() -> list[dict]:
    """UI 左栏流程卡：固定 12 项 + 模块提供的参数与就绪状态。"""
    implemented = _module_steps()
    out = []
    for base in STEP_ORDER:
        item = dict(base)
        impl = implemented.get(base["id"])
        if impl:
            item["params"] = impl["params"]
            item["requires"] = impl["requires"] or base["requires"]
            item["available"] = True
            item["module"] = impl["module"]
        else:
            item["params"] = []
            item["available"] = False
            item["module"] = None
            item["reason"] = "该步骤的模块尚未就绪"
        out.append(item)
    return out


def analyses_catalog(dataset) -> dict:
    found = _module_analyses()
    groups: list[str] = []
    items = []
    for analysis in found.values():
        ok, reason = needs_ok(dataset, analysis["needs"])
        if analysis["group"] not in groups:
            groups.append(analysis["group"])
        items.append({
            "id": analysis["id"],
            "title": analysis["title"],
            "group": analysis["group"],
            "method": analysis["method"],
            "needs": analysis["needs"],
            "params": analysis["params"],
            "available": ok,
            "reason": None if ok else reason,
            "module": analysis["module"],
        })
    # 组顺序：基本分析 → 差异分析 → 通路分析 → 免疫浸润 → 临床生存 → 高级预后 → 单细胞 → 非编码 → 其它
    preferred = ["基本分析", "差异分析", "通路分析", "免疫浸润", "临床生存", "高级预后", "单细胞", "非编码", "其它"]
    groups.sort(key=lambda g: (preferred.index(g) if g in preferred else len(preferred), g))
    items.sort(key=lambda it: (groups.index(it["group"]) if it["group"] in groups else 99, it["title"]))
    return {"analyses": items, "groups": groups}


def needs_ok(dataset, needs) -> tuple[bool, str | None]:
    if not needs:
        return True, None
    normalized = [str(n).strip() for n in needs if str(n).strip()]
    for need in normalized:
        if need.startswith("modules:"):
            kind = need.split(":", 1)[1]
            if kind not in getattr(dataset, "modules", {}):
                return False, f"需要先挂载 {kind} 数据"
            continue
        if need == "de":
            if "de" not in dataset.results and "de" not in dataset.steps:
                return False, "需要先跑「差异表达」"
            continue
        if need == "groups":
            if not getattr(dataset, "group_key", None):
                return False, "需要先指定分组列"
            continue
        if need == "clinical":
            if not dataset.has_clinical():
                return False, "样本信息里缺少生存时间/事件列，请先导入临床表"
            continue
        if need == "sc":
            if not getattr(dataset, "sc_ref", None):
                return False, "需要先指定单细胞参考（h5ad）"
            continue
        if need in ("normalize", "pca", "filter", "qc", "batch"):
            if need not in dataset.steps:
                return False, f"需要先跑「{need}」步骤"
            continue
        if need == "enrich":
            if "enrich" not in dataset.results and "enrich" not in dataset.steps:
                return False, "需要先跑「功能富集」"
            continue
        if need == "immune":
            if "immune" not in dataset.results and "immune" not in dataset.steps:
                return False, "需要先跑「免疫浸润」"
            continue
        if need == "survival":
            if "survival" not in dataset.results and "survival" not in dataset.steps:
                return False, "需要先跑「临床生存」"
            continue
        if need in ("genes", "gene_order", "file"):
            continue
        # 未知 need：不阻塞
        continue
    return True, None


# ---------------------------------------------------------------------------
# ctx
# ---------------------------------------------------------------------------

def build_ctx(dataset, outdir: str, engine_info: dict | None = None) -> dict:
    ensure_dir(outdir)
    groups = dataset.groups
    return {
        "dataset_id": dataset.id,
        "name": dataset.name,
        "genes": list(dataset.genes),
        "samples": list(dataset.samples),
        "X": dataset.X,
        "X_kind": dataset.value_kind,
        "counts": dataset.counts,
        "sample_info": dataset.sample_info,
        "group_key": dataset.group_key,
        "groups": groups,
        "species": dataset.species,
        "de": dataset.results.get("de"),
        "results": dataset.results,
        "modules": dataset.modules,
        "rng": np.random.default_rng(0),
        "workdir": os.path.dirname(outdir.rstrip(os.sep)) or outdir,
        "outdir": outdir,
        "engine": engine_info or {},
        "sc_ref": dataset.sc_ref,
        "notes": list(dataset.notes),
        "warnings": list(dataset.warnings),
    }


def _normalize_store(value):
    """把模块返回的缓存值统一成可长期保存的结构（DataFrame 留在内存）。"""
    if not isinstance(value, dict):
        return value
    out = dict(value)
    for key in ("table", "tables"):
        frame = out.get(key)
        if isinstance(frame, pd.DataFrame):
            continue
        if isinstance(frame, list) and frame and isinstance(frame[0], dict):
            try:
                out[key] = pd.DataFrame(frame)
            except Exception:  # noqa: BLE001
                continue
    if "de_table" in out and not isinstance(out["de_table"], pd.DataFrame):
        rows = out.get("de_table")
        if isinstance(rows, list) and rows and isinstance(rows[0], dict):
            out["de_table"] = pd.DataFrame(rows)
    # 纯 JSON 结构保持原样
    return out


def _apply_result(dataset, payload: dict) -> dict:
    """应用结果的副作用（replace / store / sample_columns），返回给前端的干净信封。"""
    if not isinstance(payload, dict):
        return result(title="完成")
    outgoing = dict(payload)
    replace = outgoing.pop("replace", None) or outgoing.pop("_replace", None)
    if isinstance(replace, dict) and replace:
        dataset.apply_replace(replace)
    # de_store 是「一个结果」（含 group_col/case/control/table…），不是 {store_name: value}
    de_store = outgoing.pop("de_store", None)
    if isinstance(de_store, dict):
        try:
            dataset.results["de"] = _normalize_store(de_store)
        except Exception:  # noqa: BLE001
            pass
    store: dict = {}
    for key in ("result_store", "store", "_store"):
        value = outgoing.pop(key, None)
        if isinstance(value, dict):
            store.update(value)
    for key, value in store.items():
        try:
            dataset.results[str(key)] = _normalize_store(value)
        except Exception:  # noqa: BLE001
            continue
    cols = outgoing.get("sample_columns")
    vals = outgoing.get("sample_values")
    if cols and vals:
        dataset.apply_sample_columns(cols, vals)
    # 保证信封字段齐全
    summary = outgoing.get("summary")
    if not isinstance(summary, dict):
        summary = {"title": "", "bullets": [], "metrics": {}}
    summary.setdefault("title", "")
    summary.setdefault("bullets", [])
    summary.setdefault("metrics", {})
    outgoing["summary"] = summary
    for key in ("tables", "figures", "warnings", "notes"):
        if not isinstance(outgoing.get(key), list):
            outgoing[key] = []
    return jsonable(outgoing)


# ---------------------------------------------------------------------------
# 结构化执行记录（EXPORT-01 / EXPORT-02 的数据来源）
#
# 只在 registry 层补记录，不改 Dataset 定义：每次成功/失败的步骤与分析都往
# dataset.history / dataset.analyses_run 里附加「真实执行」的结构化字段：
#   metrics      模块返回的 summary.metrics（actual_method/design/value_kind 等）
#   notes / warnings  模块返回的执行说明与降级原因
#   data_version 运行**之前**的数据指纹（矩阵/基因/样本/分组），用于判断记录
#                是否属于当前数据版本，也用于导出脚本的逐步重放校验
# 取不到的字段宁可写「未记录」，也不按按钮名补写方法。
# ---------------------------------------------------------------------------

def _sha256_text(text: str) -> str:
    return hashlib.sha256(str(text).encode("utf-8")).hexdigest()


def _column_names(frame) -> list:
    """pandas Index 的 in 判断要显式转 list（Index 的布尔值有歧义）。"""
    columns = getattr(frame, "columns", None)
    if columns is None:
        return []
    try:
        return [str(item) for item in list(columns)]
    except TypeError:
        return []


def _data_version(dataset) -> dict:
    """当前数据状态指纹：与导出脚本里的同名函数保持同一算法。"""
    X = np.ascontiguousarray(np.asarray(getattr(dataset, "X", np.zeros((0, 0))), dtype=np.float32))
    genes = [str(g) for g in (getattr(dataset, "genes", None) or [])]
    samples = [str(s) for s in (getattr(dataset, "samples", None) or [])]
    group_key = getattr(dataset, "group_key", None)
    group_sha256 = None
    info = getattr(dataset, "sample_info", None)
    if group_key and info is not None and group_key in _column_names(info):
        labels = [str(v) for v in info[group_key].astype(str).tolist()]
        group_sha256 = _sha256_text("\n".join(labels))
    return {
        "n_samples": int(X.shape[0]) if X.ndim == 2 else 0,
        "n_genes": int(X.shape[1]) if X.ndim == 2 else 0,
        "value_kind": str(getattr(dataset, "value_kind", "unknown") or "unknown"),
        "matrix_sha256": hashlib.sha256(X.data if X.size else b"").hexdigest(),
        "genes_sha256": _sha256_text("\n".join(genes)),
        "samples_sha256": _sha256_text("\n".join(samples)),
        "group_key": None if not group_key else str(group_key),
        "group_sha256": group_sha256,
        "dataset_version": getattr(dataset, "version", None),
        "dataset_revision": getattr(dataset, "revision", None),
    }


def _payload_lists(payload) -> tuple[list[str], list[str]]:
    if not isinstance(payload, dict):
        return [], []
    notes = [str(x) for x in (payload.get("notes") or []) if str(x).strip()]
    warnings = [str(x) for x in (payload.get("warnings") or []) if str(x).strip()]
    return notes, warnings


def _payload_metrics(payload) -> dict:
    if not isinstance(payload, dict):
        return {}
    summary = payload.get("summary")
    metrics = summary.get("metrics") if isinstance(summary, dict) else None
    return jsonable(metrics or {}) if isinstance(metrics, dict) else {}


# 结果快照里保留哪些字段（表格本身太大，不写入 history；de 另存前 20 行关键列）
_SNAPSHOT_FIELDS = {
    "de": ("group_col", "case", "control", "method", "method_requested", "lite", "params",
           "labels", "samples", "groups", "n_case", "n_ctrl", "n_up", "n_down", "n_sig"),
    "enrich": ("method", "databases", "n_terms", "n_significant", "actual_method",
               "software", "software_version", "permutation", "gene_set_source", "coverage"),
    "immune": ("engine", "requested", "methods_ok", "methods_failed", "primary_method",
               "primary_label", "primary_derived", "normalize", "n_samples", "n_cells"),
    "survival": ("time_col", "event_col", "cutoff", "gene", "logrank_p", "cindex"),
    "batch": ("batch_key", "method", "n_batches", "batch_r2_before", "batch_r2_after"),
}
_SNAPSHOT_FIELDS_BY_STEP = {
    "de": _SNAPSHOT_FIELDS["de"],
    "prognosis": ("time_col", "event_col", "genes", "cindex", "model", "model_type"),
}


def _payload_store(payload) -> dict:
    """模块返回里被 registry 写入 dataset.results 的部分（执行当时的结果）。"""
    out: dict = {}
    if not isinstance(payload, dict):
        return out
    de_store = payload.get("de_store")
    if isinstance(de_store, dict):
        out["de"] = de_store
    for key in ("result_store", "store", "_store"):
        value = payload.get(key)
        if isinstance(value, dict):
            out.update(value)
    return out


def _result_snapshot(payload, step: str) -> dict:
    """执行当时的结果快照：后续高级分析可能覆盖 dataset.results，Methods/重放要用这一刻的值。"""
    snapshot: dict = {}
    for key, store in _payload_store(payload).items():
        if not isinstance(store, dict) or not store:
            continue
        fields = _SNAPSHOT_FIELDS_BY_STEP.get(step) or _SNAPSHOT_FIELDS.get(key, ())
        item = {}
        for field in fields:
            if field in store:
                item[field] = jsonable(store[field])
        if key == "de":
            rows = _expected_de_rows(store)
            if rows:
                item["table_head"] = rows
        snapshot[str(key)] = item
    return snapshot


def _attach_step_record(dataset, payload, before: dict | None, error=None, step: str | None = None) -> None:
    """给最近一条 history 补结构化执行字段（不改变既有字段含义）。"""
    history = getattr(dataset, "history", None)
    if not isinstance(history, list) or not history:
        return
    entry = history[-1]
    try:
        notes, warnings = _payload_lists(payload)
        entry["metrics"] = _payload_metrics(payload)
        entry["notes"] = notes
        entry["warnings"] = warnings
        if before is not None:
            entry["data_version"] = jsonable(before)
        if payload is not None and step:
            snapshot = _result_snapshot(payload, step)
            if snapshot:
                entry["result_snapshot"] = snapshot
        if error:
            entry["error"] = str(error)
    except Exception:  # noqa: BLE001 - 记录失败不能影响主流程
        pass


def _attach_analysis_record(dataset, params, payload, before: dict | None, error=None) -> None:
    runs = getattr(dataset, "analyses_run", None)
    if not isinstance(runs, list) or not runs:
        return
    entry = runs[-1]
    try:
        notes, warnings = _payload_lists(payload)
        entry["params"] = jsonable(dict(params or {}))
        entry["metrics"] = _payload_metrics(payload)
        entry["notes"] = notes
        entry["warnings"] = warnings
        if before is not None:
            entry["data_version"] = jsonable(before)
        if error:
            entry["error"] = str(error)
    except Exception:  # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# 执行
# ---------------------------------------------------------------------------

def run_step(dataset, step_id: str, params: dict, workdir: str, engine_info: dict | None = None) -> dict:
    implemented = _module_steps()
    step = implemented.get(step_id)
    if step is None:
        raise EngineError("unsupported", f"步骤「{step_id}」尚未实现（对应模块未就绪）")
    requires = step["requires"] or next(
        (b["requires"] for b in STEP_ORDER if b["id"] == step_id), []
    )
    missing = [r for r in requires if r not in dataset.steps]
    if missing:
        titles = {b["id"]: b["title"] for b in STEP_ORDER}
        raise EngineError(
            "prerequisite",
            f"需要先完成：{'、'.join(titles.get(m, m) for m in missing)}",
        )
    outdir = ensure_dir(os.path.join(workdir, "analyses", dataset.id, f"step-{step_id}"))
    ctx = build_ctx(dataset, outdir, engine_info)
    params = dict(params or {})
    if dataset.group_key and "group_col" not in params:
        params["group_col"] = dataset.group_key
    t0 = time.time()
    dataset_id_before = dataset.id
    before = _data_version(dataset)  # 步骤输入状态：记录到 history，供 EXPORT-01/02 使用
    try:
        payload = step["run"](ctx, params, outdir)
    except EngineError as error:
        dataset.add_history(step_id, params, False, int((time.time() - t0) * 1000), summary=error.message, engine=engine_info)
        _attach_step_record(dataset, None, before, error=error.message, step=step_id)
        raise
    elapsed = int((time.time() - t0) * 1000)
    clean = _apply_result(dataset, payload)
    dataset.mark_step(step_id)
    bullets = (clean.get("summary") or {}).get("bullets") or []
    dataset.add_history(
        step_id,
        params,
        True,
        elapsed,
        summary=str(bullets[0]) if bullets else "",
        engine=engine_info,
    )
    _attach_step_record(dataset, payload, before, step=step_id)
    out = dict(clean)
    out.update({
        "step": step_id,
        "steps": list(dataset.steps),
        "current": dataset.current,
        "dataset": dataset.summary(),
        "ms": elapsed,
    })
    return out


def run_analysis(dataset, analysis_id: str, params: dict, workdir: str, engine_info: dict | None = None) -> dict:
    analyses = _module_analyses()
    analysis = analyses.get(analysis_id)
    if analysis is None:
        raise EngineError("not_found", f"分析「{analysis_id}」不存在（模块可能未就绪）")
    ok, reason = needs_ok(dataset, analysis["needs"])
    if not ok:
        raise EngineError("prerequisite", reason or "缺少前置数据")
    outdir = ensure_dir(os.path.join(workdir, "analyses", dataset.id, f"an-{analysis_id}"))
    ctx = build_ctx(dataset, outdir, engine_info)
    params = dict(params or {})
    if dataset.group_key and "group_col" not in params:
        params["group_col"] = dataset.group_key
    started = time.time()
    before = _data_version(dataset)
    try:
        payload = analysis["run"](ctx, params, outdir)
    except EngineError as error:
        dataset.add_analysis_run(analysis_id, False, int((time.time() - started) * 1000), error.message)
        _attach_analysis_record(dataset, params, None, before, error=error.message)
        raise
    elapsed = int((time.time() - started) * 1000)
    clean = _apply_result(dataset, payload)
    bullets = (clean.get("summary") or {}).get("bullets") or []
    dataset.add_analysis_run(analysis_id, True, elapsed, str(bullets[0]) if bullets else "")
    _attach_analysis_record(dataset, params, payload, before)
    out = dict(clean)
    out.update({
        "analysis_id": analysis_id,
        "title": analysis["title"],
        "group": analysis["group"],
        "method": analysis["method"],
        "dataset": dataset.summary(),
        "ms": elapsed,
    })
    return out


def run_plot(dataset, kind: str, params: dict, workdir: str, engine_info: dict | None = None) -> dict:
    import bg_plots

    outdir = ensure_dir(os.path.join(workdir, "analyses", dataset.id, "plots"))
    ctx = build_ctx(dataset, outdir, engine_info)
    return bg_plots.render(ctx, kind, dict(params or {}), outdir)


# ---------------------------------------------------------------------------
# 样本表 / 模块 / 预览
# ---------------------------------------------------------------------------

def handle_samples(dataset, args: dict) -> dict:
    op = str(args.get("op") or "list")
    info = dataset.sample_info
    if op == "list":
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset)}
    if op == "update":
        samples = args.get("samples") or []
        groups = args.get("groups") or []
        if isinstance(samples, list):
            for entry in samples:
                if not isinstance(entry, dict):
                    continue
                ident = str(entry.get("id") or entry.get("sample") or "")
                if ident not in info.index:
                    continue
                if entry.get("name"):
                    info.loc[ident, "display_name"] = str(entry.get("name"))
                if entry.get("group"):
                    info.loc[ident, "group"] = str(entry.get("group"))
        if groups:
            dataset.notes = [n for n in dataset.notes if not n.startswith("分组：")]
            dataset.notes.append("分组：" + "/".join(str(g) for g in groups))
        if dataset.group_key is None and "group" in info.columns:
            dataset.set_group("group")
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset)}
    if op == "auto_group":
        pattern = str(args.get("group_regex") or args.get("pattern") or "")
        if not pattern:
            raise EngineError("bad_request", "请提供 group_regex 正则（第 1 个捕获组=大组）")
        import re

        regex = re.compile(pattern)
        values = []
        for name in info.index:
            match = regex.search(str(name))
            values.append(match.group(1) if match and match.groups() else "other")
        info["group"] = values
        dataset.set_group("group")
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset)}
    if op == "derive":
        import re

        name = str(args.get("name") or "derived")
        pattern = str(args.get("group_regex") or args.get("source") or "")
        column = args.get("column")
        regex = re.compile(pattern) if pattern else None
        values = []
        for name_value in info.index:
            if regex is not None:
                match = regex.search(str(name_value))
                values.append(match.group(1) if match and match.groups() else "")
            elif column and str(column) in info.columns:
                values.append(info.loc[name_value, str(column)])
            else:
                values.append("")
        info[name] = values
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset)}
    if op == "set_group":
        column = str(args.get("column") or "")
        dataset.set_group(column)
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset)}
    if op == "join":
        path = str(args.get("path") or "")
        if not path:
            raise EngineError("bad_request", "请提供样本信息表路径（path）")
        frame = _read_info_table(path)
        key = args.get("key")
        if key and str(key) in frame.columns:
            frame = frame.set_index(str(key))
        frame.index = [str(i) for i in frame.index]
        matched = [i for i in info.index if i in frame.index]
        if not matched:
            # 尝试大小写/空白不敏感匹配
            lower = {str(i).strip().lower(): i for i in frame.index}
            remap = {}
            for ident in info.index:
                cand = lower.get(str(ident).strip().lower())
                if cand:
                    remap[cand] = ident
            if remap:
                frame = frame.rename(index=remap)
                matched = [i for i in info.index if i in frame.index]
        if not matched:
            raise EngineError(
                "bad_request",
                f"临床表没有样本能对上（表内示例：{list(frame.index[:3])}；数据集示例：{list(info.index[:3])}）",
            )
        for col in frame.columns:
            if col in info.columns and col != "sample":
                info[f"{col}_joined"] = frame[col].reindex(info.index)
            else:
                info[col] = frame[col].reindex(info.index)
        dataset.notes.append(f"已 join 样本信息表：{os.path.basename(path)}（匹配 {len(matched)}/{len(info)} 样本）")
        if dataset.group_key is None:
            guess = dataset._guess_group()
            if guess:
                dataset.group_key = guess
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset), "matched": len(matched)}
    if op == "remove":
        names = [str(x) for x in (args.get("samples") or [])]
        keep = [i for i in info.index if i not in set(names)]
        if not keep:
            raise EngineError("bad_request", "不能移除全部样本")
        mask = np.array([s in set(keep) for s in dataset.samples])
        dataset.X = dataset.X[mask]
        if dataset.counts is not None:
            dataset.counts = dataset.counts[mask]
        dataset.samples = [s for s in dataset.samples if s in set(keep)]
        dataset.sample_info = info.loc[keep]
        dataset.steps = [s for s in dataset.steps if s in ("qc", "filter")]
        dataset.results.clear()
        return {"dataset": dataset.summary(), "sample_info": _sample_info_table(dataset), "removed": len(names)}
    raise EngineError("bad_request", f"不支持的样本操作：{op}")


def _read_info_table(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith((".xlsx", ".xls")):
        frame = pd.read_excel(path, index_col=0)
    else:
        sep = "\t" if lower.endswith((".tsv", ".txt")) else ","
        try:
            frame = pd.read_csv(path, sep=sep, index_col=0)
        except Exception:  # noqa: BLE001
            frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    frame.columns = [str(c) for c in frame.columns]
    return frame


def _sample_info_table(dataset) -> list[dict]:
    frame = dataset.sample_info.copy()
    frame.insert(0, "sample", frame.index.astype(str))
    return jsonable(frame.reset_index(drop=True).to_dict(orient="records"))[:500]


def attach_module(dataset, kind: str, loaded: dict, path: str) -> dict:
    """把 bg_io 读到的模块数据挂到 dataset.modules[kind]。"""
    kind = str(kind)
    if kind not in ("circ", "mirna", "lncrna", "fusion", "splicing"):
        raise EngineError("bad_request", f"未知模块类型：{kind}（可选 circ/mirna/lncrna/fusion/splicing）")
    payload: dict = {"kind": kind, "path": path, "format": loaded.get("format")}
    X = loaded.get("X")
    if X is not None:
        payload["X"] = np.asarray(X, dtype=np.float32)
        payload["genes"] = list(loaded.get("genes") or loaded.get("features") or [])
        payload["samples"] = list(loaded.get("samples") or [])
        payload["n_features"] = len(payload["genes"])
    frame = loaded.get("table")
    if frame is None and isinstance(loaded.get("tables"), dict):
        tables = loaded.get("tables") or {}
        frame = next(iter(tables.values()), None)
    if isinstance(frame, pd.DataFrame):
        payload["table"] = frame
        payload["n_features"] = int(len(frame))
    info_value = loaded.get("sample_info")
    if info_value is None:
        info_value = loaded.get("info")
    payload["info"] = info_value
    payload["notes"] = list(loaded.get("notes") or [])
    payload["warnings"] = list(loaded.get("warnings") or [])
    dataset.modules[kind] = payload
    dataset.warnings.extend([w for w in payload["warnings"] if w not in dataset.warnings])
    return {
        "dataset": dataset.summary(),
        "module": {
            "kind": kind,
            "path": path,
            "format": payload.get("format"),
            "n_features": payload.get("n_features"),
            "n_samples": len(payload.get("samples") or []),
        },
    }


def preview(dataset) -> dict:
    return {
        "genes": [str(g) for g in dataset.genes[:50]],
        "n_genes": dataset.n_genes,
        "n_samples": dataset.n_samples,
        "samples": list(dataset.samples[:200]),
        "sample_columns": dataset.sample_columns(),
        "group_key": dataset.group_key,
        "groups": dataset.groups,
        "steps": list(dataset.steps),
        "current": dataset.current,
        "modules": list(dataset.modules.keys()),
        "value_kind": dataset.value_kind,
        "clinical": dataset.clinical_columns(),
        "sc_ref": dataset.sc_ref,
    }


# ---------------------------------------------------------------------------
# 建议 / 自检
# ---------------------------------------------------------------------------

def suggest(dataset) -> dict:
    steps = set(dataset.steps)
    suggestions: list[dict] = []

    def add(sid, title, detail, step=None, params=None):
        suggestions.append({"id": sid, "title": title, "detail": detail, "step": step, "params": params or {}})

    if not steps:
        add("run:qc", "样本质控", "先看文库大小、检出基因数与离群样本", "qc", {"group_col": dataset.group_key})
    if "qc" in steps and "filter" not in steps:
        add("run:filter", "过滤低表达基因", "过滤掉低表达基因与（可选）离群样本", "filter", {})
    if "filter" in steps and "normalize" not in steps:
        add("run:normalize", "标准化", "log2 转换 + 分位数标准化，便于样本间比较", "normalize", {"method": "auto"})
    if dataset.group_key and "normalize" in steps and "de" not in steps:
        add("run:de", "差异表达", f"按「{dataset.group_key}」分组做 limma 差异分析", "de",
            {"group_col": dataset.group_key})
    if "de" in steps and "enrich" not in steps:
        add("run:enrich", "功能富集", "对差异基因做 GO/KEGG/Hallmark 富集", "enrich", {"db": ["GO_BP", "KEGG"]})
    if "normalize" in steps and "immune" not in steps:
        add("run:immune", "免疫浸润", "CIBERSORT/EPIC/MCP-counter/xCell/quanTIseq/ESTIMATE/ssGSEA/TIMER", "immune", {})
    if dataset.has_clinical() and "survival" not in steps:
        add("run:survival", "临床生存", f"用 {dataset.clinical_columns()['time']} / {dataset.clinical_columns()['event']} 做 KM 与 Cox", "survival", {})
    if "de" in steps and "scmap" not in steps:
        add("run:scmap", "映射到单细胞", "把 bulk 差异签名映射到单细胞亚群，看是哪些细胞在表达", "scmap", {})
    if "normalize" in steps and not dataset.modules:
        add("analysis:ncrna", "非编码/融合/剪接", "挂载 circRNA/miRNA/lncRNA/fusion/rMATS 结果后分析", None, {})
    if "de" in steps and "enrich" in steps and "immune" in steps:
        add("export:report", "导出分析报告", "把 Methods/表格/图整理成 Markdown 报告", None, {})
    return {"suggestions": suggestions[:8]}


AUDIT_CHECKS = [
    ("qc_pass", "先质控再分析"),
    ("normalized", "标准化状态"),
    ("group_imbalanced", "组间样本平衡"),
    ("batch_confounded", "批次与分组混杂"),
    ("low_samples", "样本量"),
    ("de_missing", "差异分析"),
    ("enrich_missing", "功能富集"),
    ("immune_missing", "免疫浸润"),
    ("survival_missing", "临床生存"),
    ("sc_missing", "单细胞映射"),
    ("modules_missing", "非编码模块"),
]


def audit(dataset) -> dict:
    steps = set(dataset.steps)
    checks: list[dict] = []

    def add(cid, level, title, detail, hint=""):
        checks.append({"id": cid, "level": level, "title": title, "detail": detail, "hint": hint})

    # qc
    if "qc" in steps:
        add("qc_pass", "ok", "已做样本质控", "样本层面指标已计算")
    else:
        add("qc_pass", "warn", "还没做样本质控", "建议先看文库大小与离群样本", "运行「样本质控」")

    # normalized
    if dataset.value_kind in ("lognorm", "tpm", "fpkm"):
        add("normalized", "ok", "矩阵已适合作差异分析", f"当前值类型：{dataset.value_kind}")
    elif "normalize" in steps:
        add("normalized", "ok", "已标准化", "标准化步骤已完成")
    else:
        add("normalized", "fail", "矩阵尚未标准化", f"当前值类型：{dataset.value_kind}，直接做 DE 会有偏", "运行「标准化」")

    # group balance
    groups = dataset.groups or []
    if groups:
        counts = {g: int((dataset.group_labels() == g).sum()) for g in groups}
        small = [g for g, n in counts.items() if n < 3]
        if small:
            add("group_imbalanced", "fail", "存在样本量 <3 的分组", f"{small}", "补充样本或合并分组")
        elif len(set(counts.values())) > 1 and max(counts.values()) - min(counts.values()) >= 3:
            add("group_imbalanced", "warn", "组间样本数不均衡", f"{counts}", "注意差异分析效能")
        else:
            add("group_imbalanced", "ok", "组间样本量均衡", f"{counts}")
    else:
        add("group_imbalanced", "warn", "还没有分组列", "差异分析与生存分析需要分组/临床列", "在「样本」里设置分组列")

    # batch confounded
    info = dataset.sample_info
    batch_col = next((c for c in info.columns if str(c).lower() in ("batch", "batch_id", "plate", "run", "cohort")), None)
    if batch_col and dataset.group_key:
        try:
            ct = pd.crosstab(info[batch_col].astype(str), info[dataset.group_key].astype(str))
            if (ct.values > 0).sum(axis=1).max() == 1 and ct.shape[0] > 1:
                add("batch_confounded", "fail", "批次与分组完全混杂", f"批次列 {batch_col} 与分组 {dataset.group_key} 一一对应", "无法区分批次与生物学差异")
            else:
                add("batch_confounded", "ok", "批次与分组未完全混杂", f"批次列 {batch_col}")
        except Exception:  # noqa: BLE001
            pass
    if "batch" not in steps and batch_col and ct_multi(info, batch_col):
        add("batch_confounded", "warn", "存在多批次但未做去批次", f"批次列 {batch_col}", "运行「去批次」")

    # low samples
    if dataset.n_samples < 6:
        add("low_samples", "warn", "样本量偏少", f"当前 {dataset.n_samples} 个样本，统计效能有限")
    else:
        add("low_samples", "ok", "样本量可用", f"当前 {dataset.n_samples} 个样本")

    # downstream
    add("de_missing", "ok" if "de" in steps else "warn", "差异分析", "已完成" if "de" in steps else "尚未做差异表达", "" if "de" in steps else "运行「差异表达」")
    add("enrich_missing", "ok" if "enrich" in steps else "warn", "功能富集", "已完成" if "enrich" in steps else "尚未做功能富集", "" if "enrich" in steps else "运行「功能富集」")
    add("immune_missing", "ok" if "immune" in steps else "warn", "免疫浸润", "已完成" if "immune" in steps else "尚未做免疫浸润", "" if "immune" in steps else "运行「免疫浸润」")
    if dataset.has_clinical():
        add("survival_missing", "ok" if "survival" in steps else "warn", "临床生存", "已完成" if "survival" in steps else "检测到临床列但尚未做生存分析", "" if "survival" in steps else "运行「临床生存」")
    add("sc_missing", "ok" if "scmap" in steps else "warn", "单细胞映射", "已完成" if "scmap" in steps else "尚未映射到单细胞", "" if "scmap" in steps else "在「单细胞」tab 指定参考 h5ad")
    add("modules_missing", "ok" if dataset.modules else "warn", "非编码模块", "已挂载：" + ",".join(dataset.modules.keys()) if dataset.modules else "尚未挂载 circRNA/miRNA/lncRNA/fusion/剪接数据", "" if dataset.modules else "在「非编码」tab 挂载数据")
    checks = [c for c in checks if c["id"] in dict(AUDIT_CHECKS)]
    fails = sum(1 for c in checks if c["level"] == "fail")
    warns = sum(1 for c in checks if c["level"] == "warn")
    score = max(0, min(100, 100 - 30 * fails - 12 * warns))
    level = "fail" if fails else ("warn" if warns else "ok")
    return {"score": score, "level": level, "checks": checks}


def ct_multi(info, col) -> bool:
    try:
        return int(info[col].nunique(dropna=True)) > 1
    except Exception:  # noqa: BLE001
        return False


RECIPES = [
    {
        "id": "standard",
        "title": "标准分析流程",
        "detail": "质控→过滤→标准化→降维→差异→富集→免疫",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {"method": "auto"}},
            {"step": "pca", "params": {"n_pcs": 20}},
            {"step": "de", "params": {"method": "limma"}},
            {"step": "enrich", "params": {"db": ["GO_BP", "KEGG", "HALLMARK"], "method": "ora"}},
            {"step": "immune", "params": {}},
        ],
    },
    {
        "id": "geo",
        "title": "GEO 多数据集合并",
        "detail": "质控→过滤→标准化→去批次→降维→差异",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {"method": "auto"}},
            {"step": "batch", "params": {"method": "combat"}},
            {"step": "pca", "params": {"n_pcs": 20}},
            {"step": "de", "params": {"method": "limma"}},
        ],
    },
    {
        "id": "quick",
        "title": "快速差异",
        "detail": "质控→过滤→标准化→差异",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {"method": "auto"}},
            {"step": "de", "params": {"method": "limma"}},
        ],
    },
    {
        "id": "clinical",
        "title": "临床预后",
        "detail": "标准化→差异→生存→预后",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {"method": "auto"}},
            {"step": "de", "params": {"method": "limma"}},
            {"step": "survival", "params": {}},
            {"step": "prognosis", "params": {"model": "cox"}},
        ],
    },
]


def recipe(recipe_id: str | None = None):
    if not recipe_id:
        return {"recipes": RECIPES}
    for item in RECIPES:
        if item["id"] == recipe_id:
            return item
    raise EngineError("not_found", f"预设流程不存在：{recipe_id}")


# ---------------------------------------------------------------------------
# 导出
# ---------------------------------------------------------------------------

def export_dataset(dataset, what: str, out_path: str | None, workdir: str, engine_info: dict | None = None) -> dict:
    what = str(what or "summary_json")
    export_dir = ensure_dir(os.path.join(workdir, "exports"))
    stamp = time.strftime("%Y%m%d-%H%M%S")
    base = f"{dataset.id}-{what}-{stamp}"

    def target(ext: str) -> str:
        if out_path:
            return out_path
        return os.path.join(export_dir, f"{base}.{ext}")

    if what == "expr_csv":
        path = target("csv")
        frame = pd.DataFrame(dataset.X, index=dataset.samples, columns=dataset.genes)
        frame.to_csv(path, encoding="utf-8")
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what == "sample_info_csv":
        path = target("csv")
        frame = dataset.sample_info.copy()
        frame.insert(0, "sample", frame.index.astype(str))
        frame.to_csv(path, index=False, encoding="utf-8")
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what in ("de_csv", "enrich_csv", "immune_csv", "survival_csv"):
        key = {"de_csv": "de", "enrich_csv": "enrich", "immune_csv": "immune", "survival_csv": "survival"}[what]
        stored = dataset.results.get(key)
        if not stored:
            raise EngineError("prerequisite", f"还没有 {key} 结果，先运行对应步骤")
        path = target("csv")
        frame = _store_to_frame(stored)
        frame.to_csv(path, index=index_wanted(key), encoding="utf-8")
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what == "summary_json":
        path = target("json")
        payload = {
            "dataset": dataset.summary(history_limit=50),
            "steps": list(dataset.steps),
            "history": dataset.history,
            "analyses_run": dataset.analyses_run,
            "engine": jsonable(engine_info or {}),
            "results_available": sorted(dataset.results.keys()),
            "modules": list(dataset.modules.keys()),
            "exported_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(jsonable(payload), handle, ensure_ascii=False, indent=2)
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what == "protocol":
        path = target("json")
        records = _execution_records(dataset, engine_info)
        payload = {
            "plugin": "dsh-geobulk",
            "dataset": dataset.summary(history_limit=50),
            "source": {"path": dataset.source_path, "format": dataset.source_format},
            "engine": jsonable(engine_info or {}),
            "steps": [{"step": h.get("step"), "params": h.get("params"), "ok": h.get("ok"), "ms": h.get("ms")}
                      for h in dataset.history],
            "execution_record": jsonable({
                "generated_at": records["generated_at"],
                "data_version": records["data_version"],
                "dataset_version": records["dataset_version"],
                "dataset_revision": records["dataset_revision"],
                "steps": [
                    {
                        "id": record["id"],
                        "params": record["params"],
                        "requested_method": record["requested_method"],
                        "actual_method": record["actual_method"],
                        "design": record["design"],
                        "contrast": record["contrast"],
                        "fallback": record["fallback"],
                        "data_version": record["data_version"],
                        "notes": record["notes"],
                        "warnings": record["warnings"],
                    }
                    for record in records["steps"]
                ],
                "excluded": records["excluded"],
                "analyses": records["analyses"],
                "unrecorded": records["unrecorded"],
            }),
            "seed": 0,
            "exported_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(jsonable(payload), handle, ensure_ascii=False, indent=2)
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what == "methods":
        path = target("md")
        text = _methods_markdown(dataset, engine_info)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what == "report":
        path = target("md")
        text = _report_markdown(dataset, engine_info)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

    if what == "script":
        path = target("py")
        text, manifest = _script_python(dataset, engine_info, filename=os.path.basename(path))
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        manifest_path = os.path.splitext(path)[0] + ".manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as handle:
            json.dump(jsonable(manifest), handle, ensure_ascii=False, indent=2)
        label = manifest.get("label") or {}
        return {
            "file": path,
            "bytes": os.path.getsize(path),
            "what": what,
            "manifest": manifest_path,
            "manifest_bytes": os.path.getsize(manifest_path),
            "replay_status": str(label.get("status") or ""),
            "replayable": bool(label.get("replayable")),
            "label": str(label.get("text") or ""),
            "replay_reasons": list(label.get("blockers") or []) + list(label.get("notes") or []),
            "script_sha256": _file_sha256(path),
        }

    raise EngineError("bad_request", f"不支持的导出类型：{what}")


def index_wanted(key: str) -> bool:
    return True


def _store_to_frame(stored) -> pd.DataFrame:
    if isinstance(stored, pd.DataFrame):
        return stored
    if not isinstance(stored, dict):
        return pd.DataFrame()
    for key in ("table", "de_table"):
        frame = stored.get(key)
        if isinstance(frame, pd.DataFrame):
            return frame
        if isinstance(frame, list) and frame and isinstance(frame[0], dict):
            return pd.DataFrame(frame)
    # 免疫/生存：dict of 列
    candidate = {k: v for k, v in stored.items() if isinstance(v, (list, tuple, np.ndarray))}
    if candidate:
        try:
            return pd.DataFrame(candidate)
        except Exception:  # noqa: BLE001
            pass
    try:
        return pd.DataFrame([{k: v for k, v in stored.items() if not isinstance(v, (dict, list))}])
    except Exception:  # noqa: BLE001
        return pd.DataFrame()


def _engine_line(engine_info: dict | None) -> str:
    engine_info = engine_info or {}
    py = (engine_info.get("python") or "").strip()
    r = (engine_info.get("r") or {}).get("version") if isinstance(engine_info.get("r"), dict) else None
    parts = [f"dsh-geobulk 工作台（Python {py or 'unknown'}"]
    if r:
        parts.append(f"R {r}")
    return ", ".join(parts) + ")"


# ---------------------------------------------------------------------------
# 执行记录 → Methods（EXPORT-01）
#
# 原则：只描述 dataset.history / dataset.analyses_run / dataset.results 里真实记录到
# 的、成功执行的、且与当前数据/分组一致的步骤；取不到的字段写「未记录」，绝不按
# 按钮名或"通常做法"补写。
# ---------------------------------------------------------------------------

_NOT_RECORDED = "未记录"

_DE_METHOD_TEXT = {
    "limma": "limma moderated t-statistics (limma::lmFit + eBayes + topTable)",
    "deseq2": "DESeq2 negative binomial generalized linear model with Wald tests (DESeq2::results)",
    "welch": "Welch's two-sample t-test (unequal variances; Python implementation)",
}

_NORMALIZE_TEXT = {
    "cpm": "Counts were normalized to log2 counts per million (log2(CPM+1)), implemented in the workbench",
    "tmm": "Counts were normalized with the TMM procedure (log2 TMM-normalized CPM), implemented in the workbench",
    "quantile": "Expression values were quantile-normalized",
    "median": "Expression values were median-centered",
    "log": "Expression values were log2-transformed without further normalization",
    "none": "No normalization was applied",
}

_GENERIC_STEP_TEXT = {
    "scmap": "Bulk signatures were projected onto the single-cell reference recorded for this dataset",
    "ncrna": "Non-coding RNA / fusion / alternative-splicing modules were analysed on the attached module data",
}


def _fmt_value(value, digits: int = 4) -> str | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number in (float("inf"), float("-inf")):
        return None
    if abs(number - round(number)) < 1e-9 and abs(number) < 1e7:
        return f"{int(round(number))}"
    return f"{number:.{digits}g}"


def _r_version(engine: dict | None) -> str | None:
    info = (engine or {}).get("r")
    if isinstance(info, dict) and info.get("version"):
        return str(info["version"])
    return None


def _last_success(dataset, step: str) -> dict | None:
    for entry in reversed(list(getattr(dataset, "history", None) or [])):
        if str(entry.get("step")) == step and entry.get("ok"):
            return entry
    return None


def _result_stores(dataset) -> dict:
    """dataset.results 里的结构化结果（兼容 {key: {...}} 嵌套写法）。"""
    out: dict = {}
    results = getattr(dataset, "results", None) or {}
    if not isinstance(results, dict):
        return out
    for key in ("de", "enrich", "immune", "survival", "prognosis", "batch"):
        value = results.get(key)
        if not isinstance(value, dict) or not value:
            continue
        inner = value.get(key)
        if isinstance(inner, dict) and len(value) == 1:
            value = inner
        out[key] = value
    return out


def _de_consistency(dataset, store: dict | None = None) -> dict:
    """差异结果是否仍对应当前分组/样本（分组改过就必须重新判定）。

    store 为空时用 dataset.results["de"]；应优先传执行当时的 result_snapshot。
    """
    if not isinstance(store, dict) or not store:
        results = getattr(dataset, "results", None) or {}
        store = results.get("de") if isinstance(results, dict) else None
    if not isinstance(store, dict) or not store:
        return {"status": "missing", "detail": "没有 de 结果缓存"}
    group_col = str(store.get("group_col") or "")
    samples = [str(s) for s in (store.get("samples") or [])]
    labels = [str(v) for v in (store.get("labels") or [])]
    if not group_col or not samples or len(samples) != len(labels):
        return {"status": "incomplete", "detail": "de 记录缺少 group_col/samples/labels"}
    info = getattr(dataset, "sample_info", None)
    if info is None or group_col not in _column_names(info):
        return {"status": "stale", "detail": f"样本信息里已没有分组列「{group_col}」"}
    index = [str(i) for i in list(getattr(info, "index", None) if getattr(info, "index", None) is not None else [])]
    current = []
    for sample in samples:
        if sample not in index:
            return {"status": "stale", "detail": f"样本 {sample} 已不在当前数据集"}
        current.append(str(info.loc[sample, group_col]))
    if current != labels:
        return {"status": "stale", "detail": "group labels differ from the execution record"}
    if sorted(str(s) for s in (getattr(dataset, "samples", None) or [])) != sorted(samples):
        return {"status": "stale", "detail": "de 记录的样本集合与当前数据集不一致"}
    return {"status": "ok", "detail": None}


def _fallback_reason(record: dict) -> str | None:
    for item in list(record.get("warnings") or []) + list(record.get("notes") or []):
        text_value = str(item)
        if any(token in text_value for token in ("回退", "不可用", "unavailable", "fallback", "改换")):
            return text_value.strip()
    return None


# -- 各步骤的描述器（只写记录里有的内容） -----------------------------------

def _describe_qc(record: dict) -> None:
    metrics = record["metrics"]
    text = "Sample-level quality control was performed"
    count = _fmt_value(metrics.get("n_samples"))
    if count:
        text += f" on {count} samples"
    details = []
    library = _fmt_value(metrics.get("library_median"))
    if library:
        details.append(f"median library size {library}")
    detected = _fmt_value(metrics.get("detect_median"))
    if detected:
        details.append(f"median detected genes {detected}")
    outliers = _fmt_value(metrics.get("n_outliers"))
    if outliers is not None:
        details.append(f"{outliers} outlier sample(s) flagged")
    if details:
        text += " (" + "; ".join(details) + ")"
    text += "."
    if metrics.get("value_kind"):
        record["input_value_kind"] = str(metrics["value_kind"])
        text += f" Input matrix scale: `{metrics['value_kind']}`."
    else:
        record["unrecorded"].append("质控输入矩阵尺度")
    if not metrics:
        record["unrecorded"].append("质控指标")
    record["description"].append(text)


def _describe_filter(record: dict) -> None:
    metrics = record["metrics"]
    params = record["params"]
    text = "Low-expression genes were filtered"
    criteria = []
    min_expr = _fmt_value(params.get("min_expr"))
    if min_expr is not None:
        frac = params.get("min_samples_frac")
        try:
            pct = f"{float(frac) * 100:g}%"
        except (TypeError, ValueError):
            pct = None
        criteria.append(f"expression ≥ {min_expr}" + (f" in ≥ {pct} of samples" if pct else ""))
    if criteria:
        text += " (" + "; ".join(criteria) + ")"
    if metrics:
        before = _fmt_value(metrics.get("n_genes_before"))
        after = _fmt_value(metrics.get("n_genes_after"))
        removed = _fmt_value(metrics.get("removed_genes"))
        if before and after:
            text += f", retaining {after} of {before} genes"
            if removed is not None:
                text += f" ({removed} removed)"
        removed_samples = _fmt_value(metrics.get("removed_samples"))
        if removed_samples and removed_samples != "0":
            text += f"; {removed_samples} sample(s) removed"
    else:
        record["unrecorded"].append("过滤前后基因/样本数")
    record["description"].append(text + ".")


def _describe_normalize(record: dict) -> None:
    params = record["params"]
    metrics = record["metrics"]
    requested = str(params.get("method") or "") or None
    actual = str(metrics.get("method") or "") or None
    if not actual:
        match = re.search(r"标准化方法[：:]\s*([^\s（(]+)", str(record.get("summary") or ""))
        if match:
            actual = match.group(1)
    if not actual:
        record["description"].append(
            "A normalization step was executed, but the applied method is not recorded (未记录).")
        record["unrecorded"].append("标准化实际方法")
        return
    record["requested_method"] = requested or actual
    record["actual_method"] = actual
    log_transformed = metrics.get("log_transformed")
    before_kind = metrics.get("value_kind_before")
    after_kind = metrics.get("value_kind_after")
    if actual == "none":
        kind = after_kind or before_kind or _NOT_RECORDED
        text = f"No normalization was applied; the matrix remained on the original `{kind}` scale"
    else:
        text = _NORMALIZE_TEXT.get(actual) or f"Normalization was performed with method `{actual}`"
    detail = []
    if requested:
        detail.append(f"requested method: `{requested}`")
    log_request = params.get("log_transform")
    if log_request not in (None, ""):
        detail.append(f"requested log transform: `{log_request}`")
    if log_transformed is not None:
        detail.append(f"log2 applied: {'yes' if log_transformed else 'no'}")
    if detail:
        text += " (" + "; ".join(detail) + ")"
    text += "."
    if before_kind or after_kind:
        text += f" Matrix scale: `{before_kind or _NOT_RECORDED}` → `{after_kind or _NOT_RECORDED}`."
        record["input_value_kind"] = before_kind or record.get("input_value_kind")
        record["output_value_kind"] = after_kind or record.get("output_value_kind")
    else:
        record["unrecorded"].append("标准化前后矩阵尺度")
    if actual == "quantile":
        text += " Ties were averaged against the reference distribution."
    if actual == "none":
        text += (" Fold changes are reported on the original scale and must not be interpreted "
                 "as log2 fold changes.")
    record["description"].append(text)


def _describe_batch(record: dict) -> None:
    metrics = record["metrics"]
    params = record["params"]
    method = str(metrics.get("method") or params.get("method") or "")
    lowered = method.lower()
    if "combat" in lowered:
        text = "Batch effects were corrected with ComBat"
    elif "removebatcheffect" in lowered:
        text = "Batch effects were regressed out with limma::removeBatchEffect"
    elif method in ("none", ""):
        text = "No batch correction was applied"
    else:
        text = f"Batch correction was applied with method `{method}`"
    details = []
    batch_key = params.get("batch_key") or metrics.get("batch_key")
    if batch_key:
        details.append(f"batch column `{batch_key}`")
    n_batches = _fmt_value(metrics.get("n_batches"))
    if n_batches:
        details.append(f"{n_batches} batches")
    if details:
        text += " (" + ", ".join(details) + ")"
    record["description"].append(text + ".")
    record["actual_method"] = method or None


def _describe_pca(record: dict) -> None:
    params = record["params"]
    metrics = record["metrics"]
    n_pcs = _fmt_value(params.get("n_pcs") or metrics.get("n_pcs"))
    text = "Principal component analysis was performed for sample-level exploration"
    if n_pcs:
        text += f" ({n_pcs} components)"
    text += "; PCA results were used for exploration only and not for statistical inference."
    record["description"].append(text)


def _describe_de(record: dict, dataset, stores: dict) -> None:
    metrics = record["metrics"]
    params = record["params"]
    store = stores.get("de") or {}
    consistency = _de_consistency(dataset, store)
    if consistency.get("status") == "stale":
        record["excluded"] = True
        record["exclude_reason"] = (
            "the differential-expression result does not match the current sample groups "
            f"({consistency.get('detail')})"
        )
        return
    case = str(store.get("case") or metrics.get("case") or "")
    control = str(store.get("control") or metrics.get("control") or "")
    requested = metrics.get("method_requested") or store.get("method_requested")
    actual = metrics.get("actual_method") or metrics.get("method") or store.get("method")
    if not actual:
        record["description"].append(
            "A differential-expression step was executed, but the applied method is not recorded (未记录).")
        record["unrecorded"].append("差异分析实际方法")
        return
    record["requested_method"] = requested or actual
    record["actual_method"] = actual
    design = metrics.get("design")
    record["design"] = design
    record["contrast"] = metrics.get("contrast") or (
        f"{case}-{control}" if case and control else None)
    n_case = _fmt_value(store.get("n_case") or metrics.get("n_case"))
    n_ctrl = _fmt_value(store.get("n_ctrl") or metrics.get("n_control"))
    text = "Differential expression"
    if case and control:
        text += f" between `{case}`"
        if n_case:
            text += f" (n={n_case})"
        text += f" and `{control}`"
        if n_ctrl:
            text += f" (n={n_ctrl})"
    method_text = _DE_METHOD_TEXT.get(actual)
    r_version = _r_version(record.get("engine"))
    if method_text is None:
        method_text = f"method `{actual}`"
    elif actual in ("limma", "deseq2"):
        if r_version:
            method_text = method_text[:-1] + f", R {r_version})"
    text += f" was assessed with {method_text}."
    if case and control:
        text += f" Contrast: `{case} - {control}` (positive log fold change = higher expression in `{case}`)."
    else:
        record["unrecorded"].append("差异分析对比方向")
    if design:
        text += f" Design: `{design}`."
    else:
        record["unrecorded"].append("差异分析设计公式")
    paired = metrics.get("paired", (store.get("params") or {}).get("paired"))
    pair_col = metrics.get("pair_col") or (store.get("params") or {}).get("pair_col")
    if paired:
        text += f" Paired design by `{pair_col}`." if pair_col else " Paired design was recorded."
    covars = metrics.get("covars") or (store.get("params") or {}).get("covars") or []
    if covars:
        text += " Covariates included: " + ", ".join(f"`{c}`" for c in covars) + "."
    input_kind = record.get("input_value_kind")
    if input_kind:
        text += f" Input matrix scale: `{input_kind}`."
    else:
        record["unrecorded"].append("差异分析输入矩阵尺度")
    padj = _fmt_value(metrics.get("padj", params.get("padj")))
    lfc = _fmt_value(metrics.get("lfc", params.get("lfc")))
    if padj is not None or lfc is not None:
        parts = []
        if padj is not None:
            parts.append(f"FDR threshold: {padj}")
        if lfc is not None:
            parts.append(f"|logFC| threshold: {lfc}")
        text += " " + "; ".join(parts) + "."
    n_genes = _fmt_value(store.get("n_genes") or metrics.get("n_genes"))
    n_sig = _fmt_value(store.get("n_sig") or metrics.get("n_sig"))
    if n_genes:
        text += f" {n_genes} genes were tested"
        if n_sig is not None:
            up = _fmt_value(store.get("n_up") or metrics.get("n_up"))
            down = _fmt_value(store.get("n_down") or metrics.get("n_down"))
            text += f"; {n_sig} significant"
            if up is not None and down is not None:
                text += f" ({up} up / {down} down)"
        text += "."
    for note in record.get("notes") or []:
        if "log2(x+1)" in str(note):
            text += " The input matrix was log2(x+1)-transformed before testing (recorded note)."
            break
    lite = metrics.get("lite", store.get("lite"))
    if lite and requested and requested != actual:
        reason = _fallback_reason(record)
        record["fallback"] = {"requested": requested, "actual": actual, "reason": reason}
        text += (f" The requested method `{requested}` is recorded as unavailable; the recorded "
                 f"fallback is `{actual}`.")
        if reason:
            text += f" Reason recorded: {reason}."
        text += " Fallback results are exploratory (lite implementation)."
    record["description"].append(text)


def _describe_enrich(record: dict, stores: dict) -> None:
    metrics = record["metrics"]
    store = stores.get("enrich") or {}
    method = str(metrics.get("method") or store.get("method") or "")
    actual = metrics.get("actual_method")
    databases = metrics.get("databases") or store.get("databases") or []
    if isinstance(databases, str):
        databases = [databases]
    if method == "gsea" or actual:
        text = "Gene-set enrichment was assessed with a weighted Kolmogorov-Smirnov statistic"
        details = []
        if actual:
            details.append(f"recorded implementation: {actual}")
        software = metrics.get("software")
        version = metrics.get("software_version")
        if software:
            details.append(f"software: {software}" + (f" {version}" if version else ""))
        if metrics.get("permutation"):
            details.append(f"permutation scheme: {metrics['permutation']}")
        if metrics.get("gene_set_source"):
            details.append(f"gene sets: {metrics['gene_set_source']}")
        if details:
            text += " (" + "; ".join(details) + ")"
        text += "."
        record["requested_method"] = method or "gsea"
        record["actual_method"] = actual or method
    else:
        text = "Over-representation analysis was performed with hypergeometric tests"
        if databases:
            text += " against: " + ", ".join(str(item) for item in databases)
        text += "."
        record["requested_method"] = method or "ora"
        record["actual_method"] = metrics.get("actual_method") or "ora"
    n_terms = _fmt_value(store.get("n_terms") or metrics.get("n_terms"))
    n_sig = _fmt_value(store.get("n_significant") or metrics.get("n_significant"))
    if n_terms:
        text += f" {n_terms} terms were tested"
        if n_sig is not None:
            text += f"; {n_sig} significant"
        text += "."
    record["description"].append(text)


def _describe_immune(record: dict, stores: dict) -> None:
    metrics = record["metrics"]
    store = stores.get("immune") or {}
    methods_ok = metrics.get("methods_ok") or store.get("methods_ok") or []
    failed = metrics.get("methods_failed") or store.get("methods_failed") or []
    engine_label = metrics.get("engine") or store.get("engine")
    primary = metrics.get("primary_method") or store.get("primary_method")
    derived = metrics.get("derived", store.get("primary_derived"))
    if methods_ok:
        text = "Immune cell infiltration was estimated with " + ", ".join(str(m) for m in methods_ok)
        record["actual_method"] = ", ".join(str(m) for m in methods_ok)
    else:
        text = "An immune-infiltration step was executed, but the successful methods are not recorded (未记录)"
        record["unrecorded"].append("免疫浸润成功方法")
    if engine_label:
        text += f" ({engine_label})"
    text += "."
    if primary and derived:
        text += f" Cell proportions were derived from {primary} and are exploratory."
    if failed:
        text += " Requested but unavailable: " + ", ".join(str(m) for m in failed) + "."
    record["description"].append(text)


def _describe_survival(record: dict, stores: dict) -> None:
    metrics = record["metrics"]
    store = stores.get("survival") or {}
    time_col = metrics.get("time_col") or store.get("time_col")
    event_col = metrics.get("event_col") or store.get("event_col")
    gene = metrics.get("gene") or store.get("gene")
    cutoff = metrics.get("cutoff") or store.get("cutoff")
    logrank = metrics.get("logrank_p", store.get("logrank_p"))
    hr = metrics.get("HR")
    multi = metrics.get("cox_multi_cindex")
    cindex = metrics.get("cindex", store.get("cindex"))
    if not time_col and not event_col:
        record["description"].append(
            "A survival analysis was executed, but the time/event columns are not recorded (未记录).")
        record["unrecorded"].append("生存分析时间/事件列")
        return
    text = "Kaplan-Meier estimates and log-rank tests were used for survival analysis"
    details = []
    if time_col:
        details.append(f"time column `{time_col}`")
    if event_col:
        details.append(f"event column `{event_col}`")
    if details:
        text += " (" + ", ".join(details) + ")"
    text += "."
    if gene:
        text += f" Patients were split by `{gene}` (`{cutoff or 'median'}` cut-off)."
    logrank_value = _fmt_value(logrank, 3)
    if logrank_value is not None:
        text += f" The recorded log-rank P value was {logrank_value}."
    record["actual_method"] = "Kaplan-Meier + log-rank"
    if hr is not None:
        text += (" A univariable Cox proportional hazards model was fitted with "
                 f"`{gene}` as the predictor (hazard ratio and P value recorded).")
        record["actual_method"] += " + Cox"
    if multi is not None:
        text += " A multivariable Cox proportional hazards model was also fitted."
    cindex_value = _fmt_value(cindex)
    if cindex_value is not None:
        text += f" The concordance index was recorded ({cindex_value})."
    record["description"].append(text)


def _describe_prognosis(record: dict, stores: dict) -> None:
    metrics = record["metrics"]
    store = stores.get("prognosis") or stores.get("survival") or {}
    model = metrics.get("model") or metrics.get("method") or store.get("model")
    model_kind = metrics.get("model_kind") or store.get("model_type")
    n_genes = _fmt_value(metrics.get("n_genes") or (len(store.get("genes") or []) or None))
    scope = metrics.get("evaluation_scope")
    text = "A prognostic risk model was built"
    if model:
        text += f" with the recorded model `{model}`"
    if model_kind:
        text += f" (recorded model kind: `{model_kind}`)"
    if n_genes:
        text += f" from {n_genes} genes"
    text += "."
    if scope:
        text += f" Evaluation scope: {scope}."
    else:
        text += " The evaluation scope is not recorded (未记录)."
        record["unrecorded"].append("预后模型评估口径")
    record["description"].append(text)


def _describe_generic(record: dict, step: str) -> None:
    base = _GENERIC_STEP_TEXT.get(step, f"The `{step}` step was executed")
    record["description"].append(base + " (the detailed method record is unavailable: 未记录).")
    record["unrecorded"].append(f"{step} 方法细节")


def _new_step_record(dataset, step: str, hist: dict, engine_info: dict | None) -> dict:
    """只放原始记录字段，描述文本稍后（链式尺度补全之后）再生成。"""
    metrics = dict(hist.get("metrics") or {})
    return {
        "id": step,
        "at": hist.get("at"),
        "ms": hist.get("ms"),
        "params": dict(hist.get("params") or {}),
        "metrics": metrics,
        "notes": [str(x) for x in (hist.get("notes") or [])],
        "warnings": [str(x) for x in (hist.get("warnings") or [])],
        "summary": str(hist.get("summary") or ""),
        "engine": dict(hist.get("engine") or engine_info or {}),
        "data_version": dict(hist.get("data_version") or {}),
        "result_snapshot": dict(hist.get("result_snapshot") or {}),
        "description": [],
        "unrecorded": [],
        "requested_method": None,
        "actual_method": None,
        "input_value_kind": metrics.get("value_kind_before") or metrics.get("value_kind"),
        "output_value_kind": metrics.get("value_kind_after"),
        "design": None,
        "contrast": None,
        "fallback": None,
        "excluded": False,
        "exclude_reason": None,
    }


def _record_store_view(record: dict, fallback: dict) -> dict:
    """本步骤执行当时的结果（快照优先）——高级分析之后覆盖 dataset.results 也不影响它。"""
    view = dict(fallback or {})
    for key, value in (record.get("result_snapshot") or {}).items():
        if isinstance(value, dict) and value:
            view[str(key)] = value
    return view


def _describe_record(record: dict, dataset, stores: dict) -> None:
    step = record["id"]
    stores = _record_store_view(record, stores)
    if step == "qc":
        _describe_qc(record)
    elif step == "filter":
        _describe_filter(record)
    elif step == "normalize":
        _describe_normalize(record)
    elif step == "batch":
        _describe_batch(record)
    elif step == "pca":
        _describe_pca(record)
    elif step == "de":
        _describe_de(record, dataset, stores)
    elif step == "enrich":
        _describe_enrich(record, stores)
    elif step == "immune":
        _describe_immune(record, stores)
    elif step == "survival":
        _describe_survival(record, stores)
    elif step == "prognosis":
        _describe_prognosis(record, stores)
    else:
        _describe_generic(record, step)


def _chain_matrix_scale(records: list[dict]) -> None:
    """把矩阵尺度按执行顺序串起来（记录值链，不改变任何一步的实际结论）。"""
    current = None
    for record in records:
        if record.get("input_value_kind") is None and current:
            record["input_value_kind"] = current
        if record.get("output_value_kind") is None:
            # qc/filter/de 等步骤不改变矩阵尺度：沿用本步骤输入尺度
            record["output_value_kind"] = record.get("input_value_kind")
        if record.get("output_value_kind"):
            current = record["output_value_kind"]
    # de 的输入尺度：模块 notes 里会写「矩阵类型为 X」（比链式更直接）
    for record in records:
        if record["id"] == "de" and record.get("input_value_kind") is None:
            for note in record.get("notes") or []:
                match = re.search(r"矩阵类型为\s*([A-Za-z_]+)", str(note))
                if match:
                    record["input_value_kind"] = match.group(1)
                    break


def _execution_records(dataset, engine_info: dict | None = None) -> dict:
    """Methods / 导出脚本共用的结构化执行记录。"""
    stores = _result_stores(dataset)
    records: list[dict] = []
    excluded: list[dict] = []
    for step in [str(s) for s in (getattr(dataset, "steps", None) or [])]:
        hist = _last_success(dataset, step)
        if hist is None:
            excluded.append({"id": step, "reason": "没有成功的执行记录（history 中只有失败运行）"})
            continue
        records.append(_new_step_record(dataset, step, hist, engine_info))
    _chain_matrix_scale(records)
    kept: list[dict] = []
    for record in records:
        _describe_record(record, dataset, stores)
        if record.get("excluded"):
            excluded.append({"id": record["id"],
                             "reason": record.get("exclude_reason") or "记录不属于当前数据版本"})
            continue
        kept.append(record)
    records = kept
    de_record = next((record for record in records if record["id"] == "de"), None)
    de_store = (de_record or {}).get("result_snapshot", {}).get("de") if de_record else None
    analyses = []
    for entry in (getattr(dataset, "analyses_run", None) or []):
        if not entry.get("ok"):
            continue
        analyses.append({
            "id": str(entry.get("id")),
            "at": entry.get("at"),
            "summary": str(entry.get("summary") or ""),
            "params": dict(entry.get("params") or {}),
            "params_recorded": bool(entry.get("params")),
            "metrics": dict(entry.get("metrics") or {}),
        })
    unrecorded = sorted({str(item) for record in records for item in (record.get("unrecorded") or [])})
    result_steps = {"de", "enrich", "immune", "survival", "prognosis"}
    legacy = [record["id"] for record in records
              if record["id"] in result_steps and not (record.get("result_snapshot") or {})]
    if legacy:
        unrecorded.append("执行期结果快照（旧记录，使用当前结果缓存：" + ", ".join(legacy) + "）")
    return {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "data_version": _data_version(dataset),
        "dataset_version": getattr(dataset, "version", None),
        "dataset_revision": getattr(dataset, "revision", None),
        "steps": records,
        "excluded": excluded,
        "analyses": analyses,
        "stores": stores,
        "de_consistency": _de_consistency(dataset, de_store if isinstance(de_store, dict) else None),
        "unrecorded": unrecorded,
    }


def _statistics_lines(records: list[dict]) -> list[str]:
    de = next((record for record in records if record["id"] == "de"), None)
    if not de or not de.get("actual_method"):
        return ["- No differential-expression result is recorded for the current dataset version; "
                "no statistical test is claimed."]
    actual = str(de["actual_method"])
    method_text = _DE_METHOD_TEXT.get(actual, f"`{actual}`")
    metrics = de.get("metrics") or {}
    params = de.get("params") or {}
    padj = _fmt_value(metrics.get("padj", params.get("padj")))
    lfc = _fmt_value(metrics.get("lfc", params.get("lfc")))
    text = f"- Differential expression used {method_text}"
    text += "; P values were adjusted with the Benjamini-Hochberg procedure; tests were two-sided"
    thresholds = []
    if padj is not None:
        thresholds.append(f"FDR < {padj}")
    if lfc is not None:
        thresholds.append(f"|logFC| ≥ {lfc}")
    if thresholds:
        text += " and significance was defined as " + " and ".join(thresholds)
    lines = [text + "."]
    fallback = de.get("fallback")
    if isinstance(fallback, dict) and fallback.get("requested"):
        lines.append(f"- The requested method `{fallback['requested']}` was unavailable; the reported "
                     f"results come from `{fallback['actual']}` and are exploratory.")
    return lines


def _provenance_lines(dataset, records: dict, engine_info: dict | None) -> list[str]:
    dv = records.get("data_version") or {}
    lines = ["## Record provenance", ""]
    stamp = [
        f"matrix `{str(dv.get('matrix_sha256') or '')[:12]}`",
        f"genes `{str(dv.get('genes_sha256') or '')[:12]}`",
        f"samples `{str(dv.get('samples_sha256') or '')[:12]}`",
        f"value kind `{dv.get('value_kind') or 'unknown'}`",
    ]
    if dv.get("group_sha256"):
        stamp.append(f"group labels `{str(dv['group_sha256'])[:12]}`")
    else:
        stamp.append(f"group labels {_NOT_RECORDED}")
    lines.append(f"- Records were read from the dataset execution history at {records['generated_at']}; "
                 "data version stamp: " + ", ".join(stamp) + ".")
    version = records.get("dataset_version")
    revision = records.get("dataset_revision")
    if version is None and revision is None:
        lines.append(f"- Dataset version / revision: {_NOT_RECORDED}.")
    else:
        lines.append(f"- Dataset version / revision: version "
                     f"{version if version is not None else _NOT_RECORDED}, "
                     f"revision {revision if revision is not None else _NOT_RECORDED}.")
    described = [record["id"] for record in records.get("steps") or []]
    lines.append("- Steps described (successful, current data version): "
                 + (", ".join(described) if described else "none") + ".")
    excluded = records.get("excluded") or []
    if excluded:
        lines.append("- Excluded records: " + "; ".join(
            f"`{item['id']}`: {item['reason']}" for item in excluded) + ".")
    else:
        lines.append("- Excluded records: none.")
    unrecorded = list(records.get("unrecorded") or [])
    if unrecorded:
        lines.append(f"- Fields not recorded: {'; '.join(unrecorded)}.")
    else:
        lines.append("- Fields not recorded: none.")
    lines.append(f"- R package versions: {_NOT_RECORDED}（执行记录只保存 R 版本与包可用性）; "
                 "module-internal default parameters are applied identically on replay but are "
                 f"{_NOT_RECORDED} as explicit values.")
    if getattr(dataset, "source_path", None):
        lines.append(f"- Input source: `{os.path.basename(str(dataset.source_path))}` "
                     f"(format `{dataset.source_format}`); import options beyond the format "
                     f"(transpose / sheet / member / gene_col / delimiter) are {_NOT_RECORDED}.")
    else:
        lines.append(f"- Input source: {_NOT_RECORDED}（没有原始来源路径）.")
    return lines


def _methods_markdown(dataset, engine_info: dict | None = None) -> str:
    records = _execution_records(dataset, engine_info)
    lines = [
        "# Methods",
        "",
        f"Bulk transcriptome analysis was performed with the {_engine_line(engine_info)}.",
        f"The dataset contained {dataset.n_samples} samples and {dataset.n_genes} genes"
        + (f" (source: `{os.path.basename(dataset.source_path)}`)." if dataset.source_path else "."),
        "Every statement below is generated from the structured execution record of the current "
        "dataset version: steps without a successful, current record are not described, and "
        f"unavailable fields are marked {_NOT_RECORDED} rather than inferred.",
        "",
        "## Data processing",
        "",
    ]
    if not records["steps"]:
        lines.append("No analysis step was recorded for the current dataset version, so no "
                     "processing or statistical method is claimed.")
    for record in records["steps"]:
        description = " ".join(str(part).strip() for part in (record.get("description") or [])
                               if str(part).strip())
        if not description:
            description = f"executed (method description not recorded: {_NOT_RECORDED})."
        lines.append(f"- **{record['id']}**: {description}")
    lines += ["", "## Statistics", ""]
    lines.extend(_statistics_lines(records["steps"]))
    lines += ["", *_provenance_lines(dataset, records, engine_info), "",
              f"_Generated at {records['generated_at']}._"]
    return "\n".join(lines)


def _report_markdown(dataset, engine_info: dict | None = None) -> str:
    lines = [
        f"# Bulk 转录组分析报告：{dataset.name}",
        "",
        "## 数据概览",
        "",
        f"- 样本数：**{dataset.n_samples}**；基因数：**{dataset.n_genes}**",
        f"- 数据来源：`{dataset.source_path or '内置示例'}`（格式 {dataset.source_format}）",
        f"- 值类型：{dataset.value_kind}；物种：{dataset.species}",
        f"- 分组列：`{dataset.group_key or '未设置'}`；分组：{', '.join(dataset.groups or []) or '—'}",
        f"- 已完成步骤：{', '.join(dataset.steps) or '—'}",
        f"- 引擎：{_engine_line(engine_info)}",
        "",
        "## Methods",
        "",
        _methods_markdown(dataset, engine_info).split("\n", 1)[1] if dataset.steps else "(尚未运行分析步骤)",
        "",
        "## 表达矩阵（前 10 基因 × 前 10 样本）",
        "",
        _matrix_markdown(dataset),
        "",
        "## 差异表达（前 30）",
        "",
        _de_markdown(dataset),
        "",
        "## 功能富集（前 20）",
        "",
        _enrich_markdown(dataset),
        "",
        "## 免疫浸润（前 10 样本）",
        "",
        _immune_markdown(dataset),
        "",
        "## 分析记录",
        "",
        *[f"- {h.get('at')} `{h.get('step')}` {'✓' if h.get('ok') else '✗'} {h.get('summary') or ''}" for h in dataset.history[-20:]],
        "",
        f"_报告生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}_",
    ]
    return "\n".join(lines)


def _matrix_markdown(dataset) -> str:
    if not dataset.genes or not dataset.samples:
        return "—"
    genes = dataset.genes[:10]
    samples = dataset.samples[:10]
    idx = [dataset.genes.index(g) for g in genes]
    sidx = [dataset.samples.index(s) for s in samples]
    header = "| gene | " + " | ".join(samples) + " |"
    sep = "| --- | " + " | ".join(["---"] * len(samples)) + " |"
    rows = [header, sep]
    for gi, gene in zip(idx, genes):
        vals = [f"{float(dataset.X[si, gi]):.2f}" for si in sidx]
        rows.append(f"| {gene} | " + " | ".join(vals) + " |")
    return "\n".join(rows)


def _de_markdown(dataset) -> str:
    stored = dataset.results.get("de")
    if not stored:
        return "（尚未做差异表达）"
    frame = _store_to_frame(stored)
    if frame.empty:
        return "（无结果）"
    cols = [c for c in ("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val", "regulation") if c in frame.columns]
    frame = frame[cols].head(30)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = [header, sep]
    for _, row in frame.iterrows():
        rows.append("| " + " | ".join(_fmt_cell(row[c]) for c in cols) + " |")
    return "\n".join(rows)


def _enrich_markdown(dataset) -> str:
    stored = dataset.results.get("enrich")
    if not stored:
        return "（尚未做功能富集）"
    frame = _store_to_frame(stored)
    if frame.empty:
        return "（无显著条目）"
    cols = [c for c in ("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "Count") if c in frame.columns]
    if not cols:
        cols = list(frame.columns[:6])
    frame = frame[cols].head(20)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = [header, sep]
    for _, row in frame.iterrows():
        rows.append("| " + " | ".join(_fmt_cell(row[c]) for c in cols) + " |")
    return "\n".join(rows)


def _immune_markdown(dataset) -> str:
    stored = dataset.results.get("immune")
    if not stored:
        return "（尚未做免疫浸润）"
    frame = _store_to_frame(stored)
    if frame.empty:
        return "（无结果）"
    frame = frame.head(10)
    cols = list(frame.columns[:8])
    header = "| " + " | ".join([str(c) for c in cols]) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = [header, sep]
    for _, row in frame.iterrows():
        rows.append("| " + " | ".join(_fmt_cell(row[c]) for c in cols) + " |")
    return "\n".join(rows)


def _fmt_cell(value) -> str:
    try:
        if isinstance(value, float):
            if abs(value) < 1e-3 and value != 0:
                return f"{value:.2e}"
            return f"{value:.4g}"
        return str(value).replace("|", "\\|")
    except Exception:  # noqa: BLE001
        return str(value)


REPLAY_MANIFEST_SCHEMA = "dsh-geobulk/replay-manifest/1"
PLUGIN_ROOT = os.path.dirname(HERE)


def _file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _plugin_meta() -> dict:
    meta = {"name": "dsh-geobulk", "version": "unknown", "lib_dir": HERE}
    try:
        with open(os.path.join(PLUGIN_ROOT, "package.json"), "r", encoding="utf-8") as handle:
            meta["version"] = str(json.load(handle).get("version") or "unknown")
    except Exception:  # noqa: BLE001
        pass
    return meta


def _sample_info_csv(dataset) -> str:
    info = getattr(dataset, "sample_info", None)
    if info is None:
        return ""
    try:
        frame = info.copy()
        frame.index = [str(item) for item in frame.index]
        frame.index.name = str(frame.index.name or "sample")
        return frame.to_csv(sep="\t")
    except Exception:  # noqa: BLE001
        return ""


def _expected_de_rows(store: dict, limit: int = 20) -> list[dict]:
    if not isinstance(store, dict):
        return []
    try:
        frame = _store_to_frame(store)
    except Exception:  # noqa: BLE001
        return []
    if frame is None or frame.empty or "gene" not in frame.columns:
        return []
    if "P.Value" in frame.columns:
        try:
            frame = frame.sort_values("P.Value")
        except Exception:  # noqa: BLE001
            pass
    rows = []
    for _, row in frame.head(limit).iterrows():
        item = {"gene": str(row.get("gene"))}
        for column in ("logFC", "P.Value", "adj.P.Val"):
            if column in frame.columns:
                try:
                    item[column] = float(row[column])
                except (TypeError, ValueError):
                    continue
        rows.append(item)
    return rows


def _step_snapshot_store(record: dict | None, fallback: dict, key: str) -> dict:
    """优先用该步骤执行当时的结果快照（避免被后续高级分析覆盖）。"""
    snapshot = (record or {}).get("result_snapshot") or {}
    value = snapshot.get(key)
    if isinstance(value, dict) and value:
        return value
    value = (fallback or {}).get(key)
    return value if isinstance(value, dict) else {}


def _expected_results(records: dict) -> dict:
    stores = records["stores"]
    steps = records["steps"]
    expected: dict = {}
    de_record = next((record for record in steps if record["id"] == "de"), None)
    store = _step_snapshot_store(de_record, stores, "de")
    if de_record or store:
        metrics = (de_record or {}).get("metrics") or {}
        params = (store.get("params") or {})
        payload = {
            "method": metrics.get("actual_method") or metrics.get("method") or store.get("method"),
            "method_requested": metrics.get("method_requested") or store.get("method_requested"),
            "lite": bool(metrics.get("lite", store.get("lite"))),
            "case": store.get("case") or metrics.get("case"),
            "control": store.get("control") or metrics.get("control"),
            "contrast": (de_record or {}).get("contrast"),
            "design": metrics.get("design"),
            "paired": bool(params.get("paired")),
            "pair_col": params.get("pair_col"),
            "covars": list(params.get("covars") or []),
            "n_up": store.get("n_up"),
            "n_down": store.get("n_down"),
            "n_sig": store.get("n_sig"),
            "rows": store.get("table_head") or _expected_de_rows(store),
        }
        if any(value is not None for value in payload.values()):
            expected["de"] = payload
    enrich_record = next((record for record in steps if record["id"] == "enrich"), None)
    enrich_store = _step_snapshot_store(enrich_record, stores, "enrich")
    if enrich_store or enrich_record:
        metrics = (enrich_record or {}).get("metrics") or {}
        expected["enrich"] = {
            "method": metrics.get("method") or enrich_store.get("method"),
            "actual_method": metrics.get("actual_method"),
            "databases": enrich_store.get("databases"),
            "n_terms": enrich_store.get("n_terms"),
            "n_significant": enrich_store.get("n_significant"),
        }
    immune_record = next((record for record in steps if record["id"] == "immune"), None)
    immune_store = _step_snapshot_store(immune_record, stores, "immune")
    if immune_store or immune_record:
        metrics = (immune_record or {}).get("metrics") or {}
        expected["immune"] = {
            "methods_ok": metrics.get("methods_ok") or immune_store.get("methods_ok"),
            "primary_method": metrics.get("primary_method") or immune_store.get("primary_method"),
            "engine": metrics.get("engine") or immune_store.get("engine"),
        }
    survival_record = next((record for record in steps if record["id"] == "survival"), None)
    survival_store = _step_snapshot_store(survival_record, stores, "survival")
    if survival_store or survival_record:
        metrics = (survival_record or {}).get("metrics") or {}
        expected["survival"] = {
            "time_col": metrics.get("time_col") or survival_store.get("time_col"),
            "event_col": metrics.get("event_col") or survival_store.get("event_col"),
            "gene": metrics.get("gene") or survival_store.get("gene"),
            "logrank_p": metrics.get("logrank_p", survival_store.get("logrank_p")),
            "HR": metrics.get("HR"),
            "cindex": metrics.get("cindex", survival_store.get("cindex")),
        }
    return jsonable(expected)


def _replay_plan(dataset, engine_info: dict | None = None) -> dict:
    """判断导出物能否从同一输入重放，并返回 manifest（EXPORT-02）。"""
    records = _execution_records(dataset, engine_info)
    blockers: list[str] = []
    notes: list[str] = []
    source_path = getattr(dataset, "source_path", None)
    source_text = str(source_path) if source_path else ""
    is_url = bool(re.match(r"^(?:https?|ftp)://", source_text, re.I))
    source_sha = None
    source_bytes = None
    if not source_text:
        blockers.append("数据集没有原始来源路径（demo/挂载数据），无法从同一输入重放")
    elif is_url:
        blockers.append("输入来自 URL，导出物没有保存可校验的本地副本")
    elif not os.path.exists(source_text):
        blockers.append(f"原输入文件不存在：{source_text}")
    else:
        try:
            source_sha = _file_sha256(source_text)
            source_bytes = os.path.getsize(source_text)
        except OSError as error:
            blockers.append(f"原输入文件无法读取：{error}")
    consistency = records.get("de_consistency") or {}
    if consistency.get("status") == "stale":
        blockers.append("差异分析记录与当前分组不一致（记录已过期，不能当作当前结论）")
    modules = sorted((getattr(dataset, "modules", None) or {}).keys())
    if modules:
        blockers.append("已挂载外部模块数据（{}），重放未覆盖".format(", ".join(modules)))
    module_steps = [record["id"] for record in records["steps"] if record["id"] in ("scmap", "ncrna")]
    if module_steps:
        blockers.append("步骤 {} 依赖外部挂载数据/单细胞参考，重放未覆盖".format(", ".join(module_steps)))
    if getattr(dataset, "sc_ref", None):
        blockers.append("单细胞参考（sc_ref）未纳入重放")
    analyses = records.get("analyses") or []
    if analyses:
        notes.append("{} 个高级分析未纳入重放：{}".format(
            len(analyses), ", ".join(item["id"] for item in analyses)))
    if not records["steps"]:
        notes.append("没有已执行步骤：重放只做导入与终态校验")
    import_fingerprint = None
    if records["steps"]:
        import_fingerprint = dict(records["steps"][0].get("data_version") or {})
    if not import_fingerprint:
        notes.append("导入边界指纹未记录（旧记录）：重放只核对终态与关键结果，不校验导入边界")
    final = _data_version(dataset)
    info_csv = _sample_info_csv(dataset)
    operations = []
    for index, record in enumerate(records["steps"]):
        engine = record.get("engine") or {}
        operations.append({
            "seq": index + 1,
            "kind": "step",
            "id": record["id"],
            "params": jsonable(record.get("params") or {}),
            "summary": record.get("summary") or "",
            "at": record.get("at"),
            "ms": record.get("ms"),
            "input_version": jsonable(record.get("data_version") or {}),
            "engine": {
                "python": (engine.get("python") if isinstance(engine, dict) else None),
                "r": _r_version(engine),
            },
        })
    if blockers:
        status, label_text, replayable = "record", "记录（不可独立复现）", False
    elif analyses:
        status, label_text, replayable = "partial", "部分记录（步骤可重放；高级分析未重放）", False
    else:
        status, label_text, replayable = "replayable", "可复现脚本（replay）", True
    engine_info = engine_info or {}
    manifest = {
        "schema": REPLAY_MANIFEST_SCHEMA,
        "generated_at": records["generated_at"],
        "plugin": _plugin_meta(),
        "label": {
            "status": status,
            "text": label_text,
            "replayable": replayable,
            "blockers": list(blockers),
            "notes": list(notes),
        },
        "dataset": {
            "id": getattr(dataset, "id", None),
            "name": getattr(dataset, "name", None),
            "n_samples": int(getattr(dataset, "n_samples", 0) or 0),
            "n_genes": int(getattr(dataset, "n_genes", 0) or 0),
            "species": getattr(dataset, "species", None),
            "value_kind": getattr(dataset, "value_kind", None),
            "group_key": getattr(dataset, "group_key", None),
            "groups": getattr(dataset, "groups", None),
            "version": records.get("dataset_version"),
            "revision": records.get("dataset_revision"),
            "data_version": final,
        },
        "input": {
            "source_path": source_text or None,
            "source_format": getattr(dataset, "source_format", None),
            "exists": bool(source_text) and os.path.exists(source_text) if source_text else False,
            "sha256": source_sha,
            "bytes": source_bytes,
            "sample_info_csv": info_csv,
            "sample_info_sha256": _sha256_text(info_csv) if info_csv else None,
            "import_fingerprint": jsonable(import_fingerprint),
            "genes_sha256": final["genes_sha256"],
            "samples_sha256": final["samples_sha256"],
            "genes": [str(g) for g in (getattr(dataset, "genes", None) or [])],
            "samples": [str(s) for s in (getattr(dataset, "samples", None) or [])],
            "format_options_not_recorded": [
                "transpose", "sheet", "member", "gene_col", "delimiter", "object",
            ],
        },
        "seed": 0,
        "operations": operations,
        "analyses": jsonable(records.get("analyses") or []),
        "expected": dict(
            {
                "n_samples": final["n_samples"],
                "n_genes": final["n_genes"],
                "value_kind": final["value_kind"],
                "matrix_sha256": final["matrix_sha256"],
                "genes_sha256": final["genes_sha256"],
                "samples_sha256": final["samples_sha256"],
                "group_key": final["group_key"],
                "group_sha256": final["group_sha256"],
                "steps": [record["id"] for record in records["steps"]],
                "step_fingerprints": [
                    {"id": record["id"], "input": jsonable(record.get("data_version") or {})}
                    for record in records["steps"]
                ],
            },
            **_expected_results(records),
        ),
        "versions": {
            "python": engine_info.get("python"),
            "platform": engine_info.get("platform"),
            "packages": jsonable(engine_info.get("packages") or {}),
            "r": jsonable(engine_info.get("r") or {}),
        },
        "replay": {
            "scope": "import -> steps (dataset.steps order) -> verification",
            "tolerance": 1e-8,
            "not_replayed": list(blockers) + list(notes),
        },
    }
    return jsonable(manifest)


def _compact_manifest(manifest: dict) -> dict:
    """内嵌到脚本里的清单：去掉超长的基因/样本全量列表，保留校验值与头部清单。"""
    compact = json.loads(json.dumps(jsonable(manifest), ensure_ascii=False))
    inputs = compact.get("input") or {}
    genes = inputs.pop("genes", []) or []
    samples = inputs.pop("samples", []) or []
    inputs["genes_head"] = genes[:200]
    inputs["samples_head"] = samples[:200]
    compact["input"] = inputs
    return compact


def _embed_manifest(manifest: dict) -> str:
    payload = json.dumps(manifest, ensure_ascii=False)
    if '"""' in payload:
        payload = json.dumps(manifest, ensure_ascii=True)
    return payload


_REPLAY_SCRIPT = r'''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
''' + r'''"""
dsh-geobulk 导出物 — __LABEL_TEXT__

生成时间：__GENERATED_AT__（dsh-geobulk __TOOL_VERSION__）
数据集：__DATASET_NAME__（样本 __N_SAMPLES__ / 基因 __N_GENES__，矩阵尺度 __VALUE_KIND__）
重放范围：从同一输入导入 → 按记录顺序重跑 __OPERATION_COUNT__ 个步骤 → 与执行记录逐项核对
不可重放项：__BLOCKER_TEXT__
未纳入重放：__NOTE_TEXT__

本脚本不是通用示例模板：它按该数据集的执行记录重放同一流程（同一导入入口、同一模块、
同一方法参数），并在结束时逐项核对样本/基因/矩阵尺度/对比方向/实际方法/关键结果。
需要：python3 + numpy/pandas + 生成结果时的插件源码目录（manifest.plugin.lib_dir）+ R 环境。
默认在脚本所在目录下写 replay_out/，核对报告为 replay_out/replay_report.json。

用法：
    python3 __SCRIPT_FILENAME__
    python3 __SCRIPT_FILENAME__ --source /path/to/input --lib <plugin>/lib --workdir ./out
    python3 __SCRIPT_FILENAME__ --options-json '{"transpose": true}'   # 导入选项未记录时手工指定
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
import time

MANIFEST = json.loads(r"""__MANIFEST_JSON__""")
SIDECAR_NAME = "__SIDECAR_NAME__"
CORE_KEYS = (
    "n_samples", "n_genes", "value_kind", "matrix_sha256",
    "genes_sha256", "samples_sha256", "group_key", "group_sha256",
)
IMPORT_KEYS = (
    "n_samples", "n_genes", "value_kind", "matrix_sha256", "genes_sha256", "samples_sha256",
)


def _sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _column_names(frame):
    # pandas Index 的 in 判断要显式转 list（Index 的布尔值有歧义）
    columns = getattr(frame, "columns", None)
    if columns is None:
        return []
    try:
        return [str(item) for item in list(columns)]
    except TypeError:
        return []


def _fingerprint(dataset, numpy):
    X = numpy.ascontiguousarray(numpy.asarray(getattr(dataset, "X", numpy.zeros((0, 0))), dtype=numpy.float32))
    genes = [str(g) for g in (getattr(dataset, "genes", None) or [])]
    samples = [str(s) for s in (getattr(dataset, "samples", None) or [])]
    group_key = getattr(dataset, "group_key", None)
    group_sha256 = None
    info = getattr(dataset, "sample_info", None)
    if group_key and info is not None and group_key in _column_names(info):
        labels = [str(v) for v in info[group_key].astype(str).tolist()]
        group_sha256 = _sha256_bytes("\n".join(labels).encode("utf-8"))
    return {
        "n_samples": int(X.shape[0]) if X.ndim == 2 else 0,
        "n_genes": int(X.shape[1]) if X.ndim == 2 else 0,
        "value_kind": str(getattr(dataset, "value_kind", "unknown") or "unknown"),
        "matrix_sha256": _sha256_bytes(X.data if X.size else b""),
        "genes_sha256": _sha256_bytes("\n".join(genes).encode("utf-8")),
        "samples_sha256": _sha256_bytes("\n".join(samples).encode("utf-8")),
        "group_key": None if not group_key else str(group_key),
        "group_sha256": group_sha256,
    }


def _load_manifest():
    here = os.path.dirname(os.path.abspath(__file__))
    side = os.path.join(here, SIDECAR_NAME)
    if os.path.exists(side):
        try:
            with open(side, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            if isinstance(data, dict) and data.get("schema") == MANIFEST.get("schema"):
                return data, side
        except Exception as error:  # noqa: BLE001
            print("[replay] WARN manifest 读取失败，改用内嵌副本：%s" % (error,))
    return MANIFEST, "embedded"


def _check(checks, name, ok, detail=""):
    checks.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    suffix = "" if ok or not detail else "  " + str(detail)
    print("[replay] %s %s%s" % ("OK  " if ok else "FAIL", name, suffix))
    return bool(ok)


def _compare_fingerprint(checks, name, recorded, actual, keys):
    recorded = recorded if isinstance(recorded, dict) else {}
    bad = []
    for key in keys:
        if recorded.get(key) != actual.get(key):
            bad.append("%s: recorded=%r replay=%r" % (key, recorded.get(key), actual.get(key)))
    return _check(checks, name, not bad, "; ".join(bad[:4]))


def _close(left, right, tol=1e-8):
    try:
        a = float(left)
        b = float(right)
    except (TypeError, ValueError):
        return False
    if a != a or b != b:
        return (a != a) and (b != b)
    return abs(a - b) <= tol * max(1.0, abs(b))


def _verify_de(checks, dataset, expected_de, tolerance):
    results = getattr(dataset, "results", None) or {}
    store = results.get("de") if isinstance(results, dict) else None
    if not isinstance(store, dict) or not store:
        _check(checks, "de result", False, "重放后没有 de 结果缓存")
        return None
    settled = _check(
        checks, "de method/contrast",
        all(store.get(key) == expected_de.get(key)
            for key in ("method", "method_requested", "lite", "case", "control", "n_up", "n_down", "n_sig")
            if key in expected_de),
        "recorded=%s replay=%s" % (
            {key: expected_de.get(key) for key in ("method", "case", "control")},
            {key: store.get(key) for key in ("method", "case", "control")},
        ),
    )
    frame = store.get("table")
    index = None
    if frame is not None and "gene" in _column_names(frame):
        try:
            index = frame.set_index("gene")
        except Exception:  # noqa: BLE001
            index = None
    rows = expected_de.get("rows") or []
    actual_rows = []
    missing = []
    worst = 0.0
    for row in rows:
        gene = str(row.get("gene"))
        if index is None or gene not in index.index:
            missing.append(gene)
            continue
        actual = {"gene": gene}
        for column in ("logFC", "P.Value", "adj.P.Val"):
            if column not in row:
                continue
            try:
                value = float(index.loc[gene, column])
            except (TypeError, ValueError):
                continue
            actual[column] = value
            reference = float(row[column])
            worst = max(worst, abs(value - reference) / max(1.0, abs(reference)))
        actual_rows.append(actual)
    rows_ok = bool(rows) and not missing and worst <= float(tolerance)
    _check(checks, "de rows (%d recorded)" % len(rows), rows_ok,
           "missing=%s max_rel=%.2e" % (missing[:3], worst))
    return {
        "ok": bool(settled and rows_ok),
        "method": store.get("method"),
        "method_requested": store.get("method_requested"),
        "lite": store.get("lite"),
        "case": store.get("case"),
        "control": store.get("control"),
        "n_up": store.get("n_up"),
        "n_down": store.get("n_down"),
        "n_sig": store.get("n_sig"),
        "rows": actual_rows,
    }


def _replay(manifest, args):
    started_at = time.time()
    checks = []
    label = manifest.get("label") or {}
    plugin = manifest.get("plugin") or {}
    lib_dir = args.lib or plugin.get("lib_dir") or ""
    if not lib_dir or not os.path.isdir(lib_dir):
        _check(checks, "plugin lib_dir", False,
               "找不到插件源码目录：%r（用 --lib 指定）" % (lib_dir,))
        return 1, checks, None
    if lib_dir not in sys.path:
        sys.path.insert(0, lib_dir)
    import numpy as np
    import pandas as pd
    import bg_io
    import bg_registry
    import bg_state

    workdir = os.path.abspath(args.workdir or os.path.join(os.getcwd(), "replay_out"))
    os.makedirs(workdir, exist_ok=True)
    inputs = manifest.get("input") or {}
    expected = manifest.get("expected") or {}
    operations = [op for op in (manifest.get("operations") or []) if op.get("kind") == "step"]

    source = args.source or inputs.get("source_path") or ""
    if not source:
        _check(checks, "input source", False, "没有记录输入路径（可用 --source 指定）")
        return 1, checks, None
    if not re_match_url(source):
        source = os.path.abspath(source)
    if not os.path.exists(source):
        _check(checks, "input exists", False, "输入文件不存在：%s" % (source,))
        return 1, checks, None
    _check(checks, "input exists", True, source)
    recorded_sha = inputs.get("sha256")
    if recorded_sha:
        actual_sha = _sha256_file(source)
        if not _check(checks, "input sha256", actual_sha == recorded_sha,
                      "sha256 mismatch: recorded=%s actual=%s" % (recorded_sha, actual_sha)):
            return 1, checks, None
    else:
        _check(checks, "input sha256", False, "manifest 未记录输入 sha256")

    info_csv = inputs.get("sample_info_csv") or ""
    sample_info = None
    if info_csv:
        expected_info_sha = inputs.get("sample_info_sha256")
        _check(checks, "sample_info sha256",
               not expected_info_sha or expected_info_sha == _sha256_bytes(info_csv.encode("utf-8")),
               "manifest 内记录的样本表校验值不一致")
        sample_info = pd.read_csv(io.StringIO(info_csv), sep="\t", index_col=0)
        sample_info.index = [str(item) for item in sample_info.index]
        sample_info.index.name = "sample"
    else:
        _check(checks, "sample_info", False, "manifest 未保存样本信息表（分组/协变量无法恢复）")

    recorded_format = str(inputs.get("source_format") or "")
    options = {}
    if recorded_format and recorded_format not in ("unknown", "auto", "None"):
        options["data_format"] = recorded_format
    if args.options_json:
        try:
            options.update(json.loads(args.options_json))
        except Exception as error:  # noqa: BLE001
            _check(checks, "options_json", False, str(error))
            return 1, checks, None
    loaded = bg_io.load_source(source, workdir=workdir, options=options)
    value_kind = loaded.get("value_kind") or (loaded.get("meta") or {}).get("value_kind") or "unknown"
    genes = [str(g) for g in loaded.get("genes") or []]
    samples = [str(s) for s in loaded.get("samples") or []]
    X = np.asarray(loaded["X"])
    counts = loaded.get("counts")
    if counts is None and str(value_kind) == "counts":
        counts = X
    if sample_info is not None:
        sample_info = sample_info.reindex(pd.Index(samples))
    meta_dataset = manifest.get("dataset") or {}
    dataset = bg_state.Dataset(
        "replay", str(meta_dataset.get("name") or "replay"), X=X, genes=genes,
        samples=samples, counts=counts, sample_info=sample_info, source_path=source,
        source_format=str(inputs.get("source_format") or "unknown"),
        species=str(meta_dataset.get("species") or "unknown"),
        value_kind=str(value_kind),
    )
    group_key = meta_dataset.get("group_key")
    if group_key and group_key in dataset.sample_info.columns:
        dataset.group_key = str(group_key)
    _check(checks, "import", True,
           "%d samples x %d genes, value_kind=%s" % (dataset.n_samples, dataset.n_genes, dataset.value_kind))

    import_fingerprint = operations[0].get("input_version") if operations else None
    if isinstance(import_fingerprint, dict) and import_fingerprint:
        _compare_fingerprint(checks, "import fingerprint", import_fingerprint,
                             _fingerprint(dataset, np), IMPORT_KEYS)
    else:
        _check(checks, "import fingerprint", True, "没有导入边界记录，跳过（只核对终态）")

    engine = None
    try:
        import bg_engine
        engine = bg_engine.Worker(workdir).engine_info()
    except Exception:  # noqa: BLE001
        engine = None

    for index, operation in enumerate(operations, start=1):
        step_id = str(operation.get("id") or "")
        params = dict(operation.get("params") or {})
        step_started = time.time()
        try:
            bg_registry.run_step(dataset, step_id, params, workdir, engine)
        except Exception as error:  # noqa: BLE001
            _check(checks, "step %d %s" % (index, step_id), False,
                   "%s: %s" % (type(error).__name__, error))
            return 1, checks, None
        _check(checks, "step %d %s" % (index, step_id), True, "%.1fs" % (time.time() - step_started))

    fingerprint = _fingerprint(dataset, np)
    _compare_fingerprint(checks, "final state", expected, fingerprint, CORE_KEYS)
    if expected.get("steps") is not None:
        _check(checks, "step list", list(expected.get("steps") or []) == list(dataset.steps),
               "recorded=%s replay=%s" % (expected.get("steps"), list(dataset.steps)))

    tolerance = (manifest.get("replay") or {}).get("tolerance") or 1e-8
    de_report = None
    expected_de = expected.get("de") or None
    if isinstance(expected_de, dict) and expected_de:
        de_report = _verify_de(checks, dataset, expected_de, tolerance)

    status = "verified" if all(item["ok"] for item in checks) else "failed"
    report = {
        "status": status,
        "label": label.get("text"),
        "manifest_source": manifest.get("_manifest_source"),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "workdir": workdir,
        "checks": checks,
        "dataset": {
            "n_samples": int(dataset.n_samples),
            "n_genes": int(dataset.n_genes),
            "value_kind": str(dataset.value_kind),
            "matrix_sha256": fingerprint["matrix_sha256"],
            "group_key": dataset.group_key,
            "steps": list(dataset.steps),
        },
        "de": de_report,
        "recorded": {
            "steps": expected.get("steps"),
            "de": {key: (expected_de or {}).get(key) for key in ("method", "case", "control")}
            if isinstance(expected_de, dict) else None,
        },
        "elapsed_s": round(time.time() - started_at, 2),
    }
    report_path = os.path.join(workdir, "replay_report.json")
    with open(report_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print("[replay] report: %s" % (report_path,))
    if status != "verified":
        print("REPLAY STATUS: failed")
        return 1, checks, report
    if (label or {}).get("status") == "partial":
        print("[replay] 注意：%s" % (label.get("text"),))
        print("REPLAY STATUS: steps-verified")
        return 0, checks, report
    print("REPLAY STATUS: verified")
    return 0, checks, report


def re_match_url(text):
    lowered = str(text or "").lower()
    return lowered.startswith(("http://", "https://", "ftp://"))


def main():
    parser = argparse.ArgumentParser(description="dsh-geobulk 重放脚本（按执行记录恢复同一分析流程）")
    parser.add_argument("--source", default=None, help="输入文件路径（默认用 manifest 记录值）")
    parser.add_argument("--lib", default=None, help="插件 lib 目录（默认用 manifest 记录值）")
    parser.add_argument("--workdir", default=None, help="工作目录（默认 ./replay_out）")
    parser.add_argument("--options-json", default=None, help="额外的 bg_io.load_source 导入选项（JSON）")
    args = parser.parse_args()
    manifest, manifest_source = _load_manifest()
    manifest["_manifest_source"] = manifest_source
    label = manifest.get("label") or {}
    print("[replay] 产物标签：%s" % (label.get("text") or "未标注",))
    print("[replay] manifest：%s" % (manifest_source,))
    for item in label.get("blockers") or []:
        print("[replay] 不可重放项：%s" % (item,))
    for item in label.get("notes") or []:
        print("[replay] 未纳入重放：%s" % (item,))
    try:
        code, _checks, _report = _replay(manifest, args)
    except Exception as error:  # noqa: BLE001
        print("[replay] ERROR %s: %s" % (type(error).__name__, error))
        print("REPLAY STATUS: failed")
        return 1
    return code


if __name__ == "__main__":
    sys.exit(main())
'''


def _script_python(dataset, engine_info: dict | None = None,
                   filename: str | None = None) -> tuple[str, dict]:
    """导出「可重放脚本」+ manifest（EXPORT-02）。

    脚本按执行记录调用同一导入入口（bg_io.load_source 与 bg_registry.run_step）重跑
    同一流程，并逐项核对；同目录同时写出 <name>.manifest.json（输入清单与校验值、
    样本/基因映射、按顺序的操作与实际参数、依赖版本、随机种子、期望终态）。
    """
    manifest = _replay_plan(dataset, engine_info)
    label = manifest.get("label") or {}
    blocked = label.get("blockers") or []
    notes = label.get("notes") or []
    tool = manifest.get("plugin") or {}
    script_name = str(filename or f"{getattr(dataset, 'id', 'dataset')}-replay.py")
    compact = _compact_manifest(manifest)
    text = _REPLAY_SCRIPT
    replacements = {
        "__LABEL_TEXT__": str(label.get("text") or "记录（不可独立复现）"),
        "__GENERATED_AT__": str(manifest.get("generated_at") or ""),
        "__TOOL_VERSION__": str(tool.get("version") or "unknown"),
        "__DATASET_NAME__": str((manifest.get("dataset") or {}).get("name") or ""),
        "__N_SAMPLES__": str((manifest.get("dataset") or {}).get("n_samples") or 0),
        "__N_GENES__": str((manifest.get("dataset") or {}).get("n_genes") or 0),
        "__VALUE_KIND__": str((manifest.get("dataset") or {}).get("value_kind") or "unknown"),
        "__OPERATION_COUNT__": str(len(manifest.get("operations") or [])),
        "__BLOCKER_TEXT__": "；".join(blocked) if blocked else "无",
        "__NOTE_TEXT__": "；".join(notes) if notes else "无",
        "__SCRIPT_FILENAME__": script_name,
        "__SIDECAR_NAME__": os.path.splitext(script_name)[0] + ".manifest.json",
        "__MANIFEST_JSON__": _embed_manifest(compact),
    }
    for token, value in replacements.items():
        if token == "__MANIFEST_JSON__":
            text = text.replace(token, str(value))
            continue
        # 头部注入的是自由文本：去掉换行/三引号，避免破坏生成脚本的 docstring
        safe = str(value).replace("\n", " ").replace("\r", " ").replace('\"\"\"', "'''")
        text = text.replace(token, safe)
    return text, manifest
