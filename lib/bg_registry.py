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

import importlib
import json
import os
import pkgutil
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
    try:
        payload = step["run"](ctx, params, outdir)
    except EngineError as error:
        dataset.add_history(step_id, params, False, int((time.time() - t0) * 1000), summary=error.message, engine=engine_info)
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
    try:
        payload = analysis["run"](ctx, params, outdir)
    except EngineError as error:
        dataset.add_analysis_run(analysis_id, False, int((time.time() - started) * 1000), error.message)
        raise
    elapsed = int((time.time() - started) * 1000)
    clean = _apply_result(dataset, payload)
    bullets = (clean.get("summary") or {}).get("bullets") or []
    dataset.add_analysis_run(analysis_id, True, elapsed, str(bullets[0]) if bullets else "")
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
        payload = {
            "plugin": "dsh-geobulk",
            "dataset": dataset.summary(history_limit=50),
            "source": {"path": dataset.source_path, "format": dataset.source_format},
            "engine": jsonable(engine_info or {}),
            "steps": [{"step": h.get("step"), "params": h.get("params"), "ok": h.get("ok"), "ms": h.get("ms")}
                      for h in dataset.history],
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
        text = _script_python(dataset)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        return {"file": path, "bytes": os.path.getsize(path), "what": what}

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


def _methods_markdown(dataset, engine_info: dict | None = None) -> str:
    lines = [
        "# Methods",
        "",
        f"Bulk transcriptome analysis was performed with the {_engine_line(engine_info)}.",
        f"The dataset contained {dataset.n_samples} samples and {dataset.n_genes} genes"
        + (f" (source: `{os.path.basename(dataset.source_path)}`)." if dataset.source_path else "."),
        "",
        "## Data processing",
    ]
    step_text = {
        "qc": "Sample-level quality control assessed library size, detected genes and outlier samples.",
        "filter": "Low-expression genes were removed.",
        "normalize": "Expression values were log-transformed and quantile-normalized.",
        "batch": "Batch effects were corrected when multiple batches were present.",
        "pca": "Principal component analysis was used for sample-level dimensionality reduction.",
        "de": "Differential expression analysis was performed.",
        "enrich": "Functional enrichment analysis was performed on the differential gene list.",
        "immune": "Immune cell infiltration was estimated with deconvolution and signature-based methods.",
        "survival": "Kaplan-Meier curves and Cox proportional hazards models were used for survival analysis.",
        "prognosis": "A risk model was built and evaluated with time-dependent ROC and C-index.",
        "scmap": "Bulk signatures were projected onto a single-cell reference dataset to identify expressing cell subsets.",
        "ncrna": "Non-coding RNA, fusion gene and alternative splicing analyses were performed.",
    }
    for step in dataset.steps:
        if step in step_text:
            params = next((h.get("params") for h in reversed(dataset.history) if h.get("step") == step), {}) or {}
            lines.append(f"- **{step}**: {step_text[step]}" + (f" Parameters: `{json.dumps(jsonable(params), ensure_ascii=False)}`." if params else ""))
    lines += [
        "",
        "## Statistics",
        "Differential expression used moderated t-statistics with Benjamini-Hochberg FDR control. "
        "Enrichment used hypergeometric tests (ORA) or weighted Kolmogorov-Smirnov statistics (GSEA) against public gene sets. "
        "Survival comparisons used the log-rank test and Cox proportional hazards models. "
        "All tests were two-sided and P values below 0.05 were considered significant unless otherwise stated.",
        "",
        f"_Generated at {time.strftime('%Y-%m-%d %H:%M:%S')}._",
    ]
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


def _script_python(dataset) -> str:
    steps = list(dataset.steps)
    lines = [
        "# 由 dsh-geobulk 生成的分析脚本（复现记录）",
        "# 依赖：pandas / numpy / scipy / lifelines；R 官方实现（limma/clusterProfiler/IOBR）请见 methods。",
        "import numpy as np",
        "import pandas as pd",
        "",
        f"# 数据集：{dataset.name}；样本 {dataset.n_samples}；基因 {dataset.n_genes}",
        f"# 来源：{dataset.source_path or 'demo'}",
        f"# 分组列：{dataset.group_key}；分组：{dataset.groups}",
        "",
        "# 1) 读取表达矩阵（此处用 CSV 示例，按实际来源替换）",
        "# expr = pd.read_csv('expr.csv', index_col=0)",
        "",
        "# 2) 已执行的步骤：",
    ]
    for step in steps:
        params = next((h.get("params") for h in reversed(dataset.history) if h.get("step") == step), {}) or {}
        lines.append(f"#    - {step}: {json.dumps(jsonable(params), ensure_ascii=False)}")
    lines += [
        "",
        "if __name__ == '__main__':",
        f"    print('dsh-geobulk 复现脚本：{dataset.name}')",
    ]
    return "\n".join(lines)
