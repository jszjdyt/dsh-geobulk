"""dsh-geobulk — 非编码 RNA / 融合基因 / 可变剪接（bulk-ncrna 负责）

步骤 `ncrna`：按 `kind` 分发到 5 个模态，数据来自 `ctx["modules"][kind]`
（矩阵型给 `X`(样本×特征)+`genes`(特征名)+`samples`；表格型给 `table`），
也可用 `module_path` 直接读一个表格文件（pandas 读，不依赖 bg_io）。

模态实现（CONTRACT §8.5）：
  * circRNA  —— ID 坐标解析（chr1:12345-12789:+ / chr1_12345_12789_+ / hsa_circ_0001234）；
                Welch t + BH 差异（共享 bg_common.welch_de 统计核心）；linear 列 → back-spliced/total
                比例；matplotlib 极坐标圈图。
  * miRNA    —— 名称规范化（hsa-mir-21 → hsa-miR-21-5p，5p/3p 臂、let-7 家族）；差异分析；
                内置「常见 miRNA → 靶基因」简化表（≥40 个家族）+ miRNA-mRNA 负相关 → 靶基因表
                （诚实标注 lite：不做全基因组靶点预测/种子匹配）；ceRNA 三元组候选。
  * lncRNA   —— lncrna_class（LINC… / …-AS1 / …-IT1 / ENSG / 内置已知名）；差异分析；
                与 mRNA 共表达网络（|r|>0.6 & FDR<0.05，networkx spring 布局）。
  * fusion   —— 列识别（STAR-Fusion `#FusionName`/`LeftGene`、Arriba `gene1`+`gene2`、通用
                `fusion`/`gene1`+`gene2`）→ 标准化表；top 融合对条形图 + 样本×融合矩阵；
                内置癌基因（≥120）+ 激酶（≥90）过滤，标注 in-frame / known。
  * splicing —— rMATS（IJC_SAMPLE_1/SJC_SAMPLE_1/IncLevel1/2，含重复 ID 列）→ dPSI + Welch t + BH；
                SUPPA/通用 PSI 矩阵 → dPSI + Welch t + BH；火山图 / 事件类型计数 / top 事件热图。

图：一律 `bg_common.new_fig` / `save_fig`（白底、英文文字、dpi≈140）。
约定：只依赖标准库 + numpy/pandas/scipy（函数内惰性 import）+ networkx（惰性 import）
      + bg_common（契约 §3.3）；只读 ctx，不原地修改输入。
"""

from __future__ import annotations

import gzip
import math
import os
import re
import tempfile

import numpy as np
import pandas as pd

import bg_common as bc
from bg_common import EngineError, export_csv, palette, result, save_fig, table

__all__ = [
    "STEPS",
    "ANALYSES",
    "run_ncrna",
    "parse_circ_id",
    "normalize_mirna",
    "classify_lncrna",
]


KINDS = ("circ", "mirna", "lncrna", "fusion", "splicing")
KIND_LABELS = {
    "circ": "circRNA",
    "mirna": "miRNA",
    "lncrna": "lncRNA",
    "fusion": "融合基因",
    "splicing": "可变剪接",
}
KIND_TITLES = {
    "circ": "circRNA 分析",
    "mirna": "miRNA 分析",
    "lncrna": "lncRNA 分析",
    "fusion": "融合基因分析",
    "splicing": "可变剪接分析",
}

# 图内一律英文（契约 §2.5）；中文只出现在表格/摘要/notes
_PADJ_DEFAULT = 0.05
_LFC_DEFAULT = 0.0


# ---------------------------------------------------------------------------
# 通用：outdir / 模块取用 / 文件读取
# ---------------------------------------------------------------------------

def _outdir(ctx, outdir) -> str:
    path = outdir or (ctx or {}).get("outdir")
    if not path:
        path = tempfile.mkdtemp(prefix="bg-ncrna-")
    os.makedirs(path, exist_ok=True)
    return path


def _params(params) -> dict:
    return dict(params or {})


def _module(ctx, kind: str):
    """取 ctx["modules"][kind]（不存在/结构不对 → None）。"""
    modules = (ctx or {}).get("modules") or {}
    mod = modules.get(kind) if isinstance(modules, dict) else None
    return mod if isinstance(mod, dict) else None


def _module_ready(ctx, kind: str) -> bool:
    mod = _module(ctx, kind)
    if not mod:
        return False
    return isinstance(mod.get("X"), (np.ndarray, list)) or isinstance(mod.get("table"), pd.DataFrame)


def _missing_module_error(kind: str, ctx=None) -> EngineError:
    mod = _module(ctx, kind)
    hint = ""
    if mod:
        hint = "（模块已挂载但数据为空：请确认 attach 的文件里含表达矩阵/表格）"
    return EngineError(
        "prerequisite",
        f"还没有可用的 {KIND_LABELS[kind]} 数据：请先在「非编码」tab 挂载（attach module={kind}），"
        f"或用参数 module_path 指定表格文件{hint}",
    )


def _read_head(path, nbytes: int = 2 * 1024 * 1024) -> str:
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read(nbytes)
    except OSError as exc:
        raise EngineError("not_found", f"文件读取失败：{path}（{exc}）") from exc


def _sniff_delimiter(head: str) -> str:
    line = next((ln for ln in str(head).splitlines() if ln.strip()), "")
    counts = {
        "\t": line.count("\t"),
        ",": line.count(","),
        ";": line.count(";"),
        "|": line.count("|"),
    }
    best = max(counts, key=counts.get)
    return best if counts[best] > 0 else "whitespace"


def _read_table_file(path, warnings=None) -> pd.DataFrame:
    """用 pandas 读模块数据表格（csv/tsv/txt/gz/xlsx），不依赖 bg_io。"""
    raw = str(path or "").strip()
    if not raw:
        raise EngineError("bad_request", "module_path 为空")
    full = os.path.abspath(os.path.expanduser(raw))
    if not os.path.exists(full):
        raise EngineError("not_found", f"模块数据文件不存在：{full}")
    if os.path.isdir(full):
        raise EngineError("bad_request", f"module_path 需要数据文件而不是目录：{full}")
    name = os.path.basename(full)
    lower = full.lower()
    if lower.endswith((".xlsx", ".xls", ".xlsm")):
        try:
            frame = pd.read_excel(full)
        except Exception as exc:  # noqa: BLE001
            raise EngineError("bad_request", f"Excel 读取失败（{name}）：{exc}") from exc
    else:
        head = _read_head(full)
        if not str(head).strip():
            raise EngineError("bad_request", f"文件是空的：{name}")
        delimiter = _sniff_delimiter(head)
        try:
            if delimiter == "whitespace":
                frame = pd.read_csv(full, sep=r"\s+", engine="python", low_memory=False)
            else:
                frame = pd.read_csv(full, sep=delimiter, low_memory=False)
        except Exception as exc:  # noqa: BLE001
            raise EngineError("bad_request", f"表格读取失败（{name}）：{exc}") from exc
    if frame is None or frame.empty:
        raise EngineError("bad_request", f"表格没有数据行：{name}")
    if warnings is not None:
        warnings.append(f"已用 module_path 读取表格：{name}（{len(frame)} 行 × {len(frame.columns)} 列）")
    return frame


def _is_numeric_text(value) -> bool:
    try:
        float(str(value).replace(",", "").strip())
        return True
    except (TypeError, ValueError):
        return False


def _dedupe(names, prefix: str = "feature") -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for item in names:
        base = str(item).strip() or prefix
        count = seen.get(base, 0)
        if count == 0:
            out.append(base)
        else:
            out.append(f"{base}#{count + 1}")
        seen[base] = count + 1
    return out


_ID_CANDIDATES = {
    "circ": ["circ_id", "circrna", "circ", "junction", "id", "feature", "gene"],
    "mirna": ["mirna", "mirna_id", "mature_id", "name", "id", "feature", "gene"],
    "lncrna": ["lncrna", "lncrna_id", "gene", "gene_id", "id", "name", "feature"],
    "splicing": ["event_id", "id", "event", "feature", "name", "gene"],
}

# circRNA 表里常见的注释数值列（不是样本列）
_CIRC_ANNOT_NUMERIC = ("linear", "linear_expr", "linear_counts", "linear_expression", "circ_total", "total")


def _matrix_from_frame(frame: pd.DataFrame, kind: str, warnings=None) -> dict:
    """宽表 → (X 样本×特征, genes, samples, info)。也接受长表（feature/sample/value）。"""
    warnings = warnings if warnings is not None else []
    frame = frame.copy()
    frame.columns = [str(c) for c in frame.columns]
    lower = {str(c).strip().lower(): c for c in frame.columns}

    # 0) 长表透视为宽表
    feat_col = next((lower[k] for k in ("feature", "circ_id", "circrna", "mirna", "lncrna",
                                        "gene", "event_id", "event", "junction", "id") if k in lower), None)
    sample_col = next((lower[k] for k in ("sample", "sample_id", "samples", "gsm") if k in lower), None)
    value_col = next((lower[k] for k in ("value", "count", "counts", "expression", "expr",
                                         "tpm", "rpkm", "psi") if k in lower), None)
    if feat_col and sample_col and value_col:
        try:
            wide = frame.pivot_table(index=feat_col, columns=sample_col, values=value_col, aggfunc="mean")
            frame = wide.reset_index()
            lower = {str(c).strip().lower(): c for c in frame.columns}
            warnings.append("检测到长表（feature/sample/value），已透视为 特征×样本 宽表")
        except Exception:  # noqa: BLE001
            pass

    # 1) 特征 ID 列
    id_col = None
    for candidate in _ID_CANDIDATES.get(kind, ["id", "feature", "gene"]):
        if candidate in lower:
            id_col = lower[candidate]
            break
    if id_col is None:
        first = frame.columns[0]
        ratio = sum(1 for v in frame[first].astype(str) if _is_numeric_text(v)) / max(1, len(frame))
        if ratio < 0.5:
            id_col = first
    if id_col is None:
        raise EngineError("bad_request", f"表格里没有特征 ID 列（{KIND_LABELS[kind]} 名称/ID 列）")

    # 2) 样本列 / 注释列
    numeric_cols: list[str] = []
    annot_cols: list[str] = []
    linear = None
    for col in frame.columns:
        if col == id_col:
            continue
        series = frame[col]
        ratio = sum(1 for v in series.astype(str) if _is_numeric_text(v)) / max(1, len(series))
        is_numeric = ratio >= 0.8
        name_low = str(col).strip().lower()
        if kind == "circ" and is_numeric and name_low in _CIRC_ANNOT_NUMERIC:
            linear = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
            continue
        if is_numeric:
            numeric_cols.append(col)
        else:
            annot_cols.append(col)
    if not numeric_cols:
        raise EngineError("bad_request", "表格里没有数值表达列（样本列），无法做矩阵分析")

    genes = _dedupe([str(v) for v in frame[id_col]])
    samples = _dedupe([str(c) for c in numeric_cols])
    values = frame[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
    info = frame[annot_cols].copy().reset_index(drop=True)
    if linear is not None and len(linear) == len(frame):
        info["linear"] = linear
    X = values.T.astype(np.float64)
    if X.shape[0] != len(samples) or X.shape[1] != len(genes):
        raise EngineError("internal", "矩阵解析后维度与特征/样本数不一致")
    info.index = pd.Index(genes, name="feature")
    return {"X": X, "genes": genes, "samples": samples, "info": info}


def _matrix_inputs(ctx, params, kind: str, *, warnings=None, notes=None) -> dict:
    """取某矩阵型模态的输入（优先 module_path，其次 ctx["modules"][kind]）。"""
    params = _params(params)
    warnings = list(warnings or [])
    notes = list(notes or [])
    mod = _module(ctx, kind)
    path = str(params.get("module_path") or "").strip()
    X = genes = samples = info = None
    source = ""
    if path:
        frame = _read_table_file(path, warnings)
        parsed = _matrix_from_frame(frame, kind, warnings)
        X, genes, samples, info = parsed["X"], parsed["genes"], parsed["samples"], parsed["info"]
        source = os.path.basename(path)
        notes.append(f"模块数据来自 module_path：{path}")
    elif mod:
        raw_x = mod.get("X")
        if isinstance(raw_x, (np.ndarray, list)) and len(raw_x):
            X = np.asarray(raw_x, dtype=np.float64)
            genes = [str(g) for g in (mod.get("genes") or [])]
            samples = [str(s) for s in (mod.get("samples") or [])]
            raw_info = mod.get("info")
            info = raw_info if isinstance(raw_info, pd.DataFrame) else None
        elif isinstance(mod.get("table"), pd.DataFrame):
            parsed = _matrix_from_frame(mod["table"], kind, warnings)
            X, genes, samples, info = parsed["X"], parsed["genes"], parsed["samples"], parsed["info"]
            notes.append("已挂载模块只有表格，按宽表解析出表达矩阵")
        source = str(mod.get("path") or "已挂载模块")
        notes.extend([str(n) for n in (mod.get("notes") or [])][:3])
        warnings.extend([str(w) for w in (mod.get("warnings") or [])][:3])
    else:
        raise _missing_module_error(kind, ctx)

    if X is None:
        raise EngineError("prerequisite", f"{KIND_LABELS[kind]} 模块里没有表达矩阵")
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise EngineError("bad_request", f"{KIND_LABELS[kind]} 表达矩阵维度不正确（应为 样本×特征）")
    genes = [str(g) for g in (genes or [])]
    samples = [str(s) for s in (samples or [])]
    if len(genes) != X.shape[1] and len(genes) == X.shape[0] and len(samples) == X.shape[1]:
        X = X.T
        warnings.append("表达矩阵为 特征×样本 方向，已转置为 样本×特征")
    if len(genes) != X.shape[1]:
        if not genes:
            genes = [f"F{i + 1}" for i in range(X.shape[1])]
            warnings.append("模块没有特征名，已用 F1… 编号")
        else:
            raise EngineError("bad_request", f"表达矩阵列数（{X.shape[1]}）与特征数（{len(genes)}）不一致")
    if not samples:
        samples = [f"S{i + 1}" for i in range(X.shape[0])]
        warnings.append("模块没有样本名，已用 S1… 编号")
    if len(samples) != X.shape[0]:
        raise EngineError("bad_request", f"表达矩阵行数（{X.shape[0]}）与样本数（{len(samples)}）不一致")
    if info is None or len(info) != len(genes):
        info = pd.DataFrame(index=pd.Index(genes, name="feature"))
    else:
        info = info.copy()
        info.index = pd.Index(genes, name="feature")

    groups = _resolve_groups(ctx, mod, samples, params)
    if groups:
        notes.append(f"分组来自 {groups['source']}：{', '.join(groups['groups'])}")
    return {
        "kind": kind,
        "X": X,
        "genes": genes,
        "samples": samples,
        "info": info,
        "linear": _linear_column(info, genes),
        "groups": groups,
        "source": source,
        "warnings": warnings,
        "notes": notes,
    }


def _linear_column(info: pd.DataFrame, genes: list[str]):
    """circRNA 模块里的 linear（线性转录本）列：返回 (values[n_features], column_name) 或 (None, None)。"""
    if not isinstance(info, pd.DataFrame) or info.empty:
        return None, None
    for col in info.columns:
        name_low = str(col).strip().lower()
        if name_low in _CIRC_ANNOT_NUMERIC:
            series = pd.to_numeric(info[col], errors="coerce")
            if len(series) == len(genes) and series.notna().any():
                return series.to_numpy(dtype=np.float64), str(col)
    return None, None


# ---------------------------------------------------------------------------
# 通用：分组解析 / 差异统计
# ---------------------------------------------------------------------------

_GROUP_COL_PREFERRED = (
    "group", "groups", "condition", "group_name", "groupname", "class", "type",
    "status", "phenotype", "arm", "cohort", "subtype", "disease", "treatment",
)
_CASE_HINTS = (
    "case", "tumor", "tumour", "cancer", "disease", "treatment", "treated", "treat",
    "resist", "resistant", "high", "mut", "mutant", "ko", "knockout", "oe",
    "overexpress", "drug", "td", "experiment",
)
_CONTROL_HINTS = (
    "control", "ctrl", "normal", "healthy", "vehicle", "wt", "wildtype",
    "untreated", "nc", "negative", "low", "sensitive", "adjacent",
)


def _group_column(info: pd.DataFrame, requested=None):
    """在样本信息表里挑分组列：优先显式指定，其次常见命名，最后第一个 2–10 取值的非数值列。"""
    if not isinstance(info, pd.DataFrame) or len(info.columns) == 0:
        return None
    columns = [str(c) for c in info.columns]
    if requested:
        for col in columns:
            if col == str(requested):
                return col
        for col in columns:
            if col.lower() == str(requested).lower():
                return col
        return None
    lower = {col.lower(): col for col in columns}
    for name in _GROUP_COL_PREFERRED:
        if name in lower:
            return lower[name]
    total = len(info)
    for col in columns:
        series = info[col]
        if pd.api.types.is_numeric_dtype(series):
            continue
        n_unique = int(series.nunique(dropna=True))
        if 2 <= n_unique <= min(10, max(2, total - 1)):
            return col
    return None


def _resolve_groups(ctx, mod, samples, params):
    """给模块样本求分组标签：params.sample_groups → module.sample_info → ctx.sample_info → rMATS 命名。"""
    params = _params(params)
    samples = [str(s) for s in (samples or [])]
    if not samples:
        return None
    mapping: dict[str, str] = {}
    source = ""

    explicit = params.get("sample_groups")
    if isinstance(explicit, dict) and explicit:
        mapping = {str(k): str(v) for k, v in explicit.items()}
        source = "params.sample_groups"
    elif isinstance(explicit, (list, tuple)) and len(explicit) == len(samples):
        mapping = dict(zip(samples, [str(v) for v in explicit]))
        source = "params.sample_groups"

    if not mapping and isinstance(mod, dict):
        sinfo = mod.get("sample_info")
        if isinstance(sinfo, pd.DataFrame) and len(sinfo.columns):
            gcol = _group_column(sinfo, "group")
            if gcol:
                mapping = {str(k): str(v) for k, v in sinfo[gcol].items()}
                source = f"module.sample_info.{gcol}"

    if not mapping:
        info = (ctx or {}).get("sample_info")
        if isinstance(info, pd.DataFrame) and len(info.columns):
            wanted = params.get("group_col") or (ctx or {}).get("group_key")
            gcol = _group_column(info, wanted)
            if gcol:
                lookup = {str(k): str(v) for k, v in info[gcol].items()}
                lower_lookup = {k.lower(): k for k in lookup}
                for sample in samples:
                    key = sample if sample in lookup else lower_lookup.get(sample.lower())
                    if key is not None:
                        mapping[sample] = str(lookup[key])
                source = f"sample_info.{gcol}"
                covered = sum(1 for s in samples if s in mapping)
                if covered < max(2, int(0.6 * len(samples))):
                    mapping = {}

    if not mapping:
        if all(re.match(r"^group[12][_\-]?\d*$", str(s), re.I) for s in samples):
            mapping = {
                s: ("group1" if re.match(r"^group1", str(s), re.I) else "group2")
                for s in samples
            }
            source = "rMATS IncLevel 组列"

    if not mapping:
        return None
    labels = [mapping.get(s, "") for s in samples]
    if not all(labels):
        covered = sum(1 for value in labels if value)
        if covered < max(2, int(0.6 * len(samples))):
            return None
        labels = [value or "unknown" for value in labels]
    groups = [g for g in dict.fromkeys(labels) if g and g.strip() and g.lower() not in ("nan", "none", "na")]
    if len(groups) < 2:
        return None
    return {"labels": labels, "groups": groups, "source": source}


def _match_group(value, groups):
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    if text in groups:
        return text
    for group in groups:
        if group.lower() == text.lower():
            return group
    raise EngineError("bad_request", f"分组「{text}」不在样本分组里（可选：{', '.join(groups)}）")


def _case_control(groups, params):
    """实验组/对照组：显式参数 → 常见命名启发 → groups 首尾。"""
    params = _params(params)
    groups = [str(g) for g in groups]
    if len(groups) < 2:
        raise EngineError("prerequisite", f"分组数不足 2（当前：{', '.join(groups)}），无法做组间差异检验")
    case = _match_group(params.get("case"), groups)
    control = _match_group(params.get("control"), groups)
    lower = {g.lower(): g for g in groups}
    if not case or not control:
        guess_case = next((lower[h] for h in _CASE_HINTS if h in lower), None)
        guess_ctrl = next((lower[h] for h in _CONTROL_HINTS if h in lower), None)
        if guess_case and guess_ctrl and guess_case != guess_ctrl:
            case = case or guess_case
            control = control or guess_ctrl
        elif guess_case:
            case = case or guess_case
            control = control or next((g for g in groups if g != case), groups[0])
        elif guess_ctrl:
            control = control or guess_ctrl
            case = case or next((g for g in groups if g != control), groups[-1])
        else:
            control = control or groups[0]
            case = case or next((g for g in reversed(groups) if g != control), groups[-1])
    if case == control:
        raise EngineError("bad_request", f"实验组与对照组相同（{case}），请在参数里指定 case/control")
    return case, control


def _log2_if_needed(X, warnings, label: str) -> tuple[np.ndarray, bool]:
    arr = np.asarray(X, dtype=np.float64)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0 or finite.min() < 0:
        return arr, False
    integer_like = bool(np.allclose(finite, np.round(finite)))
    if integer_like or finite.max() > 50:
        warnings.append(f"{label} 表达量按 log2(x+1) 转换后做差异检验（logFC 为 log2 尺度差异）")
        return np.log2(np.clip(arr, 0, None) + 1.0), True
    return arr, False


def _diff_frame(inputs, params, *, label: str, warnings, notes):
    """组间 Welch t + BH（共享 bg_common.welch_de 统计核心）。"""
    groups = inputs.get("groups")
    if not groups:
        raise EngineError(
            "prerequisite",
            f"{label} 模块样本与样本表对不上，无法分组做差异分析："
            "请在「样本」里设置分组列，或给 ncrna 显式传 sample_groups / case / control",
        )
    case, control = _case_control(groups["groups"], params)
    labels = groups["labels"]
    Xv, transformed = _log2_if_needed(inputs["X"], warnings, label)
    if transformed:
        notes.append(f"{label} 差异分析在 log2 尺度上进行（Welch t + BH）")
    padj_cut = float(params.get("padj") or _PADJ_DEFAULT)
    lfc_cut = float(params.get("lfc") or _LFC_DEFAULT)
    de = bc.welch_de(
        Xv,
        labels,
        case,
        control,
        genes=inputs["genes"],
        lfc=lfc_cut,
        padj_cut=padj_cut,
    )
    return de, case, control


_DE_COLUMNS = [
    {"key": "gene", "label": "特征", "type": "string"},
    {"key": "logFC", "label": "logFC", "type": "number"},
    {"key": "AveExpr", "label": "平均表达", "type": "number"},
    {"key": "t", "label": "t 统计量", "type": "number"},
    {"key": "P.Value", "label": "P 值", "type": "number"},
    {"key": "adj.P.Val", "label": "FDR", "type": "number"},
    {"key": "regulation", "label": "方向", "type": "string"},
]


def _de_table_obj(frame: pd.DataFrame, tid: str, title: str, note: str | None = None,
                  extra_columns=None, lead_columns=None) -> dict:
    columns = list(lead_columns or []) + list(_DE_COLUMNS)
    if extra_columns:
        existing = {c["key"] for c in columns}
        for item in extra_columns:
            if item["key"] not in existing:
                columns.append(item)
    return table(tid, title, columns, frame, note=note)


def _de_summary(de: pd.DataFrame, label: str, case: str, control: str, padj_cut: float) -> tuple[str, dict]:
    sig = de[de["regulation"].isin(["up", "down"])]
    n_up = int((de["regulation"] == "up").sum())
    n_down = int((de["regulation"] == "down").sum())
    bullet = (
        f"{label}：{case} vs {control}，FDR<{padj_cut:g} 的差异特征 {len(sig)} 个"
        f"（上调 {n_up} / 下调 {n_down}）"
    )
    metrics = {
        "n_features": int(len(de)),
        "n_de": int(len(sig)),
        "n_up": n_up,
        "n_down": n_down,
    }
    return bullet, metrics


# ---------------------------------------------------------------------------
# 通用：绘图（白底、英文文字）
# ---------------------------------------------------------------------------

_CJK_RE = re.compile(r"[\u2e80-\u9fff\u3000-\u303f\uff00-\uffef]")


def _ascii_text(text, fallback: str = "") -> str:
    value = str(text)
    return value if not _CJK_RE.search(value) else fallback


def _fig_volcano(frame: pd.DataFrame, outdir: str, name: str, title: str, *,
                 xcol: str, pcol: str, labelcol: str, padj_cut: float = 0.05,
                 x_cut: float = 0.0, x_label: str = "log2 fold change",
                 y_label: str = "-log10 FDR", top_label: int = 8,
                 up_label: str = "Up", down_label: str = "Down"):
    x = pd.to_numeric(frame[xcol], errors="coerce")
    p = pd.to_numeric(frame[pcol], errors="coerce")
    x = x.fillna(0.0)
    p = p.fillna(1.0).clip(lower=1e-300)
    y = -np.log10(p)
    significant = p <= padj_cut
    if x_cut > 0:
        sig_up = significant & (x >= x_cut)
        sig_down = significant & (x <= -x_cut)
    else:
        sig_up = significant & (x > 0)
        sig_down = significant & (x < 0)
    neutral = ~(sig_up | sig_down)
    fig, ax = bc.new_fig(7.6, 5.4)
    ax.scatter(x[neutral], y[neutral], s=13, color="#cfd4dc", alpha=0.8, label="NS", edgecolor="none")
    ax.scatter(x[sig_down], y[sig_down], s=20, color="#4e79a7", alpha=0.85, label=down_label, edgecolor="none")
    ax.scatter(x[sig_up], y[sig_up], s=20, color="#e15759", alpha=0.85, label=up_label, edgecolor="none")
    if 0 < padj_cut < 1:
        ax.axhline(-math.log10(padj_cut), color="#8c8c8c", lw=0.8, ls="--")
    if x_cut > 0:
        ax.axvline(x_cut, color="#8c8c8c", lw=0.8, ls="--")
        ax.axvline(-x_cut, color="#8c8c8c", lw=0.8, ls="--")
    order = np.argsort(-(y.to_numpy() * (significant.to_numpy().astype(float) + 1e-9)))[:max(0, top_label)]
    for idx in order:
        if idx < 0 or idx >= len(frame):
            continue
        if not np.isfinite(float(x.iloc[idx])) or float(y.iloc[idx]) <= 0:
            continue
        label = str(frame[labelcol].iloc[idx])
        ax.annotate(label[:20], (float(x.iloc[idx]), float(y.iloc[idx])),
                    fontsize=6.5, color="#575d66", xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(loc="best")
    return save_fig(fig, outdir, name, title=None)


def _fig_barh(labels, values, outdir: str, name: str, title: str, *,
              x_label: str = "Count", colors=None, max_items: int = 15,
              annotate: bool = True):
    pairs = [(str(l), float(v)) for l, v in zip(labels, values) if v is not None and np.isfinite(float(v))]
    pairs = pairs[:max_items]
    if not pairs:
        return None
    pairs = pairs[::-1]
    names = [p[0][:34] for p in pairs]
    vals = [p[1] for p in pairs]
    fig, ax = bc.new_fig(7.8, max(2.8, 0.32 * len(vals) + 1.3))
    bar_colors = colors if isinstance(colors, (list, tuple)) and len(colors) == len(vals) else palette(len(vals))
    ax.barh(names, vals, color=bar_colors)
    if annotate:
        span = max(abs(v) for v in vals) or 1.0
        for i, value in enumerate(vals):
            ax.text(value + span * 0.01, i, f"{value:.3g}", va="center", fontsize=7.5, color="#575d66")
    ax.set_xlabel(x_label)
    ax.set_title(title)
    ax.margins(x=0.14)
    return save_fig(fig, outdir, name, title=None)


def _fig_bar(labels, values, outdir: str, name: str, title: str, *,
             y_label: str = "Count", colors=None, rotate: int = 30, annotate: bool = True):
    pairs = [(str(l), float(v)) for l, v in zip(labels, values) if v is not None and np.isfinite(float(v))]
    if not pairs:
        return None
    names = [p[0][:24] for p in pairs]
    vals = [p[1] for p in pairs]
    fig, ax = bc.new_fig(max(6.2, 0.62 * len(vals) + 3.2), 4.2)
    bar_colors = colors if isinstance(colors, (list, tuple)) and len(colors) == len(vals) else palette(len(vals))
    ax.bar(np.arange(len(vals)), vals, color=bar_colors, width=0.7)
    if annotate:
        span = max(abs(v) for v in vals) or 1.0
        for i, value in enumerate(vals):
            ax.text(i, value + span * 0.015, f"{value:.3g}", ha="center", fontsize=7.5, color="#575d66")
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(names, rotation=rotate, ha="right")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.margins(y=0.14)
    return save_fig(fig, outdir, name, title=None)


def _fig_heat(frame: pd.DataFrame, outdir: str, name: str, title: str, *,
              x_label: str = "Samples", y_label: str = "Features",
              cbar_label: str = "row z-score", max_rows: int = 40):
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return None
    mat = frame.to_numpy(dtype=np.float64)
    mat = np.nan_to_num(mat, nan=0.0)
    frame = frame.head(max_rows)
    mat = mat[:max_rows]
    z = bc.zscore_rows(mat)
    height = max(3.2, 0.24 * mat.shape[0] + 1.6)
    width = max(6.0, 0.34 * mat.shape[1] + 2.6)
    fig, ax = bc.new_fig(width, height)
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.grid(False)
    ax.set_xticks(np.arange(mat.shape[1]))
    rot = 90 if mat.shape[1] > 8 else 30
    ax.set_xticklabels([str(c)[:16] for c in frame.columns], rotation=rot, fontsize=7)
    step = 1 if mat.shape[0] <= 30 else int(math.ceil(mat.shape[0] / 30))
    ax.set_yticks(np.arange(0, mat.shape[0], step))
    ax.set_yticklabels([str(i)[:30] for i in frame.index[::step]], fontsize=7)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(cbar_label)
    return save_fig(fig, outdir, name, title=None)


def _fig_scatter(x, y, outdir: str, name: str, title: str, *,
                 x_label: str = "x", y_label: str = "y", colors=None,
                 labels=None, top_labels: int = 0, identity: bool = False,
                 log_axes: bool = False):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if not mask.any():
        return None
    xs, ys = x[mask], y[mask]
    fig, ax = bc.new_fig(6.6, 5.2)
    if isinstance(colors, (list, tuple)) and len(colors) == len(x):
        colors_arr = np.asarray([str(c) for c in colors], dtype=object)[mask]
        for color in dict.fromkeys(colors_arr.tolist()):
            sel = colors_arr == color
            ax.scatter(xs[sel], ys[sel], s=18, alpha=0.8, color=color, label=str(color), edgecolor="none")
        ax.legend(loc="best")
    else:
        ax.scatter(xs, ys, s=18, alpha=0.8, color="#4e79a7", edgecolor="none")
    if identity:
        lo = float(min(xs.min(), ys.min()))
        hi = float(max(xs.max(), ys.max()))
        ax.plot([lo, hi], [lo, hi], color="#b0b6bf", lw=0.9, ls="--")
    if top_labels > 0 and labels is not None:
        order = np.argsort(-(np.abs(ys - np.median(ys))))
        for idx in order[:top_labels]:
            ax.annotate(str(labels[idx])[:20], (xs[idx], ys[idx]), fontsize=6.5, color="#575d66",
                        xytext=(3, 3), textcoords="offset points")
    if log_axes:
        ax.set_xscale("symlog")
        ax.set_yscale("symlog")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    return save_fig(fig, outdir, name, title=None)


def _fig_network(edges: pd.DataFrame, outdir: str, name: str, title: str, *,
                 top_nodes: int = 60, max_edges: int = 120,
                 left_label: str = "lncRNA", right_label: str = "mRNA"):
    """共表达/ceRNA 网络图（networkx 不可用时退化为圆环布局）。"""
    if not isinstance(edges, pd.DataFrame) or edges.empty:
        return None
    frame = edges.head(max_edges).copy()
    graph = {}
    for _, row in frame.iterrows():
        left = str(row.get("lncrna") or row.get("source") or row.get("mirna") or "")
        right = str(row.get("mrna") or row.get("target") or row.get("target_mrna") or "")
        if not left or not right:
            continue
        graph.setdefault(left, set()).add(right)
        graph.setdefault(right, set()).add(left)
    if not graph:
        return None
    # 节点裁剪：按度取前 top_nodes
    degree = {node: len(neigh) for node, neigh in graph.items()}
    keep = set(sorted(degree, key=lambda n: -degree[n])[:max(1, top_nodes)])
    nodes = [n for n in graph if n in keep]
    left_nodes = [n for n in nodes if n in (edges.get("lncrna", pd.Series(dtype=str)).astype(str).tolist()
                                            + edges.get("mirna", pd.Series(dtype=str)).astype(str).tolist()
                                            + edges.get("source", pd.Series(dtype=str)).astype(str).tolist())]
    left_set = set(left_nodes)
    pos = None
    try:
        import networkx as nx  # 惰性 import（§3.3）

        graph_obj = nx.Graph()
        for _, row in frame.iterrows():
            left = str(row.get("lncrna") or row.get("source") or row.get("mirna") or "")
            right = str(row.get("mrna") or row.get("target") or row.get("target_mrna") or "")
            if not left or not right or left not in keep or right not in keep:
                continue
            weight = abs(float(row.get("r") or row.get("weight") or 0.0))
            graph_obj.add_edge(left, right, weight=weight,
                               sign=("pos" if float(row.get("r") or 0.0) >= 0 else "neg"))
        if graph_obj.number_of_edges() == 0:
            return None
        pos = nx.spring_layout(graph_obj, seed=0, k=1.4 / math.sqrt(max(1, graph_obj.number_of_nodes())))
        pos = {node: (float(xy[0]), float(xy[1])) for node, xy in pos.items()}
        graph = {node: set(graph_obj.neighbors(node)) for node in graph_obj.nodes}
    except Exception:  # noqa: BLE001 - 无 networkx / 布局失败 → 圆环
        pos = None
    if pos is None:
        nodes = list(graph)[:max_nodes_default(top_nodes)]
        pos = {}
        for i, node in enumerate(nodes):
            theta = 2 * math.pi * i / max(1, len(nodes))
            pos[node] = (math.cos(theta), math.sin(theta))
    fig, ax = bc.new_fig(8.0, 6.4)
    ax.set_axis_off()
    for _, row in frame.iterrows():
        left = str(row.get("lncrna") or row.get("source") or row.get("mirna") or "")
        right = str(row.get("mrna") or row.get("target") or row.get("target_mrna") or "")
        if left not in pos or right not in pos:
            continue
        r = float(row.get("r") or row.get("weight") or 0.0)
        ax.plot([pos[left][0], pos[right][0]], [pos[left][1], pos[right][1]],
                color=("#d95252" if r >= 0 else "#4e79a7"),
                lw=max(0.5, min(2.4, abs(r) * 2.0)), alpha=0.55, zorder=1)
    xs = [pos[n][0] for n in pos]
    ys = [pos[n][1] for n in pos]
    colors = ["#f28e2b" if node in left_set else "#4e79a7" for node in pos]
    sizes = [60 + 26 * len(graph.get(node, ())) for node in pos]
    ax.scatter(xs, ys, s=sizes, color=colors, alpha=0.9, edgecolor="white", linewidth=0.8, zorder=2)
    top_degree = sorted(pos, key=lambda n: -len(graph.get(n, ())))[:18]
    for node in top_degree:
        ax.annotate(node[:22], pos[node], fontsize=6.8, color="#3c424a",
                    xytext=(4, 4), textcoords="offset points", zorder=3)
    ax.set_title(title)
    ax.text(0.99, 0.01, f"nodes={len(pos)} edges={len(frame)}", transform=ax.transAxes,
            ha="right", fontsize=7.5, color="#7a828c")
    return save_fig(fig, outdir, name, title=None)


def max_nodes_default(value: int) -> int:
    return max(8, int(value or 60))


def _safe_csv(frame: pd.DataFrame, outdir: str, name: str, label: str, data_out: list):
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return
    try:
        path = export_csv(frame, outdir, name)
        data_out.append(bc.data_file("table_csv", path, label))
    except Exception:  # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# 统计工具：相关矩阵（靶基因/共表达/ceRNA 共用）
# ---------------------------------------------------------------------------

def _zscore_columns(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=np.float64)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    mean = np.nanmean(arr, axis=0)
    std = np.nanstd(arr, axis=0)
    std = np.where(np.isfinite(std) & (std > 0), std, 1.0)
    return np.nan_to_num((arr - mean) / std, nan=0.0, posinf=0.0, neginf=0.0)


def _pearson_matrix(A, B) -> tuple[np.ndarray, np.ndarray]:
    """A(n×p) 与 B(n×q) 的逐列 Pearson 相关 → (r[p×q], p 值[p×q])。"""
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    if A.ndim == 1:
        A = A.reshape(-1, 1)
    if B.ndim == 1:
        B = B.reshape(-1, 1)
    n = A.shape[0]
    Az = _zscore_columns(A)
    Bz = _zscore_columns(B)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (Az.T @ Bz) / max(n, 1)
    r = np.clip(np.nan_to_num(r, nan=0.0), -1.0, 1.0)
    if n > 2:
        with np.errstate(divide="ignore", invalid="ignore"):
            t = r * np.sqrt((n - 2) / np.clip(1.0 - r * r, 1e-12, None))
        try:
            from scipy import stats

            p = 2.0 * stats.t.sf(np.abs(t), n - 2)
        except Exception:  # noqa: BLE001 - scipy 不可用时用正态近似
            p = np.vectorize(lambda value: math.erfc(abs(float(value)) / math.sqrt(2.0)))(t)
    else:
        p = np.ones_like(r)
    p = np.clip(np.nan_to_num(p, nan=1.0), 0.0, 1.0)
    return r, p


def _bh_matrix(pvals: np.ndarray) -> np.ndarray:
    """对整张 p 值矩阵做一次 BH，返回同形状 q 值。"""
    arr = np.asarray(pvals, dtype=np.float64)
    flat = arr.reshape(-1)
    return bc.bh_fdr(flat).reshape(arr.shape)


def _align_samples(target_samples, source_samples):
    """返回 (target 行下标, source 行下标, 共同样本名) 按 target 顺序。"""
    lookup = {str(s): i for i, s in enumerate(source_samples)}
    tgt_idx, src_idx, names = [], [], []
    for i, sample in enumerate(target_samples):
        key = lookup.get(str(sample))
        if key is not None:
            tgt_idx.append(i)
            src_idx.append(key)
            names.append(str(sample))
    return np.asarray(tgt_idx, dtype=int), np.asarray(src_idx, dtype=int), names


def _top_variable(matrix: np.ndarray, limit: int) -> np.ndarray:
    arr = np.asarray(matrix, dtype=np.float64)
    if arr.shape[1] <= limit:
        return np.arange(arr.shape[1])
    var = np.nanvar(arr, axis=0)
    return np.argsort(-np.nan_to_num(var))[:limit]


def _is_mrna_like(name: str) -> bool:
    upper = str(name).upper()
    if not upper or upper.startswith(("HSA-", "HSA_CIRC_")):
        return False
    if re.match(r"^LINC\d+", upper) or re.search(r"-AS\d*$", upper) or re.search(r"-IT\d*$", upper):
        return False
    if upper in _LNCRNA_KNOWN_UPPER:
        return False
    return True


# ---------------------------------------------------------------------------
# circRNA（CONTRACT §8.5：ID 解析 / 差异 / back-spliced 比例 / 圈图）
# ---------------------------------------------------------------------------

_CIRC_COORD_RE = re.compile(
    r"^(chr[0-9XYM]+)[:\-_](\d+)[\-_|:](\d+)(?:[:\-_]?([+-]))?$", re.I
)
_CIRC_JOIN_RE = re.compile(r"^(chr[0-9XYM]+)[:\-_](\d+)[:\-_](\d+)[:\-_]?([+-])$", re.I)


def parse_circ_id(circ_id) -> dict:
    """解析 circRNA ID。

    支持：`chr1:12345-12789:+`、`chr1_12345_12789_+`、`chr1:12345|12789:+`、
    `chr1:12345-12789`（无链）、`hsa_circ_0001234`、`circ_000123`、Ensembl ID。
    """
    raw = str(circ_id).strip()
    out = {
        "circ_id": raw,
        "chrom": "",
        "start": None,
        "end": None,
        "strand": "",
        "length": None,
        "coord_ok": False,
        "naming": "unknown",
    }
    match = _CIRC_COORD_RE.match(raw) or _CIRC_JOIN_RE.match(raw)
    if match:
        chrom, start, end, strand = match.groups()
        start_i, end_i = int(start), int(end)
        if end_i >= start_i:
            out.update(
                chrom=chrom,
                start=start_i,
                end=end_i,
                strand=strand or "",
                length=end_i - start_i,
                coord_ok=True,
                naming="genomic_coordinate",
            )
            return out
    upper = raw.upper()
    if re.match(r"^HSA_CIRC_?\d+", upper):
        out["naming"] = "hsa_circ_id"
    elif re.match(r"^CIRC_?\d+", upper):
        out["naming"] = "circ_id"
    elif re.match(r"^ENS[A-Z]*G\d+", upper):
        out["naming"] = "ensembl_id"
    elif re.match(r"^CHR[0-9XYM]+[:\-_]", upper) or ":" in raw:
        out["naming"] = "region_like"
    elif "CIRC" in upper:
        out["naming"] = "circ_name"
    return out


def _circ_annotation(genes) -> pd.DataFrame:
    frame = pd.DataFrame([parse_circ_id(gene) for gene in genes])
    if frame.empty:
        frame = pd.DataFrame(columns=["circ_id", "chrom", "start", "end", "strand", "length",
                                      "coord_ok", "naming"])
    return frame


_CIRC_ANNOT_COLUMNS = [
    {"key": "circ_id", "label": "circRNA ID", "type": "string"},
    {"key": "chrom", "label": "染色体", "type": "string"},
    {"key": "start", "label": "起点", "type": "number"},
    {"key": "end", "label": "终点", "type": "number"},
    {"key": "strand", "label": "链", "type": "string"},
    {"key": "length", "label": "长度(bp)", "type": "number"},
    {"key": "naming", "label": "ID 命名类型", "type": "string"},
    {"key": "coord_ok", "label": "坐标可解析", "type": "boolean"},
]


def _circ_de_output(ctx, params, outdir, inputs):
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    tables: list[dict] = []
    figures: list[dict] = []
    de, case, control = _diff_frame(inputs, params, label="circRNA", warnings=warnings, notes=notes)
    annotation = _circ_annotation(inputs["genes"])
    frame = de.copy()
    frame = frame.join(annotation.set_index("circ_id")[["chrom", "start", "end", "strand",
                                                         "length", "naming"]], on="gene")
    table_obj = _de_table_obj(
        frame,
        "circ_de",
        "circRNA 差异分析",
        note="Welch t + BH（bg_common.welch_de 统计核心）；regulation 依据 FDR 与 |logFC| 阈值",
        lead_columns=[
            {"key": "chrom", "label": "染色体", "type": "string"},
            {"key": "start", "label": "起点", "type": "number"},
            {"key": "end", "label": "终点", "type": "number"},
            {"key": "strand", "label": "链", "type": "string"},
        ],
    )
    tables.append(table_obj)
    figures.append(_fig_volcano(
        frame, outdir, "circ_volcano", "circRNA differential expression",
        xcol="logFC", pcol="adj.P.Val", labelcol="gene",
        padj_cut=float(params.get("padj") or _PADJ_DEFAULT),
        x_cut=float(params.get("lfc") or _LFC_DEFAULT),
        up_label="Up in case", down_label="Down in case",
    ))
    top = frame[frame["regulation"].isin(["up", "down"])].head(20)
    if len(top) >= 3:
        idx = [inputs["genes"].index(gene) for gene in top["gene"] if gene in set(inputs["genes"])]
        if idx:
            matrix = pd.DataFrame(
                np.asarray(inputs["X"])[:, idx].T,
                index=[inputs["genes"][i] for i in idx],
                columns=inputs["samples"],
            )
            heat = _fig_heat(matrix, outdir, "circ_de_heat", "Top circRNAs (row z-score)",
                             x_label="Samples", y_label="circRNA")
            if heat:
                figures.append(heat)
    bullet, metrics = _de_summary(frame, "circRNA", case, control,
                                  float(params.get("padj") or _PADJ_DEFAULT))
    _safe_csv(frame, outdir, "circ_de", "circRNA 差异表", data_out)
    notes.append("差异检验为 Welch t + BH（与 bulk DE 共用 bg_common 统计核心，未使用 limma/DESeq2 模型）")
    return {
        "bullets": [bullet],
        "metrics": metrics,
        "tables": tables,
        "figures": [f for f in figures if f],
        "warnings": warnings,
        "notes": notes,
        "data_out": data_out,
    }


def _circ_ratio_output(ctx, params, outdir, inputs):
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    linear = inputs.get("linear")
    if linear is None:
        raise EngineError(
            "prerequisite",
            "circRNA 数据里没有 linear（线性转录本）列，无法计算 back-spliced/total 比例："
            "请提供含 `linear` 列的表格（circ 表达 + 对应线性表达）",
        )
    genes = inputs["genes"]
    X = np.asarray(inputs["X"], dtype=np.float64)
    lin = np.asarray(linear, dtype=np.float64).reshape(-1)
    if lin.size != len(genes):
        raise EngineError("bad_request", f"linear 列长度（{lin.size}）与 circRNA 数（{len(genes)}）不一致")
    back_mean = np.nanmean(X, axis=0) if X.size else np.zeros(len(genes))
    total = back_mean + lin
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(total > 0, back_mean / np.where(total > 0, total, 1.0), np.nan)
    annotation = _circ_annotation(genes)
    frame = pd.DataFrame({
        "circ_id": genes,
        "chrom": annotation["chrom"].to_numpy(),
        "back_spliced_mean": back_mean,
        "linear": lin,
        "total": total,
        "circ_ratio": ratio,
    }).sort_values("circ_ratio", ascending=False, na_position="last").reset_index(drop=True)
    tables = [table(
        "circ_ratio",
        "circRNA back-spliced / total 比例",
        [
            {"key": "circ_id", "label": "circRNA ID", "type": "string"},
            {"key": "chrom", "label": "染色体", "type": "string"},
            {"key": "back_spliced_mean", "label": "平均 back-spliced 表达", "type": "number"},
            {"key": "linear", "label": "linear 表达", "type": "number"},
            {"key": "total", "label": "总表达", "type": "number"},
            {"key": "circ_ratio", "label": "circ/linear 比例", "type": "number"},
        ],
        frame,
        note="circ_ratio = mean(back-spliced) / (mean(back-spliced) + linear)；仅对同时有 circ 与 linear 定量的特征计算",
    )]
    figures = []
    if len(frame) >= 3:
        figures.append(_fig_scatter(
            np.log1p(lin), np.log1p(back_mean), outdir, "circ_ratio_scatter",
            "Back-spliced vs linear expression",
            x_label="log1p(linear)", y_label="log1p(back-spliced)",
            labels=genes, top_labels=6,
        ))
        top = frame.dropna(subset=["circ_ratio"]).head(12)
        if len(top):
            bar = _fig_barh(top["circ_id"], top["circ_ratio"], outdir, "circ_ratio_top",
                            "Top circRNA back-spliced ratio", x_label="circ / (circ + linear)")
            if bar:
                figures.append(bar)
    valid = frame["circ_ratio"].dropna()
    bullet = (
        f"circ/linear 比例：{len(valid)} 个特征可计算，中位数 "
        f"{float(valid.median()) if len(valid) else float('nan'):.2f}，>0.5 的有 {int((valid > 0.5).sum())} 个"
    ) if len(valid) else "circ/linear 比例：没有可计算的特征（linear 全为 0 或缺失）"
    _safe_csv(frame, outdir, "circ_ratio", "circRNA back-spliced 比例表", data_out)
    return {
        "bullets": [bullet],
        "metrics": {"n_ratio": int(len(valid)),
                    "median_ratio": round(float(valid.median()), 4) if len(valid) else None},
        "tables": tables,
        "figures": [f for f in figures if f],
        "warnings": warnings,
        "notes": notes + ["circ_ratio 按特征均值计算；若需要逐样本比例，请提供样本级 linear 矩阵"],
        "data_out": data_out,
    }


def _circ_circos_output(ctx, params, outdir, inputs):
    warnings: list[str] = []
    notes: list[str] = []
    annotation = _circ_annotation(inputs["genes"])
    usable = annotation[annotation["coord_ok"]].copy()
    if len(usable) < 3:
        warnings.append(
            f"只有 {len(usable)} 个 circRNA ID 能解析出染色体坐标，圈图至少需要 3 个（如 chr1:1000-2000:+）"
        )
        return {"bullets": [], "metrics": {"n_circos": 0}, "tables": [], "figures": [],
                "warnings": warnings, "notes": notes, "data_out": []}
    mean_expr = pd.Series(np.nanmean(np.asarray(inputs["X"], dtype=np.float64), axis=0),
                          index=inputs["genes"])
    usable["expr"] = usable["circ_id"].map(mean_expr).fillna(0.0)
    chrom_order = [f"chr{i}" for i in list(range(1, 23))] + ["chrX", "chrY", "chrM"]
    chroms = sorted(usable["chrom"].unique().tolist(),
                    key=lambda c: (chrom_order.index(c) if c in chrom_order else 99, c))
    if not chroms:
        return {"bullets": [], "metrics": {"n_circos": 0}, "tables": [], "figures": [],
                "warnings": warnings, "notes": notes, "data_out": []}
    expr_span = float(usable["expr"].max()) or 1.0
    fig, ax = bc.new_fig(7.4, 7.4, subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    try:
        ax.spines["polar"].set_visible(False)
    except Exception:  # noqa: BLE001
        pass
    n_chr = len(chroms)
    step = 2 * math.pi / max(1, n_chr)
    gap = min(0.08, step * 0.12)
    colors = palette(n_chr)
    for i, chrom in enumerate(chroms):
        theta0 = i * step + gap / 2.0
        theta1 = (i + 1) * step - gap / 2.0
        ax.plot([theta0, theta1], [1.0, 1.0], color="#c8ccd4", lw=1.6, solid_capstyle="round")
        mid = (theta0 + theta1) / 2.0
        ax.text(mid, 1.16, _ascii_text(chrom, f"chr{i + 1}"), ha="center", va="center",
                fontsize=7.5, color="#3c424a")
        sub = usable[usable["chrom"] == chrom]
        span = float(sub["end"].max()) if len(sub) else 1.0
        span = span if span > 0 else 1.0
        for _, row in sub.iterrows():
            frac = float(np.clip(((row["start"] + row["end"]) / 2.0) / span, 0.0, 1.0))
            theta = theta0 + frac * (theta1 - theta0)
            radius = 0.94 + 0.28 * float(np.log1p(row["expr"]) / math.log1p(expr_span))
            color = "#e15759" if str(row["strand"]) == "+" else ("#4e79a7" if str(row["strand"]) == "-"
                                                                 else "#8c8c8c")
            ax.scatter([theta], [radius], s=26, color=color, alpha=0.8,
                       edgecolor="white", linewidth=0.4, zorder=3)
        _ = colors[i]
    ax.set_title("circRNA genomic landscape (polar)", pad=18)
    ax.text(0.5, -0.04, "radius = log1p(mean expression) · color = strand",
            transform=ax.transAxes, ha="center", fontsize=7.5, color="#7a828c")
    figure = save_fig(fig, outdir, "circos", title=None)
    return {
        "bullets": [f"圈图：{len(usable)} 个可定位 circRNA 分布在 {n_chr} 条染色体上"],
        "metrics": {"n_circos": int(len(usable)), "n_chrom": int(n_chr)},
        "tables": [],
        "figures": [figure],
        "warnings": warnings,
        "notes": notes + ["圈图按染色体等分角度、用 ID 内坐标定位；无坐标的 hsa_circ_… ID 只出现在注释表"],
        "data_out": [],
    }


def _circ_expression_output(inputs, outdir, top_n: int = 12):
    """无坐标/无分组时的兜底图：平均表达最高的 circRNA。"""
    mean_expr = np.nanmean(np.asarray(inputs["X"], dtype=np.float64), axis=0)
    order = np.argsort(-np.nan_to_num(mean_expr))[:max(3, top_n)]
    labels = [inputs["genes"][i] for i in order]
    values = [float(mean_expr[i]) for i in order]
    figure = _fig_barh(labels, values, outdir, "circ_expression",
                       "Top expressed circRNAs", x_label="Mean expression")
    return {"bullets": [], "metrics": {}, "tables": [], "figures": [figure] if figure else [],
            "warnings": [], "notes": [], "data_out": []}


def _step_circ(ctx, params, outdir):
    warnings: list[str] = []
    notes: list[str] = []
    tables: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []
    metrics: dict = {}
    sections: list[dict] = []

    inputs = _circ_inputs(ctx, params, warnings=warnings, notes=notes)
    annotation = _circ_annotation(inputs["genes"])
    n_coord = int(annotation["coord_ok"].sum())
    n_named = int((annotation["naming"].isin(["hsa_circ_id", "circ_id", "circ_name"])).sum())
    tables.append(table(
        "circ_annotation",
        "circRNA ID 解析",
        _CIRC_ANNOT_COLUMNS,
        annotation,
        note="坐标型 ID 可解析染色体位置（用于圈图）；hsa_circ_… 等命名型 ID 只保留名称",
    ))
    bullets.append(
        f"circRNA：{len(inputs['genes'])} 个 circ × {len(inputs['samples'])} 样本"
        f"（来源 {inputs['source']}）；ID 解析出坐标 {n_coord} 个、circ 命名型 {n_named} 个"
    )
    metrics.update({"n_circ": int(len(inputs["genes"])), "n_samples": int(len(inputs["samples"])),
                    "n_coord": n_coord})

    if inputs["groups"] is None:
        warnings.append("circRNA 样本与样本表的分组列对不上，跳过了差异分析（可用 sample_groups / case / control 指定）")
    else:
        try:
            section = _circ_de_output(ctx, params, outdir, inputs)
            sections.append(section)
            metrics.update(section.get("metrics") or {})
        except EngineError as error:
            warnings.append(f"circRNA 差异分析跳过：{error.message}")

    linear, linear_col = inputs.get("linear"), None
    if linear is not None:
        linear_col = next((str(c) for c in inputs["info"].columns
                           if str(c).strip().lower() in _CIRC_ANNOT_NUMERIC), None)
        try:
            section = _circ_ratio_output(ctx, params, outdir, inputs)
            sections.append(section)
            metrics.update(section.get("metrics") or {})
        except EngineError as error:
            warnings.append(f"circ/linear 比例跳过：{error.message}")
    else:
        warnings.append("表里没有 linear（线性转录本）列，跳过 back-spliced/total 比例（circ_ratio 分析需要它）")

    sections.append(_circ_circos_output(ctx, params, outdir, inputs))

    figures = [f for section in sections for f in (section.get("figures") or []) if f]
    if not figures:
        sections.append(_circ_expression_output(inputs, outdir))
    for section in sections:
        tables.extend(section.get("tables") or [])
        data_out.extend(section.get("data_out") or [])
        bullets.extend(section.get("bullets") or [])
        metrics.update(section.get("metrics") or {})
        warnings.extend(section.get("warnings") or [])
        notes.extend(section.get("notes") or [])
    payload = _finish(KIND_TITLES["circ"], sections, bullets=bullets, tables=tables, data_out=data_out)
    payload["warnings"] = _dedupe_list(warnings + (payload.get("warnings") or []))
    payload["notes"] = _dedupe_list(notes + (payload.get("notes") or []))
    payload["summary"]["metrics"].update(metrics)
    return payload


# ---------------------------------------------------------------------------
# miRNA（名称规范化 / 差异 / 靶基因 / ceRNA）
# ---------------------------------------------------------------------------

# 常见 precursor → mature 映射（≥60 条；未收录时按规则猜测并标注）
_MIRNA_PRECURSOR_MAP = {
    "hsa-mir-21": "hsa-miR-21-5p", "hsa-mir-155": "hsa-miR-155-5p",
    "hsa-mir-146a": "hsa-miR-146a-5p", "hsa-mir-146b": "hsa-miR-146b-5p",
    "hsa-mir-34a": "hsa-miR-34a-5p", "hsa-mir-34b": "hsa-miR-34b-5p",
    "hsa-mir-34c": "hsa-miR-34c-5p",
    "hsa-mir-200a": "hsa-miR-200a-3p", "hsa-mir-200b": "hsa-miR-200b-3p",
    "hsa-mir-200c": "hsa-miR-200c-3p",
    "hsa-mir-9-1": "hsa-miR-9-5p", "hsa-mir-9-2": "hsa-miR-9-5p",
    "hsa-mir-9-3": "hsa-miR-9-5p",
    "hsa-mir-let-7a-1": "hsa-let-7a-5p", "hsa-mir-let-7a-2": "hsa-let-7a-5p",
    "hsa-mir-let-7a-3": "hsa-let-7a-5p", "hsa-mir-let-7b": "hsa-let-7b-5p",
    "hsa-mir-let-7c": "hsa-let-7c-5p", "hsa-mir-let-7d": "hsa-let-7d-5p",
    "hsa-mir-let-7e": "hsa-let-7e-5p", "hsa-mir-let-7f-1": "hsa-let-7f-5p",
    "hsa-mir-let-7f-2": "hsa-let-7f-5p", "hsa-mir-let-7g": "hsa-let-7g-5p",
    "hsa-mir-let-7i": "hsa-let-7i-5p",
    "hsa-mir-10b": "hsa-miR-10b-5p", "hsa-mir-15a": "hsa-miR-15a-5p",
    "hsa-mir-16-1": "hsa-miR-16-5p", "hsa-mir-16-2": "hsa-miR-16-5p",
    "hsa-mir-17": "hsa-miR-17-5p", "hsa-mir-18a": "hsa-miR-18a-5p",
    "hsa-mir-19a": "hsa-miR-19a-3p", "hsa-mir-19b-1": "hsa-miR-19b-3p",
    "hsa-mir-19b-2": "hsa-miR-19b-3p", "hsa-mir-20a": "hsa-miR-20a-5p",
    "hsa-mir-20b": "hsa-miR-20b-5p", "hsa-mir-22": "hsa-miR-22-3p",
    "hsa-mir-23a": "hsa-miR-23a-3p", "hsa-mir-23b": "hsa-miR-23b-3p",
    "hsa-mir-24-1": "hsa-miR-24-3p", "hsa-mir-24-2": "hsa-miR-24-3p",
    "hsa-mir-25": "hsa-miR-25-3p", "hsa-mir-26a-1": "hsa-miR-26a-5p",
    "hsa-mir-26a-2": "hsa-miR-26a-5p", "hsa-mir-26b": "hsa-miR-26b-5p",
    "hsa-mir-27a": "hsa-miR-27a-3p", "hsa-mir-27b": "hsa-miR-27b-3p",
    "hsa-mir-29a": "hsa-miR-29a-3p", "hsa-mir-29b-1": "hsa-miR-29b-3p",
    "hsa-mir-29b-2": "hsa-miR-29b-3p", "hsa-mir-29c": "hsa-miR-29c-3p",
    "hsa-mir-30a": "hsa-miR-30a-5p", "hsa-mir-30b": "hsa-miR-30b-5p",
    "hsa-mir-30c-1": "hsa-miR-30c-5p", "hsa-mir-31": "hsa-miR-31-5p",
    "hsa-mir-92a-1": "hsa-miR-92a-3p", "hsa-mir-92a-2": "hsa-miR-92a-3p",
    "hsa-mir-93": "hsa-miR-93-5p", "hsa-mir-96": "hsa-miR-96-5p",
    "hsa-mir-101-1": "hsa-miR-101-3p", "hsa-mir-101-2": "hsa-miR-101-3p",
    "hsa-mir-124-1": "hsa-miR-124-3p", "hsa-mir-124-2": "hsa-miR-124-3p",
    "hsa-mir-124-3": "hsa-miR-124-3p",
    "hsa-mir-125a": "hsa-miR-125a-5p", "hsa-mir-125b-1": "hsa-miR-125b-5p",
    "hsa-mir-125b-2": "hsa-miR-125b-5p", "hsa-mir-126": "hsa-miR-126-3p",
    "hsa-mir-133a-1": "hsa-miR-133a-3p", "hsa-mir-133a-2": "hsa-miR-133a-3p",
    "hsa-mir-141": "hsa-miR-141-3p", "hsa-mir-143": "hsa-miR-143-3p",
    "hsa-mir-145": "hsa-miR-145-5p", "hsa-mir-148a": "hsa-miR-148a-3p",
    "hsa-mir-150": "hsa-miR-150-5p", "hsa-mir-181a-1": "hsa-miR-181a-5p",
    "hsa-mir-181a-2": "hsa-miR-181a-5p", "hsa-mir-181b-1": "hsa-miR-181b-5p",
    "hsa-mir-182": "hsa-miR-182-5p", "hsa-mir-191": "hsa-miR-191-5p",
    "hsa-mir-203": "hsa-miR-203a-3p", "hsa-mir-205": "hsa-miR-205-5p",
    "hsa-mir-210": "hsa-miR-210-3p", "hsa-mir-214": "hsa-miR-214-3p",
    "hsa-mir-221": "hsa-miR-221-3p", "hsa-mir-222": "hsa-miR-222-3p",
    "hsa-mir-223": "hsa-miR-223-3p", "hsa-mir-296": "hsa-miR-296-5p",
    "hsa-mir-335": "hsa-miR-335-5p", "hsa-mir-375": "hsa-miR-375-3p",
    "hsa-mir-451a": "hsa-miR-451a", "hsa-mir-486-1": "hsa-miR-486-5p",
    "hsa-mir-486-2": "hsa-miR-486-5p", "hsa-mir-493": "hsa-miR-493-5p",
    "hsa-mir-497": "hsa-miR-497-5p", "hsa-mir-503": "hsa-miR-503-5p",
}

# 内置「常见 miRNA → 靶基因」简化表（≥40 个家族；lite：非全基因组靶点预测，未做种子序列匹配）
_MIRNA_TARGETS = {
    "miR-21": ["PTEN", "PDCD4", "SPRY2", "RECK", "TIMP3", "BTG2"],
    "miR-155": ["SOCS1", "INPP5D", "TP53INP1", "FOXO3", "CEBPB", "BACH1"],
    "let-7a": ["KRAS", "HMGA2", "MYC", "LIN28B", "IGF2BP1", "CCND1"],
    "let-7b": ["KRAS", "HMGA2", "MYC", "LIN28B", "CCND1"],
    "let-7c": ["KRAS", "HMGA2", "MYC", "LIN28B"],
    "let-7d": ["KRAS", "HMGA2", "MYC"],
    "let-7e": ["KRAS", "HMGA2"],
    "let-7f": ["KRAS", "HMGA2", "MYC"],
    "let-7g": ["KRAS", "HMGA2", "MYC"],
    "miR-145": ["SOX2", "KLF4", "MUC1", "MYC", "FSCN1", "YES1"],
    "miR-34a": ["SIRT1", "CD44", "NOTCH1", "BCL2", "MET"],
    "miR-34b": ["CD44", "NOTCH1", "BCL2"],
    "miR-34c": ["CD44", "NOTCH1", "BCL2", "MET"],
    "miR-200a": ["ZEB1", "ZEB2", "CTNNB1"],
    "miR-200b": ["ZEB1", "ZEB2", "VEGFA", "DNMT1"],
    "miR-200c": ["ZEB1", "ZEB2", "VIM", "FN1", "CDH11"],
    "miR-9": ["CDH1", "FOXO1", "SOX2", "REST"],
    "miR-10b": ["HOXD10", "RHOC", "PTEN"],
    "miR-17": ["PTEN", "E2F1", "BCL2L11", "CDKN1A"],
    "miR-18a": ["ESR1", "ATM", "HIF1A"],
    "miR-19a": ["PTEN", "SOCS1", "ATXN1"],
    "miR-19b": ["PTEN", "SOCS1"],
    "miR-20a": ["E2F1", "PTEN", "BMPR2", "CDKN1A"],
    "miR-92a": ["PTEN", "CDH1", "ITGA5"],
    "miR-93": ["PTEN", "E2F1", "CDKN1A"],
    "miR-106a": ["PTEN", "E2F1", "CDKN1A"],
    "miR-221": ["PTEN", "CDKN1B", "BBC3", "TIMP3", "KIT"],
    "miR-222": ["PTEN", "CDKN1B", "TIMP3", "KIT"],
    "miR-126": ["VEGFA", "PIK3R2", "SPRED1", "VCAM1"],
    "miR-143": ["KRAS", "MYO6", "BCL2", "MAPK7"],
    "miR-146a": ["TRAF6", "IRAK1", "CXCR4", "STAT1"],
    "miR-146b": ["TRAF6", "IRAK1", "CXCR4"],
    "miR-148a": ["MET", "DNMT1", "ROCK1", "NR1I2"],
    "miR-101": ["EZH2", "PTGS2", "MCL1", "ROCK1"],
    "miR-125a": ["ERBB2", "ERBB3", "LIN28B", "VEGFA"],
    "miR-125b": ["BAK1", "TP53", "ERBB2", "ERBB3", "LIN28B"],
    "miR-133a": ["FGFR1", "EGFR", "TAGLN2"],
    "miR-141": ["ZEB1", "ZEB2", "TGFB2"],
    "miR-26a": ["EZH2", "CCND2", "CCNE2", "MTDH"],
    "miR-26b": ["EZH2", "CCND2", "CCNE2"],
    "miR-27a": ["ZEB1", "EGFR", "MET"],
    "miR-29a": ["DNMT3A", "DNMT3B", "COL1A1", "BCL2", "MCL1"],
    "miR-29b": ["DNMT3A", "DNMT3B", "COL1A1", "MCL1"],
    "miR-29c": ["DNMT3A", "DNMT3B", "COL1A1"],
    "miR-30a": ["SNAI1", "BECN1", "VIM"],
    "miR-30b": ["SNAI1", "BECN1"],
    "miR-31": ["SATB2", "MET", "RHOA"],
    "miR-203": ["SNAI2", "ABL1", "TP63"],
    "miR-205": ["ZEB1", "ZEB2", "E2F1"],
    "miR-210": ["GPD1L", "SDHD", "E2F3", "ISCU"],
    "miR-214": ["PTEN", "EZH2", "TWIST1"],
    "miR-375": ["PDPK1", "JAK2", "YAP1", "CTNNB1"],
    "miR-451a": ["CAB39", "YWHAZ", "MIF"],
    "miR-486": ["PTEN", "FOXO1", "SIRT1"],
    "miR-15a": ["BCL2", "CCND1", "MCL1", "WNT3A"],
    "miR-16": ["BCL2", "CCND1", "MCL1", "WNT3A"],
    "miR-497": ["BCL2", "CCND1", "WEE1"],
    "miR-503": ["CCND1", "CHEK1", "VEGFA"],
    "miR-22": ["PTEN", "HDAC4", "ESR1", "SIRT1"],
    "miR-181a": ["PTEN", "ATM", "BCL2"],
    "miR-181b": ["PTEN", "ATM", "CYLD"],
    "miR-182": ["FOXO3", "MITF", "ADCY6"],
    "miR-124": ["CDK4", "CDK6", "SNAI2", "STAT3"],
    "miR-137": ["CDK6", "EZH2", "MITF"],
    "miR-138": ["ZEB2", "EZH2", "HIF1A", "TWIST1"],
    "miR-139": ["HOXB5", "CXCR4", "IGF1R"],
    "miR-140": ["HDAC4", "VEGFA"],
    "miR-150": ["MYB", "BCL2", "CXCR4"],
    "miR-152": ["DNMT1", "MET", "PIK3CA"],
    "miR-153": ["SNAI1", "ZEB2"],
    "miR-194": ["RHOB", "CDH1"],
    "miR-199a": ["HIF1A", "MET"],
    "miR-204": ["SOX4", "MCL1"],
    "miR-223": ["FBXW7", "MEF2C", "NLRP3"],
    "miR-25": ["PTEN", "E2F1"],
    "miR-335": ["SOX4"],
    "miR-24": ["MYC", "BCL2L11", "CDKN1B"],
    "miR-23a": ["PTEN", "SMAD3"],
    "miR-196a": ["HOXB8", "ANXA1"],
    "miR-196b": ["HOXB8", "HOXA9"],
    "miR-7": ["EGFR", "RAF1", "IRS1", "IRS2"],
    "miR-128": ["BMI1", "E2F3", "EGFR"],
    "miR-424": ["CCND1", "CHEK1", "VEGFA"],
}


def _mirna_family_body(body: str) -> str:
    text = str(body).strip().lower()
    return re.sub(r"-\d+$", "", text)


def normalize_mirna(name) -> dict:
    """miRNA 名称规范化：`hsa-mir-21` → `hsa-miR-21-5p`；识别 5p/3p 臂、let-7 家族、precursor。"""
    raw = str(name).strip()
    text = re.sub(r"\s+", "", raw.replace("_", "-"))
    text = re.sub(r"-{2,}", "-", text)
    low = text.lower()
    species = "hsa"
    match = re.match(r"^(hsa|mmu|rno|gga|dre|cel|dme|xtr|bta|ssc|ptr|ocu|mdo|efu|mne)-", low)
    if match:
        species = match.group(1)
        low = low[match.end():]
    arm = ""
    match = re.search(r"-(5p|3p)$", low)
    if match:
        arm = match.group(1)
        low = low[:-3]
    precursor = False
    guessed = False
    form = ""
    if low.startswith("mir-let-7"):
        low = low[4:]
    if low.startswith("let-7") or low.startswith("let7"):
        body = low.replace("let7", "let-7", 1) if low.startswith("let7") else low
        family_body = _mirna_family_body(body)
        family_key = family_body if family_body.startswith("let-") else f"let-{family_body}"
        if not arm:
            arm = "5p"
            guessed = True
        canonical = f"{species}-{family_key}-{arm}"
        form = "mature (let-7 family)"
    elif low.startswith("mir-") and re.match(r"^mir-\d", low):
        precursor = True
        body = low[4:]
        family_body = _mirna_family_body(body)
        family_key = f"miR-{family_body}"
        mapped = _MIRNA_PRECURSOR_MAP.get(f"{species}-mir-{body}")
        if mapped:
            canonical = mapped
            arm_match = re.search(r"-(5p|3p)$", canonical)
            if arm_match:
                arm = arm_match.group(1)
            form = "precursor→mature（内置映射）"
        else:
            guess = arm or "5p"
            arm = guess
            canonical = f"{species}-{family_key}-{guess}"
            guessed = True
            form = "precursor→mature（规则猜测 5p）"
    else:
        body = low
        for prefix in ("microrna-", "mirna-", "mir-"):
            if body.startswith(prefix):
                body = body[len(prefix):]
                break
        family_body = _mirna_family_body(body)
        family_key = f"miR-{family_body}"
        canonical = f"{species}-{family_key}-{arm}" if arm else f"{species}-{family_key}"
        form = "mature" if arm else "mature（无臂注释）"
    return {
        "original": raw,
        "canonical": canonical,
        "species": species,
        "family": family_key,
        "arm": arm,
        "precursor": precursor,
        "guessed": guessed,
        "form": form,
    }


def _mirna_targets_for(family: str) -> list[str]:
    fam = str(family)
    if fam in _MIRNA_TARGETS:
        return _MIRNA_TARGETS[fam]
    lower_map = {key.lower(): key for key in _MIRNA_TARGETS}
    if fam.lower() in lower_map:
        return _MIRNA_TARGETS[lower_map[fam.lower()]]
    match = re.match(r"^(miR-\d+)[a-z]$", fam)
    if match and match.group(1) in _MIRNA_TARGETS:
        return _MIRNA_TARGETS[match.group(1)]
    if fam.lower().startswith("let-7"):
        return _MIRNA_TARGETS.get("let-7a", [])
    return []


def _mirna_annotation(genes) -> pd.DataFrame:
    rows = [normalize_mirna(gene) for gene in genes]
    frame = pd.DataFrame(rows)
    if frame.empty:
        frame = pd.DataFrame(columns=["original", "canonical", "species", "family", "arm",
                                      "precursor", "guessed", "form"])
    return frame


_MIRNA_ANNOT_COLUMNS = [
    {"key": "original", "label": "原始名称", "type": "string"},
    {"key": "canonical", "label": "规范化名称", "type": "string"},
    {"key": "species", "label": "物种前缀", "type": "string"},
    {"key": "family", "label": "家族", "type": "string"},
    {"key": "arm", "label": "臂", "type": "string"},
    {"key": "form", "label": "形式", "type": "string"},
    {"key": "precursor", "label": "是否 precursor", "type": "boolean"},
    {"key": "guessed", "label": "是否规则猜测", "type": "boolean"},
]


def _mirna_inputs(ctx, params, *, warnings=None, notes=None) -> dict:
    inputs = _matrix_inputs(ctx, params, "mirna", warnings=warnings, notes=notes)
    inputs["annotation"] = _mirna_annotation(inputs["genes"])
    return inputs


def _mirna_de_output(ctx, params, outdir, inputs):
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    de, case, control = _diff_frame(inputs, params, label="miRNA", warnings=warnings, notes=notes)
    annotation = inputs["annotation"].set_index("original")
    frame = de.copy()
    frame = frame.join(annotation[["canonical", "family", "arm"]], on="gene")
    tables = [_de_table_obj(
        frame,
        "mirna_de",
        "miRNA 差异分析",
        note="名称先规范化（hsa-mir-X → hsa-miR-X-5p）；Welch t + BH，regulation 依据 FDR 与 |logFC| 阈值",
        lead_columns=[
            {"key": "canonical", "label": "规范化名称", "type": "string"},
            {"key": "family", "label": "家族", "type": "string"},
            {"key": "arm", "label": "臂", "type": "string"},
        ],
    )]
    figures = [_fig_volcano(
        frame, outdir, "mirna_volcano", "miRNA differential expression",
        xcol="logFC", pcol="adj.P.Val", labelcol="gene",
        padj_cut=float(params.get("padj") or _PADJ_DEFAULT),
        x_cut=float(params.get("lfc") or _LFC_DEFAULT),
        up_label="Up in case", down_label="Down in case",
    )]
    top = frame[frame["regulation"].isin(["up", "down"])].head(20)
    if len(top) >= 3:
        genes = inputs["genes"]
        index = {gene: i for i, gene in enumerate(genes)}
        idx = [index[gene] for gene in top["gene"] if gene in index]
        if idx:
            matrix = pd.DataFrame(
                np.asarray(inputs["X"])[:, idx].T,
                index=[genes[i] for i in idx],
                columns=inputs["samples"],
            )
            heat = _fig_heat(matrix, outdir, "mirna_de_heat", "Top miRNAs (row z-score)",
                             x_label="Samples", y_label="miRNA")
            if heat:
                figures.append(heat)
    bullet, metrics = _de_summary(frame, "miRNA", case, control,
                                  float(params.get("padj") or _PADJ_DEFAULT))
    _safe_csv(frame, outdir, "mirna_de", "miRNA 差异表", data_out)
    notes.append("差异检验为 Welch t + BH（与 bulk DE 共用 bg_common 统计核心）")
    return {"bullets": [bullet], "metrics": metrics, "tables": tables,
            "figures": [f for f in figures if f], "warnings": warnings, "notes": notes,
            "data_out": data_out}


def _mirna_target_output(ctx, params, outdir, inputs):
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    annotation = inputs["annotation"]
    ctx_genes = [str(g) for g in (ctx.get("genes") or [])]
    ctx_X = ctx.get("X")
    ctx_samples = [str(s) for s in (ctx.get("samples") or [])]
    gene_index = {gene: i for i, gene in enumerate(ctx_genes)}
    rows: list[dict] = []
    missing_targets: set[str] = set()
    common: list[str] = []
    r_mat = p_mat = None
    if isinstance(ctx_X, (np.ndarray, list)) and len(ctx_genes) and len(ctx_samples):
        mir_rows, mrna_rows, common = _align_samples(inputs["samples"], ctx_samples)
        if len(common) >= 4:
            A = np.asarray(inputs["X"], dtype=np.float64)[mir_rows, :]
            B = np.asarray(ctx_X, dtype=np.float64)[mrna_rows, :]
            r_mat, p_mat = _pearson_matrix(A, B)
            notes.append(f"miRNA-mRNA 相关用 {len(common)} 个共同样本计算（Pearson）")
        else:
            warnings.append(f"miRNA 样本与 bulk 样本只有 {len(common)} 个交集（<4），跳过相关计算，只给内置靶基因表")
    else:
        warnings.append("ctx 里没有 mRNA 表达矩阵（X/genes），只给内置靶基因注释、不给相关")

    for j, original in enumerate(inputs["genes"]):
        info = annotation.iloc[j]
        targets = _mirna_targets_for(info["family"])
        if not targets:
            continue
        for target in targets:
            target_idx = gene_index.get(target)
            if target_idx is None:
                missing_targets.add(target)
                continue
            r = float(r_mat[j, target_idx]) if r_mat is not None else None
            p = float(p_mat[j, target_idx]) if p_mat is not None else None
            rows.append({
                "mirna": str(original),
                "canonical": info["canonical"],
                "family": info["family"],
                "target": target,
                "r": r,
                "p": p,
            })
    if not rows:
        warnings.append(
            "内置靶基因表里的靶点在这份 bulk mRNA 矩阵里一个都没找到："
            "请确认 bulk 数据与 miRNA 来自同一物种、基因用 symbol 命名"
        )
        return {"bullets": [], "metrics": {"n_pairs": 0}, "tables": [], "figures": [],
                "warnings": warnings, "notes": notes, "data_out": []}
    frame = pd.DataFrame(rows)
    if frame["p"].notna().any():
        frame["padj"] = np.nan
        valid = frame["p"].notna()
        frame.loc[valid, "padj"] = bc.bh_fdr(frame.loc[valid, "p"].to_numpy())
    else:
        frame["padj"] = np.nan
    frame["neg_corr"] = frame["r"].apply(lambda value: bool(np.isfinite(value) and value < 0)
                                         if value is not None else False)
    frame["support"] = np.where(
        frame["neg_corr"] & (pd.to_numeric(frame["padj"], errors="coerce") <= 0.05),
        "内置靶点 + 显著负相关",
        np.where(frame["neg_corr"], "内置靶点 + 负相关（未达显著）", "仅内置靶点（无相关证据）"),
    )
    frame = frame.sort_values(["neg_corr", "padj", "r"], ascending=[False, True, True]).reset_index(drop=True)
    tables = [table(
        "mirna_target",
        "miRNA → 靶基因（lite）",
        [
            {"key": "mirna", "label": "miRNA", "type": "string"},
            {"key": "canonical", "label": "规范化名称", "type": "string"},
            {"key": "family", "label": "家族", "type": "string"},
            {"key": "target", "label": "靶基因", "type": "string"},
            {"key": "r", "label": "相关系数 r", "type": "number"},
            {"key": "p", "label": "P 值", "type": "number"},
            {"key": "padj", "label": "FDR", "type": "number"},
            {"key": "support", "label": "证据", "type": "string"},
        ],
        frame,
        note="lite 实现：内置常见 miRNA→靶基因 curated 表（≥40 家族）；相关为 bulk miRNA 与 mRNA 的 Pearson；"
             "未做全基因组靶点预测/种子序列匹配",
    )]
    figures = []
    neg = frame[frame["neg_corr"] & frame["padj"].notna()].head(15)
    if len(neg):
        figures.append(_fig_barh(
            (neg["mirna"].astype(str) + " → " + neg["target"].astype(str)).tolist(),
            (-np.log10(np.clip(neg["padj"].to_numpy(dtype=float), 1e-300, None))).tolist(),
            outdir, "mirna_target_bar", "Top miRNA-target pairs (negative correlation)",
            x_label="-log10 FDR", max_items=15,
        ))
    if r_mat is not None and len(common):
        figs = _fig_scatter(
            frame["r"].fillna(0.0).to_numpy(dtype=float),
            (-np.log10(np.clip(frame["p"].fillna(1.0).to_numpy(dtype=float), 1e-300, None))),
            outdir, "mirna_target_volcano", "miRNA-target relationship",
            x_label="Pearson r (miRNA vs mRNA)", y_label="-log10 P",
        )
        figures.append(figs)
    metrics = {
        "n_pairs": int(len(frame)),
        "n_mirna_with_target": int(frame["mirna"].nunique()),
        "n_negative": int(frame["neg_corr"].sum()),
        "n_negative_significant": int((frame["neg_corr"] & (frame["padj"] <= 0.05)).sum()),
        "n_missing_targets": int(len(missing_targets)),
    }
    bullets = [
        f"内置靶基因表命中 {metrics['n_pairs']} 对（{metrics['n_mirna_with_target']} 个 miRNA、"
        f"{int(frame['target'].nunique())} 个靶基因）；负相关 {metrics['n_negative']} 对，"
        f"其中显著 {metrics['n_negative_significant']} 对"
    ]
    if missing_targets:
        warnings.append(f"有 {len(missing_targets)} 个内置靶基因不在 bulk mRNA 矩阵里（如 "
                        f"{', '.join(sorted(missing_targets)[:6])}）")
    _safe_csv(frame, outdir, "mirna_target", "miRNA 靶基因表", data_out)
    notes.append("lite 标注：内置靶基因表是常见文献靶点的简化集合，不能替代 miRDB/TargetScan 全基因组预测")
    return {"bullets": bullets, "metrics": metrics, "tables": tables,
            "figures": [f for f in figures if f], "warnings": warnings, "notes": notes,
            "data_out": data_out}


def _cerna_output(ctx, params, outdir, mirna_inputs, *, strict: bool = False):
    """ceRNA 三元组候选：lncRNA ↑ 与 miRNA ↓，且与 miRNA 靶 mRNA ↑。"""
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    min_r = abs(float(params.get("min_abs_r") or 0.4))
    ctx_genes = [str(g) for g in (ctx.get("genes") or [])]
    ctx_samples = [str(s) for s in (ctx.get("samples") or [])]
    ctx_X = ctx.get("X")
    if not isinstance(ctx_X, (np.ndarray, list)) or not ctx_genes:
        if strict:
            raise EngineError("prerequisite", "ctx 里没有 bulk mRNA 矩阵（X/genes），无法构建 ceRNA 网络")
        warnings.append("没有 bulk mRNA 矩阵，跳过 ceRNA 候选")
        return {"bullets": [], "metrics": {}, "tables": [], "figures": [], "warnings": warnings,
                "notes": notes, "data_out": []}
    ctx_X = np.asarray(ctx_X, dtype=np.float64)

    # lncRNA 来源：优先已挂载的 lncRNA 模块，其次 bulk 里命名像 lncRNA 的基因
    lnc_mod = _module(ctx, "lncrna")
    lnc_names: list[str] = []
    lnc_X = None
    if lnc_mod and isinstance(lnc_mod.get("X"), (np.ndarray, list)) and len(lnc_mod.get("X")):
        lnc_names = [str(g) for g in (lnc_mod.get("genes") or [])]
        lnc_X = np.asarray(lnc_mod["X"], dtype=np.float64)
        lnc_samples = [str(s) for s in (lnc_mod.get("samples") or [])]
    else:
        lnc_idx = [i for i, gene in enumerate(ctx_genes) if _looks_lncrna(gene)]
        if not lnc_idx:
            if strict:
                raise EngineError(
                    "prerequisite",
                    "没有找到 lncRNA 数据：请先挂载 lncRNA 模块，或确保 bulk 基因里含 LINC…/-AS1 等 lncRNA 名称",
                )
            warnings.append("没有 lncRNA 数据（未挂载模块、bulk 里也没有 lncRNA 名称），跳过 ceRNA 候选")
            return {"bullets": [], "metrics": {}, "tables": [], "figures": [], "warnings": warnings,
                    "notes": notes, "data_out": []}
        lnc_names = [ctx_genes[i] for i in lnc_idx]
        lnc_X = ctx_X[:, lnc_idx]
        lnc_samples = ctx_samples

    mir_samples = mirna_inputs["samples"]
    mir_rows, mrna_rows, common = _align_samples(mir_samples, ctx_samples)
    if len(common) < 4:
        if strict:
            raise EngineError("prerequisite", f"miRNA 与 bulk 样本只有 {len(common)} 个交集（<4），无法算 ceRNA 相关")
        warnings.append(f"miRNA 与 bulk 样本交集只有 {len(common)} 个（<4），跳过 ceRNA 候选")
        return {"bullets": [], "metrics": {}, "tables": [], "figures": [], "warnings": warnings,
                "notes": notes, "data_out": []}
    sampled = set(common)
    lnc_sample_idx = [i for i, s in enumerate(lnc_samples) if str(s) in sampled]
    lnc_common_names = [str(lnc_samples[i]) for i in lnc_sample_idx]
    # 只保留三套矩阵共同覆盖的样本（保持顺序一致）
    order = [s for s in common if s in set(lnc_common_names)]
    if len(order) < 4:
        warnings.append(f"miRNA/bulk/lncRNA 三者交集只有 {len(order)} 个样本（<4），跳过 ceRNA 候选")
        return {"bullets": [], "metrics": {}, "tables": [], "figures": [], "warnings": warnings,
                "notes": notes, "data_out": []}
    mir_lookup = {str(s): i for i, s in enumerate(mir_samples)}
    ctx_lookup = {str(s): i for i, s in enumerate(ctx_samples)}
    lnc_lookup = {str(s): i for i, s in enumerate(lnc_samples)}
    mir_rows = [mir_lookup[s] for s in order]
    ctx_rows = [ctx_lookup[s] for s in order]
    lnc_rows = [lnc_lookup[s] for s in order]

    mir_sel = _top_variable(np.asarray(mirna_inputs["X"], dtype=np.float64)[mir_rows, :], 20)
    mrna_sel = _top_variable(ctx_X[ctx_rows, :], 150)
    lnc_sel = _top_variable(lnc_X[lnc_rows, :], 100)
    A = np.asarray(mirna_inputs["X"], dtype=np.float64)[mir_rows, :][:, mir_sel]
    B = ctx_X[ctx_rows, :][:, mrna_sel]
    C = lnc_X[lnc_rows, :][:, lnc_sel]
    r_mp, p_mp = _pearson_matrix(A, B)
    r_lm, p_lm = _pearson_matrix(C, B)
    r_lmi, _p_lmi = _pearson_matrix(C, A)
    mir_names = [mirna_inputs["genes"][i] for i in mir_sel]
    mrna_names = [ctx_genes[i] for i in mrna_sel]
    lnc_names_sel = [lnc_names[i] for i in lnc_sel]
    rows = []
    for a, mir_name in enumerate(mir_names):
        for b, mrna_name in enumerate(mrna_names):
            r_mp_val = float(r_mp[a, b])
            if r_mp_val > -min_r:
                continue
            for c, lnc_name in enumerate(lnc_names_sel):
                r_lm_val = float(r_lm[c, b])
                r_lmi_val = float(r_lmi[c, a])
                if r_lm_val < min_r or r_lmi_val > -0.2:
                    continue
                score = (-r_mp_val) * r_lm_val
                rows.append({
                    "lncrna": lnc_name,
                    "mirna": mir_name,
                    "mrna": mrna_name,
                    "r_lncrna_mrna": round(r_lm_val, 4),
                    "r_mirna_mrna": round(r_mp_val, 4),
                    "r_lncrna_mirna": round(r_lmi_val, 4),
                    "score": round(score, 4),
                })
    if not rows:
        warnings.append(f"没有满足条件的 ceRNA 三元组（|r|>{min_r:g}）；样本数少或信号弱时常见")
        return {"bullets": [], "metrics": {"n_cerna": 0}, "tables": [], "figures": [],
                "warnings": warnings, "notes": notes, "data_out": []}
    frame = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    if len(frame) > 500:
        frame = frame.head(500)
    tables = [table(
        "cerna_candidates",
        "ceRNA 三元组候选",
        [
            {"key": "lncrna", "label": "lncRNA", "type": "string"},
            {"key": "mirna", "label": "miRNA", "type": "string"},
            {"key": "mrna", "label": "靶 mRNA", "type": "string"},
            {"key": "r_lncrna_mrna", "label": "lncRNA–mRNA r", "type": "number"},
            {"key": "r_mirna_mrna", "label": "miRNA–mRNA r", "type": "number"},
            {"key": "r_lncrna_mirna", "label": "lncRNA–miRNA r", "type": "number"},
            {"key": "score", "label": "得分", "type": "number"},
        ],
        frame,
        note="候选规则：miRNA–mRNA r ≤ -min_abs_r、lncRNA–mRNA r ≥ min_abs_r、lncRNA–miRNA r ≤ -0.2；"
             "score = (-r_miRNA_mRNA)·r_lncRNA_mRNA。仅相关证据，非因果/非海绵实验验证",
    )]
    top = frame.head(12)
    figures = []
    if len(top):
        figures.append(_fig_barh(
            (top["lncrna"].astype(str) + " – " + top["mirna"].astype(str) + " – " + top["mrna"].astype(str)).tolist(),
            top["score"].tolist(), outdir, "cerna_top", "Top ceRNA triplets (correlation score)",
            x_label="score", max_items=12,
        ))
        network = _fig_network(frame.head(80), outdir, "cerna_network", "ceRNA candidate network")
        if network:
            figures.append(network)
    metrics = {"n_cerna": int(len(frame)),
               "n_lncrna": int(frame["lncrna"].nunique()),
               "n_mirna": int(frame["mirna"].nunique()),
               "n_mrna": int(frame["mrna"].nunique())}
    bullets = [f"ceRNA 候选 {len(frame)} 个三元组（{metrics['n_lncrna']} lncRNA / "
               f"{metrics['n_mirna']} miRNA / {metrics['n_mrna']} mRNA），共同样本 {len(order)} 个"]
    _safe_csv(frame, outdir, "cerna_candidates", "ceRNA 候选表", data_out)
    notes.append("ceRNA 只是相关性候选（海绵假说），需要实验验证；未做 miRNA 结合位点预测")
    return {"bullets": bullets, "metrics": metrics, "tables": tables,
            "figures": [f for f in figures if f], "warnings": warnings, "notes": notes,
            "data_out": data_out}


def _step_mirna(ctx, params, outdir):
    warnings: list[str] = []
    notes: list[str] = []
    tables: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []
    metrics: dict = {}
    sections: list[dict] = []

    inputs = _mirna_inputs(ctx, params, warnings=warnings, notes=notes)
    annotation = inputs["annotation"]
    tables.append(table(
        "mirna_annotation",
        "miRNA 名称规范化",
        _MIRNA_ANNOT_COLUMNS,
        annotation,
        note="规范化规则：hsa-mir-21 → hsa-miR-21-5p（内置 precursor 映射，未收录时按规则猜 5p 并标记 guessed）；"
             "识别 5p/3p 臂与 let-7 家族；不改动原始行，差异仍按原始特征计算",
    ))
    n_precursor = int(annotation["precursor"].sum())
    n_guessed = int(annotation["guessed"].sum())
    dup_canonical = int(annotation["canonical"].duplicated().sum())
    bullets.append(
        f"miRNA：{len(inputs['genes'])} 个特征 × {len(inputs['samples'])} 样本（来源 {inputs['source']}）；"
        f"规范化后 {int(annotation['canonical'].nunique())} 个唯一 mature 名称，precursor {n_precursor} 个"
        + (f"、规则猜测臂 {n_guessed} 个" if n_guessed else "")
    )
    if dup_canonical:
        warnings.append(
            f"有 {dup_canonical} 行规范化后与其它行同名（例如 hsa-miR-let-7a-5p → hsa-let-7a-5p），"
            "差异分析仍按原始行计算，命名表里已保留别名"
        )
    metrics.update({"n_mirna": int(len(inputs["genes"])), "n_samples": int(len(inputs["samples"])),
                    "n_precursor": n_precursor, "n_guessed": n_guessed})

    if inputs["groups"] is None:
        warnings.append("miRNA 样本与样本表的分组列对不上，跳过差异分析")
    else:
        try:
            section = _mirna_de_output(ctx, params, outdir, inputs)
            sections.append(section)
            metrics.update(section.get("metrics") or {})
        except EngineError as error:
            warnings.append(f"miRNA 差异分析跳过：{error.message}")

    try:
        section = _mirna_target_output(ctx, params, outdir, inputs)
        sections.append(section)
        metrics.update(section.get("metrics") or {})
    except EngineError as error:
        warnings.append(f"miRNA 靶基因分析跳过：{error.message}")

    try:
        section = _cerna_output(ctx, params, outdir, inputs, strict=False)
        sections.append(section)
        metrics.update(section.get("metrics") or {})
    except EngineError as error:
        warnings.append(f"ceRNA 候选跳过：{error.message}")

    sections.append(_circ_expression_output(
        {"X": inputs["X"], "genes": inputs["genes"], "samples": inputs["samples"]}, outdir,
    ))
    for section in sections:
        tables.extend(section.get("tables") or [])
        data_out.extend(section.get("data_out") or [])
        bullets.extend(section.get("bullets") or [])
        metrics.update(section.get("metrics") or {})
        warnings.extend(section.get("warnings") or [])
        notes.extend(section.get("notes") or [])
    payload = _finish(KIND_TITLES["mirna"], sections, bullets=bullets, tables=tables, data_out=data_out)
    payload["warnings"] = _dedupe_list(warnings + (payload.get("warnings") or []))
    payload["notes"] = _dedupe_list(notes + (payload.get("notes") or []))
    payload["summary"]["metrics"].update(metrics)
    return payload


# ---------------------------------------------------------------------------
# lncRNA（分类 / 差异 / 共表达网络）
# ---------------------------------------------------------------------------

# 内置已知 lncRNA → 类别（无 GTF 时的启发式注释，来源：常见 lncRNA 文献命名）
_LNCRNA_KNOWN = {
    "MALAT1": "lincRNA", "NEAT1": "lincRNA", "XIST": "lincRNA", "HOTAIR": "lincRNA",
    "PVT1": "lincRNA", "MEG3": "lincRNA", "H19": "lincRNA", "TUG1": "lincRNA",
    "GAS5": "lincRNA", "DANCR": "lincRNA", "NORAD": "lincRNA", "HULC": "lincRNA",
    "UCA1": "lincRNA", "CCAT1": "lincRNA", "CCAT2": "lincRNA", "HOTTIP": "lincRNA",
    "BANCR": "lincRNA", "TINCR": "lincRNA", "HAR1A": "lincRNA", "TP53TG1": "lincRNA",
    "PRNCR1": "lincRNA", "PCAT1": "lincRNA", "PCAT6": "lincRNA", "PCAT7": "lincRNA",
    "ANRIL": "antisense", "KCNQ1OT1": "antisense", "SNHG1": "snoRNA_host",
    "SNHG3": "snoRNA_host", "SNHG5": "snoRNA_host", "SNHG6": "snoRNA_host",
    "SNHG7": "snoRNA_host", "SNHG8": "snoRNA_host", "SNHG12": "snoRNA_host",
    "SNHG15": "snoRNA_host", "SNHG16": "snoRNA_host", "SNHG17": "snoRNA_host",
    "SNHG20": "snoRNA_host", "SNORD3A": "snoRNA_host", "SNORD3B": "snoRNA_host",
    "ZFAS1": "antisense", "ZEB2-AS1": "antisense", "ZEB1-AS1": "antisense",
    "DLX6-AS1": "antisense", "FOXD3-AS1": "antisense", "PAX8-AS1": "antisense",
    "ZNFX1-AS1": "antisense", "HOXA11-AS": "antisense", "HOXD-AS1": "antisense",
    "CDKN2B-AS1": "antisense", "SPRY4-IT1": "intronic", "LINC00152": "lincRNA",
    "LINC00261": "lincRNA", "LINC00312": "lincRNA", "LINC00472": "lincRNA",
    "LINC00511": "lincRNA", "LINC00662": "lincRNA", "LINC00673": "lincRNA",
    "LINC00707": "lincRNA", "LINC00958": "lincRNA", "LINC01006": "lincRNA",
    "LINC01116": "lincRNA", "LINC01234": "lincRNA", "LINC01314": "lincRNA",
    "LINC01410": "lincRNA", "LINC01503": "lincRNA", "LINC01614": "lincRNA",
    "LINC01705": "lincRNA", "LINC01857": "lincRNA", "LINC01963": "lincRNA",
    "LINC02035": "lincRNA", "LINC02159": "lincRNA", "LINC02268": "lincRNA",
    "LINC02381": "lincRNA", "LINC02446": "lincRNA", "LINC02532": "lincRNA",
    "LINC02605": "lincRNA", "LINC02731": "lincRNA", "LINC02807": "lincRNA",
    "MIR31HG": "lincRNA", "PTCSC3": "lincRNA", "BISPR": "lincRNA",
}
_LNCRNA_KNOWN_UPPER = {str(k).upper() for k in _LNCRNA_KNOWN}


def classify_lncrna(name) -> tuple[str, str]:
    """按 ID 命名规则 + 内置已知名单给出 lncRNA 类别与依据（无 GTF 时是启发式）。"""
    raw = str(name).strip()
    upper = raw.upper()
    if not upper:
        return "unknown", "空名称"
    if upper.startswith("LINC"):
        return "lincRNA", "LINC 前缀"
    if re.search(r"-AS\d*$", upper):
        return "antisense", "-AS 后缀"
    if re.search(r"-IT\d*$", upper):
        return "intronic", "-IT 后缀"
    if re.search(r"-DT\d*$", upper):
        return "divergent", "-DT 后缀"
    if upper.startswith("SNHG"):
        return "snoRNA_host", "SNHG 前缀"
    if upper.startswith("SNORD") or upper.startswith("SNORA"):
        return "snoRNA_host", "snoRNA 前缀"
    if upper.startswith("PCAT"):
        return "lincRNA", "PCAT 前缀（前列腺癌相关转录本）"
    if re.search(r"-AS\d*$", upper) is None and upper in _LNCRNA_KNOWN_UPPER:
        return _LNCRNA_KNOWN_CLASS.get(upper, "lincRNA"), "内置已知 lncRNA 名单"
    if upper.startswith(("ENSG", "ENSMUSG")):
        return "ensembl_id", "Ensembl ID（需 GTF 判定 biotype）"
    if re.match(r"^CHR[0-9XYM]+[:\-_]", upper) or ":" in raw:
        return "genomic_region", "坐标型 ID"
    region = _parse_region(raw)
    if region is not None:
        return "genomic_region", "坐标型 ID"
    return "unknown", "无注释依据（需 GTF/参考注释）"


def _looks_lncrna(name) -> bool:
    class_name, _basis = classify_lncrna(name)
    return class_name not in ("unknown", "ensembl_id")


def _lncrna_annotation(genes, info=None) -> pd.DataFrame:
    rows = []
    lengths = None
    if isinstance(info, pd.DataFrame) and "length" in [str(c).lower() for c in info.columns]:
        col = next(c for c in info.columns if str(c).lower() == "length")
        lengths = pd.to_numeric(info[col], errors="coerce").tolist()
    for i, gene in enumerate(genes):
        class_name, basis = classify_lncrna(gene)
        region = _parse_region(str(gene))
        length = None
        if lengths is not None and i < len(lengths) and np.isfinite(lengths[i]):
            length = int(lengths[i])
        elif region is not None and region["end"] is not None:
            length = int(region["end"] - region["start"])
        rows.append({
            "gene": str(gene),
            "class": class_name,
            "basis": basis,
            "chrom": region["chrom"] if region else "",
            "start": region["start"] if region else None,
            "end": region["end"] if region else None,
            "length": length,
        })
    return pd.DataFrame(rows)


_LNCRNA_ANNOT_COLUMNS = [
    {"key": "gene", "label": "lncRNA", "type": "string"},
    {"key": "class", "label": "类别", "type": "string"},
    {"key": "basis", "label": "判定依据", "type": "string"},
    {"key": "chrom", "label": "染色体", "type": "string"},
    {"key": "start", "label": "起点", "type": "number"},
    {"key": "end", "label": "终点", "type": "number"},
    {"key": "length", "label": "长度(bp)", "type": "number"},
]


def _lncrna_inputs(ctx, params, *, warnings=None, notes=None) -> dict:
    inputs = _matrix_inputs(ctx, params, "lncrna", warnings=warnings, notes=notes)
    inputs["annotation"] = _lncrna_annotation(inputs["genes"], inputs.get("info"))
    return inputs


def _lncrna_class_output(ctx, params, outdir, inputs):
    annotation = inputs["annotation"]
    class_counts = annotation["class"].value_counts().to_dict()
    frame = annotation.sort_values(["class", "gene"]).reset_index(drop=True)
    tables = [table(
        "lncrna_class",
        "lncRNA 类别注释",
        _LNCRNA_ANNOT_COLUMNS,
        frame,
        note="按 ID 命名规则（LINC…/-AS1/-IT1/SNHG…）+ 内置已知 lncRNA 名单做启发式分类；"
             "提供 GTF 时可做精确 biotype 注释",
    )]
    sorted_counts = sorted(class_counts.items(), key=lambda item: -item[1])
    figures = [_fig_bar(
        [item[0] for item in sorted_counts],
        [item[1] for item in sorted_counts],
        outdir, "lncrna_class_bar", "lncRNA class distribution",
        y_label="Count", colors=palette(len(sorted_counts)),
    )]
    bullets = [
        f"lncRNA 分类：{len(frame)} 个特征 → " +
        "、".join(f"{key} {value}" for key, value in sorted_counts[:6])
    ]
    n_unknown = int((frame["class"] == "unknown").sum())
    warnings = []
    if n_unknown:
        warnings.append(f"{n_unknown} 个 lncRNA 名称无法按命名规则分类（类别 unknown），需要 GTF 才能精确注释")
    notes = ["无 GTF 时长度/biotype 未知；分类是命名启发式，不等同于 GENCODE biotype"]
    _safe_csv(frame, outdir, "lncrna_class", "lncRNA 类别表", [])
    return {"bullets": bullets, "metrics": {"n_lncrna": int(len(frame)),
                                            "classes": class_counts,
                                            "n_unknown": n_unknown},
            "tables": tables, "figures": [f for f in figures if f],
            "warnings": warnings, "notes": notes, "data_out": []}


def _lncrna_de_output(ctx, params, outdir, inputs):
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    de, case, control = _diff_frame(inputs, params, label="lncRNA", warnings=warnings, notes=notes)
    annotation = inputs["annotation"].set_index("gene")
    frame = de.copy()
    frame = frame.join(annotation[["class", "basis"]], on="gene")
    tables = [_de_table_obj(
        frame,
        "lncrna_de",
        "lncRNA 差异分析",
        note="Welch t + BH；class 来自命名启发式注释（无 GTF 时）",
        lead_columns=[{"key": "class", "label": "类别", "type": "string"}],
    )]
    figures = [_fig_volcano(
        frame, outdir, "lncrna_volcano", "lncRNA differential expression",
        xcol="logFC", pcol="adj.P.Val", labelcol="gene",
        padj_cut=float(params.get("padj") or _PADJ_DEFAULT),
        x_cut=float(params.get("lfc") or _LFC_DEFAULT),
        up_label="Up in case", down_label="Down in case",
    )]
    top = frame[frame["regulation"].isin(["up", "down"])].head(20)
    if len(top) >= 3:
        index = {gene: i for i, gene in enumerate(inputs["genes"])}
        idx = [index[gene] for gene in top["gene"] if gene in index]
        if idx:
            matrix = pd.DataFrame(
                np.asarray(inputs["X"])[:, idx].T,
                index=[inputs["genes"][i] for i in idx],
                columns=inputs["samples"],
            )
            heat = _fig_heat(matrix, outdir, "lncrna_de_heat", "Top lncRNAs (row z-score)",
                             x_label="Samples", y_label="lncRNA")
            if heat:
                figures.append(heat)
    bullet, metrics = _de_summary(frame, "lncRNA", case, control,
                                  float(params.get("padj") or _PADJ_DEFAULT))
    _safe_csv(frame, outdir, "lncrna_de", "lncRNA 差异表", data_out)
    notes.append("差异检验为 Welch t + BH（与 bulk DE 共用 bg_common 统计核心）")
    return {"bullets": [bullet], "metrics": metrics, "tables": tables,
            "figures": [f for f in figures if f], "warnings": warnings, "notes": notes,
            "data_out": data_out}


def _lncrna_network_output(ctx, params, outdir, inputs):
    """lncRNA-mRNA 共表达网络：|r|>0.6 且 BH FDR<0.05。"""
    warnings: list[str] = []
    notes: list[str] = []
    data_out: list[dict] = []
    min_r = abs(float(params.get("min_abs_r") or 0.6))
    fdr_cut = float(params.get("fdr") or 0.05)
    ctx_genes = [str(g) for g in (ctx.get("genes") or [])]
    ctx_samples = [str(s) for s in (ctx.get("samples") or [])]
    ctx_X = ctx.get("X")
    if not isinstance(ctx_X, (np.ndarray, list)) or not ctx_genes:
        raise EngineError("prerequisite", "ctx 里没有 bulk mRNA 矩阵（X/genes），无法做 lncRNA–mRNA 共表达")
    ctx_X = np.asarray(ctx_X, dtype=np.float64)
    lnc_rows, ctx_rows, common = _align_samples(inputs["samples"], ctx_samples)
    if len(common) < 5:
        raise EngineError(
            "prerequisite",
            f"lncRNA 与 bulk 样本只有 {len(common)} 个交集（共表达至少需要 5 个）",
        )
    A = np.asarray(inputs["X"], dtype=np.float64)[lnc_rows, :]
    B = ctx_X[ctx_rows, :]
    lnc_sel = _top_variable(A, 150)
    mrna_sel = _top_variable(B, 300)
    A = A[:, lnc_sel]
    B = B[:, mrna_sel]
    r, p = _pearson_matrix(A, B)
    q = _bh_matrix(p)
    mask = (np.abs(r) >= min_r) & (q <= fdr_cut)
    lnc_names = [inputs["genes"][i] for i in lnc_sel]
    mrna_names = [ctx_genes[i] for i in mrna_sel]
    annotation = inputs["annotation"].set_index("gene")
    la, lb = np.where(mask)
    rows = []
    for a, b in zip(la.tolist(), lb.tolist()):
        lnc = lnc_names[a]
        rows.append({
            "lncrna": lnc,
            "mrna": mrna_names[b],
            "r": round(float(r[a, b]), 4),
            "p": float(p[a, b]),
            "padj": float(q[a, b]),
            "direction": "positive" if r[a, b] >= 0 else "negative",
            "lncrna_class": str(annotation["class"].get(lnc, "unknown")),
        })
    if not rows:
        warnings.append(
            f"没有 |r|>{min_r:g} 且 FDR<{fdr_cut:g} 的 lncRNA–mRNA 共表达边"
            f"（共同样本 {len(common)} 个）；样本少或信号弱时常见"
        )
        return {"bullets": [], "metrics": {"n_edges": 0, "n_common_samples": len(common)},
                "tables": [], "figures": [], "warnings": warnings, "notes": notes,
                "data_out": []}
    frame = pd.DataFrame(rows).sort_values(["padj", "r"], key=lambda s: s.abs() if s.name == "r" else s,
                                           ascending=[True, False]).reset_index(drop=True)
    tables = [table(
        "lncrna_network_edges",
        "lncRNA–mRNA 共表达边",
        [
            {"key": "lncrna", "label": "lncRNA", "type": "string"},
            {"key": "mrna", "label": "mRNA", "type": "string"},
            {"key": "r", "label": "相关系数 r", "type": "number"},
            {"key": "p", "label": "P 值", "type": "number"},
            {"key": "padj", "label": "FDR", "type": "number"},
            {"key": "direction", "label": "方向", "type": "string"},
            {"key": "lncrna_class", "label": "lncRNA 类别", "type": "string"},
        ],
        frame,
        note=f"共表达筛选：|r|≥{min_r:g} 且 BH FDR≤{fdr_cut:g}；Pearson 相关，跨 {len(common)} 个共同样本",
    )]
    figures = []
    network = _fig_network(frame.head(120), outdir, "lncrna_network", "lncRNA-mRNA co-expression network")
    if network:
        figures.append(network)
    top = frame.head(15)
    bar = _fig_barh(
        (top["lncrna"].astype(str) + " – " + top["mrna"].astype(str)).tolist(),
        top["r"].abs().tolist(), outdir, "lncrna_network_top",
        "Top lncRNA-mRNA co-expression pairs", x_label="|r|", max_items=15,
    )
    if bar:
        figures.append(bar)
    metrics = {
        "n_edges": int(len(frame)),
        "n_lncrna": int(frame["lncrna"].nunique()),
        "n_mrna": int(frame["mrna"].nunique()),
        "n_common_samples": int(len(common)),
        "n_positive": int((frame["r"] > 0).sum()),
    }
    bullets = [
        f"共表达网络：{metrics['n_edges']} 条边（{metrics['n_lncrna']} lncRNA / {metrics['n_mrna']} mRNA，"
        f"正相关 {metrics['n_positive']}），共同样本 {len(common)} 个"
    ]
    _safe_csv(frame, outdir, "lncrna_network_edges", "lncRNA 共表达边表", data_out)
    notes.append("共表达是相关性网络（非因果/非调控方向），需实验验证；网络图节点数/边数做了上限裁剪")
    return {"bullets": bullets, "metrics": metrics, "tables": tables,
            "figures": [f for f in figures if f], "warnings": warnings, "notes": notes,
            "data_out": data_out}


def _step_lncrna(ctx, params, outdir):
    warnings: list[str] = []
    notes: list[str] = []
    tables: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []
    metrics: dict = {}
    sections: list[dict] = []

    inputs = _lncrna_inputs(ctx, params, warnings=warnings, notes=notes)
    bullets.append(
        f"lncRNA：{len(inputs['genes'])} 个 lncRNA × {len(inputs['samples'])} 样本（来源 {inputs['source']}）"
    )

    class_section = _lncrna_class_output(ctx, params, outdir, inputs)
    sections.append(class_section)
    metrics.update(class_section.get("metrics") or {})

    if inputs["groups"] is None:
        warnings.append("lncRNA 样本与样本表的分组列对不上，跳过差异分析")
    else:
        try:
            section = _lncrna_de_output(ctx, params, outdir, inputs)
            sections.append(section)
            metrics.update(section.get("metrics") or {})
        except EngineError as error:
            warnings.append(f"lncRNA 差异分析跳过：{error.message}")

    try:
        section = _lncrna_network_output(ctx, params, outdir, inputs)
        sections.append(section)
        metrics.update(section.get("metrics") or {})
    except EngineError as error:
        warnings.append(f"lncRNA 共表达网络跳过：{error.message}")

    if not any(section.get("figures") for section in sections):
        mean_expr = np.nanmean(np.asarray(inputs["X"], dtype=np.float64), axis=0)
        order = np.argsort(-np.nan_to_num(mean_expr))[:12]
        figure = _fig_barh([inputs["genes"][i] for i in order], [float(mean_expr[i]) for i in order],
                           outdir, "lncrna_expression", "Top expressed lncRNAs",
                           x_label="Mean expression")
        if figure:
            sections.append({"figures": [figure], "bullets": [], "metrics": {}, "tables": [],
                             "warnings": [], "notes": [], "data_out": []})

    for section in sections:
        tables.extend(section.get("tables") or [])
        data_out.extend(section.get("data_out") or [])
        bullets.extend(section.get("bullets") or [])
        metrics.update(section.get("metrics") or {})
        warnings.extend(section.get("warnings") or [])
        notes.extend(section.get("notes") or [])
    payload = _finish(KIND_TITLES["lncrna"], sections, bullets=bullets, tables=tables, data_out=data_out)
    payload["warnings"] = _dedupe_list(warnings + (payload.get("warnings") or []))
    payload["notes"] = _dedupe_list(notes + (payload.get("notes") or []))
    payload["summary"]["metrics"].update(metrics)
    return payload


# === PART3 ===



# === PART3 ===
# 说明：本段补齐 PART1/PART2 依赖的公共件（circ 输入包装、结果汇总、注册表），
# 并实现 fusion gene 与可变剪接两个模态。

FUSION_ONCOGENES = """
ABL1 ABL2 AKT1 AKT2 AKT3 ALK APC AR ARAF ATF1 BCL2 BCL3 BCL6 BCL9 BCL10 BCL11A BCR BRAF
BRCA1 BRCA2 CCND1 CCND2 CCND3 CDK4 CDK6 CDX2 CREB1 CRKL CSF1R CTNNB1 DDIT3 DDR2 E2F1 EGFR
EIF4A2 ELL EML4 ERBB2 ERBB3 ERBB4 ERG ETV1 ETV4 ETV5 ETV6 EWSR1 EZH2 FGFR1 FGFR2 FGFR3
FGFR4 FLI1 FOS FOXA1 FOXO1 FUS GATA1 GATA2 GLI1 GLI2 GNA11 GNAQ GNAS HMGA2 HRAS IDH1 IDH2
JAK1 JAK2 JAK3 KIT KRAS LMO1 LMO2 LYN MDM2 MECOM MET MITF MLLT1 MLLT3 MLLT10 MLLT11 MYB MYC
MYCN MYD88 NCOA1 NCOA2 NCOA4 NF1 NF2 NFKB2 NKX2-1 NOTCH1 NOTCH2 NPM1 NRAS NR4A3 NTRK1 NTRK2
NTRK3 PAX3 PAX5 PAX7 PAX8 PBX1 PDGFB PDGFRA PDGFRB PIK3CA PIK3R1 PLAG1 PPARG PRCC PRKACA
PRKAR1A PTEN RAF1 RARA RB1 REL RET RHOA ROS1 RUNX1 RUNX1T1 SET SFPQ SMAD2 SMAD4 SMO SS18
STAT3 STAT5A STAT5B STAT6 TAL1 TAL2 TCF3 TCF7L2 TET1 TET2 TFEB TGFBR2 TLX1 TLX3 TMPRSS2
TP53 TPM3 TPM4 TRAF7 TSC1 TSC2 USP6 VHL WT1 YAP1 ZNF384
""".split()

FUSION_KINASES = """
ABL1 ABL2 AKT1 AKT2 AKT3 ALK ARAF AURKA AURKB AXL BRAF BTK CDK1 CDK2 CDK4 CDK6 CDK7 CDK9
CHEK1 CHEK2 CSF1R CSK DDR1 DDR2 EGFR EPHA2 EPHB2 ERBB2 ERBB3 ERBB4 FER FGFR1 FGFR2 FGFR3
FGFR4 FLT1 FLT3 FLT4 FYN GSK3B HCK IGF1R IKBKB INSR IRAK4 ITK JAK1 JAK2 JAK3 KDR KIT LCK
LYN MAP2K1 MAP2K2 MAPK1 MAPK3 MAPK8 MAPK14 MET MERTK MST1R MTOR NEK2 NTRK1 NTRK2 NTRK3 PAK1
PDGFRA PDGFRB PIK3CA PIK3CB PIK3CD PIK3CG PIM1 PIM2 PIM3 PLK1 PLK2 PLK3 PLK4 PRKACA PRKCA
PRKCB PRKCD PRKDC RAF1 RET ROCK1 ROCK2 ROS1 RPS6KA1 SGK1 SRC STK11 SYK TBK1 TEC TGFBR1
TGFBR2 TNIK TYK2 WEE1 ZAP70
""".split()

_FUSION_COLUMN_ALIASES = {
    "name": ["fusionname", "fusion_name", "fusion", "name", "fusionevent"],
    "gene1": ["gene1", "gene_1", "gene5", "gene3", "leftgene", "left_gene", "g1"],
    "gene2": ["gene2", "gene_2", "gene3p", "gene3", "rightgene", "right_gene", "g2"],
    "junction": ["junctionreadcount", "junction_reads", "junctionreads", "split_reads1", "splitreads1",
                 "read_count", "readcount", "junctionspan", "junc_read_count"],
    "spanning": ["spanningfragcount", "spanningfrags", "spanning_frag_count", "spanningreadcount",
                 "discordant_mates", "split_reads2", "spanning_reads"],
    "sample": ["sample", "sampleid", "sample_id", "tumor_sample", "tumorsample", "case"],
    "type": ["type", "fusiontype", "fusion_type", "svtype"],
    "frame": ["readingframe", "reading_frame", "fusion_cds", "frame", "inframe"],
    "breakpoint1": ["breakpoint1", "bp1", "site1", "leftbreakpoint"],
    "breakpoint2": ["breakpoint2", "bp2", "site2", "rightbreakpoint"],
}


_LNCRNA_KNOWN_CLASS = {str(key).upper(): value for key, value in _LNCRNA_KNOWN.items()}


def _parse_region(text):
    """解析 chr1:12345-12789:+ / chr1_12345_12789 / 1:12345-12789 形式的坐标 ID。"""
    raw = str(text or "")
    match = re.search(r"chr([0-9XYMxy]+)[:\-_](\d+)[-_:](\d+)", raw, re.I)
    if not match:
        match = re.search(r"(?<![0-9])([0-9]{1,2}|[XYM])[:_](\d+)[-_:](\d+)", raw, re.I)
    if not match:
        return None
    chrom = match.group(1).upper()
    if not chrom.startswith("CHR"):
        chrom = "chr" + chrom
    start = int(match.group(2))
    end = int(match.group(3))
    if start > end:
        start, end = end, start
    strand = ""
    tail = raw[match.end():match.end() + 3]
    if "+" in tail:
        strand = "+"
    elif "-" in tail:
        strand = "-"
    return {"chrom": chrom, "start": start, "end": end, "strand": strand,
            "length": end - start, "label": f"{chrom}:{start}-{end}"}


def _column_norm(name) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lstrip("#").lower())


def _pick_columns(frame: pd.DataFrame, aliases: dict) -> dict:
    mapping = {}
    normalized = {_column_norm(col): col for col in frame.columns}
    for slot, candidates in aliases.items():
        for candidate in candidates:
            key = _column_norm(candidate)
            if key in normalized:
                mapping[slot] = normalized[key]
                break
    return mapping


def _read_module_table(ctx, kind: str, params) -> pd.DataFrame:
    """fusion / splicing 的表格来源：module_path > ctx.modules[kind].table。"""
    params = _params(params)
    warnings: list[str] = []
    path = str(params.get("module_path") or "").strip()
    if path:
        frame = _read_table_file(path, warnings)
        if isinstance(frame, pd.DataFrame) and not frame.empty:
            return frame
    mod = _module(ctx, kind)
    if isinstance(mod, dict) and isinstance(mod.get("table"), pd.DataFrame) and not mod["table"].empty:
        return mod["table"].copy()
    if isinstance(mod, dict) and isinstance(mod.get("info"), pd.DataFrame) and not mod["info"].empty:
        return mod["info"].copy()
    if isinstance(mod, dict) and mod.get("X") is not None:
        # 矩阵型附件（如 rMATS PSI/IncLevel 也被解析成矩阵）：转成 行=特征 的表格
        matrix = np.asarray(mod.get("X"), dtype=float)
        genes = [str(g) for g in (mod.get("genes") or [])]
        samples = [str(s) for s in (mod.get("samples") or [])]
        if matrix.ndim == 2:
            if len(genes) == matrix.shape[1] and len(samples) == matrix.shape[0]:
                return pd.DataFrame(matrix.T, index=genes, columns=samples)
            if len(genes) == matrix.shape[0] and len(samples) == matrix.shape[1]:
                return pd.DataFrame(matrix, index=genes, columns=samples)
        raise EngineError("unsupported", f"{KIND_LABELS[kind]} 附件是矩阵但没有行列名，无法解析事件表；请直接挂载原始表格文件")
    raise _missing_module_error(kind, ctx)


def _normalize_fusion(frame: pd.DataFrame, warnings: list[str]) -> pd.DataFrame:
    frame = frame.copy()
    frame.columns = [str(c) for c in frame.columns]
    mapping = _pick_columns(frame, _FUSION_COLUMN_ALIASES)
    out = pd.DataFrame(index=frame.index)
    if "name" in mapping:
        out["fusion"] = frame[mapping["name"]].astype(str)
    if "gene1" in mapping:
        out["gene1"] = frame[mapping["gene1"]].astype(str)
    if "gene2" in mapping:
        out["gene2"] = frame[mapping["gene2"]].astype(str)
    if "fusion" in out.columns and ("gene1" not in out.columns or "gene2" not in out.columns):
        parts = out["fusion"].astype(str).str.replace("--", "::", regex=False).str.replace("_", "::", regex=False).str.split("::")
        out["gene1"] = parts.str[0]
        out["gene2"] = parts.str[1] if parts.str.len().max() > 1 else ""
    if "gene1" not in out.columns and "fusion" not in out.columns:
        # 退化：可能是长表（sample/gene1/gene2 各一行）
        raise EngineError("unsupported", "融合表里找不到融合名/基因列（支持 STAR-Fusion、Arriba 或 gene1+gene2 列）")
    if "fusion" not in out.columns:
        out["fusion"] = out["gene1"].astype(str) + "--" + out["gene2"].astype(str)
    for slot, default in (("junction", 0), ("spanning", 0)):
        if slot in mapping:
            values = pd.to_numeric(frame[mapping[slot]].astype(str).str.split(",").str[0], errors="coerce")
            out["junction_reads" if slot == "junction" else "spanning_reads"] = values.fillna(default)
        else:
            out["junction_reads" if slot == "junction" else "spanning_reads"] = default
    if "sample" in mapping:
        out["sample"] = frame[mapping["sample"]].astype(str)
    if "type" in mapping:
        out["type"] = frame[mapping["type"]].astype(str)
    if "frame" in mapping:
        out["frame"] = frame[mapping["frame"]].astype(str)
    if "breakpoint1" in mapping:
        out["breakpoint1"] = frame[mapping["breakpoint1"]].astype(str)
    if "breakpoint2" in mapping:
        out["breakpoint2"] = frame[mapping["breakpoint2"]].astype(str)
    out["gene1"] = out["gene1"].astype(str).str.split("^").str[0].str.strip()
    out["gene2"] = out["gene2"].astype(str).str.split("^").str[0].str.strip()
    out["reads"] = pd.to_numeric(out["junction_reads"], errors="coerce").fillna(0) + \
        pd.to_numeric(out["spanning_reads"], errors="coerce").fillna(0)
    if len(frame.columns) <= 6 and "sample" not in out.columns:
        warnings.append("融合表没有样本列：只做融合概览，无法做样本×融合矩阵")
    return out


def _fusion_flags(frame: pd.DataFrame) -> pd.DataFrame:
    oncogenes = {g.upper() for g in FUSION_ONCOGENES}
    kinases = {g.upper() for g in FUSION_KINASES}
    frame = frame.copy()
    frame["oncogene"] = frame["gene1"].str.upper().isin(oncogenes) | frame["gene2"].str.upper().isin(oncogenes)
    frame["kinase"] = frame["gene1"].str.upper().isin(kinases) | frame["gene2"].str.upper().isin(kinases)
    frame["known_driver"] = frame["oncogene"] | frame["kinase"]
    return frame


def _step_fusion(ctx, params, outdir):
    warnings: list[str] = []
    notes: list[str] = []
    frame = _normalize_fusion(_read_module_table(ctx, "fusion", params), warnings)
    frame = _fusion_flags(frame)
    sections: list[dict] = []
    tables: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []
    metrics: dict = {}
    min_reads = int(_params(params).get("min_reads") or 2)
    filtered = frame[frame["reads"] >= min_reads].copy()
    if filtered.empty:
        filtered = frame.copy()
        warnings.append(f"没有融合读到 {min_reads} 条的记录，展示全部 {len(frame)} 条")
    filtered = filtered.sort_values("reads", ascending=False).reset_index(drop=True)
    tables.append(table("fusion_table", "融合基因表", list(filtered.columns)[:12], filtered.head(500),
                        note=f"共 {len(filtered)} 条；阈值 reads≥{min_reads}"))
    path = export_csv(filtered, outdir, "fusion_normalized")
    data_out.append(bc.data_file("table_csv", path, "标准化融合表"))
    n_pairs = filtered["fusion"].nunique()
    n_drivers = int(filtered.loc[filtered["known_driver"], "fusion"].nunique())
    bullets.append(f"融合基因：{len(filtered)} 条记录 / {n_pairs} 个融合事件；已知驱动型（癌基因/激酶）{n_drivers} 个")
    metrics.update({"n_fusions": int(len(filtered)), "n_pairs": int(n_pairs), "n_drivers": n_drivers})
    counts = filtered.groupby("fusion")["reads"].sum().sort_values(ascending=False)
    figure = _fig_barh(list(counts.index)[:15], list(counts.values)[:15], outdir, "fusion_top",
                       "Top fusion events", x_label="Supporting reads")
    if figure:
        sections.append({"figures": [figure]})
    if "sample" in filtered.columns and filtered["sample"].nunique() > 1:
        matrix = filtered.pivot_table(index="sample", columns="fusion", values="reads", aggfunc="sum", fill_value=0)
        matrix = matrix.loc[:, matrix.sum(axis=0).sort_values(ascending=False).index[:40]]
        tables.append(table("fusion_matrix", "样本×融合矩阵", [{"key": "sample", "label": "Sample", "type": "string"}] +
                            [{"key": str(c), "label": str(c)[:40], "type": "number"} for c in matrix.columns],
                            matrix.reset_index().rename(columns={"index": "sample"}).head(200)))
        heat = _fig_heat(matrix.T, outdir, "fusion_heat", "Fusion support (samples)")
        if heat:
            sections.append({"figures": [heat]})
        metrics["n_samples"] = int(matrix.shape[0])
        recurring = matrix.astype(bool).sum(axis=0).sort_values(ascending=False)
        bullets.append(f"复发性融合（≥2 样本）：{int((recurring >= 2).sum())} 个")
    drivers = filtered[filtered["known_driver"]]
    if not drivers.empty:
        tables.append(table("fusion_drivers", "癌基因/激酶融合", list(drivers.columns)[:11], drivers.head(200)))
        notes.append("驱动型判定基于内置癌基因/激酶名单（约 160 + 105 个），非数据库注释")
    payload = _finish(KIND_TITLES["fusion"], sections, bullets=bullets, tables=tables, data_out=data_out, metrics=metrics)
    payload["warnings"] = _dedupe_list(warnings + (payload.get("warnings") or []))
    payload["notes"] = _dedupe_list(notes + (payload.get("notes") or []))
    payload["store"] = {"ncrna_fusion": {"n_fusions": int(len(filtered)), "n_pairs": int(n_pairs)}}
    return payload


SPLICE_EVENT_ORDER = ["SE", "RI", "A5SS", "A3SS", "MXE", "AF", "AL", "AP", "AT", "OTHER"]


def _splice_event_type(row, columns) -> str:
    for column in columns:
        text = str(row.get(column) or "")
        upper = text.upper()
        for code in ("SE", "RI", "A5SS", "A3SS", "MXE", "AF", "AL", "AP", "AT"):
            if upper.startswith(code) or f":{code}:" in upper or f"_{code}_" in upper:
                return code
    return "OTHER"


def _step_splicing(ctx, params, outdir):
    warnings: list[str] = []
    notes: list[str] = []
    params = _params(params)
    frame = _read_module_table(ctx, "splicing", params)
    frame = frame.copy()
    frame.columns = [str(c) for c in frame.columns]
    normalized = {_column_norm(c): c for c in frame.columns}
    padj_cut = float(params.get("padj") or _PADJ_DEFAULT)
    psidiff_cut = float(params.get("psidiff") or params.get("delta_psi") or 0.1)
    case = _match_group(params.get("case"), list({str(v) for v in (ctx.get("groups") or [])})) if ctx.get("groups") else None
    tables: list[dict] = []
    figures: list[dict] = []
    data_out: list[dict] = []
    bullets: list[str] = []
    metrics: dict = {}

    inc1 = normalized.get("inclevel1") or normalized.get("psilevel1")
    inc2 = normalized.get("inclevel2") or normalized.get("psilevel2")
    is_rmats = bool(inc1 and inc2)

    if is_rmats:
        def _parse(value):
            text = str(value).replace(";", ",")
            parts = [p.strip() for p in text.split(",") if p.strip() != ""]
            vals = []
            for part in parts:
                try:
                    vals.append(float(part))
                except ValueError:
                    continue
            return np.array(vals, dtype=float)

        psi1 = frame[inc1].map(_parse)
        psi2 = frame[inc2].map(_parse)
        n1 = int(np.median([len(v) for v in psi1 if len(v)])) if len(psi1) else 0
        n2 = int(np.median([len(v) for v in psi2 if len(v)])) if len(psi2) else 0
        rows = []
        for idx in frame.index:
            a = psi1.loc[idx]
            b = psi2.loc[idx]
            if len(a) == 0 or len(b) == 0:
                continue
            mean1 = float(np.nanmean(a))
            mean2 = float(np.nanmean(b))
            dpsi = mean1 - mean2
            pvalue = 1.0
            if len(a) >= 2 and len(b) >= 2:
                try:
                    from scipy import stats

                    pvalue = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
                except Exception:  # noqa: BLE001
                    pvalue = 1.0
            row = {
                "event": str(frame.loc[idx, normalized.get("eventid", normalized.get("id", frame.columns[0]))])
                if normalized.get("eventid") or normalized.get("id") else f"event_{idx}",
                "type": _splice_event_type(frame.loc[idx], [c for c in frame.columns if _column_norm(c) in ("type", "eventtype", "splicetype", "id", "eventid")]),
                "gene": str(frame.loc[idx, normalized["geneid"]]) if "geneid" in normalized else "",
                "psi_group1": mean1,
                "psi_group2": mean2,
                "dPSI": dpsi,
                "pvalue": pvalue,
                "n1": int(len(a)),
                "n2": int(len(b)),
            }
            for extra in ("ijcsample1", "sjcsample1", "ijcsample2", "sjcsample2"):
                if extra in normalized:
                    row[extra.replace("ijcsample", "IJC_").replace("sjcsample", "SJC_")] = frame.loc[idx, normalized[extra]]
            rows.append(row)
        de = pd.DataFrame(rows)
        if de.empty:
            raise EngineError("prerequisite", "rMATS 表里没有可解析的 IncLevel 值")
        de["p.adjust"] = _bh_matrix(de["pvalue"].values.reshape(1, -1))[0] if len(de) else de["pvalue"]
        de = de.sort_values("p.adjust").reset_index(drop=True)
        de["significant"] = (de["p.adjust"] <= padj_cut) & (de["dPSI"].abs() >= psidiff_cut)
        notes.append(f"rMATS 解析：IncLevel1/2 分别 {n1}/{n2} 个重复；Welch t + BH（rMATS 官方用其自身模型）")
        bullets.append(f"rMATS：{len(de)} 个剪接事件，显著（FDR≤{padj_cut} 且 |dPSI|≥{psidiff_cut}）{int(de['significant'].sum())} 个")
    else:
        # PSI 矩阵：行=事件，列=样本
        sample_like = [c for c in frame.columns if pd.api.types.is_numeric_dtype(frame[c])]
        if len(sample_like) < 4:
            raise EngineError(
                "unsupported",
                "剪接表既不是 rMATS（需要 IncLevel1/IncLevel2 列），也不是 PSI 矩阵（需要 ≥4 个数值样本列）",
            )
        matrix = frame[sample_like].apply(pd.to_numeric, errors="coerce")
        groups = _resolve_groups(ctx, _module(ctx, "splicing"), sample_like, params)
        labels = groups["labels"] if groups else None
        if labels is None:
            split = max(1, len(sample_like) // 2)
            labels = ["group1"] * split + ["group2"] * (len(sample_like) - split)
            warnings.append("没有分组信息，PSI 矩阵按前后各半分组（建议在「样本」里设置分组列）")
        labels = np.asarray(labels)
        unique = [g for g in dict.fromkeys(labels.tolist())]
        if len(unique) < 2:
            raise EngineError("prerequisite", "PSI 矩阵只有一个分组，无法做组间剪接差异")
        case_label = _match_group(params.get("case"), unique) or unique[-1]
        control_label = _match_group(params.get("control"), unique) or next(g for g in unique if g != case_label)
        rows = []
        for idx, row in matrix.iterrows():
            a = row[labels == case_label].values.astype(float)
            b = row[labels == control_label].values.astype(float)
            a = a[np.isfinite(a)]
            b = b[np.isfinite(b)]
            if len(a) < 2 or len(b) < 2:
                continue
            pvalue = 1.0
            try:
                from scipy import stats

                pvalue = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
            except Exception:  # noqa: BLE001
                pvalue = 1.0
            rows.append({"event": str(idx), "type": _splice_event_type({"event": str(idx), "id": str(idx)}, ["event", "id"]),
                         "gene": "", "psi_group1": float(np.mean(a)), "psi_group2": float(np.mean(b)),
                         "dPSI": float(np.mean(a) - np.mean(b)), "pvalue": pvalue, "n1": len(a), "n2": len(b)})
        de = pd.DataFrame(rows)
        if de.empty:
            raise EngineError("prerequisite", "PSI 矩阵里没有有效事件（每行需 ≥2 个非缺失值/组）")
        de["p.adjust"] = _bh_matrix(de["pvalue"].values.reshape(1, -1))[0]
        de = de.sort_values("p.adjust").reset_index(drop=True)
        de["significant"] = (de["p.adjust"] <= padj_cut) & (de["dPSI"].abs() >= psidiff_cut)
        notes.append(f"PSI 矩阵解析：{case_label} vs {control_label}（Welch t + BH）")
        bullets.append(f"PSI：{len(de)} 个事件，显著 {int(de['significant'].sum())} 个（|dPSI|≥{psidiff_cut}）")

    tables.append(_de_table_obj(de.head(500), "splicing", "剪接事件差异", note="dPSI = group1 − group2"))
    path = export_csv(de, outdir, "splicing_dpsi")
    data_out.append(bc.data_file("table_csv", path, "剪接差异全表"))
    # 火山图
    significant = de["significant"].values
    fig, ax = bc.new_fig(7.6, 5.8)
    neglog = -np.log10(np.clip(de["p.adjust"].astype(float).values, 1e-300, 1))
    ax.scatter(de["dPSI"][~significant], neglog[~significant], s=12, color="#b9bfc8", alpha=0.6, label="ns")
    ax.scatter(de["dPSI"][significant], neglog[significant], s=14, color="#c0392b", alpha=0.85,
               label=f"significant ({int(significant.sum())})")
    ax.axvline(psidiff_cut, color="#8a919c", linestyle="--", linewidth=0.8)
    ax.axvline(-psidiff_cut, color="#8a919c", linestyle="--", linewidth=0.8)
    ax.set_xlabel("dPSI")
    ax.set_ylabel("-log10(adjusted P)")
    ax.set_title("Splicing differential events")
    ax.legend(fontsize=8)
    figures.append(save_fig(fig, outdir, "splicing_volcano.png", "Splicing volcano"))
    # 事件类型统计
    if "type" in de.columns and de["type"].astype(str).str.len().gt(0).any():
        counts = de["type"].astype(str).str.upper().value_counts()
        counts = counts.reindex([c for c in SPLICE_EVENT_ORDER if c in counts.index]).dropna()
        names = list(counts.index)
        values = [int(v) for v in counts.values]
        figure = _fig_bar(names, values, outdir, "splicing_event_type", "Splicing event types", y_label="Events")
        if figure:
            figures.append(figure)
        tables.append(table("splicing_types", "事件类型统计",
                            [{"key": "type", "label": "类型", "type": "string"},
                             {"key": "count", "label": "事件数", "type": "number"}],
                            [{"type": k, "count": v} for k, v in zip(names, values)]))
        bullets.append("事件类型：" + "、".join(f"{k}={v}" for k, v in zip(names, values)))
    # top 事件热图（若有 PSI 矩阵）
    if not is_rmats:
        top = de.head(30)["event"].astype(str).tolist()
        sub = matrix.loc[[idx for idx in matrix.index.astype(str) if idx in set(top)]]
        if not sub.empty:
            heat = _fig_heat(sub, outdir, "splicing_heat", "Top differential splicing events")
            if heat:
                figures.append(heat)
    metrics.update({"n_events": int(len(de)), "n_significant": int(de["significant"].sum()),
                    "n_sig_up": int(((de["significant"]) & (de["dPSI"] > 0)).sum()),
                    "n_sig_down": int(((de["significant"]) & (de["dPSI"] < 0)).sum())})
    payload = _finish(KIND_TITLES["splicing"], [], bullets=bullets, tables=tables, data_out=data_out, metrics=metrics)
    payload["figures"] = figures
    payload["warnings"] = _dedupe_list(warnings + (payload.get("warnings") or []))
    payload["notes"] = _dedupe_list(notes + (payload.get("notes") or []))
    return payload


# ---------------------------------------------------------------------------
# 公共件（PART1/PART2 依赖）
# ---------------------------------------------------------------------------

def _circ_inputs(ctx, params, *, warnings=None, notes=None):
    return _matrix_inputs(ctx, params, "circ", warnings=warnings, notes=notes)


def _dedupe_list(items) -> list[str]:
    seen = set()
    out = []
    for item in items or []:
        text = str(item)
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return out


def _finish(title: str, sections, *, bullets=None, tables=None, data_out=None, metrics=None) -> dict:
    """把若干 section 汇总成 §2.5 结果信封。"""
    figures: list[dict] = []
    table_list: list[dict] = list(tables or [])
    bullet_list: list[str] = list(bullets or [])
    warning_list: list[str] = []
    note_list: list[str] = []
    data_list: list[dict] = list(data_out or [])
    metric: dict = dict(metrics or {})
    for section in sections or []:
        if not isinstance(section, dict):
            continue
        for item in section.get("figures") or []:
            if item and item not in figures:
                figures.append(item)
        for item in section.get("tables") or []:
            if item and item not in table_list:
                table_list.append(item)
        for item in section.get("bullets") or []:
            if item and item not in bullet_list:
                bullet_list.append(item)
        for item in section.get("warnings") or []:
            warning_list.append(str(item))
        for item in section.get("notes") or []:
            note_list.append(str(item))
        for item in section.get("data_out") or []:
            if item and item not in data_list:
                data_list.append(item)
        try:
            metric.update(section.get("metrics") or {})
        except Exception:  # noqa: BLE001
            pass
    return result(title=title, bullets=_dedupe_list(bullet_list), metrics=metric, tables=table_list,
                  figures=figures, warnings=_dedupe_list(warning_list), notes=_dedupe_list(note_list),
                  data_out=data_list)


# ---------------------------------------------------------------------------
# 步骤调度与注册
# ---------------------------------------------------------------------------

STEP_PARAMS = [
    {"key": "kind", "label": "数据类型", "type": "select", "default": "circ",
     "options": [{"value": "circ", "label": "circRNA"}, {"value": "mirna", "label": "miRNA"},
                 {"value": "lncrna", "label": "lncRNA"}, {"value": "fusion", "label": "融合基因"},
                 {"value": "splicing", "label": "可变剪接"}]},
    {"key": "module_path", "label": "数据文件（可选，留空用已挂载模块）", "type": "file", "default": ""},
    {"key": "case", "label": "实验组", "type": "text", "default": ""},
    {"key": "control", "label": "对照组", "type": "text", "default": ""},
    {"key": "group_col", "label": "分组列", "type": "column", "default": None},
    {"key": "padj", "label": "FDR 阈值", "type": "number", "default": 0.05, "min": 0, "max": 1, "step": 0.01},
    {"key": "lfc", "label": "logFC / dPSI 阈值", "type": "number", "default": 1.0, "min": 0, "max": 10, "step": 0.1},
    {"key": "min_reads", "label": "融合最少 reads", "type": "number", "default": 2, "min": 0, "max": 1000, "step": 1},
]


def run_ncrna(ctx, params, outdir):
    kind = str(_params(params).get("kind") or "circ").lower()
    if kind == "circ":
        return _step_circ(ctx, params, outdir)
    if kind == "mirna":
        return _step_mirna(ctx, params, outdir)
    if kind == "lncrna":
        return _step_lncrna(ctx, params, outdir)
    if kind == "fusion":
        return _step_fusion(ctx, params, outdir)
    if kind == "splicing":
        return _step_splicing(ctx, params, outdir)
    raise EngineError("bad_request", f"未知非编码模块类型：{kind}（可选 circ/mirna/lncrna/fusion/splicing）")


def _analysis_for(kind: str, title: str, method: str):
    def run(ctx, params, outdir):
        merged = dict(params or {})
        merged["kind"] = kind
        return run_ncrna(ctx, merged, outdir)
    return {"id": f"{kind}_analysis", "title": title, "group": "非编码",
            "method": method, "needs": [f"modules:{kind}"], "params": [], "run": run}


def _run_circ_ratio(ctx, params, outdir):
    merged = dict(params or {})
    merged["kind"] = "circ"
    return _step_circ(ctx, merged, outdir)


def _run_fusion_filter(ctx, params, outdir):
    payload = _step_fusion(ctx, params, outdir)
    payload["summary"]["title"] = "融合基因筛选（癌基因/激酶）"
    return payload


def _run_splicing_only(ctx, params, outdir):
    merged = dict(params or {})
    merged["kind"] = "splicing"
    return _step_splicing(ctx, merged, outdir)


STEPS = [
    {
        "id": "ncrna",
        "title": "非编码/融合/剪接",
        "requires": ["normalize"],
        "params": STEP_PARAMS,
        "run": run_ncrna,
    }
]

ANALYSES = [
    {"id": "circ_de", "title": "circRNA 差异分析", "group": "非编码",
     "method": "Welch t + BH（复用 welch_de）", "needs": ["modules:circ"], "params": [], "run": _analysis_for("circ", "circRNA 差异分析", "Welch t + BH")["run"]},
    {"id": "circ_ratio", "title": "circRNA 背靠背比值", "group": "非编码",
     "method": "back-spliced / total ratio", "needs": ["modules:circ"], "params": [], "run": _run_circ_ratio},
    {"id": "circ_circos", "title": "circRNA 基因组圈图", "group": "非编码",
     "method": "解析 chr:start-end 画极坐标圈图", "needs": ["modules:circ"], "params": [], "run": _run_circ_ratio},
    {"id": "mirna_de", "title": "miRNA 差异分析", "group": "非编码",
     "method": "miRBase 名规范化 + Welch t + BH", "needs": ["modules:mirna"], "params": [], "run": _analysis_for("mirna", "miRNA 差异分析", "Welch t + BH")["run"]},
    {"id": "mirna_target", "title": "miRNA 靶基因与 ceRNA", "group": "非编码",
     "method": "内置靶基因表 + 负相关（lite）", "needs": ["modules:mirna"], "params": [], "run": _analysis_for("mirna", "miRNA 靶基因", "内置靶表 + 相关")["run"]},
    {"id": "lncrna_class", "title": "lncRNA 分类注释", "group": "非编码",
     "method": "命名规则 + 内置已知名单", "needs": ["modules:lncrna"], "params": [], "run": _analysis_for("lncrna", "lncRNA 分类", "命名规则")["run"]},
    {"id": "lncrna_de", "title": "lncRNA 差异分析", "group": "非编码",
     "method": "Welch t + BH", "needs": ["modules:lncrna"], "params": [], "run": _analysis_for("lncrna", "lncRNA 差异", "Welch t + BH")["run"]},
    {"id": "lncrna_network", "title": "lncRNA-mRNA 共表达网络", "group": "非编码",
     "method": "Pearson 相关 + FDR，networkx 布局", "needs": ["modules:lncrna"], "params": [], "run": _analysis_for("lncrna", "lncRNA 网络", "相关网络")["run"]},
    {"id": "fusion_summary", "title": "融合基因概览", "group": "非编码",
     "method": "STAR-Fusion/Arriba 表解析 + 支持 reads 排序", "needs": ["modules:fusion"], "params": [], "run": _step_fusion},
    {"id": "fusion_filter", "title": "融合基因筛选（癌基因/激酶）", "group": "非编码",
     "method": "内置癌基因/激酶名单过滤", "needs": ["modules:fusion"], "params": [], "run": _run_fusion_filter},
    {"id": "splicing_dpsi", "title": "可变剪接 dPSI", "group": "非编码",
     "method": "rMATS IncLevel 或 PSI 矩阵 + Welch t + BH", "needs": ["modules:splicing"], "params": [], "run": _run_splicing_only},
    {"id": "splicing_volcano", "title": "剪接火山图", "group": "非编码",
     "method": "dPSI vs -log10(FDR)", "needs": ["modules:splicing"], "params": [], "run": _run_splicing_only},
    {"id": "splicing_event", "title": "剪接事件类型统计", "group": "非编码",
     "method": "SE/RI/A5SS/A3SS/MXE 计数", "needs": ["modules:splicing"], "params": [], "run": _run_splicing_only},
]
