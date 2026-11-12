"""
dsh-geobulk — 数据集状态层

职责：
  * Dataset：一个 bulk 数据集的内存表示（表达矩阵 + 基因 + 样本 + 样本信息 + 步骤状态 + 历史 + 模块）
  * State：多数据集注册表（新建/查找/关闭/切换/清空）
  * demo 数据生成（§2.7）：36 样本 × 4000 基因 × 2 组 × 2 批次 × 临床列 × 免疫/增殖信号

矩阵约定：X / counts 均为 numpy float32，形状 (n_samples, n_genes)。
"""

from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from bg_common import EngineError, jsonable

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO_POOL_PATH = os.path.join(HERE, "demo_gene_pool.txt")
DEMO_N_GENES = 4000
DEMO_N_PER_GROUP = 18
DEMO_SEED = 20261112

# demo 一定包含的真实基因（免疫/增殖/干扰素/EMT），保证免疫浸润与富集有信号
CURATED_GENES = """
CD8A CD8B GZMA GZMB GZMK GZMH PRF1 NKG7 GNLY KLRD1 KLRB1 CD3D CD3E CD3G CD2 CD4 IL7R
FOXP3 CTLA4 PDCD1 LAG3 HAVCR2 TIGIT TOX ICOS CD27 CD28
MS4A1 CD79A CD79B CD19 IGHM IGHD IGHG1 IGKC MZB1 JCHAIN
LYZ CD14 FCGR3A FCGR1A S100A8 S100A9 ITGAM ITGAX CSF1R CD68 CD163 MRC1
FCER1A CLEC9A CD1C LILRA4 IRF8 TCF4
TPSAB1 CPA3 MS4A2 KIT
HBB HBA1 HBA2 ALAS2
PECAM1 VWF CDH5 KDR ENG
COL1A1 COL1A2 COL3A1 COL5A2 ACTA2 TAGLN FN1 VIM S100A4
EPCAM KRT8 KRT18 KRT19 KRT7 MUC1
MKI67 TOP2A CCNB1 CCNB2 CDK1 AURKA AURKB PLK1 BUB1 BUB1B TYMS RRM2
PCNA MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 CCNE1 CCNE2 CDK2 E2F1 E2F3
STAT1 STAT2 IFNG IFNGR1 CXCL9 CXCL10 CXCL11 ISG15 IFI44L IFIT1 IFIT3 OAS1 OAS2 MX1 MX2
IL6 CXCL8 TNF IL1B CCL2 CCL5 CXCL12 CXCR4 CCR7 SELL
TGFB1 TGFBR1 EGF EGFR ERBB2 TP53 MYC KRAS NRAS PIK3CA PTEN AKT1 MTOR
CCND1 CDKN1A CDKN2A RB1 BCL2 BAX BAK1 CASP3 CASP8 CASP9
NFKB1 RELA RELB JAK1 JAK2 STAT3 STAT5A SRC FYN
MMP2 MMP9 MMP14 SNAI1 SNAI2 TWIST1 ZEB1 ZEB2 CDH1 CDH2
WNT5A NOTCH1 HES1 VEGFA HIF1A LDHA PKM HK2 G6PD IDH1 IDH2
GAPDH ACTB TUBB B2M HLA-A HLA-B HLA-C HLA-DRA HLA-DRB1
""".split()

# demo 差异基因分组（用于造出有生物学意义的富集结果）
UP_CELLCYCLE = """
MKI67 TOP2A CCNB1 CCNB2 CDK1 AURKA AURKB PLK1 BUB1 BUB1B TYMS RRM2 PCNA
MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 CCNE1 CCNE2 CDK2 E2F1 E2F3
""".split()
UP_INTERFERON = """
STAT1 STAT2 IFNG CXCL9 CXCL10 CXCL11 ISG15 IFI44L IFIT1 IFIT3 OAS1 OAS2 MX1 MX2
""".split()
UP_REPAIR = """
BRCA1 BRCA2 RAD51 ATM ATR CHEK1 CHEK2 PALB2 FANCA FANCD2 MLH1 MSH2 MSH6
""".split()
DOWN_IMMUNE = """
CD8A GZMB GZMA PRF1 NKG7 GNLY CD3D CD3E CD3G CD2 IL7R CD27 CD28 ICOS
MS4A1 CD79A CD79B CD19 IGHM IGHG1 IGKC MZB1 JCHAIN
""".split()
DOWN_ADHESION = """
CDH1 EPCAM KRT8 KRT18 KRT19 COL1A1 COL1A2 FN1 VIM
""".split()

TIME_PATTERNS = ("time", "times", "os_time", "os.time", "survival_time", "dfs_time", "rfs_time", "pfs_time", "month", "months", "days", "follow")
EVENT_PATTERNS = ("event", "status", "dead", "death", "censor", "os_event", "vital")


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def _dedupe_columns(names) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for raw in names:
        name = str(raw).strip() or "col"
        if name in seen:
            seen[name] += 1
            out.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            out.append(name)
    return out


def dedupe_samples(names) -> list[str]:
    return _dedupe_columns(names)


def column_kind(series: pd.Series) -> str:
    """判断样本信息列的种类：categorical / numeric / time / event。"""
    key = str(series.name or "").lower()
    if series.dropna().empty:
        return "categorical"
    as_str = series.astype(str).str.lower().str.strip()
    uniq = set(as_str.dropna().unique())
    if uniq and uniq <= {"0", "1", "0.0", "1.0"}:
        if any(pat in key for pat in EVENT_PATTERNS):
            return "event"
    if any(pat in key for pat in EVENT_PATTERNS) and pd.api.types.is_numeric_dtype(series):
        return "event"
    if any(pat in key for pat in TIME_PATTERNS) and pd.api.types.is_numeric_dtype(series):
        return "time"
    if pd.api.types.is_numeric_dtype(series) and series.nunique(dropna=True) > 12:
        return "numeric"
    return "categorical"


class Dataset:
    """一个 bulk 数据集。"""

    def __init__(
        self,
        ds_id: str,
        name: str,
        X,
        genes,
        samples,
        *,
        counts=None,
        sample_info: pd.DataFrame | None = None,
        source_path: str | None = None,
        source_format: str | None = None,
        species: str = "human",
        value_kind: str = "unknown",
        notes=None,
        warnings=None,
        group_key: str | None = None,
    ):
        arr = np.asarray(X, dtype=np.float32)
        if arr.ndim != 2:
            raise EngineError("bad_request", "表达矩阵必须是二维（样本 × 基因）")
        self.id = ds_id
        self.name = name
        self.X = arr
        self.counts = None if counts is None else np.asarray(counts, dtype=np.float32)
        self.genes = [str(g) for g in genes]
        self.samples = dedupe_samples([str(s) for s in samples])
        if len(self.genes) != arr.shape[1]:
            raise EngineError("bad_request", f"基因数（{len(self.genes)}）与矩阵列数（{arr.shape[1]}）不一致")
        if len(self.samples) != arr.shape[0]:
            raise EngineError("bad_request", f"样本数（{len(self.samples)}）与矩阵行数（{arr.shape[0]}）不一致")
        if sample_info is None:
            sample_info = pd.DataFrame(index=pd.Index(self.samples, name="sample"))
        else:
            info = sample_info.copy()
            info.index = [str(i) for i in info.index][: len(self.samples)]
            if len(info) != len(self.samples):
                # 对齐到样本顺序，缺失补 NA
                info = info.reindex(pd.Index(self.samples, name=info.index.name or "sample"))
            info.index = pd.Index(self.samples, name=info.index.name or "sample")
            self.sample_info = info
        if sample_info is None:
            self.sample_info.index = pd.Index(self.samples, name="sample")
        self.source_path = source_path
        self.source_format = source_format or "unknown"
        self.species = species or "unknown"
        self.value_kind = value_kind or "unknown"
        self.notes = list(notes or [])
        self.warnings = list(warnings or [])
        self.steps: list[str] = []
        self.current: str | None = None
        self.history: list[dict] = []
        self.analyses_run: list[dict] = []
        self.modules: dict[str, dict] = {}
        self.results: dict[str, dict] = {}
        self.sc_ref: str | None = None
        self.group_key: str | None = None
        self.created_at = now_iso()
        if group_key:
            self.set_group(group_key)
        else:
            guess = self._guess_group()
            if guess:
                self.group_key = guess

    # -- 基本属性 ---------------------------------------------------------
    @property
    def n_samples(self) -> int:
        return int(self.X.shape[0])

    @property
    def n_genes(self) -> int:
        return int(self.X.shape[1])

    def _guess_group(self) -> str | None:
        for prefer in ("group", "Group", "GROUP", "condition", "Condition", "tissue", "Tissue",
                       "disease", "Disease", "type", "Type", "category", "phenotype", "stage"):
            if prefer in self.sample_info.columns:
                col = self.sample_info[prefer]
                if 1 < col.nunique(dropna=True) <= 20:
                    return prefer
        for col in self.sample_info.columns:
            if 1 < self.sample_info[col].nunique(dropna=True) <= 10:
                if pd.api.types.is_numeric_dtype(self.sample_info[col]) and self.sample_info[col].nunique() > 5:
                    continue
                return str(col)
        return None

    def set_group(self, key: str) -> None:
        if key not in self.sample_info.columns:
            raise EngineError("bad_request", f"样本信息里没有列「{key}」")
        values = self.sample_info[key].astype(str)
        order = list(dict.fromkeys(values.tolist()))
        self.group_key = key
        self.group_levels = order

    @property
    def groups(self) -> list[str] | None:
        if not self.group_key or self.group_key not in self.sample_info.columns:
            return None
        values = self.sample_info[self.group_key].astype(str)
        return [v for v in dict.fromkeys(values.tolist()) if v not in ("nan", "None", "")]

    def group_labels(self) -> np.ndarray:
        if not self.group_key:
            return np.array(["all"] * self.n_samples)
        if self.group_key not in self.sample_info.columns:
            return np.array(["all"] * self.n_samples)
        return self.sample_info[self.group_key].astype(str).values

    def sample_columns(self) -> list[dict]:
        out = []
        for col in self.sample_info.columns:
            series = self.sample_info[col]
            out.append({
                "key": str(col),
                "kind": column_kind(series),
                "n_unique": int(series.nunique(dropna=True)),
            })
        return out

    def clinical_columns(self) -> dict:
        """识别时间/事件列。"""
        time_col = event_col = None
        for col in self.sample_info.columns:
            kind = column_kind(self.sample_info[col])
            lname = str(col).lower()
            if kind == "time" and time_col is None:
                time_col = str(col)
            if kind == "event" and event_col is None:
                event_col = str(col)
        return {"time": time_col, "event": event_col}

    def has_clinical(self) -> bool:
        cols = self.clinical_columns()
        return bool(cols["time"] and cols["event"])

    def summary(self, history_limit: int = 20) -> dict:
        out = {
            "id": self.id,
            "name": self.name,
            "path": self.source_path,
            "source_format": self.source_format,
            "n_samples": self.n_samples,
            "n_genes": self.n_genes,
            "species": self.species,
            "value_kind": self.value_kind,
            "sample_columns": self.sample_columns(),
            "group_key": self.group_key,
            "groups": self.groups,
            "clinical": self.clinical_columns(),
            "steps": list(self.steps),
            "current": self.current,
            "modules": [{"kind": k, **{kk: vv for kk, vv in v.items() if kk in ("n", "path", "format", "n_features")}} for k, v in self.modules.items()],
            "warnings": list(self.warnings),
            "notes": list(self.notes),
            "history": list(self.history[-history_limit:]),
            "analyses_run": list(self.analyses_run[-history_limit:]),
            "has_counts": self.counts is not None,
            "has_de": "de" in self.results,
            "has_immune": "immune" in self.results,
            "has_survival": "survival" in self.results,
            "sc_ref": self.sc_ref,
            "created_at": self.created_at,
        }
        return jsonable(out)

    # -- 历史 -------------------------------------------------------------
    def add_history(self, step: str, params: dict, ok: bool, ms: int, summary: str = "", engine: dict | None = None) -> None:
        self.history.append({
            "at": now_iso(),
            "step": str(step),
            "params": jsonable(params or {}),
            "ok": bool(ok),
            "ms": int(ms),
            "summary": str(summary or ""),
            "engine": jsonable(engine or {}),
        })
        if len(self.history) > 50:
            self.history = self.history[-50:]

    def add_analysis_run(self, analysis_id: str, ok: bool, ms: int, summary: str = "") -> None:
        self.analyses_run.append({"at": now_iso(), "id": str(analysis_id), "ok": bool(ok), "ms": int(ms), "summary": str(summary or "")})
        if len(self.analyses_run) > 100:
            self.analyses_run = self.analyses_run[-100:]

    def mark_step(self, step: str) -> None:
        # 重跑某步时清空其下游结果（history 保留）
        order = ["qc", "filter", "normalize", "batch", "pca", "de", "enrich", "immune",
                 "survival", "prognosis", "scmap", "ncrna"]
        if step in order:
            idx = order.index(step)
            downstream = set(order[idx + 1:])
            self.steps = [s for s in self.steps if s not in downstream and s != step]
            self.steps.append(step)
            # 下游缓存清理
            for key in list(self.results.keys()):
                if key in downstream:
                    self.results.pop(key, None)
        else:
            if step not in self.steps:
                self.steps.append(step)
        self.current = step

    def reorder_steps(self) -> None:
        order = ["qc", "filter", "normalize", "batch", "pca", "de", "enrich", "immune",
                 "survival", "prognosis", "scmap", "ncrna"]
        self.steps = [s for s in order if s in self.steps] + [s for s in self.steps if s not in order]

    def apply_replace(self, replace: dict) -> None:
        """步骤返回 replace 时更新矩阵/基因/样本（保持 shape 一致）。"""
        if not replace:
            return
        if "X" in replace:
            arr = np.asarray(replace["X"], dtype=np.float32)
            if arr.ndim != 2:
                raise EngineError("internal", "replace.X 必须是二维")
            if "genes" in replace:
                if len(replace["genes"]) != arr.shape[1]:
                    raise EngineError("internal", "replace.genes 与矩阵列数不一致")
                self.genes = [str(g) for g in replace["genes"]]
            if "samples" in replace:
                if len(replace["samples"]) != arr.shape[0]:
                    raise EngineError("internal", "replace.samples 与矩阵行数不一致")
                self.samples = [str(s) for s in replace["samples"]]
                self.sample_info = self.sample_info.reindex(pd.Index(self.samples, name=self.sample_info.index.name or "sample"))
                self.sample_info.index = pd.Index(self.samples, name=self.sample_info.index.name or "sample")
            self.X = arr
        if "counts" in replace:
            counts = replace["counts"]
            self.counts = None if counts is None else np.asarray(counts, dtype=np.float32)
        if "sample_info" in replace and isinstance(replace["sample_info"], pd.DataFrame):
            self.sample_info = replace["sample_info"]
        if "value_kind" in replace:
            self.value_kind = str(replace["value_kind"])
        if "genes" in replace and "X" not in replace:
            self.genes = [str(g) for g in replace["genes"]]
        if "samples" in replace and "X" not in replace:
            self.samples = [str(s) for s in replace["samples"]]

    def apply_sample_columns(self, columns, values) -> None:
        """分析/步骤返回的新样本列写回 sample_info。"""
        if not columns:
            return
        if isinstance(columns, dict):
            columns = [columns]
        if not values or not isinstance(values, dict):
            return
        for item in columns:
            key = str(item.get("key") if isinstance(item, dict) else item)
            if key not in values:
                continue
            vals = values[key]
            if len(vals) != self.n_samples:
                continue
            try:
                self.sample_info[key] = list(vals)
            except Exception:  # noqa: BLE001
                continue


class State:
    """多数据集注册表。"""

    def __init__(self):
        self.datasets: dict[str, Dataset] = {}
        self.active: str | None = None
        self.counter = 0

    def new(self, **kwargs) -> Dataset:
        self.counter += 1
        ds_id = f"ds{self.counter}"
        dataset = Dataset(ds_id, **kwargs)
        self.datasets[ds_id] = dataset
        if not self.active:
            self.active = ds_id
        return dataset

    def add(self, dataset: Dataset) -> Dataset:
        self.datasets[dataset.id] = dataset
        if not self.active:
            self.active = dataset.id
        return dataset

    def get(self, ds_id: str | None = None) -> Dataset:
        target = ds_id or self.active
        if not target:
            raise EngineError("bad_request", "还没有打开任何数据集，请先导入数据或载入示例")
        dataset = self.datasets.get(target)
        if dataset is None:
            raise EngineError("not_found", f"数据集 {target} 不存在（可能已关闭）")
        return dataset

    def find_by_path(self, path: str) -> Dataset | None:
        target = os.path.abspath(path)
        for dataset in self.datasets.values():
            if dataset.source_path and os.path.abspath(dataset.source_path) == target:
                return dataset
        return None

    def close(self, ds_id: str) -> None:
        if ds_id not in self.datasets:
            raise EngineError("not_found", f"数据集 {ds_id} 不存在")
        self.datasets.pop(ds_id, None)
        if self.active == ds_id:
            self.active = next(iter(self.datasets), None)

    def reset(self) -> None:
        self.datasets.clear()
        self.active = None

    def summaries(self) -> list[dict]:
        return [ds.summary() for ds in self.datasets.values()]


# ---------------------------------------------------------------------------
# demo 数据
# ---------------------------------------------------------------------------

def _load_pool() -> list[str]:
    try:
        with open(DEMO_POOL_PATH, "r", encoding="utf-8") as handle:
            genes = [line.strip() for line in handle if line.strip()]
        if genes:
            return genes
    except OSError:
        pass
    return []


def demo_gene_list(n_genes: int = DEMO_N_GENES, seed: int = DEMO_SEED) -> list[str]:
    """确定性的 demo 基因清单：精选基因在前，其余从真实符号池补充。"""
    rng = np.random.default_rng(seed)
    pool = _load_pool()
    genes: list[str] = []
    seen = set()
    for gene in CURATED_GENES:
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    if pool:
        pool = list(dict.fromkeys(pool))
        order = rng.permutation(len(pool))
        for idx in order:
            gene = pool[int(idx)]
            if gene in seen:
                continue
            genes.append(gene)
            seen.add(gene)
            if len(genes) >= n_genes:
                break
    i = 1
    while len(genes) < n_genes:
        gene = f"GENE{i:04d}"
        i += 1
        if gene in seen:
            continue
        genes.append(gene)
        seen.add(gene)
    return genes[:n_genes]


def _index_of(genes: list[str], names) -> list[int]:
    lookup = {g: i for i, g in enumerate(genes)}
    out = []
    for name in names:
        idx = lookup.get(name)
        if idx is not None and idx not in out:
            out.append(idx)
    return out


def make_demo(state: State, name: str = "demo-bulk") -> Dataset:
    """内置示例：36 样本 × 4000 基因，两组两批次，带临床与免疫/增殖信号。"""
    rng = np.random.default_rng(DEMO_SEED)
    genes = demo_gene_list()
    n_genes = len(genes)
    n_case = DEMO_N_PER_GROUP
    n_ctrl = DEMO_N_PER_GROUP
    n_samples = n_case + n_ctrl

    case_idx = np.arange(n_case)
    ctrl_idx = np.arange(n_case, n_samples)
    labels = np.array(["case"] * n_case + ["control"] * n_ctrl)
    batch = np.array(["B1", "B2"] * (n_samples // 2))[:n_samples]

    # 基础 log2 表达（基因均值 + 样本文库因子 + 噪声）
    gene_mean = rng.normal(6.3, 1.35, n_genes)
    lib_factor = np.exp(rng.normal(0.0, 0.22, n_samples))
    log_x = gene_mean[None, :] + rng.normal(0.0, 0.42, (n_samples, n_genes))

    # 分组效应：增殖/干扰素/修复上调，免疫/黏附下调
    up_idx = _index_of(genes, UP_CELLCYCLE) + _index_of(genes, UP_INTERFERON) + _index_of(genes, UP_REPAIR)
    up_idx += [i for i in rng.choice(n_genes, 90, replace=False).tolist()]
    down_idx = _index_of(genes, DOWN_IMMUNE) + _index_of(genes, DOWN_ADHESION)
    down_idx += [i for i in rng.choice(n_genes, 80, replace=False).tolist()]
    up_idx = sorted(set(up_idx))
    down_idx = sorted(set(down_idx) - set(up_idx))
    log_x[np.ix_(case_idx, up_idx)] += rng.uniform(1.3, 2.4, (n_case, len(up_idx)))
    log_x[np.ix_(ctrl_idx, down_idx)] += rng.uniform(1.3, 2.2, (n_ctrl, len(down_idx)))

    # 批次效应：批次 2 在 ~600 个基因上整体偏移
    batch_genes = rng.choice(n_genes, min(600, n_genes), replace=False)
    batch2 = batch == "B2"
    log_x[np.ix_(np.where(batch2)[0], batch_genes)] += rng.normal(0.45, 0.12, (int(batch2.sum()), len(batch_genes)))

    # 免疫潜变量：控制免疫基因共表达（供免疫浸润方法使用）
    immune_latent = rng.normal(0.0, 1.0, n_samples)
    immune_latent[case_idx] -= 0.7  # 病例组免疫浸润更低
    immune_genes = _index_of(genes, """
        CD8A CD8B GZMA GZMB GZMK PRF1 NKG7 GNLY KLRD1 CD3D CD3E CD3G CD2 IL7R
        MS4A1 CD79A CD79B CD19 IGHM IGKC MZB1 LYZ CD14 FCGR3A ITGAM CD68 CD163 MRC1
        FCER1A CD1C LILRA4 TPSAB1 CPA3 HBB PECAM1 VWF COL1A1 COL1A2 ACTA2
        STAT1 ISG15 IFIT1 CXCL9 CXCL10 HLA-A HLA-B HLA-DRA B2M
    """.split())
    loads = rng.uniform(0.35, 0.75, len(immune_genes))
    log_x[:, immune_genes] += immune_latent[:, None] * loads[None, :]

    # 风险签名（细胞周期基因均值）驱动生存
    risk_genes = _index_of(genes, UP_CELLCYCLE)
    risk_raw = log_x[:, risk_genes].mean(axis=1)
    risk_z = (risk_raw - risk_raw.mean()) / (risk_raw.std() + 1e-9)
    hazard = 0.018 * np.exp(0.55 * risk_z)
    surv_time = -np.log(rng.uniform(0.02, 0.98, n_samples)) / hazard
    censor_time = 60.0
    observed_time = np.minimum(surv_time, censor_time)
    event = (surv_time <= censor_time).astype(int)
    # 确保两组都有事件
    event[:4] = np.maximum(event[:4], 1)

    stage = rng.choice(["I", "II", "III", "IV"], n_samples, p=[0.2, 0.3, 0.3, 0.2])
    age = np.clip(rng.normal(62, 9, n_samples), 30, 88).round(0)
    sex = rng.choice(["male", "female"], n_samples)
    patient = [f"P{i + 1:03d}" for i in range(n_samples)]

    sample_names = [f"{'C' if lab == 'case' else 'N'}{i + 1:02d}" for i, lab in enumerate(labels)]
    sample_info = pd.DataFrame(
        {
            "sample": sample_names,
            "group": labels,
            "batch": batch,
            "patient_id": patient,
            "age": age,
            "sex": sex,
            "stage": stage,
            "os_time": observed_time.round(1),
            "os_event": event,
        }
    ).set_index("sample")

    # 转成 counts（负二项近似泊松），并保留 log 参考矩阵
    mu = np.power(2.0, log_x) * lib_factor[:, None]
    counts = rng.poisson(mu).astype(np.float32)
    counts[counts < 0] = 0

    dataset = Dataset(
        ds_id="",
        name=name,
        X=counts,
        genes=genes,
        samples=sample_names,
        counts=counts,
        sample_info=sample_info,
        source_path=None,
        source_format="demo",
        species="human",
        value_kind="counts",
        notes=["内置示例数据：36 样本 × 4000 基因，两组（case/control）× 两批次，带临床列"],
        warnings=[],
        group_key="group",
    )
    state.counter += 1
    dataset.id = f"ds{state.counter}"
    state.add(dataset)
    return dataset


def build_sample_info(samples, info: pd.DataFrame | None, extra: dict | None = None) -> pd.DataFrame:
    """把导入层给的样本信息整理成规整 DataFrame（index=样本名）。"""
    frame = pd.DataFrame(index=pd.Index([str(s) for s in samples], name="sample"))
    if info is not None and len(info):
        src = info.copy()
        src.columns = [str(c) for c in src.columns]
        # 尝试按索引对齐；对齐失败则按位置
        try:
            src.index = [str(i) for i in src.index]
            aligned = src.reindex(frame.index)
            if aligned.isna().all(axis=None) and len(src) == len(frame):
                aligned = src.copy()
                aligned.index = frame.index
            frame = frame.join(aligned)
        except Exception:  # noqa: BLE001
            pass
    if extra:
        for key, values in extra.items():
            if key in frame.columns:
                continue
            if isinstance(values, (list, tuple, np.ndarray)) and len(values) == len(frame):
                frame[str(key)] = list(values)
    return frame
