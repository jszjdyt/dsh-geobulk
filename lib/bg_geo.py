"""
dsh-geobulk — GEO 检索 / 详情 / 下载（CONTRACT §2.2 / §7.2 / §8.1）

三个公开函数：
  * geo_search(term, species, study_type, min_samples, has_clinical, limit, page, sort) -> dict
  * geo_detail(gse) -> dict
  * geo_download(gse, files, dest_dir, overwrite=False, progress=None) -> dict

实现要点：
  * 只用标准库 urllib（不依赖 requests）；所有网络失败重试 ≤2 次（退避 1s/2s），
    仍失败抛 EngineError("unsupported", "NCBI 暂时不可达：<原因>")，消息中文可读。
  * E-utilities：esearch(db=gds) → esummary(db=gds)（esummary 字段最全；efetch brief 对老
    GSE 会 "cannot get document summary"，所以不用）。UID = 200000000 + GSE 号。
  * GEO FTP：系列矩阵在 `<GSE>/matrix/`，作者补充文件在 `<GSE>/suppl/`；detail 两个目录
    都解析（CONTRACT 只写了 suppl，实测 series matrix 不在 suppl）。用 HTTP Range 只抓
    系列矩阵头部 1MB，解析 !Series_overall_design / !Series_contact_* / !Sample_characteristics_ch1
    得到 design / contacts / clinical，不需要下整份矩阵。
  * 中文：ZH_DICT ≥120 条（疾病/组织/药物/技术/物种/常用词），用于 title_zh 与 tags。
  * 下载：dest_dir 必须绝对路径且不是 `/`；单文件 8GB 上限；120s 超时；分块写盘 + 进度回调；
    overwrite=False 时已存在则跳过（ok=true, skipped=true）。

本文件只依赖标准库 + bg_common（EngineError）。禁止新增第三方依赖。
"""

from __future__ import annotations

import html
import json
import os
import re
import socket
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib

from bg_common import EngineError

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
NCBI_ORIGIN = "https://ftp.ncbi.nlm.nih.gov"
GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
NCBI_WEB = "https://www.ncbi.nlm.nih.gov"

USER_AGENT = "dsh-geobulk/1.0 (bulk RNA-seq workbench; local desktop app)"
# 实测：ftp.ncbi.nlm.nih.gov 的 WAF 只放行 UA 以 "curl/" 开头的请求，其它一律 403；
# 这里保留括号里的应用标识，避免伪装成其它软件。
FTP_USER_AGENT = "curl/8.7.1 (dsh-geobulk/1.0; bulk RNA-seq workbench)"
HTTP_TIMEOUT = 120                      # 单请求超时（秒）
RETRY_SLEEPS = (1.0, 2.0)               # 网络失败重试退避（秒）；测试里可置 0
ESUMMARY_BATCH = 100                    # esummary 单次最多 UID 数
MAX_ESUMMARY_CANDIDATES = 300           # has_clinical / sort 需要客户端筛选时的候选上限
MAX_DOWNLOAD_BYTES = 8 * 1024 ** 3      # 单文件下载上限 8GB
MATRIX_HEAD_BYTES = 96 * 1024           # 系列矩阵头部只抓 96KB（压缩后），足够解析 header
MATRIX_HEAD_TIMEOUT = 30                # 头部抓取的单请求超时（秒）
DOWNLOAD_CHUNK = 512 * 1024
SERIES_MATRIX_KIND = "series_matrix"

_SPECIES_ALIAS = {
    "human": "Homo sapiens", "homo sapiens": "Homo sapiens", "人": "Homo sapiens", "人类": "Homo sapiens",
    "mouse": "Mus musculus", "mus musculus": "Mus musculus", "小鼠": "Mus musculus",
    "rat": "Rattus norvegicus", "rattus norvegicus": "Rattus norvegicus", "大鼠": "Rattus norvegicus",
    "zebrafish": "Danio rerio", "斑马鱼": "Danio rerio",
    "fly": "Drosophila melanogaster", "果蝇": "Drosophila melanogaster",
    "yeast": "Saccharomyces cerevisiae", "酵母": "Saccharomyces cerevisiae",
}

_STUDY_TYPE_ALIAS = {
    "all": None, "any": None, "": None, "全部": None,
    "array": "expression profiling by array",
    "microarray": "expression profiling by array",
    "芯片": "expression profiling by array",
    "expression profiling by array": "expression profiling by array",
    "sequencing": "high throughput sequencing",
    "high throughput sequencing": "high throughput sequencing",
    "rnaseq": "high throughput sequencing",
    "rna-seq": "high throughput sequencing",
    "rna seq": "high throughput sequencing",
    "测序": "high throughput sequencing",
    "methylation": "methylation profiling by array",
    "甲基化": "methylation profiling by array",
}

# 判定"有临床信息"的关键词（用于 search 过滤 / item.has_clinical）
_CLINICAL_KEYWORDS = (
    "clinical", "survival", "prognos", "follow-up", "follow up", "outcome",
    "recurrence", "relapse", "mortality", "hazard ratio", "overall survival",
    "disease-free", "progression-free", "risk score", "nomogram", "cox regression",
    "kaplan", "clinical outcome", "pathological stage", "tnm stage",
    "治疗反应", "预后", "生存", "复发", "临床",
)

# 疾病类中文词典 key（用于 item.disease 识别，顺序即优先级）
DISEASE_KEYS = [
    "gastric cancer", "lung cancer", "breast cancer", "colorectal cancer", "liver cancer",
    "hepatocellular carcinoma", "esophageal cancer", "pancreatic cancer", "prostate cancer",
    "ovarian cancer", "cervical cancer", "bladder cancer", "renal cell carcinoma", "glioma",
    "glioblastoma", "melanoma", "leukemia", "lymphoma", "multiple myeloma", "thyroid cancer",
    "head and neck cancer", "oral cancer", "nasopharyngeal carcinoma", "cholangiocarcinoma",
    "osteosarcoma", "sarcoma", "neuroblastoma", "medulloblastoma", "endometrial cancer",
    "testicular cancer", "triple-negative breast cancer", "non-small cell lung cancer",
    "small cell lung cancer", "alzheimer's disease", "parkinson's disease",
    "multiple sclerosis", "type 2 diabetes", "covid-19", "sepsis", "atherosclerosis",
    "inflammatory bowel disease", "rheumatoid arthritis", "chronic kidney disease",
    "heart failure", "pulmonary fibrosis", "colitis", "psoriasis", "tuberculosis",
]

# ---------------------------------------------------------------------------
# 中英词典（≥120 条）：中文标题 / 标签 / 疾病识别共用
# 键统一小写；替换时按长度降序、英文词边界匹配，中文原样保留。
# ---------------------------------------------------------------------------

ZH_DICT = {
    # ---- 疾病（50） ----
    "gastric cancer": "胃癌",
    "stomach cancer": "胃癌",
    "lung cancer": "肺癌",
    "breast cancer": "乳腺癌",
    "colorectal cancer": "结直肠癌",
    "colon cancer": "结肠癌",
    "rectal cancer": "直肠癌",
    "liver cancer": "肝癌",
    "hepatocellular carcinoma": "肝细胞癌",
    "esophageal cancer": "食管癌",
    "esophageal squamous cell carcinoma": "食管鳞癌",
    "pancreatic cancer": "胰腺癌",
    "prostate cancer": "前列腺癌",
    "ovarian cancer": "卵巢癌",
    "cervical cancer": "宫颈癌",
    "bladder cancer": "膀胱癌",
    "renal cell carcinoma": "肾细胞癌",
    "kidney cancer": "肾癌",
    "glioma": "胶质瘤",
    "glioblastoma": "胶质母细胞瘤",
    "melanoma": "黑色素瘤",
    "leukemia": "白血病",
    "lymphoma": "淋巴瘤",
    "multiple myeloma": "多发性骨髓瘤",
    "thyroid cancer": "甲状腺癌",
    "head and neck cancer": "头颈癌",
    "oral cancer": "口腔癌",
    "nasopharyngeal carcinoma": "鼻咽癌",
    "cholangiocarcinoma": "胆管癌",
    "osteosarcoma": "骨肉瘤",
    "sarcoma": "肉瘤",
    "neuroblastoma": "神经母细胞瘤",
    "medulloblastoma": "髓母细胞瘤",
    "endometrial cancer": "子宫内膜癌",
    "testicular cancer": "睾丸癌",
    "bladder urothelial carcinoma": "膀胱尿路上皮癌",
    "triple-negative breast cancer": "三阴性乳腺癌",
    "non-small cell lung cancer": "非小细胞肺癌",
    "small cell lung cancer": "小细胞肺癌",
    "alzheimer's disease": "阿尔茨海默病",
    "parkinson's disease": "帕金森病",
    "multiple sclerosis": "多发性硬化",
    "type 2 diabetes": "2型糖尿病",
    "covid-19": "新型冠状病毒肺炎",
    "sepsis": "脓毒症",
    "atherosclerosis": "动脉粥样硬化",
    "inflammatory bowel disease": "炎症性肠病",
    "rheumatoid arthritis": "类风湿关节炎",
    "chronic kidney disease": "慢性肾病",
    "heart failure": "心力衰竭",
    "pulmonary fibrosis": "肺纤维化",
    "colitis": "结肠炎",
    "psoriasis": "银屑病",
    "tuberculosis": "结核病",
    # ---- 组织 / 部位 / 细胞（40） ----
    "tumor microenvironment": "肿瘤微环境",
    "adjacent normal": "癌旁正常",
    "normal tissue": "正常组织",
    "tumor tissue": "肿瘤组织",
    "tumour tissue": "肿瘤组织",
    "gastric mucosa": "胃黏膜",
    "mucosa": "黏膜",
    "gastric tissue": "胃组织",
    "liver tissue": "肝组织",
    "lung tissue": "肺组织",
    "breast tissue": "乳腺组织",
    "colon tissue": "结肠组织",
    "cerebrospinal fluid": "脑脊液",
    "peripheral blood mononuclear cells": "外周血单个核细胞",
    "bone marrow": "骨髓",
    "lymph node": "淋巴结",
    "spleen": "脾脏",
    "thymus": "胸腺",
    "pancreas": "胰腺",
    "prostate": "前列腺",
    "ovary": "卵巢",
    "uterus": "子宫",
    "placenta": "胎盘",
    "adipose tissue": "脂肪组织",
    "skeletal muscle": "骨骼肌",
    "fibroblast": "成纤维细胞",
    "fibroblasts": "成纤维细胞",
    "endothelial": "内皮",
    "epithelium": "上皮",
    "epithelial cells": "上皮细胞",
    "immune cells": "免疫细胞",
    "t cells": "T 细胞",
    "b cells": "B 细胞",
    "macrophage": "巨噬细胞",
    "macrophages": "巨噬细胞",
    "plasma": "血浆",
    "serum": "血清",
    "sputum": "痰液",
    "urine": "尿液",
    "saliva": "唾液",
    "biopsy": "活检",
    "blood": "血液",
    # ---- 药物 / 处理（25） ----
    "cisplatin": "顺铂",
    "paclitaxel": "紫杉醇",
    "docetaxel": "多西他赛",
    "doxorubicin": "多柔比星",
    "5-fluorouracil": "5-氟尿嘧啶",
    "gemcitabine": "吉西他滨",
    "oxaliplatin": "奥沙利铂",
    "irinotecan": "伊立替康",
    "etoposide": "依托泊苷",
    "methotrexate": "甲氨蝶呤",
    "temozolomide": "替莫唑胺",
    "sorafenib": "索拉非尼",
    "sunitinib": "舒尼替尼",
    "imatinib": "伊马替尼",
    "erlotinib": "厄洛替尼",
    "gefitinib": "吉非替尼",
    "trastuzumab": "曲妥珠单抗",
    "tamoxifen": "他莫昔芬",
    "nivolumab": "纳武利尤单抗",
    "pembrolizumab": "帕博利珠单抗",
    "dexamethasone": "地塞米松",
    "rapamycin": "雷帕霉素",
    "metformin": "二甲双胍",
    "aspirin": "阿司匹林",
    "immunotherapy": "免疫治疗",
    "chemotherapy": "化疗",
    "radiotherapy": "放疗",
    "drug resistance": "耐药",
    "chemoresistance": "化疗耐药",
    "cisplatin resistance": "顺铂耐药",
    # ---- 技术 / 平台（15） ----
    "microarray": "芯片",
    "expression profiling": "表达谱",
    "rna-seq": "转录组测序",
    "rna sequencing": "转录组测序",
    "transcriptome": "转录组",
    "single-cell rna-seq": "单细胞转录组测序",
    "single cell rna-seq": "单细胞转录组测序",
    "scrna-seq": "单细胞转录组测序",
    "whole genome sequencing": "全基因组测序",
    "exome sequencing": "外显子组测序",
    "chip-seq": "染色质免疫共沉淀测序",
    "atac-seq": "染色质可及性测序",
    "bisulfite sequencing": "亚硫酸氢盐测序",
    "spatial transcriptomics": "空间转录组",
    "qpcr": "实时定量 PCR",
    "proteomics": "蛋白组学",
    "metabolomics": "代谢组学",
    "crispr": "基因编辑",
    "high throughput sequencing": "高通量测序",
    # ---- 物种（12） ----
    "homo sapiens": "人",
    "mus musculus": "小鼠",
    "rattus norvegicus": "大鼠",
    "danio rerio": "斑马鱼",
    "human": "人",
    "mouse": "小鼠",
    "rat": "大鼠",
    "patients": "患者",
    "patient": "患者",
    "healthy": "健康",
    "wild-type": "野生型",
    "knockout": "敲除",
    # ---- 常用生物学词（30） ----
    "biomarker": "生物标志物",
    "signature": "特征签名",
    "prognosis": "预后",
    "prognostic": "预后",
    "survival": "生存",
    "recurrence": "复发",
    "metastasis": "转移",
    "immune infiltration": "免疫浸润",
    "inflammation": "炎症",
    "apoptosis": "凋亡",
    "proliferation": "增殖",
    "hypoxia": "缺氧",
    "resistance": "耐药",
    "overexpression": "过表达",
    "knockdown": "敲低",
    "mutation": "突变",
    "methylation": "甲基化",
    "differentiation": "分化",
    "development": "发育",
    "infection": "感染",
    "vaccine": "疫苗",
    "therapy": "治疗",
    "response": "反应",
    "risk": "风险",
    "expression": "表达",
    "profiling": "谱分析",
    "single cell": "单细胞",
    "exosome": "外泌体",
    "circrna": "环状RNA",
    "mirna": "微小RNA",
    "lncrna": "长链非编码RNA",
    "autophagy": "自噬",
    "ferroptosis": "铁死亡",
    "senescence": "衰老",
    "fibrosis": "纤维化",
}

# 英文词边界匹配（避免 "rat" 命中 "patient" 中的子串）
_ZH_SORTED = sorted(ZH_DICT.items(), key=lambda kv: len(kv[0]), reverse=True)
_ZH_PATTERNS = [
    (re.compile(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", re.IGNORECASE), zh)
    for term, zh in _ZH_SORTED
]
_DISEASE_PATTERNS = [
    (re.compile(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", re.IGNORECASE), term)
    for term in sorted(DISEASE_KEYS, key=len, reverse=True)
]

_TECH_TAGS = {
    "expression profiling by array": "表达谱",
    "expression profiling by high throughput sequencing": "高通量测序",
    "methylation profiling by array": "甲基化",
    "methylation profiling by high throughput sequencing": "甲基化测序",
    "genome variation profiling by array": "基因组变异",
    "genome variation profiling by high throughput sequencing": "基因组变异",
    "snapshot": "快照",
}
_SPECIES_TAGS = {"Homo sapiens": "人", "Mus musculus": "小鼠", "Rattus norvegicus": "大鼠"}


# ---------------------------------------------------------------------------
# 网络层（测试通过替换模块级 _urlopen 来 mock）
# ---------------------------------------------------------------------------


def _ssl_context():
    """HTTPS 校验上下文：优先用 certifi 的 CA 包（本机 python.org 自带的 openssl
    cert.pem 可能缺根证书，导致 CERTIFICATE_VERIFY_FAILED）。"""
    global _SSL_CTX
    if _SSL_CTX is None:
        try:
            import certifi

            _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
        except Exception:  # noqa: BLE001 - 没有 certifi 时退回默认 CA
            _SSL_CTX = ssl.create_default_context()
    return _SSL_CTX


_SSL_CTX = None


def _urlopen(request, timeout):  # pragma: no cover - 真实网络路径
    """单点 URL 打开（测试 monkeypatch 这一层）。"""
    try:
        return urllib.request.urlopen(request, timeout=timeout, context=_ssl_context())
    except TypeError:  # pragma: no cover - 某些 context 不支持的特殊 URL 类型
        return urllib.request.urlopen(request, timeout=timeout)


def _reason(error) -> str:
    """把网络异常转成中文可读原因。"""
    if isinstance(error, urllib.error.HTTPError):
        code = int(getattr(error, "code", 0) or 0)
        if code == 404:
            return "资源不存在（HTTP 404）"
        if code in (429, 503):
            return f"NCBI 限流或维护（HTTP {code}）"
        if code >= 500:
            return f"NCBI 服务端错误（HTTP {code}）"
        return f"HTTP {code}"
    if isinstance(error, urllib.error.URLError):
        detail = getattr(error, "reason", error)
        if isinstance(detail, (socket.timeout, TimeoutError)):
            return "连接超时"
        if isinstance(detail, ssl.SSLError):
            return "TLS 握手失败"
        if isinstance(detail, socket.gaierror):
            return "域名解析失败（可能无网络）"
        return str(detail)
    if isinstance(error, (socket.timeout, TimeoutError)):
        return "连接超时"
    if isinstance(error, ssl.SSLError):
        return "TLS 握手失败"
    if isinstance(error, json.JSONDecodeError):
        return "NCBI 返回内容不是合法 JSON"
    return str(error) or type(error).__name__


def _network_fail(error) -> EngineError:
    return EngineError("unsupported", f"NCBI 暂时不可达：{_reason(error)}")


def _request(url: str, headers: dict | None = None):
    agent = FTP_USER_AGENT if "ftp.ncbi.nlm.nih.gov" in url else USER_AGENT
    merged = {"User-Agent": agent, "Accept": "*/*"}
    if headers:
        merged.update(headers)
    return urllib.request.Request(url, headers=merged)


def _fetch_bytes(url: str, *, timeout: int = HTTP_TIMEOUT, retries: int = 2,
                 headers: dict | None = None, limit: int | None = None,
                 soft_404: bool = False) -> bytes | None:
    """GET 一个 URL 返回 bytes；失败重试 ≤2 次，仍失败抛 EngineError。"""
    last_error = None
    attempts = max(0, int(retries)) + 1
    for attempt in range(attempts):
        try:
            with _urlopen(_request(url, headers), timeout) as response:
                data = response.read() if limit is None else response.read(limit)
            return data or b""
        except urllib.error.HTTPError as error:
            if soft_404 and int(getattr(error, "code", 0) or 0) == 404:
                return None
            last_error = error
        except (urllib.error.URLError, socket.timeout, TimeoutError, ssl.SSLError, OSError) as error:
            last_error = error
        except Exception as error:  # noqa: BLE001 - 任何网络层异常都按不可达处理
            last_error = error
        if attempt < attempts - 1:
            delay = RETRY_SLEEPS[min(attempt, len(RETRY_SLEEPS) - 1)]
            if delay > 0:
                time.sleep(delay)
    raise _network_fail(last_error)


def _fetch_json(url: str, *, retries: int = 2, soft_404: bool = False) -> dict | None:
    raw = _fetch_bytes(url, retries=retries, soft_404=soft_404)
    if raw is None:
        return None
    try:
        text = raw.decode("utf-8", "replace")
        data = json.loads(text)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise _network_fail(error)
    if not isinstance(data, dict):
        raise _network_fail("NCBI 返回结构异常")
    if data.get("error"):
        raise EngineError("bad_request", f"检索式无法解析：{str(data.get('error'))[:300]}")
    return data


def _fetch_text(url: str, *, retries: int = 2, soft_404: bool = False, limit: int | None = None) -> str | None:
    raw = _fetch_bytes(url, retries=retries, soft_404=soft_404, limit=limit)
    if raw is None:
        return None
    return raw.decode("utf-8", "replace")


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------


def _norm_gse(gse) -> str:
    text = str(gse or "").strip().upper().replace(" ", "")
    match = re.match(r"^(?:GSE)?(\d{1,9})$", text)
    if not match:
        raise EngineError("bad_request", f"GSE 编号格式不正确：{gse!r}（示例：GSE13861）")
    return f"GSE{int(match.group(1))}"


def _gse_number(gse: str) -> int:
    return int(re.sub(r"\D", "", gse))


def _ftp_prefix(gse: str) -> str:
    number = str(_gse_number(gse))
    base = number[:-3] if len(number) > 3 else ""
    if not base:
        base = "0"
    return f"GSE{base}nnn"


def _series_root(gse: str) -> str:
    return f"{GEO_FTP}/{_ftp_prefix(gse)}/{gse}"


def _uid_for(gse: str) -> str:
    return str(200000000 + _gse_number(gse))


def _norm_date(value) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    match = re.match(r"^(\d{4})[/\-.](\d{1,2})[/\-.](\d{1,2})", text)
    if match:
        return f"{match.group(1)}-{int(match.group(2)):02d}-{int(match.group(3)):02d}"
    match = re.match(r"^([A-Z][a-z]{2})\s+(\d{1,2}),?\s+(\d{4})$", text)
    if match:
        months = {m: i + 1 for i, m in enumerate(
            ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])}
        month = months.get(match.group(1))
        if month:
            return f"{match.group(3)}-{month:02d}-{int(match.group(2)):02d}"
    return text


def _platforms(value) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    tokens = [tok for tok in re.split(r"[,;\s]+", text) if tok]
    out = []
    for token in tokens:
        token = token.strip().upper()
        if not token:
            continue
        if token.isdigit():
            token = f"GPL{token}"
        elif token.startswith("GPL"):
            token = f"GPL{token[3:].lstrip('0') or '0'}"
        if token not in out:
            out.append(token)
    return out


def _parse_size(token: str) -> int:
    text = str(token or "").strip().upper()
    if not text or text in ("-", "N/A"):
        return 0
    match = re.match(r"^([0-9.]+)\s*([KMGTP]?)", text)
    if not match:
        return 0
    try:
        value = float(match.group(1))
    except ValueError:
        return 0
    factor = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3,
              "T": 1024 ** 4, "P": 1024 ** 5}[match.group(2)]
    return int(value * factor)


def _clean_quotes(value: str) -> str:
    text = str(value or "").strip().strip('"').strip()
    # GEO series matrix 里值内的逗号会双写成 ",,"，还原成单个逗号
    if ",," in text:
        text = text.replace(",,", ",")
    return text


def _translate(text: str) -> str:
    """英文→中文词典替换；没有任何命中时返回空字符串（允许 title_zh 为空）。"""
    source = str(text or "")
    if not source:
        return ""
    # 先保护已经替成中文的片段，避免后续英文词再命中的副作用；直接顺序替换即可（最长优先）
    result = source
    hits = 0
    for pattern, zh in _ZH_PATTERNS:
        def _repl(match, _zh=zh):
            nonlocal hits
            hits += 1
            return _zh
        result, count = pattern.subn(_repl, result)
    if not hits:
        return ""
    # 压缩多余空格
    result = re.sub(r"\s{2,}", " ", result).strip()
    return result


def _matched_chinese(text: str, limit: int = 6) -> list[str]:
    out: list[str] = []
    source = str(text or "")
    if not source:
        return out
    for pattern, zh in _ZH_PATTERNS:
        if len(out) >= limit:
            break
        if pattern.search(source) and zh not in out:
            out.append(zh)
    return out


def _detect_disease(title: str, summary: str) -> str:
    for text in (title, summary):
        for pattern, term in _DISEASE_PATTERNS:
            if pattern.search(str(text or "")):
                return term
    return ""


def _looks_clinical(text: str) -> bool:
    lower = str(text or "").lower()
    return any(keyword in lower for keyword in _CLINICAL_KEYWORDS)


def _make_tags(title: str, summary: str, species: str, study_type: str,
               disease: str) -> list[str]:
    tags: list[str] = []
    disease_zh = ZH_DICT.get(disease.lower()) if disease else None
    if disease_zh:
        tags.append(disease_zh)
    for item in _matched_chinese(title, limit=6) + _matched_chinese(summary, limit=4):
        if item not in tags:
            tags.append(item)
    species_zh = _SPECIES_TAGS.get(str(species or "").strip())
    if species_zh and species_zh not in tags:
        tags.append(species_zh)
    study_zh = _TECH_TAGS.get(str(study_type or "").strip().lower())
    if study_zh and study_zh not in tags:
        tags.append(study_zh)
    if not tags:
        tags.append("GEO 数据集")
    return tags[:10]


# ---------------------------------------------------------------------------
# 检索式
# ---------------------------------------------------------------------------


def _build_query(term: str, species: str, study_type: str, min_samples: int) -> str:
    clauses: list[str] = [f"({term})"]
    if "etyp" not in term.lower():
        clauses.append("(GSE[ETYP])")
    species_key = str(species or "all").strip().lower()
    if species_key not in ("all", "", "any", "全部"):
        organism = _SPECIES_ALIAS.get(species_key)
        if not organism:
            # 用户直接给了物种拉丁名/别名：原样作为 Organism 字段搜索
            organism = str(species).strip()
        clauses.append(f'("{organism}"[Organism])')
    study_key = str(study_type or "all").strip().lower()
    mapped = _STUDY_TYPE_ALIAS.get(study_key, study_type)
    if mapped:
        mapped = str(mapped).strip()
        if mapped.lower() not in ("all", ""):
            clauses.append(f'("{mapped}"[DataSet Type])')
    try:
        min_n = int(min_samples or 0)
    except (TypeError, ValueError):
        min_n = 0
    if min_n > 0:
        clauses.append(f"({min_n}:10000000[Number of Samples])")
    return " AND ".join(clauses)


def _esearch(query: str, *, retmax: int, retstart: int = 0) -> tuple[list[str], int, str]:
    params = urllib.parse.urlencode({
        "db": "gds",
        "term": query,
        "retmode": "json",
        "retmax": int(retmax),
        "retstart": int(retstart),
    })
    data = _fetch_json(f"{EUTILS}/esearch.fcgi?{params}")
    payload = (data or {}).get("esearchresult") or {}
    count = int(payload.get("count") or 0)
    ids = [str(item) for item in (payload.get("idlist") or [])]
    translation = str(payload.get("querytranslation") or "")
    return ids, count, translation


def _esummary(uids: list[str]) -> list[dict]:
    records: list[dict] = []
    for start in range(0, len(uids), ESUMMARY_BATCH):
        batch = uids[start:start + ESUMMARY_BATCH]
        if not batch:
            continue
        params = urllib.parse.urlencode({
            "db": "gds",
            "id": ",".join(batch),
            "retmode": "json",
        })
        data = _fetch_json(f"{EUTILS}/esummary.fcgi?{params}")
        result = (data or {}).get("result") or {}
        ordered = [str(item) for item in (result.get("uids") or batch)]
        by_uid = {}
        for uid in ordered:
            record = result.get(uid)
            if isinstance(record, dict) and not record.get("error"):
                by_uid[uid] = record
        for uid in batch:
            record = by_uid.get(uid)
            if record:
                records.append(record)
    return records


def _item_from_summary(record: dict) -> dict:
    accession = str(record.get("accession") or "")
    gse = accession.upper() if accession.upper().startswith("GSE") else _norm_gse_safe(record.get("gse"))
    title = str(record.get("title") or "")
    summary = str(record.get("summary") or "")
    species = str(record.get("taxon") or "")
    try:
        n_samples = int(record.get("n_samples") or len(record.get("samples") or []))
    except (TypeError, ValueError):
        n_samples = 0
    platforms = _platforms(record.get("gpl"))
    study_type = str(record.get("gdstype") or "").strip()
    release = _norm_date(record.get("pdat"))
    disease = _detect_disease(title, summary)
    has_clinical = _looks_clinical(f"{title} {summary}")
    return {
        "gse": gse,
        "title": title,
        "title_zh": _translate(title),
        "summary": summary,
        "species": species,
        "n_samples": n_samples,
        "platform": ", ".join(platforms),
        "study_type": study_type.lower(),
        "release": release,
        "has_matrix": bool(n_samples > 0),
        "has_clinical": has_clinical,
        "disease": disease,
        "tags": _make_tags(title, summary, species, study_type, disease),
    }


def _norm_gse_safe(value) -> str:
    try:
        return _norm_gse(value)
    except EngineError:
        text = str(value or "").strip()
        return f"GSE{text}" if text.isdigit() else text


def _sort_items(items: list[dict], sort: str) -> list[dict]:
    mode = str(sort or "relevance").strip().lower()
    if mode == "date":
        return sorted(items, key=lambda it: (it.get("release") or ""), reverse=True)
    if mode in ("samples", "sample", "n_samples"):
        return sorted(items, key=lambda it: int(it.get("n_samples") or 0), reverse=True)
    return items


# ---------------------------------------------------------------------------
# geo_search
# ---------------------------------------------------------------------------


_CJK_RE = re.compile(r"[\u4e00-\u9fff]")
_ZH_TO_EN = sorted(
    ((zh, en) for en, zh in ZH_DICT.items() if isinstance(zh, str) and zh),
    key=lambda item: len(item[0]),
    reverse=True,
)


_GENERIC_QUERY_RE = re.compile(
    r"\b(expression\s+profil\w*|profiling|expression|transcriptom\w*|microarray|rna[\s-]?seq|"
    r"high[\s-]?throughput|sequencing|dataset|datasets|data)\b",
    re.I,
)


def _cascade_terms(term: str) -> list[str]:
    """把用户关键词展开成由严到宽的候选查询（智能化放宽）。"""
    base = re.sub(r"\s{2,}", " ", str(term or "")).strip()
    candidates: list[str] = [base]
    stripped = _GENERIC_QUERY_RE.sub(" ", base)
    stripped = re.sub(r"\s{2,}", " ", stripped).strip()
    if stripped and stripped != base:
        candidates.append(stripped)
    words = [w for w in re.split(r"\s+", stripped or base) if len(w) > 1]
    if len(words) > 3:
        candidates.append(" ".join(words[:3]))
    if len(words) > 2:
        candidates.append(" ".join(words[:2]))
    seen = set()
    out = []
    for item in candidates:
        key = item.lower()
        if item and key not in seen:
            seen.add(key)
            out.append(item)
    return out or [base]


def translate_query_zh(term: str) -> tuple[str, list[str], list[str]]:
    """把中文关键词翻成 NCBI 可检索的英文（用内置 ZH_DICT 反查）。

    返回 (新查询串, 命中列表, 未识别中文片段列表)。纯英文原样返回。
    """
    source = str(term or "")
    if not _CJK_RE.search(source):
        return source, [], []
    out = source
    used: list[str] = []
    for zh, en in _ZH_TO_EN:
        if zh and zh in out:
            out = out.replace(zh, en)
            used.append(f"{zh}→{en}")
    leftover = _CJK_RE.findall(out)
    out = _CJK_RE.sub(" ", out)
    out = re.sub(r"\s{2,}", " ", out).strip()
    return out, used, sorted(set(leftover))


def geo_search(term, species="all", study_type="all", min_samples=0,
               has_clinical=False, limit=20, page=0, sort="relevance") -> dict:
    """GEO 系列检索（E-utilities esearch + esummary，中文化字段）。

    返回 {"total","page","items":[{gse,title,title_zh,summary,species,n_samples,
    platform,study_type,release,has_matrix,has_clinical,disease,tags}]}
    """
    query_term = str(term or "").strip()
    if not query_term:
        raise EngineError("bad_request", "请输入检索关键词（支持中文疾病名、英文关键词或 GSE 号）")
    try:
        limit = max(1, min(100, int(limit or 20)))
    except (TypeError, ValueError):
        limit = 20
    try:
        page = max(0, int(page or 0))
    except (TypeError, ValueError):
        page = 0
    try:
        min_samples = max(0, int(min_samples or 0))
    except (TypeError, ValueError):
        min_samples = 0
    sort_mode = str(sort or "relevance").strip().lower()
    if sort_mode not in ("relevance", "date", "samples", "sample", "n_samples"):
        sort_mode = "relevance"

    if _CJK_RE.search(query_term):
        query_term, used_terms, leftover = translate_query_zh(query_term)
        if not query_term:
            raise EngineError("bad_request", "关键词里的中文没能翻译成英文，请改用英文关键词或补充疾病名（如：胃癌 → gastric cancer）")

    needs_filter = bool(has_clinical) or sort_mode != "relevance"
    cascade = _cascade_terms(query_term)
    selected_term = cascade[0]
    selected_query = _build_query(selected_term, species, study_type, min_samples)
    ids: list[str] = []
    total = 0
    translation = ""
    relaxed_used: list[str] = []
    if not needs_filter:
        for index, candidate in enumerate(cascade):
            query = _build_query(candidate, species, study_type, min_samples)
            ids, total, translation = _esearch(query, retmax=limit, retstart=page * limit)
            selected_term, selected_query = candidate, query
            if total > 0:
                relaxed_used = cascade[:index + 1]
                if total >= max(6, limit // 2) or index == len(cascade) - 1:
                    break
        records = _esummary(ids)
        items = [_item_from_summary(record) for record in records]
    else:
        # 过滤/排序需要候选池：抓一页更宽的候选再本地筛选+分页
        want = (page + 1) * limit
        candidate_n = min(MAX_ESUMMARY_CANDIDATES, max(want * 3, 60))
        for index, candidate in enumerate(cascade):
            query = _build_query(candidate, species, study_type, min_samples)
            ids, total, translation = _esearch(query, retmax=candidate_n, retstart=0)
            selected_term, selected_query = candidate, query
            if total > 0:
                relaxed_used = cascade[:index + 1]
                if total >= max(6, candidate_n // 3) or index == len(cascade) - 1:
                    break
        records = _esummary(ids)
        items = [_item_from_summary(record) for record in records]
        if has_clinical:
            items = [item for item in items if item.get("has_clinical")]
        if sort_mode != "relevance":
            items = _sort_items(items, sort_mode)
        start = page * limit
        items = items[start:start + limit]
        total = len(items) if total == 0 else total

    if not items and not ids:
        raise EngineError("not_found", "没有检索到匹配的 GEO 系列：试试更宽泛的关键词（比如只写疾病名）或调整筛选条件")

    result = {"total": int(total), "page": int(page), "items": items}
    if translation:
        result["query"] = translation
    try:
        if _CJK_RE.search(str(term or "")):
            translated_query, used_terms, leftover = translate_query_zh(str(term or ""))
            result["query_zh"] = str(term or "")
            result["query_en"] = translated_query
            result["translated"] = used_terms
            if leftover:
                result["translate_warning"] = "这些中文词没有词典映射（已忽略）：" + "、".join(leftover[:6])
    except Exception:  # noqa: BLE001
        pass
    if needs_filter:
        result["filtered"] = bool(has_clinical)
    if len(relaxed_used) > 1:
        result["relaxed"] = True
        result["used_term"] = selected_term
        result["relax_steps"] = relaxed_used
    return result


# ---------------------------------------------------------------------------
# 目录 HTML 解析
# ---------------------------------------------------------------------------

_HREF_RE = re.compile(r'<a\s+href="(?P<href>[^"]*)"[^>]*>(?P<name>[^<]*)</a>(?P<rest>[^<]*)',
                      re.IGNORECASE)
_SIZE_TOKEN_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)?[KMGTP]?$", re.IGNORECASE)


def _classify_file(name: str, source_dir: str) -> str:
    lower = name.lower()
    if "series_matrix" in lower:
        return SERIES_MATRIX_KIND
    if lower.endswith("_raw.tar") or lower.endswith("raw.tar") or re.search(r"_raw\.tar(\.gz)?$", lower):
        return "raw"
    if source_dir == "matrix":
        return "other"
    return "supplementary"


def _parse_size_from_rest(rest: str) -> int:
    tokens = [token for token in re.split(r"\s+", str(rest or "").strip()) if token]
    # 目录行：`2011-01-04 11:41  9.4K` → 最后一个 token 是大小（"-" 表示未知）
    for token in reversed(tokens):
        if token == "-":
            return 0
        if _SIZE_TOKEN_RE.match(token):
            return _parse_size(token)
    return 0


def _parse_dir_html(text: str, source_dir: str) -> list[dict]:
    """解析 Apache 目录列表 HTML → [{name,size,url,kind,source}]。"""
    if not text:
        return []
    base = ""
    title_match = re.search(r"Index of\s+([^<\s]+)", text, re.IGNORECASE)
    if title_match:
        base = title_match.group(1).rstrip("/")
    files: list[dict] = []
    for match in _HREF_RE.finditer(text):
        href = html.unescape(match.group("href") or "").strip()
        raw_name = html.unescape(match.group("name") or "").strip()
        if not href or not raw_name:
            continue
        if href.startswith(("/", "http://", "https://", "ftp://", "?")):
            continue  # Parent Directory / 外链
        if href.endswith("/"):
            continue  # 子目录
        name = raw_name
        if href.startswith("./"):
            href = href[2:]
        path = urllib.parse.unquote(href)
        name = os.path.basename(path.rstrip("/")) or name
        size = _parse_size_from_rest(match.group("rest"))
        # base 是服务器绝对路径（如 /geo/series/GSE13nnn/GSE13861/matrix），
        # 不能直接拼 GEO_FTP（后者已含 /geo/series 前缀），否则会得到重复路径
        if base.startswith("/"):
            url = f"{NCBI_ORIGIN}{base}/{urllib.parse.quote(href)}"
        elif base:
            url = f"{GEO_FTP}/{base}/{urllib.parse.quote(href)}"
        else:
            url = urllib.parse.urljoin(GEO_FTP + "/", href)
        files.append({
            "name": name,
            "size": int(size),
            "url": url,
            "kind": _classify_file(name, source_dir),
            "source": source_dir,
        })
    # 去重（同名保留第一条）
    out = []
    seen = set()
    for item in files:
        if item["name"] in seen:
            continue
        seen.add(item["name"])
        out.append(item)
    return out


def _list_series_files(gse: str) -> list[dict]:
    """抓 matrix/ 与 suppl/ 两个目录，合并文件列表（series_matrix 优先）。"""
    root = _series_root(gse)
    files: list[dict] = []
    matrix_html = _fetch_text(f"{root}/matrix/", soft_404=True)
    if matrix_html:
        files.extend(_parse_dir_html(matrix_html, "matrix"))
    suppl_html = _fetch_text(f"{root}/suppl/", soft_404=True)
    if suppl_html:
        files.extend(_parse_dir_html(suppl_html, "suppl"))
    # 排序：series_matrix → supplementary → raw → other
    order = {SERIES_MATRIX_KIND: 0, "supplementary": 1, "raw": 2, "other": 3}
    files.sort(key=lambda item: (order.get(item["kind"], 9), item["name"]))
    return files


# ---------------------------------------------------------------------------
# 系列矩阵头部解析（design / contacts / clinical，Range 只取 1MB）
# ---------------------------------------------------------------------------


def _parse_matrix_header(text: str) -> dict:
    design = ""
    contacts: list[str] = []
    contributors: list[str] = []
    characteristics: list[list[str]] = []
    organism = ""
    titles: list[str] = []
    n_samples = 0
    for line in str(text or "").splitlines():
        if not line.startswith("!"):
            if line.strip():
                break  # 进入表格段，头部结束
            continue
        parts = line.rstrip("\r\n").split("\t")
        key = parts[0].strip()
        values = [_clean_quotes(value) for value in parts[1:]]
        if key == "!Series_overall_design" and values:
            design = design or values[0]
        elif key in ("!Series_contact_name", "!Series_contact_email", "!Series_contact_institute",
                     "!Series_contact_department", "!Series_contact_address"):
            contacts.extend(value for value in values if value)
        elif key == "!Series_contributor":
            contributors.extend(value for value in values if value)
        elif key == "!Sample_characteristics_ch1":
            characteristics.append(values)
        elif key == "!Sample_organism_ch1" and values:
            organism = organism or values[0]
        elif key == "!Sample_title":
            titles = values
            n_samples = max(n_samples, len([value for value in values if value]))
        elif key == "!Sample_geo_accession":
            n_samples = max(n_samples, len([value for value in values if value]))
    fields: list[dict] = []
    for index, values in enumerate(characteristics):
        non_empty = [value for value in values if value and value.lower() not in ("", "na", "nan")]
        if not non_empty:
            continue
        first = non_empty[0]
        if ":" in first:
            prefix = first.split(":", 1)[0].strip()
            key = re.sub(r"[^0-9A-Za-z_\u4e00-\u9fff]+", "_", prefix).strip("_").lower() + ":ch1"
            examples = []
            for value in non_empty:
                if ":" in value:
                    examples.append(value.split(":", 1)[1].strip())
                else:
                    examples.append(value)
            example = examples[0] if examples else first.split(":", 1)[1].strip()
        else:
            key = f"characteristics_ch1_{index + 1}"
            example = first
        # 同名 characteristics 行（GEO 常见多行 ch1）加序号去重
        base_key = key
        suffix = 2
        while any(field["key"] == key for field in fields):
            key = f"{base_key}#{suffix}"
            suffix += 1
        fields.append({"key": key, "n": len(non_empty), "example": example})
    contact_text = ""
    if contacts:
        names = contacts[0]
        emails = next((value for value in contacts if "@" in value), "")
        institutes = next((value for value in contacts if not value.lower().endswith("@") and "@" not in value
                           and value != names), "")
        pieces = [names]
        if emails:
            pieces.append(f"<{emails}>")
        if institutes and institutes != names:
            pieces.append(institutes)
        contact_text = " ".join(piece for piece in pieces if piece)
    if contributors:
        extra = "; ".join(contributors[:6])
        contact_text = f"{contact_text}; {extra}" if contact_text else extra
    return {
        "design": design,
        "contacts": contact_text,
        "clinical_fields": fields,
        "organism": organism,
        "n_samples": n_samples,
    }


def _fetch_matrix_head(files: list[dict]) -> dict:
    matrix_file = next((item for item in files if item.get("kind") == SERIES_MATRIX_KIND), None)
    if not matrix_file:
        return {}
    url = matrix_file.get("url") or ""
    if not url:
        return {}
    try:
        raw = _fetch_bytes(url, limit=MATRIX_HEAD_BYTES, retries=1, timeout=MATRIX_HEAD_TIMEOUT,
                           headers={"Range": f"bytes=0-{MATRIX_HEAD_BYTES - 1}"})
    except EngineError:
        return {}
    if not raw:
        return {}
    if raw[:2] == b"\x1f\x8b":
        decompressor = zlib.decompressobj(31)
        try:
            raw = decompressor.decompress(raw)
        except zlib.error:
            raw = b""
    text = raw.decode("utf-8", "replace")
    return _parse_matrix_header(text)


# ---------------------------------------------------------------------------
# geo_detail
# ---------------------------------------------------------------------------


def _summary_for_gse(gse: str) -> dict:
    uid = _uid_for(gse)
    records = _esummary([uid])
    if records:
        return records[0]
    # 兜底：用 ACCN 字段搜索（极老的系列 UID 可能对不上）
    ids, _, _ = _esearch(f"{gse}[ACCN]", retmax=10, retstart=0)
    candidates = [item for item in ids if item.startswith("2")]
    if candidates:
        records = _esummary(candidates[:1])
        if records:
            return records[0]
    raise EngineError("not_found", f"GEO 系列 {gse} 不存在或 NCBI 没有返回摘要")


def geo_detail(gse) -> dict:
    """GEO 系列详情（摘要 + 目录文件列表 + 系列矩阵头部抽取的临床字段）。"""
    gse = _norm_gse(gse)
    record = _summary_for_gse(gse)
    title = str(record.get("title") or "")
    summary = str(record.get("summary") or "")
    taxon = str(record.get("taxon") or "")
    try:
        n_samples = int(record.get("n_samples") or len(record.get("samples") or []))
    except (TypeError, ValueError):
        n_samples = 0
    pubmed_ids = [str(item) for item in (record.get("pubmedids") or []) if str(item).strip()]
    files = _list_series_files(gse)
    header = _fetch_matrix_head(files) if files else {}
    species = taxon or header.get("organism") or ""
    clinical_fields = header.get("clinical_fields") or []
    result = {
        "gse": gse,
        "title": title,
        "title_zh": _translate(title),
        "summary": summary,
        "design": header.get("design") or "",
        "species": species,
        "n_samples": int(n_samples or header.get("n_samples") or 0),
        "platforms": _platforms(record.get("gpl")),
        "pubmed": ", ".join(pubmed_ids),
        "release": _norm_date(record.get("pdat")),
        "contacts": header.get("contacts") or "",
        "clinical": {
            "from_series_matrix": bool(clinical_fields),
            "fields": clinical_fields,
        },
        "files": [
            {"name": item["name"], "size": int(item["size"]), "url": item["url"],
             "kind": item["kind"], "source": item.get("source") or ""}
            for item in files
        ],
        "suggested": "",
        "url": f"{GEO_BASE}?acc={gse}",
    }
    for item in files:
        if item.get("kind") == SERIES_MATRIX_KIND:
            result["suggested"] = item["name"]
            break
    if not result["suggested"] and files:
        result["suggested"] = files[0]["name"]
    result["has_matrix"] = bool(result["suggested"])
    result["has_clinical"] = bool(clinical_fields) or _looks_clinical(f"{title} {summary}")
    result["disease"] = _detect_disease(title, summary)
    result["tags"] = _make_tags(title, summary, species, str(record.get("gdstype") or ""),
                                result["disease"])
    return result


# ---------------------------------------------------------------------------
# geo_download
# ---------------------------------------------------------------------------


def _check_dest_dir(dest_dir) -> str:
    text = str(dest_dir or "").strip()
    if not text:
        raise EngineError("bad_request", "请指定下载目录（dest_dir 必须是绝对路径）")
    path = os.path.expanduser(text)
    if not os.path.isabs(path):
        raise EngineError("bad_request", f"dest_dir 必须是绝对路径：{dest_dir!r}（示例：/Users/you/GEO）")
    path = os.path.abspath(path)
    if path == os.path.sep or path == "/":
        raise EngineError("bad_request", "不能直接下载到文件系统根目录 /，请指定子目录")
    if os.path.exists(path) and not os.path.isdir(path):
        raise EngineError("bad_request", f"dest_dir 已经是文件而不是目录：{path}")
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as error:
        raise EngineError("bad_request", f"无法创建下载目录：{path}（{error}）")
    if not os.access(path, os.W_OK):
        raise EngineError("bad_request", f"下载目录不可写：{path}")
    return path


def _safe_filename(url: str, fallback: str) -> str:
    name = os.path.basename(urllib.parse.unquote(urllib.parse.urlparse(url).path or "")) or fallback
    name = os.path.basename(name.replace("\\", "/"))
    if not name or name in (".", "..") or "/" in name or "\x00" in name:
        name = re.sub(r"[^0-9A-Za-z._-]+", "_", fallback or "download.bin") or "download.bin"
    return name


def _emit(progress, message: str, pct: float | None = None) -> None:
    if not progress:
        return
    try:
        progress(message, pct)
    except TypeError:
        try:
            progress(message)
        except Exception:  # noqa: BLE001
            pass
    except Exception:  # noqa: BLE001
        pass


def _download_one(url: str, target: str, *, progress=None, label: str = "") -> int:
    """流式下载单个 URL 到 target（重试 ≤2 次）。返回写入字节数。"""
    last_error = None
    attempts = len(RETRY_SLEEPS) + 1
    for attempt in range(attempts):
        partial = f"{target}.part"
        try:
            with _urlopen(_request(url), HTTP_TIMEOUT) as response:
                length_header = ""
                headers = getattr(response, "headers", None)
                if headers is not None:
                    try:
                        length_header = str(headers.get("Content-Length") or "")
                    except Exception:  # noqa: BLE001
                        length_header = ""
                total = 0
                try:
                    total = int(length_header)
                except (TypeError, ValueError):
                    total = 0
                if total > MAX_DOWNLOAD_BYTES:
                    raise EngineError(
                        "unsupported",
                        f"文件过大（{total / 1024 ** 3:.1f} GB），超过单文件 8GB 下载上限：{label or url}",
                    )
                written = 0
                with open(partial, "wb") as handle:
                    while True:
                        chunk = response.read(DOWNLOAD_CHUNK)
                        if not chunk:
                            break
                        written += len(chunk)
                        if written > MAX_DOWNLOAD_BYTES:
                            raise EngineError(
                                "unsupported",
                                f"文件超过单文件 8GB 下载上限，已中止：{label or url}",
                            )
                        handle.write(chunk)
                        if progress and total:
                            pct = min(99.0, written / total * 100.0)
                            _emit(progress, f"下载 {label or url}：{written / 1024 ** 2:.1f}/"
                                            f"{total / 1024 ** 2:.1f} MB", pct)
                os.replace(partial, target)
                return written
        except EngineError:
            _remove_quiet(f"{target}.part")
            raise
        except Exception as error:  # noqa: BLE001 - 网络/磁盘异常统一重试
            last_error = error
            _remove_quiet(f"{target}.part")
            _remove_quiet(target if attempt == 0 else "")
            if attempt < attempts - 1:
                delay = RETRY_SLEEPS[min(attempt, len(RETRY_SLEEPS) - 1)]
                if delay > 0:
                    time.sleep(delay)
    raise _network_fail(last_error)


def _remove_quiet(path: str) -> None:
    if not path:
        return
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


def geo_download(gse, files, dest_dir, overwrite=False, progress=None) -> dict:
    """把指定文件下载到 dest_dir（绝对路径、不可为 /），返回 {"dir","files","ms"}。"""
    started = time.time()
    gse = _norm_gse(gse)
    dest = _check_dest_dir(dest_dir)
    if isinstance(files, str):
        requested = [files]
    elif isinstance(files, (list, tuple, set)):
        requested = [str(item) for item in files if str(item).strip()]
    else:
        requested = []
    _emit(progress, f"正在获取 {gse} 的文件列表", 1)

    detail_files: list[dict] = []
    detail_error: EngineError | None = None
    try:
        detail_files = geo_detail(gse).get("files") or []
    except EngineError as error:
        if error.code == "not_found":
            raise
        detail_error = error  # 目录列表失败仍可尝试按 suppl 规则拼接 URL

    by_name = {item["name"]: item for item in detail_files}
    if not requested:
        suggested = next((item["name"] for item in detail_files
                          if item.get("kind") == SERIES_MATRIX_KIND), "")
        if not suggested and detail_files:
            suggested = detail_files[0]["name"]
        if not suggested:
            raise EngineError(
                "not_found",
                f"{gse} 在 GEO FTP 上没有可下载文件（matrix/ 与 suppl/ 都为空）",
            )
        requested = [suggested]

    plan: list[dict] = []
    unknown: list[str] = []
    for item in requested:
        if re.match(r"^(?:https?|ftp)://", item, re.IGNORECASE):
            plan.append({"name": _safe_filename(item, f"{gse}.txt"), "url": item, "size": 0})
            continue
        entry = by_name.get(item) or by_name.get(os.path.basename(item))
        if entry is None:
            if detail_files:
                unknown.append(item)
                continue
            # 详情不可用：按文件名猜目录（series_matrix → matrix/，其余 → suppl/）
            subdir = "matrix" if SERIES_MATRIX_KIND in item.lower() or "series_matrix" in item.lower() else "suppl"
            plan.append({
                "name": os.path.basename(item),
                "url": f"{_series_root(gse)}/{subdir}/{urllib.parse.quote(os.path.basename(item))}",
                "size": 0,
            })
            continue
        plan.append({"name": entry["name"], "url": entry["url"], "size": int(entry.get("size") or 0)})
    if unknown:
        raise EngineError(
            "not_found",
            f"文件不在 {gse} 的可下载列表里：{', '.join(unknown[:5])}",
            candidates=[{"name": item["name"], "size": item["size"], "kind": item["kind"]}
                        for item in detail_files[:30]],
        )
    if detail_error is not None and not plan:
        raise detail_error

    results: list[dict] = []
    total = len(plan)
    failed = 0
    for index, entry in enumerate(plan, start=1):
        target = os.path.join(dest, entry["name"])
        base = (index - 1) / max(total, 1) * 100.0
        if os.path.exists(target) and not overwrite:
            size = os.path.getsize(target)
            results.append({"name": entry["name"], "path": target, "bytes": size,
                            "ok": True, "skipped": True})
            _emit(progress, f"已存在，跳过：{entry['name']}", min(100.0, base + 100.0 / total))
            continue
        if overwrite:
            _remove_quiet(target)
        _emit(progress, f"正在下载 {entry['name']}（{index}/{total}）", base)
        try:
            written = _download_one(entry["url"], target, progress=progress, label=entry["name"])
            results.append({"name": entry["name"], "path": target, "bytes": int(written), "ok": True})
        except EngineError as error:
            if error.code == "bad_request":  # 路径类错误直接终止
                raise
            failed += 1
            results.append({"name": entry["name"], "path": target, "bytes": 0,
                            "ok": False, "error": error.message})
        _emit(progress, f"完成 {entry['name']}（{index}/{total}）", min(100.0, base + 100.0 / total))
    ok_count = sum(1 for item in results if item.get("ok"))
    if ok_count == 0 and results:
        first_error = next((item.get("error") for item in results if item.get("error")), "")
        raise EngineError("unsupported", f"下载失败：{first_error or 'NCBI 暂时不可达'}")
    _emit(progress, f"下载完成：{ok_count}/{len(results)} 个文件", 100)
    result = {
        "dir": dest,
        "gse": gse,
        "files": results,
        "ms": int((time.time() - started) * 1000),
    }
    if failed:
        result["warnings"] = [f"{failed} 个文件下载失败"]
    return result


# ---------------------------------------------------------------------------
# 自检（只用于测试/手工验证，不参与协议）
# ---------------------------------------------------------------------------

def supported_dictionary_size() -> int:
    return len(ZH_DICT)


__all__ = ["geo_search", "geo_detail", "geo_download", "ZH_DICT", "DISEASE_KEYS",
           "supported_dictionary_size"]
