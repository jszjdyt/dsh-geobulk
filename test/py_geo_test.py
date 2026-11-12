#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — GEO 检索/详情/下载的离线测试（CONTRACT §11）

不联网：把 bg_geo 的网络单点 `_urlopen` 换成 FakeUrlOpen，返回固定的
E-utilities JSON / FTP 目录 HTML / series matrix 头部文本 / 文件字节流，
断言 geo_search / geo_detail / geo_download 的解析、过滤、排序、路径校验、
重试与错误处理全部正确。

用法：cd plugins/dsh-geobulk && python3 test/py_geo_test.py
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lib"))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

import bg_geo  # noqa: E402

CHECKS = 0
FAILURES: list[str] = []


def check(condition, label: str, detail: str = "") -> bool:
    global CHECKS
    CHECKS += 1
    if condition:
        return True
    FAILURES.append(f"{label}{('：' + detail) if detail else ''}")
    print(f"  FAIL {label}" + (f"：{detail}" if detail else ""))
    return False


def check_equal(got, want, label: str) -> bool:
    return check(got == want, label, f"got={got!r} want={want!r}")


def check_in(needle, haystack, label: str) -> bool:
    return check(needle in haystack, label, f"{needle!r} 不在 {str(haystack)[:160]!r}")


def check_raises(fn, code: str, label: str) -> Exception | None:
    global CHECKS
    CHECKS += 1
    try:
        fn()
    except bg_geo.EngineError as error:
        if code and error.code != code:
            FAILURES.append(f"{label}：错误码 {error.code} != {code}（{error.message}）")
            print(f"  FAIL {label}：错误码 {error.code} != {code}")
            return error
        return error
    except Exception as error:  # noqa: BLE001
        FAILURES.append(f"{label}：抛了 {type(error).__name__} 而不是 EngineError：{error}")
        print(f"  FAIL {label}：抛了 {type(error).__name__} 而不是 EngineError")
        return None
    FAILURES.append(f"{label}：没有抛异常")
    print(f"  FAIL {label}：没有抛异常")
    return None


# ---------------------------------------------------------------------------
# 假网络：固定 JSON / HTML / 文件流
# ---------------------------------------------------------------------------

ESEARCH_JSON = {
    "header": {"type": "esearch"},
    "esearchresult": {
        "count": "297",
        "retmax": "2",
        "retstart": "0",
        "idlist": ["200013861", "200099999"],
        "querytranslation": "gastric cancer AND GSE[ETYP] AND \"Homo sapiens\"[Organism]",
    },
}

ESUMMARY_JSON = {
    "header": {"type": "esummary"},
    "result": {
        "uids": ["200013861", "200099999"],
        "200013861": {
            "uid": "200013861", "accession": "GSE13861", "gse": "13861",
            "title": "Gene expression signature-based novel prognostic risk score in gastric cancer",
            "summary": "Clinical survival and prognosis of gastric cancer patients treated with "
                       "cisplatin were analyzed to build a prognostic risk score.",
            "gpl": "6884", "taxon": "Homo sapiens", "entrytype": "GSE",
            "gdstype": "Expression profiling by array", "pdat": "2011/07/19",
            "suppfile": "TXT, XLS", "n_samples": 90,
            "pubmedids": ["21447720", "29725014"],
            "samples": [{"accession": "GSM348407", "title": "GAS_DIS_GC_YG62T_C2"}],
            "ftplink": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13861/",
        },
        "200099999": {
            "uid": "200099999", "accession": "GSE99999", "gse": "99999",
            "title": "Lung cancer cell line cisplatin response",
            "summary": "Cell lines were profiled after drug treatment.",
            "gpl": "24676, 6884", "taxon": "Homo sapiens", "entrytype": "GSE",
            "gdstype": "Expression profiling by high throughput sequencing",
            "pdat": "2020/02/03", "suppfile": "", "n_samples": 12,
            "pubmedids": [], "samples": [], "ftplink": "",
        },
    },
}

MATRIX_DIR_HTML = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html><head><title>Index of /geo/series/GSE13nnn/GSE13861/matrix</title></head><body>
<h1>Index of /geo/series/GSE13nnn/GSE13861/matrix</h1><pre>Name                           Last modified      Size  <hr>
<a href="/geo/series/GSE13nnn/GSE13861/">Parent Directory</a>                             -
<a href="GSE13861_series_matrix.txt.gz">GSE13861_series_matrix.txt.gz</a>  2026-10-06 07:05   14M
<hr></pre></body></html>"""

SUPPL_DIR_HTML = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html><head><title>Index of /geo/series/GSE13nnn/GSE13861/suppl</title></head><body>
<h1>Index of /geo/series/GSE13nnn/GSE13861/suppl</h1><pre>Name                                                             Last modified      Size  <hr>
<a href="/geo/series/GSE13nnn/GSE13861/">Parent Directory</a>                                                                      -
<a href="GSE13861_ClinicalInformation.xls.gz">GSE13861_ClinicalInformation.xls.gz</a> 2011-01-04 11:41  9.4K
<a href="GSE13861_RawData.txt.gz">GSE13861_RawData.txt.gz</a>                                         2008-12-10 13:51   20M
<a href="filelist.txt">filelist.txt</a>                                                     2013-01-17 14:47  172
<a href="GSE13861_RAW.tar">GSE13861_RAW.tar</a>                                                 2013-01-17 14:47  6.3M
<hr></pre></body></html>"""

MATRIX_HEAD_TEXT = """!Series_title\t"Gene expression signature-based novel prognostic risk score in gastric cancer"
!Series_geo_accession\t"GSE13861"
!Series_status\t"Public on Jul 19 2011"
!Series_pubmed_id\t"21447720"
!Series_overall_design\t"65 primary gastric adenocarcinoma, 6 GIST and 19 normal tissues were profiled."
!Series_contact_name\t"Sangbae,,Kim"
!Series_contact_email\t"sangbkim22@gmail.com"
!Series_contact_institute\t"HGSC"
!Sample_title\t"GAS_DIS_GC_YG62T_C2"\t"GAS_DIS_GC_YG23T_C1"
!Sample_geo_accession\t"GSM348407"\t"GSM348370"
!Sample_organism_ch1\t"Homo sapiens"\t"Homo sapiens"
!Sample_characteristics_ch1\t"age: 62"\t"age: 51"
!Sample_characteristics_ch1\t"sex: M"\t"sex: F"
!series_matrix_table_begin
"ID_REF"\t"GSM348407"\t"GSM348370"
"""

FILE_BYTES = b"hello geo download\n"


class FakeResponse:
    def __init__(self, data: bytes, status: int = 200, headers: dict | None = None):
        self._data = data or b""
        self._pos = 0
        self.status = status
        self.headers = headers or {"Content-Length": str(len(self._data))}

    def read(self, size: int | None = None) -> bytes:
        if size is None or size < 0:
            chunk = self._data[self._pos:]
            self._pos = len(self._data)
            return chunk
        chunk = self._data[self._pos:self._pos + size]
        self._pos += len(chunk)
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def close(self):
        return None


class FakeUrlOpen:
    """按 URL 子串路由的假 urlopen，同时统计调用次数与请求 URL。"""

    def __init__(self):
        self.calls: list[str] = []
        self.fail_times: dict[str, int] = {}     # URL 子串 → 还要失败几次
        self.fail_code = 503
        self.override: dict[str, bytes] = {}

    def fail(self, needle: str, times: int = 1) -> None:
        self.fail_times[needle] = times

    def _payload(self, url: str):
        if "esearch.fcgi" in url:
            return json.dumps(ESEARCH_JSON).encode()
        if "esummary.fcgi" in url:
            return json.dumps(ESUMMARY_JSON).encode()
        if url.rstrip("/").endswith("/matrix"):
            return MATRIX_DIR_HTML.encode()
        if url.rstrip("/").endswith("/suppl"):
            return SUPPL_DIR_HTML.encode()
        if "series_matrix.txt.gz" in url:
            return MATRIX_HEAD_TEXT.encode()
        if "filelist.txt" in url:
            return FILE_BYTES
        if "GSE13861_ClinicalInformation.xls.gz" in url:
            return b"clinical-xls-bytes"
        return None

    def __call__(self, request, timeout=None):
        url = request.full_url if hasattr(request, "full_url") else str(request)
        self.calls.append(url)
        for needle, remaining in list(self.fail_times.items()):
            if needle in url and remaining > 0:
                self.fail_times[needle] = remaining - 1
                raise urllib.error.HTTPError(url, self.fail_code, "fake failure", None, None)
        if url in self.override:
            return FakeResponse(self.override[url])
        payload = self._payload(url)
        if payload is None:
            raise urllib.error.HTTPError(url, 404, "fake not found", None, None)
        return FakeResponse(payload)


# ---------------------------------------------------------------------------
# 用例
# ---------------------------------------------------------------------------

def test_dictionary() -> None:
    print("[1] ZH_DICT / 查询式构造")
    check(len(bg_geo.ZH_DICT) >= 120, "ZH_DICT ≥120 条", str(len(bg_geo.ZH_DICT)))
    check_equal(bg_geo.ZH_DICT.get("gastric cancer"), "胃癌", "gastric cancer→胃癌")
    check_equal(bg_geo.ZH_DICT.get("cisplatin"), "顺铂", "cisplatin→顺铂")
    check_equal(bg_geo.ZH_DICT.get("homo sapiens"), "人", "homo sapiens→人")
    query = bg_geo._build_query("gastric cancer", "human", "array", 50)
    check_in("GSE[ETYP]", query, "默认加 GSE[ETYP]")
    check_in('"Homo sapiens"[Organism]', query, "human→Organism")
    check_in('"expression profiling by array"[DataSet Type]', query, "array→DataSet Type")
    check_in("50:10000000[Number of Samples]", query, "min_samples→Number of Samples")
    query_all = bg_geo._build_query("gastric cancer", "all", "all", 0)
    check('Organism' not in query_all, "species=all 不加 Organism", query_all)
    check_equal(bg_geo._norm_gse("gse 13861"), "GSE13861", "GSE 编号规范化")
    check_equal(bg_geo._ftp_prefix("GSE13861"), "GSE13nnn", "FTP 前缀规则")
    check_equal(bg_geo._uid_for("GSE13861"), "200013861", "UID 规则")


def test_search(fake: FakeUrlOpen) -> None:
    print("[2] geo_search")
    result = bg_geo.geo_search("gastric cancer", species="human", study_type="all",
                               min_samples=50, limit=20, page=0)
    check_equal(result["total"], 297, "total")
    check_equal(result["page"], 0, "page")
    check_equal(len(result["items"]), 2, "items 条数")
    item = result["items"][0]
    check_equal(item["gse"], "GSE13861", "gse")
    check_equal(item["title"], "Gene expression signature-based novel prognostic risk score in gastric cancer",
                "原文 title 不被修改")
    check_in("胃癌", item["title_zh"], "title_zh 中文化")
    check_in("胃癌", item["tags"], "tags 含胃癌")
    check_equal(item["species"], "Homo sapiens", "species")
    check_equal(item["n_samples"], 90, "n_samples")
    check_equal(item["platform"], "GPL6884", "platform 补 GPL 前缀")
    check_equal(item["study_type"], "expression profiling by array", "study_type 小写")
    check_equal(item["release"], "2011-07-19", "release 日期规范化")
    check(item["has_matrix"] is True, "has_matrix")
    check(item["has_clinical"] is True, "has_clinical（survival/prognosis 关键词）")
    check_equal(item["disease"], "gastric cancer", "disease")
    second = result["items"][1]
    check_equal(second["platform"], "GPL24676, GPL6884", "多平台")
    check(second["has_clinical"] is False, "无临床关键词 has_clinical=false")
    # 检索参数确实拼进了 URL
    search_url = next(url for url in fake.calls if "esearch.fcgi" in url)
    check_in("gastric+cancer", search_url, "检索词进入 URL")
    check_in("Homo+sapiens", search_url, "物种过滤进入 URL")
    check_in("Number+of+Samples", search_url, "样本数过滤进入 URL")
    # has_clinical 过滤 + 排序
    filtered = bg_geo.geo_search("gastric cancer", has_clinical=True, limit=20)
    check_equal(len(filtered["items"]), 1, "has_clinical 过滤后条数")
    check_equal(filtered["items"][0]["gse"], "GSE13861", "has_clinical 过滤结果")
    sorted_by_samples = bg_geo.geo_search("gastric cancer", sort="samples", limit=20)
    check_equal([it["gse"] for it in sorted_by_samples["items"]], ["GSE13861", "GSE99999"],
                "sort=samples 按样本数降序")
    sorted_by_date = bg_geo.geo_search("gastric cancer", sort="date", limit=20)
    check_equal([it["gse"] for it in sorted_by_date["items"]], ["GSE99999", "GSE13861"],
                "sort=date 按发布日期降序")
    # 参数校验
    check_raises(lambda: bg_geo.geo_search("  "), "bad_request", "空关键词→bad_request")


def test_detail(fake: FakeUrlOpen) -> None:
    print("[3] geo_detail")
    detail = bg_geo.geo_detail("gse13861")
    check_equal(detail["gse"], "GSE13861", "detail gse")
    check_equal(detail["species"], "Homo sapiens", "detail species")
    check_equal(detail["n_samples"], 90, "detail n_samples")
    check_equal(detail["platforms"], ["GPL6884"], "detail platforms")
    check_in("21447720", detail["pubmed"], "pubmed id")
    check_equal(detail["release"], "2011-07-19", "detail release")
    check_in("gastric adenocarcinoma", detail["design"], "design 来自 series matrix 头部")
    check_in("Sangbae,Kim", detail["contacts"], "contacts 解析（,, 还原成 ,）")
    check_in("sangbkim22@gmail.com", detail["contacts"], "contacts 邮箱")
    check(detail["clinical"]["from_series_matrix"] is True, "clinical from_series_matrix")
    field_keys = [field["key"] for field in detail["clinical"]["fields"]]
    check_equal(field_keys, ["age:ch1", "sex:ch1"], "clinical 字段 key")
    check_equal(detail["clinical"]["fields"][0]["n"], 2, "clinical 字段样本数")
    check_equal(detail["clinical"]["fields"][0]["example"], "62", "clinical 示例值")
    names = [item["name"] for item in detail["files"]]
    check_in("GSE13861_series_matrix.txt.gz", names, "series matrix 文件")
    check_in("GSE13861_RAW.tar", names, "RAW tar 文件")
    matrix = next(item for item in detail["files"] if item["kind"] == "series_matrix")
    check_equal(matrix["source"], "matrix", "series matrix 来自 matrix/ 目录")
    check_in("/matrix/GSE13861_series_matrix.txt.gz", matrix["url"], "series matrix 绝对 URL")
    check_equal(matrix["size"], 14 * 1024 ** 2, "14M 大小解析")
    check_equal(next(item for item in detail["files"] if item["name"].endswith("_RAW.tar"))["kind"],
                "raw", "RAW.tar → kind=raw")
    check_equal(next(item for item in detail["files"] if item["name"] == "filelist.txt")["kind"],
                "supplementary", "filelist.txt → kind=supplementary")
    check_equal(detail["suggested"], "GSE13861_series_matrix.txt.gz", "suggested 优先 series_matrix")
    check_equal(len([item for item in detail["files"] if item["kind"] == "series_matrix"]), 1,
                "detail 不重复列出 series matrix")
    check_raises(lambda: bg_geo.geo_detail("SRP000001"), "bad_request", "非 GSE 编号→bad_request")
    # 系列不存在：esummary 返回 error → esearch ACCN 兜底也空 → not_found
    original = fake._payload

    def payload_without_13861(url: str):
        if "esummary.fcgi" in url:
            data = json.loads(json.dumps(ESEARCH_JSON))
            data = {"header": {}, "result": {"uids": ["200000000"],
                                             "200000000": {"uid": "200000000", "error": "cannot get document summary"}}}
            return json.dumps(data).encode()
        if "esearch.fcgi" in url:
            data = json.loads(json.dumps(ESEARCH_JSON))
            data["esearchresult"]["idlist"] = []
            data["esearchresult"]["count"] = "0"
            return json.dumps(data).encode()
        return original(url)

    fake._payload = payload_without_13861  # type: ignore[assignment]
    try:
        check_raises(lambda: bg_geo.geo_detail("GSE424242"), "not_found", "系列不存在→not_found")
    finally:
        fake._payload = original  # type: ignore[assignment]


def test_download(fake: FakeUrlOpen) -> None:
    print("[4] geo_download")
    dest = tempfile.mkdtemp(prefix="bg-geo-test-")
    progress: list[tuple[str, float | None]] = []
    result = bg_geo.geo_download("GSE13861", ["filelist.txt"], dest,
                                 progress=lambda message, pct=None: progress.append((message, pct)))
    check_equal(os.path.abspath(dest), result["dir"], "dir 为绝对路径")
    check_equal(len(result["files"]), 1, "下载文件数")
    entry = result["files"][0]
    check(entry["ok"] is True, "ok=true")
    check_equal(entry["bytes"], len(FILE_BYTES), "bytes")
    check_equal(open(entry["path"], "rb").read(), FILE_BYTES, "落盘内容")
    check(isinstance(result["ms"], int) and result["ms"] >= 0, "ms 字段")
    check(len(progress) >= 2, "进度回调被调用", str(len(progress)))
    check(any(pct == 100 for _, pct in progress), "进度到 100")
    # 已存在 → 跳过
    again = bg_geo.geo_download("GSE13861", ["filelist.txt"], dest)
    check(again["files"][0].get("skipped") is True, "overwrite=False 已存在→skipped")
    # overwrite=True → 重下
    overwritten = bg_geo.geo_download("GSE13861", ["filelist.txt"], dest, overwrite=True)
    check(not overwritten["files"][0].get("skipped"), "overwrite=True 不跳过")
    # 指定 URL
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13861/suppl/filelist.txt"
    from_url = bg_geo.geo_download("GSE13861", [url], dest)
    check_equal(from_url["files"][0]["name"], "filelist.txt", "URL 下载文件名取 basename")
    check_equal(open(from_url["files"][0]["path"], "rb").read(), FILE_BYTES, "URL 下载内容")
    # 目录校验
    check_raises(lambda: bg_geo.geo_download("GSE13861", ["filelist.txt"], "relative/dir"),
                 "bad_request", "相对路径 dest_dir→bad_request")
    check_raises(lambda: bg_geo.geo_download("GSE13861", ["filelist.txt"], "/"),
                 "bad_request", "根目录→bad_request")
    check_raises(lambda: bg_geo.geo_download("GSE13861", ["missing.txt"], dest, overwrite=True),
                 "not_found", "不存在的文件→not_found")
    # 文件不存在（mock 404）→ 重试耗尽后 unsupported
    fake.fail("GSE13861_ClinicalInformation.xls.gz", times=99)
    try:
        error = check_raises(
            lambda: bg_geo.geo_download("GSE13861", ["GSE13861_ClinicalInformation.xls.gz"], dest),
            "unsupported", "下载持续失败→unsupported")
        if error is not None:
            check_in("NCBI 暂时不可达", error.message, "网络错误消息中文可读")
    finally:
        fake.fail_times.pop("GSE13861_ClinicalInformation.xls.gz", None)
    # 空 files → 下载 suggested（series matrix）
    empty = bg_geo.geo_download("GSE13861", [], tempfile.mkdtemp(prefix="bg-geo-test2-"))
    check_equal(empty["files"][0]["name"], "GSE13861_series_matrix.txt.gz", "files 为空→下载 suggested")
    shutil.rmtree(dest, ignore_errors=True)


def test_retry_and_soft404(fake: FakeUrlOpen, original_sleeps) -> None:
    print("[5] 重试 / 软 404 / 不可达")
    check_equal(original_sleeps, (1.0, 2.0), "默认退避 1s/2s")
    # 第一次 esearch 失败、第二次成功（RETRY_SLEEPS 已被替换为 0）
    fake.calls.clear()
    fake.fail("esearch.fcgi", times=1)
    result = bg_geo.geo_search("gastric cancer", limit=5)
    check_equal(result["total"], 297, "重试后 search 成功")
    esearch_calls = [url for url in fake.calls if "esearch.fcgi" in url]
    check_equal(len(esearch_calls), 2, "失败一次后重试一次")
    # 持续失败 → unsupported
    fake.fail("esearch.fcgi", times=99)
    try:
        error = check_raises(lambda: bg_geo.geo_search("gastric cancer"), "unsupported",
                             "esearch 持续失败→unsupported")
        if error is not None:
            check_in("NCBI 暂时不可达", error.message, "可读网络错误")
    finally:
        fake.fail_times.pop("esearch.fcgi", None)
    # suppl 目录 404（没有补充文件）→ matrix 文件仍可用
    fake.fail_code = 404
    fake.fail("suppl", times=99)
    try:
        detail = bg_geo.geo_detail("GSE13861")
        check_equal(len(detail["files"]), 1, "suppl 404 时只剩 matrix 文件")
        check_equal(detail["files"][0]["kind"], "series_matrix", "只剩 series matrix")
        check_equal(detail["suggested"], "GSE13861_series_matrix.txt.gz", "suggested 仍可用")
    finally:
        fake.fail_times.pop("suppl", None)
        fake.fail_code = 503


def test_dir_html_parser() -> None:
    print("[6] 目录 HTML 解析边界")
    files = bg_geo._parse_dir_html(MATRIX_DIR_HTML, "matrix")
    check_equal(len(files), 1, "Parent Directory 被忽略")
    check_equal(files[0]["url"],
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13861/matrix/GSE13861_series_matrix.txt.gz",
                "matrix 目录 URL 正确拼接（无重复 /geo/series）")
    files2 = bg_geo._parse_dir_html(SUPPL_DIR_HTML, "suppl")
    check_equal([item["name"] for item in files2],
                ["GSE13861_ClinicalInformation.xls.gz", "GSE13861_RawData.txt.gz",
                 "filelist.txt", "GSE13861_RAW.tar"], "suppl 文件顺序")
    check_equal(files2[0]["size"], int(9.4 * 1024), "9.4K 大小解析")
    check_equal(bg_geo._parse_size("20M"), 20 * 1024 ** 2, "20M 大小解析")
    check_equal(bg_geo._parse_size("-"), 0, "未知大小→0")
    check_equal(bg_geo._parse_dir_html("", "matrix"), [], "空 HTML→[]")


def main() -> int:
    fake = FakeUrlOpen()
    original_urlopen = bg_geo._urlopen
    original_sleeps = bg_geo.RETRY_SLEEPS
    bg_geo._urlopen = fake  # type: ignore[assignment]
    bg_geo.RETRY_SLEEPS = (0.0, 0.0)
    try:
        test_dictionary()
        test_search(fake)
        test_detail(fake)
        test_download(fake)
        test_retry_and_soft404(fake, original_sleeps)
        test_dir_html_parser()
    finally:
        bg_geo._urlopen = original_urlopen  # type: ignore[assignment]
        bg_geo.RETRY_SLEEPS = original_sleeps
    if FAILURES:
        print(f"\nFAIL {len(FAILURES)} checks\n" + "\n".join(f"  - {item}" for item in FAILURES))
        return 1
    print(f"PASS {CHECKS} checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
