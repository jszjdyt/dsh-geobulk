#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-geobulk — 下载原子替换/完整性校验 + 基因标识符映射报告（任务单第七节建议 1、2）

建议 1：`geo_download` 覆盖下载不再先删旧文件；先写 `<name>.part`，校验通过
（写入字节数 == Content-Length；缺失 Content-Length 时用 gzip/tar 可读性兜底）后
`os.replace` 原子替换；失败/中断保留旧文件并清理 `.part`；不实现断点续传（文档写明
「仅支持失败重试」）。
建议 2：独立的基因标识符 ID→symbol 映射报告（数据库/版本/物种/覆盖率/一对多与
多对一规则）；未知 ID 与「按 0 填充的缺失表达」都必须显式说明，不能当成真实零表达；
另有 `annotation_intersection` 供富集/免疫/单细胞匹配前显示 ID 类型、交集基因数与覆盖率。

离线：网络层参考 test/py_geo_test.py 的 FakeUrlOpen。R 注释用可注入的 FakeMapper 保证
确定性；文件末尾另有一条「本机 R + org.*.eg.db 可用时才跑」的真实映射集成检查。

用法：cd plugins/dsh-geobulk && python3 test/py_geo_io_regression_test.py
"""

from __future__ import annotations

import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.error
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lib"))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="bg-mpl-"))

import bg_geo  # noqa: E402
import bg_io  # noqa: E402

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
    return check(needle in haystack, label, f"{needle!r} 不在 {str(haystack)[:200]!r}")


def check_raises(fn, code: str | None, label: str):
    global CHECKS
    CHECKS += 1
    try:
        fn()
    except bg_geo.EngineError as error:
        if code and error.code != code:
            FAILURES.append(f"{label}：错误码 {error.code} != {code}（{error.message}）")
            print(f"  FAIL {label}：错误码 {error.code} != {code}")
        return error
    except Exception as error:  # noqa: BLE001
        FAILURES.append(f"{label}：抛了 {type(error).__name__} 而不是 EngineError：{error}")
        print(f"  FAIL {label}：抛了 {type(error).__name__} 而不是 EngineError")
        return None
    FAILURES.append(f"{label}：没有抛异常")
    print(f"  FAIL {label}：没有抛异常")
    return None


# ---------------------------------------------------------------------------
# 假网络（参考 test/py_geo_test.py 的 FakeUrlOpen）
# ---------------------------------------------------------------------------

ESUMMARY_JSON = {
    "header": {"type": "esummary"},
    "result": {
        "uids": ["200013861"],
        "200013861": {
            "uid": "200013861", "accession": "GSE13861", "gse": "13861",
            "title": "Gene expression signature-based novel prognostic risk score in gastric cancer",
            "summary": "Clinical survival and prognosis of gastric cancer patients were analyzed.",
            "gpl": "6884", "taxon": "Homo sapiens", "entrytype": "GSE",
            "gdstype": "Expression profiling by array", "pdat": "2011/07/19",
            "suppfile": "TXT, XLS", "n_samples": 90,
            "pubmedids": ["21447720"], "samples": [],
            "ftplink": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13861/",
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
<a href="GSE13861_RawData.txt.gz">GSE13861_RawData.txt.gz</a>                                         2008-12-10 13:51   20M
<a href="filelist.txt">filelist.txt</a>                                                     2013-01-17 14:47  172
<a href="GSE13861_RAW.tar">GSE13861_RAW.tar</a>                                                 2013-01-17 14:47  6.3M
<hr></pre></body></html>"""

MATRIX_HEAD_TEXT = """!Series_title\t"Gene expression signature-based novel prognostic risk score in gastric cancer"
!Series_geo_accession\t"GSE13861"
!Series_overall_design\t"65 primary gastric adenocarcinoma and 19 normal tissues were profiled."
!Sample_organism_ch1\t"Homo sapiens"
!Sample_characteristics_ch1\t"age: 62"
!series_matrix_table_begin
"ID_REF"\t"GSM348407"\t"GSM348370"
"""

FILE_BYTES = b"hello geo download\n"

SUPPL_ROOT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13861/suppl"
URL_FILELIST = f"{SUPPL_ROOT}/filelist.txt"
URL_RAWDATA = f"{SUPPL_ROOT}/GSE13861_RawData.txt.gz"
URL_GZ_CUSTOM = f"{SUPPL_ROOT}/custom_matrix.txt.gz"
URL_TAR_CUSTOM = f"{SUPPL_ROOT}/custom_raw.tar"


class FakeResponse:
    """最小 HTTP 响应；headers 可自定义 Content-Length，可模拟中途断网。"""

    def __init__(self, data: bytes = b"", headers: dict | None = None, fail_mid_read: bool = False):
        self._data = data or b""
        self._pos = 0
        self.headers = {"Content-Length": str(len(self._data))} if headers is None else headers
        self.status = 200
        self.fail_mid_read = fail_mid_read
        self._reads = 0

    def read(self, size: int | None = None) -> bytes:
        self._reads += 1
        if self.fail_mid_read and self._reads > 1:
            raise ConnectionResetError("fake connection reset by peer")
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
    """按 URL 子串路由的假 urlopen；未配置的 URL 一律 404。"""

    def __init__(self):
        self.calls: list[str] = []
        self.routes: dict[str, object] = {}
        self.fail_code = 404

    def route(self, needle: str, factory) -> None:
        self.routes[needle] = factory

    def _default_payload(self, url: str):
        if "esearch.fcgi" in url:
            return json.dumps({"header": {}, "esearchresult": {"idlist": [], "count": "0"}}).encode()
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
        return None

    def __call__(self, request, timeout=None):
        url = request.full_url if hasattr(request, "full_url") else str(request)
        self.calls.append(url)
        for needle, factory in self.routes.items():
            if needle in url:
                return factory()  # 工厂每次返回新响应，保证重试可重复
        payload = self._default_payload(url)
        if payload is None:
            raise urllib.error.HTTPError(url, self.fail_code, "fake not found", None, None)
        return FakeResponse(payload)


def response_factory(data: bytes, headers: dict | None = None, fail_mid_read: bool = False):
    return lambda: FakeResponse(data, headers=headers, fail_mid_read=fail_mid_read)


def stray_parts(dest: str) -> list[str]:
    try:
        return sorted(name for name in os.listdir(dest) if name.endswith(".part"))
    except OSError:
        return []


def write_old(dest: str, name: str, content: bytes) -> str:
    path = os.path.join(dest, name)
    with open(path, "wb") as handle:
        handle.write(content)
    return path


def read_bytes(path: str) -> bytes:
    with open(path, "rb") as handle:
        return handle.read()


# ---------------------------------------------------------------------------
# 建议 1：下载完整性校验 + 原子替换
# ---------------------------------------------------------------------------

def test_download_doc_contract() -> None:
    print("[S1] 下载契约文档（.part / 原子替换 / 仅支持失败重试）")
    doc = bg_geo.__doc__ or ""
    check_in(".part", doc, "S1-0a 模块文档写明 .part 临时文件")
    check_in("原子替换", doc, "S1-0b 模块文档写明原子替换")
    check_in("仅支持失败重试", doc, "S1-0c 文档写明仅支持失败重试")
    check_in("断点续传", doc, "S1-0d 文档说明未实现断点续传")


def test_download_short_response(fake: FakeUrlOpen) -> None:
    print("[S1] 短响应（实际字节 < Content-Length）保留旧文件")
    dest = tempfile.mkdtemp(prefix="bg-dl-short-")
    target = write_old(dest, "filelist.txt", b"OLD-CONTENT")
    fake.route("/suppl/filelist.txt",
               response_factory(b"short", headers={"Content-Length": str(len(FILE_BYTES) + 100)}))
    error = check_raises(
        lambda: bg_geo.geo_download("GSE13861", [URL_FILELIST], dest, overwrite=True),
        "unsupported", "S1a 短响应→下载判失败（unsupported）")
    if error is not None:
        check_in("校验", error.message, "S1a2 失败原因说明是完整性校验")
    check_equal(read_bytes(target), b"OLD-CONTENT", "S1a3 旧文件未被短响应覆盖")
    check_equal(stray_parts(dest), [], "S1a4 .part 已清理")
    check(not os.path.exists(target + ".part"), "S1a5 目标 .part 不存在")
    shutil.rmtree(dest, ignore_errors=True)


def test_download_midstream_disconnect(fake: FakeUrlOpen) -> None:
    print("[S1] 中途断网（read 抛异常）保留旧文件")
    dest = tempfile.mkdtemp(prefix="bg-dl-reset-")
    target = write_old(dest, "filelist.txt", b"OLD-CONTENT")
    fake.route("/suppl/filelist.txt",
               response_factory(FILE_BYTES, headers={"Content-Length": str(len(FILE_BYTES))},
                                fail_mid_read=True))
    error = check_raises(
        lambda: bg_geo.geo_download("GSE13861", [URL_FILELIST], dest, overwrite=True),
        "unsupported", "S1b 中途断网→下载判失败")
    if error is not None:
        check_in("NCBI", error.message, "S1b2 网络错误消息保持中文可读")
    check_equal(read_bytes(target), b"OLD-CONTENT", "S1b3 旧文件未被删除/截断")
    check_equal(stray_parts(dest), [], "S1b4 .part 已清理")
    shutil.rmtree(dest, ignore_errors=True)


def test_download_overwrite_success(fake: FakeUrlOpen) -> None:
    print("[S1] 覆盖已有文件成功 → 原子替换为新内容")
    dest = tempfile.mkdtemp(prefix="bg-dl-ok-")
    target = write_old(dest, "filelist.txt", b"OLD-CONTENT")
    new_bytes = b"new-content-v2\n" * 37
    fake.route("/suppl/filelist.txt", response_factory(new_bytes))
    result = bg_geo.geo_download("GSE13861", [URL_FILELIST], dest, overwrite=True)
    entry = result["files"][0]
    check(entry.get("ok") is True and not entry.get("skipped"), "S1c 覆盖下载 ok 且非 skipped")
    check_equal(entry["bytes"], len(new_bytes), "S1c2 返回写入字节数")
    check_equal(read_bytes(target), new_bytes, "S1c3 内容正确（旧文件被原子替换）")
    check_equal(stray_parts(dest), [], "S1c4 成功后没有残留 .part")
    # overwrite=False 且文件已存在 → 跳过，旧内容保持不变
    fake.route("/suppl/filelist.txt", response_factory(b"SHOULD-NOT-BE-WRITTEN"))
    again = bg_geo.geo_download("GSE13861", [URL_FILELIST], dest, overwrite=False)
    check(again["files"][0].get("skipped") is True, "S1c5 overwrite=False 已存在→skipped")
    check_equal(read_bytes(target), new_bytes, "S1c6 skipped 时旧文件未被改写")
    shutil.rmtree(dest, ignore_errors=True)


def test_download_mixed_batch_keeps_old(fake: FakeUrlOpen) -> None:
    print("[S1] 批量下载中一个失败 → 失败文件保留旧文件并回报")
    dest = tempfile.mkdtemp(prefix="bg-dl-mixed-")
    target_ok = write_old(dest, "filelist.txt", b"OLD-FILELIST")
    target_bad = write_old(dest, "GSE13861_RawData.txt.gz", b"OLD-RAWDATA")
    good = b"good-new-bytes\n" * 11
    fake.route("/suppl/filelist.txt", response_factory(good))
    fake.route("/suppl/GSE13861_RawData.txt.gz",
               response_factory(b"truncated-rawdata",
                                headers={"Content-Length": str(20 * 1024 ** 2)}))
    result = bg_geo.geo_download("GSE13861", [URL_FILELIST, URL_RAWDATA], dest, overwrite=True)
    entries = {item["name"]: item for item in result["files"]}
    check(entries["filelist.txt"].get("ok") is True, "S1d 第一个文件成功")
    check(entries["GSE13861_RawData.txt.gz"].get("ok") is False, "S1d2 短响应文件回报 ok=false")
    check(entries["GSE13861_RawData.txt.gz"].get("bytes") == 0, "S1d3 失败文件 bytes=0")
    check(result.get("warnings"), "S1d4 部分失败给出 warnings")
    check_equal(read_bytes(target_ok), good, "S1d5 成功文件内容正确")
    check_equal(read_bytes(target_bad), b"OLD-RAWDATA", "S1d6 失败文件旧内容完整保留")
    check_equal(stray_parts(dest), [], "S1d7 .part 已清理")
    shutil.rmtree(dest, ignore_errors=True)


def test_download_gzip_tar_fallback(fake: FakeUrlOpen) -> None:
    print("[S1] 缺 Content-Length 时用 gzip/tar 可读性兜底校验")
    good_gz = gzip.compress(b"gene\tS1\nTP53\t1\n" * 64)
    dest = tempfile.mkdtemp(prefix="bg-dl-gz-")
    target = write_old(dest, "custom_matrix.txt.gz", b"OLD-GZ")
    fake.route("custom_matrix.txt.gz", response_factory(good_gz, headers={}))
    result = bg_geo.geo_download("GSE13861", [URL_GZ_CUSTOM], dest, overwrite=True)
    check(result["files"][0].get("ok") is True, "S1e 无 Content-Length 的完整 gzip 被接受")
    check_equal(read_bytes(target), good_gz, "S1e2 gzip 内容正确")
    check_equal(stray_parts(dest), [], "S1e3 .part 已清理")

    # 截断 gzip：无 Content-Length，必须靠 gzip 可读性兜底发现
    dest2 = tempfile.mkdtemp(prefix="bg-dl-gzbad-")
    target2 = write_old(dest2, "custom_matrix.txt.gz", b"OLD-GZ-BAD")
    truncated_gz = good_gz[:max(1, len(good_gz) - 40)]
    fake.route("custom_matrix.txt.gz", response_factory(truncated_gz, headers={}))
    error = check_raises(
        lambda: bg_geo.geo_download("GSE13861", [URL_GZ_CUSTOM], dest2, overwrite=True),
        "unsupported", "S1e4 截断 gzip（无 Content-Length）→判失败")
    if error is not None:
        check_in("校验", error.message, "S1e5 失败原因说明 gzip 校验")
    check_equal(read_bytes(target2), b"OLD-GZ-BAD", "S1e6 截断 gzip 不覆盖旧文件")
    check_equal(stray_parts(dest2), [], "S1e7 .part 已清理")

    # 截断 tar：无 Content-Length
    tar_bytes = _make_tar_bytes()
    dest3 = tempfile.mkdtemp(prefix="bg-dl-tar-")
    target3 = write_old(dest3, "custom_raw.tar", b"OLD-TAR")
    fake.route("custom_raw.tar", response_factory(tar_bytes[:max(1, len(tar_bytes) // 2)], headers={}))
    error = check_raises(
        lambda: bg_geo.geo_download("GSE13861", [URL_TAR_CUSTOM], dest3, overwrite=True),
        "unsupported", "S1e8 截断 tar（无 Content-Length）→判失败")
    if error is not None:
        check_in("校验", error.message, "S1e9 失败原因说明 tar 校验")
    check_equal(read_bytes(target3), b"OLD-TAR", "S1e10 截断 tar 不覆盖旧文件")
    check_equal(stray_parts(dest3), [], "S1e11 .part 已清理")

    # 完整 tar（无 Content-Length）→ 接受
    dest4 = tempfile.mkdtemp(prefix="bg-dl-tarok-")
    target4 = write_old(dest4, "custom_raw.tar", b"OLD-TAR-OK")
    fake.route("custom_raw.tar", response_factory(tar_bytes, headers={}))
    result4 = bg_geo.geo_download("GSE13861", [URL_TAR_CUSTOM], dest4, overwrite=True)
    check(result4["files"][0].get("ok") is True, "S1e12 无 Content-Length 的完整 tar 被接受")
    check_equal(read_bytes(target4), tar_bytes, "S1e13 tar 内容正确")
    check_equal(stray_parts(dest4), [], "S1e14 .part 已清理")
    for path in (dest, dest2, dest3, dest4):
        shutil.rmtree(path, ignore_errors=True)


def _make_tar_bytes() -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        payload = b"hello tar member\n" * 500
        info = tarfile.TarInfo("member.txt")
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# 建议 2：基因标识符与物种确认
# ---------------------------------------------------------------------------

REQUIRED_REPORT_FIELDS = {
    "id_type", "id_type_label", "id_type_counts",
    "species", "species_source", "species_consistent",
    "database", "database_version", "database_keytype",
    "mapping_available", "mapping_note",
    "n_genes", "n_recognized", "n_mapped", "n_unknown", "coverage",
    "rules", "one_to_many", "many_to_one",
    "mapped_examples", "unknown_examples", "zero_imputation", "warnings",
}

HUMAN_SYMBOLS = ["TP53", "EGFR", "MYC", "BRCA1", "NOTCH1", "UNKNOWNXYZ1"]


class FakeMapper:
    """可注入的 ID→symbol 映射器（模拟 R / AnnotationDbi 的返回结构）。"""

    def __init__(self, table: dict[str, str], species: str = "human", available: bool = True):
        self.table = dict(table)
        self.species = species
        self.available = available
        self.calls: list[dict] = []

    def __call__(self, genes, species: str = "human", id_type: str = "symbol"):
        self.calls.append({"genes": list(genes), "species": species, "id_type": id_type})
        database = "org.Mm.eg.db" if species == "mouse" else "org.Hs.eg.db"
        if not self.available:
            return {"available": False, "database": database, "database_version": "",
                    "keytype": "", "mapping": {}, "one_to_many": 0, "many_to_one": 0,
                    "reason": "fake mapper unavailable"}
        mapping = {gene: self.table[gene] for gene in genes if gene in self.table}
        return {
            "available": True,
            "database": database,
            "database_version": "3.19.1",
            "keytype": {"symbol": "SYMBOL", "ensembl_gene": "ENSEMBL"}.get(id_type, "SYMBOL"),
            "mapping": mapping,
            "one_to_many": 0,
            "many_to_one": 0,
            "reason": "",
        }


HUMAN_TABLE = {symbol: symbol for symbol in HUMAN_SYMBOLS if symbol != "UNKNOWNXYZ1"}


def _missing_report_api() -> bool:
    return not hasattr(bg_io, "gene_annotation_report")


def test_annotation_report_symbol() -> None:
    print("[S2] 已是 gene symbol → 映射报告")
    if _missing_report_api():
        check(False, "S2a 存在 gene_annotation_report", "缺少 bg_io.gene_annotation_report")
        return
    mapper = FakeMapper(HUMAN_TABLE)
    report = bg_io.gene_annotation_report(HUMAN_SYMBOLS, mapper=mapper)
    check_equal(REQUIRED_REPORT_FIELDS - set(report), set(), "S2a 报告字段齐全")
    check_equal(report.get("id_type"), "symbol", "S2b id_type=symbol")
    check_equal(report.get("species"), "human", "S2c 物种由标识符推断为 human")
    check_equal(report.get("species_source"), "inferred", "S2d species_source=inferred")
    check_equal(report.get("database"), "org.Hs.eg.db", "S2e 记录映射数据库")
    check(report.get("database_version") in ("3.19.1",), "S2f 记录数据库版本",
          str(report.get("database_version")))
    check_equal(report.get("database_keytype"), "SYMBOL", "S2g 记录 keytype")
    check(report.get("mapping_available") is True, "S2h mapping_available=true")
    check_equal(report.get("n_genes"), 6, "S2i n_genes")
    check_equal(report.get("n_mapped"), 5, "S2j n_mapped")
    check_equal(report.get("n_unknown"), 1, "S2k n_unknown")
    check_equal(report.get("coverage"), round(5 / 6, 4), "S2l 覆盖率 = 5/6")
    check_equal(report.get("unknown_examples"), ["UNKNOWNXYZ1"], "S2m 未知 ID 示例明确")
    check_equal(report.get("mapped_examples", {}).get("TP53"), "TP53", "S2n mapped_examples 示例")
    rules = report.get("rules") or {}
    check(isinstance(rules, dict) and {"one_to_many", "many_to_one", "unmapped"} <= set(rules),
          "S2o 一对多/多对一/未映射规则齐全", str(rules)[:120])
    check(all(isinstance(value, str) and value.strip() for value in rules.values()),
          "S2p 规则文字非空")
    check_in("first", json.dumps(rules, ensure_ascii=False), "S2q 一对多规则说明取值策略（multiVals=first）")
    zero = report.get("zero_imputation") or {}
    check_in("真实零表达", json.dumps(zero, ensure_ascii=False), "S2r 报告区分「0 填充 ≠ 真实零表达」")
    check(isinstance(report.get("warnings"), list) and report["warnings"],
          "S2s 未映射 ID 有 warnings")
    check_in("UNKNOWNXYZ1", " ".join(report["warnings"]), "S2t warnings 点名未知 ID")
    check_equal(len(mapper.calls), 1, "S2u 只调用一次映射")
    check_equal(mapper.calls[0]["species"], "human", "S2v 映射时使用确认后的物种")
    check_equal(bg_io._classify_gene_id("X00351"), "accession", "S2w GenBank 登录号→accession")
    check_equal(bg_io._classify_gene_id("ENST00000269305"), "ensembl_transcript",
                "S2x Ensembl 转录本→ensembl_transcript")
    for symbol in ("C11orf21", "CXorf38", "LINC01234", "DKFZp451B082", "HBB-3'HS1"):
        check_equal(bg_io._classify_gene_id(symbol), "symbol", f"S2x2 合法 gene symbol→symbol（{symbol}）")
    for label in ("adenylate", "1-Sep", "ALOXs"):
        check_equal(bg_io._classify_gene_id(label), "unknown", f"S2x3 非 ID 标签→unknown（{label}）")
    acc_mapper = FakeMapper({"X00351": "ACTB", "M17885": "RPLP0"})
    acc_report = bg_io.gene_annotation_report(["X00351", "M17885", "not_an_accession"],
                                             species="human", mapper=acc_mapper)
    check_equal(acc_report.get("id_type"), "accession", "S2y GenBank 登录号报告 id_type")
    check_equal(acc_report.get("coverage"), round(2 / 3, 4), "S2z accession 覆盖率 2/3")


def test_annotation_report_probe_mixed() -> None:
    print("[S2] GPL probe ID 与 symbol 混合 → 不猜测 probe 映射")
    if _missing_report_api():
        check(False, "S2A 存在 gene_annotation_report", "缺少 bg_io.gene_annotation_report")
        return
    genes = ["1007_s_at", "1053_at", "117_at", "AFFX-BioB-5_at", "TP53", "EGFR"]
    mapper = FakeMapper({"TP53": "TP53", "EGFR": "EGFR"})
    report = bg_io.gene_annotation_report(genes, platform_species="Homo sapiens", mapper=mapper)
    check_equal(report.get("id_type"), "mixed", "S2A1 混合标识符→id_type=mixed")
    counts = report.get("id_type_counts") or {}
    check_equal(counts.get("probe_id"), 4, "S2A2 识别 4 个 probe ID")
    check_equal(counts.get("symbol"), 2, "S2A3 识别 2 个 symbol")
    check_equal(report.get("species"), "human", "S2A4 物种来自 GEO 平台元数据（platform_species）")
    check_equal(report.get("species_source"), "platform", "S2A5 species_source=platform")
    check_equal(report.get("n_mapped"), 2, "S2A6 只映射可识别的 symbol 子集")
    check_equal(report.get("n_unknown"), 4, "S2A7 probe ID 计为未映射")
    check_equal(report.get("coverage"), round(2 / 6, 4), "S2A8 覆盖率 = 2/6")
    joined = " ".join(report.get("warnings") or []) + " " + str(report.get("mapping_note") or "")
    check_in("GPL", joined, "S2A9 说明 probe 需要 GPL 注释表")
    check_in("1007_s_at", " ".join(report.get("unknown_examples") or []), "S2A10 未映射示例含 probe ID")
    check_equal(len(mapper.calls), 1, "S2A11 映射只调用一次")
    check(not (set(mapper.calls[0]["genes"]) & set(genes[:4])), "S2A12 probe ID 未送去做猜测映射")

    # 主导类型（≥90%）可映射：真实矩阵里少量 GenBank 登录号不应拖垮整表映射
    dominant = ["TP53", "EGFR", "MYC", "BRCA1", "NOTCH1", "KIT", "MET", "RB1", "PTEN",
                "KRAS", "AF107846"]
    dominant_mapper = FakeMapper({gene: gene for gene in dominant[:10]})
    report2 = bg_io.gene_annotation_report(dominant, mapper=dominant_mapper)
    check_equal(report2.get("id_type"), "mixed", "S2A13 少量登录号→整体记为 mixed")
    check_equal(report2.get("dominant_id_type"), "symbol", "S2A14 记录主导 ID 类型 symbol")
    check_equal(report2.get("mapping_available"), True, "S2A15 主导类型≥90% 时仍执行映射")
    check_equal(report2.get("n_mapped"), 10, "S2A16 只映射主导类型子集")
    check_equal(report2.get("coverage"), round(10 / 11, 4), "S2A17 覆盖率 = 10/11")
    check_equal(dominant_mapper.calls[0]["id_type"], "symbol", "S2A18 映射使用主导 keytype")
    check_in("主导类型", report2.get("mapping_note") or "", "S2A19 note 说明只映射主导类型")
    # 两种可映射类型势均力敌 → 不猜 keytype
    ambiguous = ["TP53", "EGFR", "AF107846", "AI582773"]
    ambiguous_mapper = FakeMapper({gene: gene for gene in ambiguous})
    report3 = bg_io.gene_annotation_report(ambiguous, mapper=ambiguous_mapper)
    check_equal(report3.get("mapping_available"), False, "S2A20 无主导类型时不猜映射")
    check_equal(report3.get("coverage"), 0.0, "S2A21 无主导类型覆盖率 0.0")
    check_equal(ambiguous_mapper.calls, [], "S2A22 无主导类型不调用映射器")
    check_in("主导类型", " ".join(report3.get("warnings") or []), "S2A23 warnings 说明缺少主导类型")


def test_annotation_report_unknown() -> None:
    print("[S2] 未知 ID → 退化报告（不硬造映射）")
    if _missing_report_api():
        check(False, "S2B 存在 gene_annotation_report", "缺少 bg_io.gene_annotation_report")
        return
    genes = ["gene1", "abc", "x_y"]
    mapper = FakeMapper({})
    report = bg_io.gene_annotation_report(genes, mapper=mapper)
    check_equal(report.get("id_type"), "unknown", "S2B1 未知 ID→id_type=unknown")
    check_equal(report.get("mapping_available"), False, "S2B2 mapping_available=false")
    check_equal(report.get("database"), None, "S2B3 未映射时不伪造数据库")
    check_equal(report.get("n_mapped"), 0, "S2B4 n_mapped=0")
    check_equal(report.get("coverage"), 0.0, "S2B5 覆盖率 0.0")
    check_equal(mapper.calls, [], "S2B6 未知 ID 不调用映射器")
    joined = " ".join(report.get("warnings") or []) + " " + str(report.get("mapping_note") or "")
    check("无法识别" in joined or "未执行" in joined, "S2B7 warnings 说明 ID 类型无法识别", joined[:160])
    check(REQUIRED_REPORT_FIELDS - set(report) == set(), "S2B8 退化路径字段仍齐全")


def test_annotation_species_confirmation() -> None:
    print("[S2] 物种确认（用户指定 vs 标识符推断）")
    if _missing_report_api():
        check(False, "S2C 存在 gene_annotation_report", "缺少 bg_io.gene_annotation_report")
        return
    mapper = FakeMapper(HUMAN_TABLE, species="mouse")
    report = bg_io.gene_annotation_report(["TP53", "EGFR"], species="mouse", mapper=mapper)
    check_equal(report.get("species"), "mouse", "S2C1 用户指定物种优先")
    check_equal(report.get("species_source"), "user", "S2C2 species_source=user")
    check_equal(report.get("species_consistent"), False, "S2C3 与推断物种冲突时 species_consistent=false")
    check_in("不一致", " ".join(report.get("warnings") or []), "S2C4 warnings 提示物种冲突")
    check_equal(report.get("database"), "org.Mm.eg.db", "S2C5 用指定物种的注释库")
    check_equal(mapper.calls[0]["species"], "mouse", "S2C6 映射器收到确认后的物种")


def test_annotation_intersection() -> None:
    print("[S2] 进入富集/免疫/单细胞匹配前的交集报告")
    if not hasattr(bg_io, "annotation_intersection"):
        check(False, "S2D 存在 annotation_intersection", "缺少 bg_io.annotation_intersection")
        return
    mapper = FakeMapper(HUMAN_TABLE)
    out = bg_io.annotation_intersection(["TP53", "EGFR", "UNKNOWNXYZ1"],
                                        ["TP53", "MYC", "EGFR", "KRAS"], mapper=mapper)
    check_equal(out.get("id_type"), "symbol", "S2D1 id_type=symbol")
    check_equal(out.get("n_input"), 3, "S2D2 n_input")
    check_equal(out.get("n_mapped"), 2, "S2D3 n_mapped")
    check_equal(out.get("coverage"), round(2 / 3, 4), "S2D4 覆盖率 = 可映射/输入")
    check_equal(out.get("reference_size"), 4, "S2D5 reference_size")
    check_equal(out.get("n_intersection"), 2, "S2D6 交集基因数")
    check_equal(out.get("intersection_coverage"), round(2 / 2, 4), "S2D7 交集占可映射基因")
    check_equal(out.get("reference_coverage"), round(2 / 4, 4), "S2D8 交集占参考基因集")
    check_in("TP53", json.dumps(out, ensure_ascii=False), "S2D9 报告含交集示例")


def test_annotation_load_integration() -> None:
    print("[S2] load_source 集成：meta.gene_annotation + notes/warnings")
    if _missing_report_api():
        check(False, "S2E 存在 gene_annotation_report", "缺少 bg_io.gene_annotation_report")
        return
    base = tempfile.mkdtemp(prefix="bg-io-annot-")
    symbol_path = os.path.join(base, "symbols.tsv")
    with open(symbol_path, "w", encoding="utf-8") as handle:
        handle.write("gene\tS1\tS2\tS3\n")
        handle.write("TP53\t1\t2\tNA\n")
        handle.write("EGFR\t4\t5\t6\n")
        handle.write("MYC\t7\t8\t9\n")
        handle.write("BRCA1\t1\t1\t1\n")
        handle.write("NOTCH1\t2\t2\t2\n")
        handle.write("UNKNOWNXYZ1\t0\t0\t0\n")
    probe_path = os.path.join(base, "probes.tsv")
    with open(probe_path, "w", encoding="utf-8") as handle:
        handle.write("probe\tS1\tS2\tS3\n")
        handle.write("1007_s_at\t1\t2\t3\n")
        handle.write("1053_at\t4\t5\t6\n")
        handle.write("117_at\t7\t8\t9\n")
    unknown_path = os.path.join(base, "unknown.tsv")
    with open(unknown_path, "w", encoding="utf-8") as handle:
        handle.write("gene_id\tS1\tS2\tS3\n")
        handle.write("gene1\t1\t2\t3\n")
        handle.write("abc\t4\t5\t6\n")
        handle.write("x_y\t7\t8\t9\n")

    original = bg_io._r_map_gene_ids
    mapper = FakeMapper(HUMAN_TABLE)
    bg_io._r_map_gene_ids = mapper  # type: ignore[assignment]
    try:
        payload = bg_io.load_source(symbol_path, workdir=base,
                                    options={"transpose": True})
        annotation = (payload.get("meta") or {}).get("gene_annotation") or {}
        check_equal(annotation.get("id_type"), "symbol", "S2E1 导入结果 meta.gene_annotation 记录 ID 类型")
        check_equal(annotation.get("coverage"), round(5 / 6, 4), "S2E2 导入结果覆盖率 5/6")
        note_text = " ".join(payload.get("notes") or [])
        check_in("覆盖", note_text, "S2E3 notes 显示 ID 类型与覆盖率")
        warnings_text = " ".join(payload.get("warnings") or [])
        check_in("真实零表达", warnings_text, "S2E4 warnings 说明 0 填充不代表真实零表达")
        zero = annotation.get("zero_imputation") or {}
        check(zero.get("n_missing") == 1, "S2E5 记录按 0 填充的缺失值数量",
              json.dumps(zero, ensure_ascii=False)[:120])
        check_equal(read_bytes(symbol_path).count(b"NA"), 1, "S2E6 输入确有 1 个 NA")

        mapper.calls.clear()
        payload_probe = bg_io.load_source(probe_path, workdir=base, options={"transpose": True})
        annotation_probe = (payload_probe.get("meta") or {}).get("gene_annotation") or {}
        check_equal(annotation_probe.get("id_type"), "probe_id", "S2E7 probe 矩阵识别为 probe_id")
        check_equal(annotation_probe.get("mapping_available"), False,
                    "S2E8 probe 矩阵不硬造映射")
        check_equal(annotation_probe.get("coverage"), 0.0, "S2E9 probe 覆盖率 0.0")
        check_equal(mapper.calls, [], "S2E10 probe-only 矩阵未调用 symbol 映射")
        probe_text = " ".join(payload_probe.get("warnings") or []) + " " + \
                     " ".join(payload_probe.get("notes") or [])
        check_in("GPL", probe_text, "S2E11 probe 提示需要 GPL 注释表")

        payload_unknown = bg_io.load_source(unknown_path, workdir=base, options={"transpose": True})
        annotation_unknown = (payload_unknown.get("meta") or {}).get("gene_annotation") or {}
        check_equal(annotation_unknown.get("id_type"), "unknown", "S2E12 未知 ID 导入 id_type=unknown")
        check_equal(annotation_unknown.get("coverage"), 0.0, "S2E13 未知 ID 覆盖率 0.0")
        unknown_text = " ".join(payload_unknown.get("warnings") or [])
        check("无法识别" in unknown_text or "未映射" in unknown_text or "未执行" in unknown_text,
              "S2E14 未知 ID 有明确说明", unknown_text[:160])
    finally:
        bg_io._r_map_gene_ids = original  # type: ignore[assignment]
        shutil.rmtree(base, ignore_errors=True)


def test_annotation_inspect() -> None:
    print("[S2] inspect 阶段显示 ID 类型（不建矩阵、不跑 R）")
    base = tempfile.mkdtemp(prefix="bg-io-inspect-")
    probe_path = os.path.join(base, "probes.tsv")
    with open(probe_path, "w", encoding="utf-8") as handle:
        handle.write("probe\tS1\tS2\n")
        handle.write("1007_s_at\t1\t2\n")
        handle.write("1053_at\t3\t4\n")
        handle.write("117_at\t5\t6\n")
    try:
        info = bg_io.inspect_source(probe_path, workdir=base, options={"transpose": True})
        check_equal(info.get("gene_id_type"), "probe_id", "S2F1 inspect 输出 ID 类型")
        check_in("GPL", " ".join(info.get("warnings") or []), "S2F2 inspect 提示 probe 需要 GPL")
    except Exception as error:  # noqa: BLE001
        check(False, "S2F inspect 可用", f"{type(error).__name__}: {error}")
    finally:
        shutil.rmtree(base, ignore_errors=True)


def test_annotation_r_unavailable_fallback() -> None:
    print("[S2] R 注释库不可用 → 退化为纯报告（不硬造映射）")
    if _missing_report_api():
        check(False, "S2H 存在 gene_annotation_report", "缺少 bg_io.gene_annotation_report")
        return
    original = bg_io._rscript
    bg_io._rscript = lambda: None  # type: ignore[assignment]
    try:
        report = bg_io.gene_annotation_report(["TP53", "EGFR", "MYC"], species="human")
    finally:
        bg_io._rscript = original  # type: ignore[assignment]
    check_equal(report.get("id_type"), "symbol", "S2H1 仍识别 ID 类型")
    check_equal(report.get("mapping_available"), False, "S2H2 R 不可用时 mapping_available=false")
    check_equal(report.get("coverage"), 0.0, "S2H3 覆盖率 0.0（不伪造映射）")
    check_equal(report.get("mapped_examples"), {}, "S2H4 没有伪造的映射示例")
    joined = " ".join(report.get("warnings") or []) + " " + str(report.get("mapping_note") or "")
    check_in("Rscript", joined, "S2H5 warnings 说明缺少 Rscript")


# ---------------------------------------------------------------------------
# 真实 R 映射（本机 R + org.*.eg.db 可用时才跑）
# ---------------------------------------------------------------------------

def _r_annotation_ready() -> bool:
    rscript = getattr(bg_io, "_rscript", lambda: None)()
    if not rscript:
        return False
    code = ('cat(all(vapply(c("AnnotationDbi","org.Hs.eg.db"), function(p) '
            'requireNamespace(p, quietly=TRUE), logical(1))))')
    try:
        completed = subprocess.run([rscript, "--vanilla", "-e", code],
                                   capture_output=True, text=True, timeout=120)
    except Exception:  # noqa: BLE001
        return False
    return completed.stdout.strip().upper() == "TRUE"


def test_annotation_real_r() -> None:
    print("[S2] 真实 R / org.Hs.eg.db 映射（无环境则 SKIP）")
    if not _r_annotation_ready():
        print("  SKIP 本机没有 Rscript 或 AnnotationDbi/org.Hs.eg.db")
        return
    if not hasattr(bg_io, "_r_map_gene_ids"):
        check(False, "S2G 存在 _r_map_gene_ids", "缺少 bg_io._r_map_gene_ids")
        return
    result = bg_io._r_map_gene_ids(["TP53", "EGFR", "MYC", "NOTAREALGENE_XYZ"], "human", "symbol")
    check_equal(result.get("available"), True, "S2G1 真实 R 映射可用")
    check_equal((result.get("mapping") or {}).get("TP53"), "TP53", "S2G2 TP53→TP53")
    check_equal((result.get("mapping") or {}).get("EGFR"), "EGFR", "S2G3 EGFR→EGFR")
    check("NOTAREALGENE_XYZ" not in (result.get("mapping") or {}), "S2G4 假基因未被映射")
    check(result.get("database") == "org.Hs.eg.db", "S2G5 记录数据库 org.Hs.eg.db",
          str(result.get("database")))
    check(bool(result.get("database_version")), "S2G6 记录数据库版本",
          str(result.get("database_version")))
    report = bg_io.gene_annotation_report(["TP53", "EGFR", "MYC", "NOTAREALGENE_XYZ"],
                                          species="human", mapper=bg_io._r_map_gene_ids)
    check_equal(report.get("coverage"), round(3 / 4, 4), "S2G7 真实 R 覆盖率 3/4")
    accession = bg_io._r_map_gene_ids(["X00351"], "human", "accession")
    check_equal(accession.get("available"), True, "S2G8 ACCNUM 映射可用")
    check_equal((accession.get("mapping") or {}).get("X00351"), "ACTB", "S2G9 X00351→ACTB")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    fake = FakeUrlOpen()
    original_urlopen = bg_geo._urlopen
    original_sleeps = bg_geo.RETRY_SLEEPS
    bg_geo._urlopen = fake  # type: ignore[assignment]
    bg_geo.RETRY_SLEEPS = (0.0, 0.0)
    # 每个下载用例前清空路由，避免用例间串味
    reset_routes = lambda: fake.routes.clear()  # noqa: E731
    try:
        print("== 建议 1：下载完整性校验 + 原子替换 ==")
        test_download_doc_contract()
        for case in (test_download_short_response, test_download_midstream_disconnect,
                     test_download_overwrite_success, test_download_mixed_batch_keeps_old,
                     test_download_gzip_tar_fallback):
            reset_routes()
            try:
                case(fake)
            except Exception as error:  # noqa: BLE001
                check(False, f"{case.__name__} 未预期异常",
                      f"{type(error).__name__}: {error}")
        print("== 建议 2：基因标识符与物种确认 ==")
        for case in (test_annotation_report_symbol, test_annotation_report_probe_mixed,
                     test_annotation_report_unknown, test_annotation_species_confirmation,
                     test_annotation_intersection, test_annotation_load_integration,
                     test_annotation_inspect, test_annotation_r_unavailable_fallback,
                     test_annotation_real_r):
            try:
                case()
            except Exception as error:  # noqa: BLE001
                check(False, f"{case.__name__} 未预期异常",
                      f"{type(error).__name__}: {error}")
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
