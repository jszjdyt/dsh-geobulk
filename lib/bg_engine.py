#!/usr/bin/env python3
"""
dsh-geobulk — Python worker（长驻子进程）

协议（CONTRACT.md §2）：
  stdin  : 一行一个 JSON 请求 {"id":"..","cmd":"..","args":{..}}
  stdout : 一行一个 JSON 响应 {"id":"..","ok":true,"result":{..}} /
           {"id":"..","ok":false,"error":{"code":"..","message":".."}}
           可选 {"id":"..","event":"progress","message":"..","pct":0-100}
  stdout 只允许协议行；日志一律 stderr。

启动：python3 -u lib/bg_engine.py --workdir <abs_dir>
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def _prepare_env(workdir: str) -> None:
    os.environ["MPLBACKEND"] = "Agg"
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(workdir, "mpl-cache"))
    os.environ["DS_BG_WORKDIR"] = workdir
    try:
        os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    except OSError:
        os.environ["MPLCONFIGDIR"] = "/tmp"
    for sub in ("figures", "exports", "uploads", "downloads", "analyses", "genesets", "sc-ref", "unpack"):
        try:
            os.makedirs(os.path.join(workdir, sub), exist_ok=True)
        except OSError:
            pass


class Worker:
    def __init__(self, workdir: str):
        self.workdir = workdir
        self.started = time.time()
        self.state = None
        self._engine_info_cache = None

    # -- 基础设施 ---------------------------------------------------------
    def load_state(self):
        if self.state is None:
            import bg_state

            self.state = bg_state.State()
        return self.state

    def engine_info(self) -> dict:
        if self._engine_info_cache is not None:
            return self._engine_info_cache
        import platform

        import numpy
        import pandas
        import scipy
        import sklearn

        packages = {
            "numpy": getattr(numpy, "__version__", "?"),
            "pandas": getattr(pandas, "__version__", "?"),
            "scipy": getattr(scipy, "__version__", "?"),
            "sklearn": getattr(sklearn, "__version__", "?"),
        }
        for name in ("matplotlib", "anndata", "lifelines", "statsmodels", "openpyxl", "h5py", "PIL"):
            try:
                mod = __import__(name)
                packages[name] = getattr(mod, "__version__", "?")
            except Exception:  # noqa: BLE001
                packages[name] = None
        r_info = {}
        try:
            import rbridge

            r_info = rbridge.r_info()
        except Exception as error:  # noqa: BLE001
            r_info = {"available": False, "error": str(error), "packages": {}}
        info = {
            "python": platform.python_version(),
            "platform": sys.platform,
            "mode": "builtin",
            "workdir": self.workdir,
            "packages": packages,
            "r": r_info,
            "started_at": self.started,
        }
        self._engine_info_cache = info
        return info

    # -- 命令 -------------------------------------------------------------
    def cmd_ping(self, args):
        import bg_registry

        info = dict(self.engine_info())
        module_status = bg_registry.module_status()
        info["modules"] = module_status
        info["capabilities"] = {
            "steps": bg_registry.steps_catalog(),
            "analyses": [a["id"] for a in bg_registry.analyses_catalog(self._empty_dataset_view())["analyses"]],
            "plots": ["qc_box", "qc_density", "pca_scatter", "sample_corr", "sample_dendro", "volcano",
                      "ma_plot", "de_heatmap", "gene_box", "enrich_bar", "enrich_dot", "immune_stack",
                      "immune_corr", "km_curve", "forest", "risk_plot", "roc_time", "batch_pca",
                      "circos", "splicing_heat"],
        }
        try:
            import bg_io

            info["capabilities"]["formats"] = list(getattr(bg_io, "SUPPORTED_FORMATS", []))
        except Exception:  # noqa: BLE001
            info["capabilities"]["formats"] = []
        info["geo"] = True
        return info

    def _empty_dataset_view(self):
        """供 ping 时计算分析可用性（无数据集 → 全部 needs 不满足）。"""

        class _Empty:
            id = ""
            name = ""
            modules = {}
            results = {}
            steps = []
            group_key = None
            sc_ref = None

            def has_clinical(self):
                return False

        return _Empty()

    def cmd_open(self, args):
        import bg_registry
        import bg_state

        state = self.load_state()
        if args.get("demo"):
            dataset = bg_state.make_demo(state, name=args.get("name") or "demo-bulk")
            return {"dataset": dataset.summary(), "preview": bg_registry.preview(dataset)}

        path = args.get("path") or args.get("url")
        if not path:
            raise _err("bad_request", "请提供 path（本地文件/目录）或 url")
        if not str(path).startswith(("http://", "https://", "ftp://")):
            abspath = os.path.abspath(str(path))
            existing = state.find_by_path(abspath)
            if existing is not None:
                return {"dataset": existing.summary(), "reused": True}
            path = abspath
            if not os.path.exists(path):
                raise _err("not_found", f"路径不存在：{path}")

        import bg_io

        options = {
            "data_format": args.get("data_format") or args.get("format") or "auto",
            "member": args.get("member"),
            "object": args.get("object"),
            "transpose": args.get("transpose"),
            "delimiter": args.get("delimiter"),
            "sheet": args.get("sheet"),
            "gene_col": args.get("gene_col"),
            "sample_sheet": args.get("sample_sheet"),
            "clinical": args.get("clinical"),
            "species": args.get("species") or "auto",
        }
        options = {k: v for k, v in options.items() if v is not None}
        loaded = bg_io.load_source(path, workdir=self.workdir, options=options)
        return self._make_dataset(loaded, name=args.get("name"), source=path)

    def _make_dataset(self, loaded, name=None, source=None):
        import bg_state

        state = self.load_state()
        meta = loaded.get("meta") or {}
        value_kind = loaded.get("value_kind") or meta.get("value_kind") or "unknown"
        sample_info = loaded.get("sample_info")
        X = loaded.get("X")
        counts = loaded.get("counts")
        if counts is None and str(value_kind) == "counts":
            counts = X
        ds_name = name or loaded.get("name") or (os.path.basename(str(source).rstrip(os.sep)) if source else "dataset")
        dataset = state.new(
            name=str(ds_name),
            X=X,
            genes=loaded.get("genes") or [],
            samples=loaded.get("samples") or [],
            counts=counts,
            sample_info=sample_info,
            source_path=str(source) if source else None,
            source_format=loaded.get("format") or meta.get("format") or "unknown",
            species=(meta.get("species") or "unknown"),
            value_kind=str(value_kind),
            notes=list(loaded.get("notes") or []),
            warnings=list(loaded.get("warnings") or []),
        )
        import bg_registry

        return {"dataset": dataset.summary(), "preview": bg_registry.preview(dataset)}

    def cmd_list(self, args):
        state = self.load_state()
        return {"datasets": state.summaries(), "active": state.active}

    def cmd_close(self, args):
        state = self.load_state()
        target = str(args.get("dataset") or state.active or "")
        if not target or target not in state.datasets:
            return {"closed": False, "datasets": state.summaries(), "active": state.active}
        state.close(target)
        return {"closed": True, "datasets": state.summaries(), "active": state.active}

    def cmd_select(self, args):
        state = self.load_state()
        target = str(args.get("dataset") or "")
        if target not in state.datasets:
            raise _err("not_found", f"数据集 {target} 不存在")
        state.active = target
        return {"active": state.active, "dataset": state.datasets[target].summary()}

    def cmd_reset(self, args):
        state = self.load_state()
        state.reset()
        return {"reset": True, "datasets": [], "active": None}

    def cmd_preview(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        return bg_registry.preview(dataset)

    def cmd_samples(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        payload = dict(args)
        payload.pop("dataset", None)
        return bg_registry.handle_samples(dataset, payload)

    def cmd_run(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        step = str(args.get("step") or "")
        if not step:
            raise _err("bad_request", "请提供 step")
        return bg_registry.run_step(dataset, step, args.get("params") or {}, self.workdir, self.engine_info())

    def cmd_plot(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        kind = str(args.get("kind") or "")
        return bg_registry.run_plot(dataset, kind, args.get("params") or {}, self.workdir, self.engine_info())

    def cmd_analyses(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        return bg_registry.analyses_catalog(dataset)

    def cmd_analysis(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        analysis_id = str(args.get("id") or "")
        if not analysis_id:
            raise _err("bad_request", "请提供分析 id")
        return bg_registry.run_analysis(dataset, analysis_id, args.get("params") or {}, self.workdir, self.engine_info())

    def cmd_attach(self, args):
        import bg_io
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        kind = str(args.get("module") or args.get("kind") or "")
        path = args.get("path") or args.get("url")
        if not path:
            raise _err("bad_request", "请提供模块数据路径（path）")
        if not str(path).startswith(("http://", "https://", "ftp://")):
            path = os.path.abspath(str(path))
            if not os.path.exists(path):
                raise _err("not_found", f"路径不存在：{path}")
        options = {"data_format": kind, "member": args.get("member"), "delimiter": args.get("delimiter"),
                   "sheet": args.get("sheet"), "gene_col": args.get("gene_col")}
        options = {k: v for k, v in options.items() if v is not None}
        loaded = bg_io.load_module(path, kind, workdir=self.workdir, options=options)
        return bg_registry.attach_module(dataset, kind, loaded, path)

    def cmd_export(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        return bg_registry.export_dataset(
            dataset,
            str(args.get("what") or "summary_json"),
            args.get("path"),
            self.workdir,
            self.engine_info(),
        )

    def cmd_audit(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        return bg_registry.audit(dataset)

    def cmd_suggest(self, args):
        import bg_registry

        dataset = self.load_state().get(args.get("dataset"))
        return bg_registry.suggest(dataset)

    def cmd_recipe(self, args):
        import bg_registry

        if not args.get("id"):
            return bg_registry.recipe(None)
        return bg_registry.recipe(str(args.get("id")))

    def cmd_steps(self, args):
        import bg_registry

        return {"steps": bg_registry.steps_catalog(), "modules": bg_registry.module_status()}

    def cmd_inspect(self, args):
        import bg_io

        path = args.get("path") or args.get("url")
        if not path:
            raise _err("bad_request", "请提供 path 或 url")
        if not str(path).startswith(("http://", "https://", "ftp://")):
            path = os.path.abspath(str(path))
        options = {
            "data_format": args.get("data_format") or args.get("format") or "auto",
            "member": args.get("member"),
            "object": args.get("object"),
            "delimiter": args.get("delimiter"),
            "sheet": args.get("sheet"),
        }
        options = {k: v for k, v in options.items() if v is not None}
        return bg_io.inspect_source(path, workdir=self.workdir, options=options)

    def cmd_formats(self, args):
        try:
            import bg_io

            return {"formats": getattr(bg_io, "FORMAT_INFO", []), "ids": list(getattr(bg_io, "SUPPORTED_FORMATS", []))}
        except Exception:  # noqa: BLE001
            return {"formats": [], "ids": []}

    # -- GEO --------------------------------------------------------------
    def cmd_geo_search(self, args):
        import bg_geo

        search = getattr(bg_geo, "geo_search")
        return search(
            term=str(args.get("term") or ""),
            species=str(args.get("species") or "all"),
            study_type=str(args.get("study_type") or "all"),
            min_samples=int(args.get("min_samples") or 0),
            has_clinical=bool(args.get("has_clinical")),
            limit=int(args.get("limit") or 20),
            page=int(args.get("page") or 0),
            sort=str(args.get("sort") or "relevance"),
        )

    def cmd_geo_detail(self, args):
        import bg_geo

        return bg_geo.geo_detail(str(args.get("gse") or ""))

    def cmd_geo_download(self, args):
        import bg_geo

        request_id = args.get("__request_id__")

        def progress(message, pct=None):
            if request_id:
                payload = {"id": request_id, "event": "progress", "message": str(message)}
                if pct is not None:
                    payload["pct"] = float(pct)
                _emit(payload)

        gse = str(args.get("gse") or "")
        files = args.get("files") or []
        dest = str(args.get("dest_dir") or os.path.join(self.workdir, "downloads", "GEO"))
        return bg_geo.geo_download(gse, files, dest, overwrite=bool(args.get("overwrite")), progress=progress)

    # -- 单细胞 -----------------------------------------------------------
    def cmd_sc_refs(self, args):
        path = args.get("path")
        if path:
            return self._sc_info(str(path))
        import glob

        roots = [
            os.path.join(os.getcwd(), "*.h5ad"),
            os.path.join(os.path.dirname(self.workdir.rstrip(os.sep)), "*.h5ad"),
            os.path.join(os.getcwd(), "plugins", "dsh-singlecell", "exports", "*.h5ad"),
        ]
        found = []
        for pattern in roots:
            for item in glob.glob(pattern):
                try:
                    size = os.path.getsize(item)
                except OSError:
                    continue
                found.append({"path": os.path.abspath(item), "name": os.path.basename(item), "bytes": size})
        dedup = {}
        for item in found:
            dedup[item["path"]] = item
        return {"refs": list(dedup.values())}

    def _sc_info(self, path: str):
        path = os.path.abspath(path)
        if not os.path.exists(path):
            raise _err("not_found", f"单细胞参考文件不存在：{path}")
        try:
            import anndata

            adata = anndata.read_h5ad(path, backed="r")
        except Exception as error:  # noqa: BLE001
            raise _err("bad_request", f"读取 h5ad 失败：{error}")
        obs_cols = []
        for col in list(adata.obs.columns):
            try:
                obs_cols.append({"key": str(col), "kind": "categorical" if adata.obs[col].dtype.name in ("category", "object", "bool") else "numeric",
                                 "n_unique": int(adata.obs[col].nunique())})
            except Exception:  # noqa: BLE001
                continue
        cell_types = None
        for key in ("cell_type", "celltype", "cell_type_l1", "clusters", "leiden", "seurat_clusters", "louvain"):
            if key in adata.obs.columns:
                try:
                    cell_types = [str(v) for v in list(dict.fromkeys(adata.obs[key].astype(str).tolist()))][:50]
                except Exception:  # noqa: BLE001
                    cell_types = None
                break
        result = {
            "path": path,
            "name": os.path.basename(path),
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "obs_columns": obs_cols,
            "cell_types": cell_types,
            "species": self._guess_species(list(adata.var_names[:200])),
        }
        try:
            adata.file.close()
        except Exception:  # noqa: BLE001
            pass
        return result

    @staticmethod
    def _guess_species(genes):
        upper = [g.upper() for g in genes]
        human = sum(1 for g in upper if g in ("CD3D", "CD8A", "MS4A1", "CD14", "EPCAM", "GAPDH", "ACTB"))
        mouse = sum(1 for g in upper if g in ("CD3E", "CD8A", "MS4A1", "CD14", "EPCAM", "GAPDH", "ACTB"))
        # 小鼠常见全大写基因风格与人类难区分：用 6 基因名启发（例如 Cd3d 之类混合大小写）
        mixed = sum(1 for g in genes if g[:1].isupper() and g[1:].islower())
        if mixed > len(genes) * 0.6:
            return "mouse"
        return "human" if human >= mouse else "unknown"

    def cmd_sc_set(self, args):
        dataset = self.load_state().get(args.get("dataset"))
        ref = str(args.get("ref_path") or "")
        if not ref:
            raise _err("bad_request", "请提供 ref_path")
        ref = os.path.abspath(ref)
        if not os.path.exists(ref):
            raise _err("not_found", f"单细胞参考文件不存在：{ref}")
        dataset.sc_ref = ref
        return {"dataset": dataset.summary(), "sc_ref": ref}

    def cmd_shutdown(self, args):
        return {"bye": True}

    # -- 分发 -------------------------------------------------------------
    def dispatch(self, cmd: str, args: dict):
        handler = getattr(self, f"cmd_{cmd}", None)
        if handler is None:
            raise _err("bad_request", f"未知命令：{cmd}")
        return handler(args or {})


class _EngineFailure(Exception):
    def __init__(self, code: str, message: str, **details):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


def _err(code: str, message: str, **details):
    return _EngineFailure(code, message, **details)


def _emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload, ensure_ascii=False) + "\n")
    sys.stdout.flush()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", default=os.path.join(os.path.dirname(HERE), "work"))
    parsed, _ = parser.parse_known_args()
    workdir = os.path.abspath(parsed.workdir)
    _prepare_env(workdir)

    import bg_common

    worker = Worker(workdir)
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        request_id = None
        try:
            payload = json.loads(line)
            request_id = payload.get("id")
            cmd = str(payload.get("cmd") or "")
            args = payload.get("args") or {}
            if cmd == "geo_download":
                args = dict(args)
                args["__request_id__"] = request_id
            result = worker.dispatch(cmd, args)
            _emit({"id": request_id, "ok": True, "result": bg_common.jsonable(result)})
            if cmd == "shutdown":
                return 0
        except _EngineFailure as error:
            _emit({"id": request_id, "ok": False,
                   "error": {"code": error.code, "message": error.message, **bg_common.jsonable(error.details)}})
        except bg_common.EngineError as error:
            _emit({"id": request_id, "ok": False, "error": error.to_dict()})
        except Exception as error:  # noqa: BLE001
            trace = traceback.format_exc()
            sys.stderr.write(trace + "\n")
            sys.stderr.flush()
            _emit({"id": request_id, "ok": False,
                   "error": {"code": "internal", "message": f"{type(error).__name__}: {error}",
                             "trace": trace[-2000:]}})
    return 0


if __name__ == "__main__":
    sys.exit(main())
