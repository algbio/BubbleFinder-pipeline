import os
from pathlib import Path
from snakemake.shell import shell

shell.executable("/bin/bash")

BASEDIR = Path(workflow.basedir) if "workflow" in globals() else Path(os.getcwd())
def resolve_env(p):
    p = str(p)
    return p if os.path.isabs(p) else str(BASEDIR / p)

DEFAULTS = config.get("defaults", {}) or {}
DATA_DIR = DEFAULTS.get("data_dir", "data")
OUT_DIR  = DEFAULTS.get("out_dir", "results")

ENV_PYTOOLS = resolve_env(DEFAULTS.get("envs", {}).get("pytools", "config/pytools.yml"))

SG_CONF       = DEFAULTS.get("sgraph", {}) or {}
SG_ENABLED    = bool(SG_CONF.get("enabled", False))  # OFF par défaut
SG_GZIP       = bool(SG_CONF.get("gzip", True))
SG_EDGES_ONLY = bool(SG_CONF.get("edges_only", True))
SG_ONE_BASED  = bool(SG_CONF.get("one_based", False))
SG_ADD_RC     = bool(SG_CONF.get("add_rc_edges", False))
SG_DEDUP      = bool(SG_CONF.get("dedup_edges", True))

SGRAPH_SUFFIX      = ".sgraph.gz" if SG_GZIP else ".sgraph"
SGRAPH_MAP_SUFFIX  = ".sgraph.map.tsv.gz" if SG_GZIP else ".sgraph.map.tsv"

def sgraph_path(dataset):
    return os.path.join(DATA_DIR, dataset, f"{dataset}{SGRAPH_SUFFIX}")

def sgraph_map_path(dataset):
    return os.path.join(DATA_DIR, dataset, f"{dataset}{SGRAPH_MAP_SUFFIX}")

def ds(name):
    for d in config.get("datasets", []):
        if d.get("name") == name:
            return d
    raise ValueError(f"Dataset not found in config: {name}")

def sgraph_enabled_for_dataset(name):
    d = ds(name)
    # Respect per-dataset override; fallback on global SG_ENABLED
    return bool((d.get("sgraph", {}) or {}).get("enabled", SG_ENABLED))

rule gfa_to_sgraph:
    message: "Convert cleaned GFA -> sgraph pour {wildcards.dataset}"
    input:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.cleaned.gfa")
    output:
        sgraph=os.path.join(DATA_DIR, "{dataset}", "{dataset}" + SGRAPH_SUFFIX),
        mapping=os.path.join(DATA_DIR, "{dataset}", "{dataset}" + SGRAPH_MAP_SUFFIX)
    conda: ENV_PYTOOLS
    log:
        os.path.join(OUT_DIR, "logs", "sgraph", "{dataset}.log")
    params:
        script=lambda wc: resolve_env("scripts/gfa_to_sgraph.py"),
        flags=lambda wc: " ".join(filter(None, [
            "--edges-only"   if ((ds(wc.dataset).get("sgraph", {}) or {}).get("edges_only",  SG_EDGES_ONLY)) else "",
            "--one-based"    if ((ds(wc.dataset).get("sgraph", {}) or {}).get("one_based",   SG_ONE_BASED))  else "",
            "--add-rc-edges" if ((ds(wc.dataset).get("sgraph", {}) or {}).get("add_rc_edges",SG_ADD_RC))     else "",
            "--dedup-edges"  if ((ds(wc.dataset).get("sgraph", {}) or {}).get("dedup_edges", SG_DEDUP))      else "",
        ]))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.sgraph})" "$(dirname {log})"
        python "{params.script}" "{input.gfa}" "{output.sgraph}" \
          --map-out "{output.mapping}" {params.flags} \
          2> >(tee -a "{log}" >&2)
        """

def _all_enabled_dataset_names():
    if "all_enabled_dataset_names" in globals() and callable(all_enabled_dataset_names):
        return list(all_enabled_dataset_names())
    return [d["name"] for d in (config.get("datasets", []) or []) if "name" in d]

def ALL_SGRAPHS():
    names = [n for n in _all_enabled_dataset_names() if sgraph_enabled_for_dataset(n)]
    return [sgraph_path(n) for n in names]

rule sgraphs:
    message: "Build all .sgraph pour datasets activés"
    input:
        ALL_SGRAPHS()