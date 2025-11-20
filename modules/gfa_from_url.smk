# modules/gfa_from_url.smk
# Download + (de)compress a remote GFA (supports .zst and .gz) into
# data/<dataset>/<dataset>.raw.gfa
#
# Usage: add a dataset with builder: gfa_from_url and either urls.url or urls.files

import os
import shlex
import shutil
import subprocess
import gzip
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
from snakemake.exceptions import WorkflowError
from snakemake.shell import shell
from pathlib import Path

shell.executable("/bin/bash")

BASEDIR = Path(workflow.basedir) if "workflow" in globals() else Path(os.getcwd())
def resolve_env(p):
    p = str(p)
    return p if os.path.isabs(p) else str(BASEDIR / p)

DEFAULTS = config.get("defaults", {}) or {}
DATA_DIR = DEFAULTS.get("data_dir", "data")
OUT_DIR = DEFAULTS.get("out_dir", "results")

def ds(name):
    for d in config.get("datasets", []):
        if d.get("name") == name:
            return d
    raise WorkflowError(f"Dataset not found in config: {name}")

def primary_url(name):
    ucfg = ds(name).get("urls", {}) or {}
    if ucfg.get("url"):
        return ucfg.get("url")
    files = ucfg.get("files") or []
    for f in files:
        if f:
            return f
    raise WorkflowError(f"{name}: urls.url or urls.files must be set for gfa_from_url builder")

def path_of(obj):
    """
    Convert Snakemake File/Log wrapper (or list/tuple of such) to a str path.
    """
    if isinstance(obj, (list, tuple)):
        obj = obj[0] if len(obj) > 0 else obj
    try:
        return os.fspath(obj)
    except Exception:
        if hasattr(obj, "path"):
            return str(obj.path)
        return str(obj)

rule gfa_from_url_download:
    message: "Download raw GFA (from URL) for {wildcards.dataset}"
    output:
        raw=os.path.join(DATA_DIR, "{dataset}", "{dataset}.raw.gfa")
    log:
        os.path.join(OUT_DIR, "logs", "gfa", "{dataset}.download.log")
    run:
        dataset = wildcards.dataset
        url = primary_url(dataset)
        if not url:
            raise WorkflowError(f"{dataset}: no URL found for gfa_from_url")

        dest = path_of(output.raw if hasattr(output, "raw") else output[0] if isinstance(output, (list, tuple)) else output)
        logp = path_of(log)

        os.makedirs(os.path.dirname(dest), exist_ok=True)
        os.makedirs(os.path.dirname(logp), exist_ok=True)
        part = dest + ".part"

        def log_text(msg):
            with open(logp, "a", encoding="utf-8") as lf:
                lf.write(msg + "\n")

        # 1) Try to download with urllib (follows redirects)
        log_text(f"[INFO] Starting download: {url} -> {part}")
        downloaded = False
        try:
            req = Request(url, headers={"User-Agent": "snakemake/1.0"})
            with urlopen(req, timeout=120) as r, open(part, "wb") as fh:
                shutil.copyfileobj(r, fh)
            if os.path.exists(part) and os.path.getsize(part) > 0:
                downloaded = True
                log_text(f"[INFO] Python download succeeded, size={os.path.getsize(part)}")
            else:
                log_text("[WARN] Python download produced empty or missing file")
        except Exception as e:
            log_text(f"[WARN] Python download failed: {e}")

        # 2) Fallback to wget/curl if needed
        if not downloaded:
            wget = shutil.which("wget")
            curl = shutil.which("curl")
            if wget:
                log_text(f"[INFO] Trying wget fallback")
                with open(logp, "ab") as logb:
                    rc = subprocess.run([wget, "-q", "-O", part, url], stdout=logb, stderr=logb)
                if rc.returncode == 0 and os.path.exists(part) and os.path.getsize(part) > 0:
                    downloaded = True
                    log_text(f"[INFO] wget succeeded, size={os.path.getsize(part)}")
                else:
                    log_text(f"[WARN] wget failed (rc={rc.returncode})")
            if not downloaded and curl:
                log_text(f"[INFO] Trying curl fallback")
                with open(logp, "ab") as logb:
                    rc = subprocess.run([curl, "-fL", "-o", part, url], stdout=logb, stderr=logb)
                if rc.returncode == 0 and os.path.exists(part) and os.path.getsize(part) > 0:
                    downloaded = True
                    log_text(f"[INFO] curl succeeded, size={os.path.getsize(part)}")
                else:
                    log_text(f"[WARN] curl failed (rc={rc.returncode})")

        if not downloaded:
            raise WorkflowError(f"Failed to download {url} (tried urllib, wget, curl). See {logp}")

        # 3) Decompress/move according to extension
        ubase = url.split("?", 1)[0].lower()
        try:
            if ubase.endswith(".zst"):
                zstd_bin = shutil.which("zstd") or shutil.which("zstdcat")
                if not zstd_bin:
                    raise WorkflowError("zstd not found on PATH; please install zstd or pre-decompress the file.")
                log_text(f"[INFO] Decompressing .zst with {zstd_bin}")
                with open(logp, "ab") as logb:
                    if os.path.basename(zstd_bin).endswith("zstdcat"):
                        with open(dest, "wb") as out:
                            rc = subprocess.run([zstd_bin, part], stdout=out, stderr=logb)
                            if rc.returncode != 0:
                                raise WorkflowError(f"zstdcat failed (rc={rc.returncode}); see {logp}")
                    else:
                        with open(dest, "wb") as out:
                            rc = subprocess.run([zstd_bin, "-d", "-c", part], stdout=out, stderr=logb)
                            if rc.returncode != 0:
                                raise WorkflowError(f"zstd -d failed (rc={rc.returncode}); see {logp}")
            elif ubase.endswith(".gz") or ubase.endswith(".tgz"):
                log_text("[INFO] Decompressing .gz with gzip module")
                try:
                    with gzip.open(part, "rb") as fin, open(dest, "wb") as fout:
                        shutil.copyfileobj(fin, fout)
                except Exception as e:
                    raise WorkflowError(f"gzip decompression failed: {e}")
            else:
                log_text("[INFO] Moving downloaded file to final destination")
                shutil.move(part, dest)
        finally:
            # cleanup .part if still present
            if os.path.exists(part):
                try:
                    os.remove(part)
                except Exception:
                    pass

        if not os.path.exists(dest) or os.path.getsize(dest) == 0:
            raise WorkflowError(f"Produced file {dest} missing or empty; check {logp}")

        log_text(f"[INFO] Produced {dest} (size={os.path.getsize(dest)})")
