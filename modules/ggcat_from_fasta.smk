import os, shlex, shutil, re
from urllib.request import urlopen
from urllib.error import URLError, HTTPError
from urllib.parse import urlparse, unquote
from snakemake.shell import shell
from snakemake.exceptions import WorkflowError
from pathlib import Path

shell.executable("/bin/bash")

BASEDIR = Path(workflow.basedir) if "workflow" in globals() else Path(os.getcwd())
def resolve_env(p):
    p = str(p)
    return p if os.path.isabs(p) else str(BASEDIR / p)

DEFAULTS = config.get("defaults", {})
DATA_DIR = DEFAULTS.get("data_dir", "data")
OUT_DIR = DEFAULTS.get("out_dir", "results")
ENV_GG = resolve_env(DEFAULTS.get("envs", {}).get("ggcat", "config/ggcat.yml"))
THREADS_GGCAT = int(DEFAULTS.get("threads", {}).get("ggcat", 8))
GGCAT_COMMON = DEFAULTS.get("ggcat", {}).get("common_opts", "-e -s 1")

def ds(name):
    for d in config.get("datasets", []):
        if d.get("name") == name:
            return d
    raise WorkflowError(f"Dataset not found in config: {name}")

def fna_url(dataset):
    u = (ds(dataset).get("urls", {}) or {}).get("fna_gz")
    return u

def fasta_dir(dataset):
    return os.path.join(DATA_DIR, dataset, "fasta")

def discover_fna_urls(dataset):
    """
    Discover source URLs for a ggcat_from_fasta dataset.

    Order of precedence:
      1) explicit urls.files (list of URLs) -> return sorted unique list
      2) a tarball URL urls.tar_gz (or urls.tarball) -> return [tar_url]
      3) listings: parse HTML listing(s) and select hrefs matching pattern
    """
    ucfg = ds(dataset).get("urls", {}) or {}
    # 1) explicit files, if provided
    files = [u for u in (ucfg.get("files") or []) if u]
    if files:
        return sorted(set(files))

    # 1b) support a tarball containing many FASTA files
    tar = ucfg.get("tar_gz") or ucfg.get("tarball")
    if tar:
        return [tar]

    # 2) listings (one or many)
    listings = ucfg.get("listings")
    if isinstance(listings, str) and listings.strip():
        listings = [listings]
    listings = [l for l in (listings or []) if l]

    pattern = ucfg.get("pattern", r".*\.fna(\.gz)?$")
    origin = ucfg.get("origin", "")
    limit = int(ucfg.get("limit", 0))
    timeout = int(ucfg.get("discovery_timeout", 20))

    if not listings:
        # No listings and no explicit files/tar -> no discovery in this mode
        return []

    rx = re.compile(pattern)
    urls = []
    for listing in listings:
        try:
            with urlopen(listing, timeout=timeout) as r:
                html = r.read().decode("utf-8", errors="replace")
            hrefs = re.findall(r'href="([^"]+)"', html)
        except (URLError, HTTPError, Exception) as e:
            raise WorkflowError(f"{dataset}: failed to read listing {listing}: {e}") from e

        for h in hrefs:
            if not rx.search(h or ""):
                continue
            if re.match(r"^https?://", h):
                u = h
            elif h.startswith("/"):
                u = (origin.rstrip("/") + h) if origin else listing.rstrip("/") + h
            else:
                u = listing.rstrip("/") + "/" + h
            urls.append(u)

    urls = sorted(set(urls))
    if limit and limit > 0:
        urls = urls[:limit]
    if not urls:
        raise WorkflowError(f"{dataset}: discovery returned no URLs matching {pattern}.")
    return urls

# MODE 1: multi-files (listings/files/tarball) → download → concat/decompress → combined.fa
rule ggcat_from_fasta_download_fna_multi:
    message: "ggcat_from_fasta: discover+download all FNA for {wildcards.dataset}"
    output:
        marker=os.path.join(DATA_DIR, "{dataset}", "fasta", ".downloaded.ok"),
        urlfile=os.path.join(DATA_DIR, "{dataset}", "fasta", "urls.txt")
    log:
        dl=os.path.join(OUT_DIR, "logs", "fasta", "{dataset}.download.log")
    run:
        dataset = wildcards.dataset
        urls = discover_fna_urls(dataset)
        if not urls:
            # No multi-mode URLs found -> let other branch apply (single fna_gz or tar handled elsewhere)
            # If the user explicitly configured listings but discovery returns none, this already raises upstream.
            raise WorkflowError(f"{dataset}: no .fna(.gz) discovered; check urls.listings/files/pattern.")
        os.makedirs(os.path.dirname(output.urlfile), exist_ok=True)
        os.makedirs(os.path.dirname(log.dl), exist_ok=True)
        os.makedirs(fasta_dir(dataset), exist_ok=True)

        with open(output.urlfile, "w", encoding="utf-8") as f:
            for u in urls:
                f.write(u.strip() + "\n")

        print(f"[{dataset}] Downloading {len(urls)} FNA files...")
        for u in urls:
            # Derive a safe basename from the URL path (drop query part if any)
            parsed = urlparse(u)
            base = os.path.basename(parsed.path) or os.path.basename(unquote(u)) or os.path.basename(u)
            base = unquote(base)
            dest = os.path.join(fasta_dir(dataset), base)
            shell(f"""
                set -euo pipefail
                echo "[INFO] Downloading: {shlex.quote(base)}" >> {shlex.quote(log.dl)}
                part={shlex.quote(dest)}.part
                trap 'rm -f "$part"' INT TERM EXIT
                ( wget -q -c -O "$part" {shlex.quote(u)} >> {shlex.quote(log.dl)} 2>&1 \
                  || curl -fsSL --retry 5 --retry-delay 2 -C - -o "$part" {shlex.quote(u)} >> {shlex.quote(log.dl)} 2>&1 ) \
                  || ( echo "[ERROR] Download failed: {u}; see {shlex.quote(log.dl)}" >&2; exit 1 )
                mv -f "$part" {shlex.quote(dest)}
                trap - INT TERM EXIT
            """)
        open(output.marker, "w").close()

rule ggcat_from_fasta_concat_combined:
    message: "ggcat_from_fasta: concat/decompress all FNA for {wildcards.dataset}"
    input:
        marker=rules.ggcat_from_fasta_download_fna_multi.output.marker,
        urlfile=rules.ggcat_from_fasta_download_fna_multi.output.urlfile
    output:
        fa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.combined.fa")
    log:
        concat=os.path.join(OUT_DIR, "logs", "fasta", "{dataset}.concat.log")
    run:
        dataset = wildcards.dataset
        os.makedirs(os.path.dirname(output.fa), exist_ok=True)
        os.makedirs(os.path.dirname(log.concat), exist_ok=True)

        with open(input.urlfile, "r", encoding="utf-8") as f:
            urls = [l.strip() for l in f if l.strip()]

        # Map URLs -> local downloaded filenames (use same basename derivation as downloader)
        files = []
        for u in urls:
            parsed = urlparse(u)
            base = os.path.basename(parsed.path) or os.path.basename(unquote(u)) or os.path.basename(u)
            base = unquote(base)
            files.append(os.path.join(fasta_dir(dataset), base))
        files = sorted(files, key=lambda p: os.path.basename(p))

        # If any tarball were downloaded, extract them and collect contained FASTA files
        tar_files = [fp for fp in files if re.search(r'\.tar(\.(gz|bz2|xz))?$', fp, re.I)]
        if tar_files:
            import tarfile
            extract_dir = os.path.join(fasta_dir(dataset), ".tar_extract")
            os.makedirs(extract_dir, exist_ok=True)
            for tf in tar_files:
                if not os.path.exists(tf):
                    raise WorkflowError(f"{dataset}: expected tarball {tf} not found (download failed?)")
                if not tarfile.is_tarfile(tf):
                    raise WorkflowError(f"{dataset}: file {tf} does not look like a tar archive")
                with tarfile.open(tf, 'r:*') as tarf:
                    # safety checks to prevent path traversal
                    for member in tarf.getmembers():
                        mname = member.name
                        if os.path.isabs(mname):
                            raise WorkflowError(f"{dataset}: unsafe member in tar (absolute path): {mname}")
                        # disallow '..' in path components
                        if '..' in Path(mname).parts:
                            raise WorkflowError(f"{dataset}: unsafe member in tar (parent-ref): {mname}")
                    tarf.extractall(extract_dir)

            # discover extracted fasta files
            extra = []
            for root, _, filenames in os.walk(extract_dir):
                for fn in filenames:
                    if re.search(r'\.(fna|fa|fasta)(\.gz)?$', fn, re.I):
                        extra.append(os.path.join(root, fn))
            if not extra:
                raise WorkflowError(f"{dataset}: tarball did not contain any fasta files")
            # remove tar files from list and append extracted FASTA files
            files = [f for f in files if f not in tar_files] + sorted(extra, key=lambda p: os.path.basename(p))
            files = sorted(files, key=lambda p: os.path.basename(p))

        parts = []
        for fp in files:
            if fp.endswith(".gz"):
                parts.append(f"gunzip -c {shlex.quote(fp)}")
            else:
                parts.append(f"cat {shlex.quote(fp)}")

        print(f"[{dataset}] Concatenating and decompressing {len(files)} FNA files...")
        cmd = (
            "set -euo pipefail; "
            + f'echo "[INFO] Concatenating {len(files)} files into {shlex.quote(output.fa)}" >> {shlex.quote(log.concat)}; '
            + "mkdir -p " + shlex.quote(os.path.dirname(output.fa)) + "; "
            + "( " + " && ".join(parts) + " )"
            + " > " + shlex.quote(output.fa)
            + " 2>> " + shlex.quote(log.concat)
            + " || ( echo '[ERROR] Concatenation/decompression failed; see "
            + shlex.quote(log.concat)
            + "' >&2; exit 1 ); "
            + "test -s " + shlex.quote(output.fa)
        )
        shell(cmd)

# MODE 2: a single .fna.gz explicit -> download .fna.gz
rule ggcat_from_fasta_download_fna:
    message: "ggcat_from_fasta: download FNA for {wildcards.dataset}"
    output:
        fna_gz=os.path.join(DATA_DIR, "{dataset}", "{dataset}.fna.gz")
    params:
        url=lambda wc: fna_url(wc.dataset)
    run:
        if not params.url:
            # No explicit fna_gz -> let other branch apply
            raise WorkflowError(f"{wildcards.dataset}: urls.fna_gz is required (or provide urls.listings/files for multi).")
        shell(r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.fna_gz})"
            part="{output.fna_gz}.part"
            trap 'rm -f "$part"' INT TERM EXIT
            if ! wget -c -O "$part" {params.url}; then
                rm -f "$part"
                curl -fSL -o "$part" {params.url}
            fi
            mv -f "$part" {output.fna_gz}
            trap - INT TERM EXIT
        """)

rule ggcat_from_fasta_decompress:
    message: "ggcat_from_fasta: decompress FNA for {wildcards.dataset}"
    input:
        fna_gz=rules.ggcat_from_fasta_download_fna.output.fna_gz
    output:
        fa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.combined.fa")
    run:
        shell(r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.fa})"
            sig=$(head -c2 {input.fna_gz} | od -An -tx1 | tr -d ' \n' || true)
            if [ "$sig" = "1f8b" ]; then
                if gzip -dc {input.fna_gz} > {output.fa}; then
                    :
                else
                    rc=$?
                    if [ $rc -eq 2 ]; then
                        gzip -dc {input.fna_gz} > {output.fa} || true
                    else
                        echo "[ERR] gzip failed (rc=$rc)" >&2
                        exit $rc
                    fi
                fi
            else
                cp -f {input.fna_gz} {output.fa}
            fi
            test -s {output.fa}
        """)

# ggcat build (consumes {dataset}.combined.fa regardless of branch)
rule ggcat_from_fasta_raw_gfa:
    message: "ggcat_from_fasta: ggcat build for {wildcards.dataset}"
    input:
        fa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.combined.fa")
    output:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.ggcat.fasta.gfa")
    conda: ENV_GG
    threads: THREADS_GGCAT
    log:
        gg=os.path.join(OUT_DIR, "logs", "ggcat", "{dataset}.log")
    params:
        k=lambda wc: int(ds(wc.dataset).get("ggcat", {}).get("k", 31)),
        common=GGCAT_COMMON
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gfa})" "$(dirname {log.gg})"
        ggcat build --gfa-v1 {params.common} -k {params.k} -j {threads} {input.fa} -o {output.gfa} -t {output.gfa}.temp 2>&1 | tee -a {log.gg}
        """

# Prefer the multi/concat branch over the single-file decompress branch
ruleorder: ggcat_from_fasta_concat_combined > ggcat_from_fasta_decompress