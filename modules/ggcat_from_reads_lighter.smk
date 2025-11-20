import os, re, shlex, shutil
from urllib.request import urlopen
from urllib.error import URLError, HTTPError
from snakemake.shell import shell
from snakemake.exceptions import WorkflowError
from pathlib import Path

BASEDIR = Path(workflow.basedir) if "workflow" in globals() else Path(os.getcwd())
def resolve_env(p):
    p = str(p)
    return p if os.path.isabs(p) else str(BASEDIR / p)

shell.executable("/bin/bash")

DEFAULTS = config.get("defaults", {})
DATA_DIR = DEFAULTS.get("data_dir", "data")
OUT_DIR = DEFAULTS.get("out_dir", "results")
ENV_GG = resolve_env(DEFAULTS.get("envs", {}).get("ggcat", "config/ggcat.yml"))
THREADS_GGCAT = int(DEFAULTS.get("threads", {}).get("ggcat", 8))
THREADS_LIGHTER_DEF = int(DEFAULTS.get("threads", {}).get("lighter", 8))
GGCAT_COMMON = DEFAULTS.get("ggcat", {}).get("common_opts", "-e -s 1")
LIGHTER_BIN_DEFAULT = DEFAULTS.get("tools", {}).get("lighter_bin", "lighter")

def ds(name):
    for d in config.get("datasets", []):
        if d.get("name") == name:
            return d
    raise WorkflowError(f"Dataset not found in config: {name}")

def reads_dir(dataset):
    return os.path.join(DATA_DIR, dataset, "reads")

def discover_urls(dataset):
    ucfg = ds(dataset).get("urls", {}) or {}
    files = [u for u in (ucfg.get("files") or []) if u]
    if files:
        return files

    listing = ucfg.get("listing")
    pattern = ucfg.get("pattern", r".*\.fastq(\.gz)?$")
    paired_only = bool(ucfg.get("paired_only", False))
    limit = int(ucfg.get("limit", 0))
    origin = ucfg.get("origin", "")
    timeout = int(ucfg.get("discovery_timeout", 15))

    if not listing:
        raise WorkflowError(f"{dataset}: urls.files is empty and urls.listing is not set")

    try:
        with urlopen(listing, timeout=timeout) as r:
            html = r.read().decode("utf-8", errors="replace")
        hrefs = re.findall(r'href="([^"]+)"', html)
    except (URLError, HTTPError, Exception) as e:
        raise WorkflowError(f"{dataset}: failed to read listing {listing}: {e}") from e

    rx = re.compile(pattern)
    urls = []
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

    if paired_only:
        urls = [u for u in urls if re.search(r'_(?:1|2)\.filt\.fastq\.gz$', u)]

    def sort_key(u):
        b = os.path.basename(u)
        m = re.match(r'(SRR\d+)(?:_(\d))?\.filt\.fastq\.gz$', b)
        srr = m.group(1) if m else b
        tag = int(m.group(2)) if (m and m.group(2)) else 3
        return (srr, tag, b)
    urls = sorted(set(urls), key=sort_key)

    if limit and limit > 0:
        urls = urls[:limit]

    if not urls:
        raise WorkflowError(f"{dataset}: discovery returned no URLs. Provide urls.files or adjust discovery settings.")
    return urls

rule ggcat_lighter_download_reads:
    message: "ggcat_from_reads_lighter: download reads for {wildcards.dataset}"
    output:
        marker=os.path.join(DATA_DIR, "{dataset}", "reads", ".downloaded.ok"),
        urlfile=os.path.join(DATA_DIR, "{dataset}", "reads", "urls.txt")
    log:
        dl=os.path.join(OUT_DIR, "logs", "reads", "{dataset}.download.log")
    run:
        dataset = wildcards.dataset
        urls = discover_urls(dataset)
        rdir = reads_dir(dataset)
        os.makedirs(rdir, exist_ok=True)
        os.makedirs(os.path.dirname(log.dl), exist_ok=True)

        # Minimal console feedback
        print(f"[{dataset}] Downloading {len(urls)} read files...")

        with open(output.urlfile, "w", encoding="utf-8") as f:
            for u in urls:
                f.write(u.strip() + "\n")

        for u in urls:
            dest = os.path.join(rdir, os.path.basename(u))
            shell(f"""
                set -euo pipefail
                echo "[INFO] Downloading: {os.path.basename(u)}" >> {shlex.quote(log.dl)}
                part={shlex.quote(dest)}.part
                trap 'rm -f "$part"' INT TERM EXIT
                # Quiet download; try wget first, fallback to curl; log everything.
                ( wget -q -c -O "$part" {shlex.quote(u)} >> {shlex.quote(log.dl)} 2>&1 \
                  || curl -fsSL --retry 5 --retry-delay 2 -C - -o "$part" {shlex.quote(u)} >> {shlex.quote(log.dl)} 2>&1 ) \
                  || ( echo "[ERROR] Download failed: {u}; see {log.dl}" >&2; exit 1 )
                mv -f "$part" {shlex.quote(dest)}
                trap - INT TERM EXIT
            """)
        # Mark completion
        open(output.marker, "w").close()

rule ggcat_lighter_concat_combined:
    message: "ggcat_from_reads_lighter: concat/decompress reads for {wildcards.dataset}"
    input:
        marker=rules.ggcat_lighter_download_reads.output.marker,
        urlfile=rules.ggcat_lighter_download_reads.output.urlfile
    output:
        combined=os.path.join(DATA_DIR, "{dataset}", "{dataset}.combined.fastq")
    log:
        concat=os.path.join(OUT_DIR, "logs", "reads", "{dataset}.concat.log")
    run:
        import re, shlex
        dataset = wildcards.dataset
        rdir = reads_dir(dataset)
        os.makedirs(os.path.dirname(log.concat), exist_ok=True)

        with open(input.urlfile, "r", encoding="utf-8") as f:
            urls = [l.strip() for l in f if l.strip()]
        files = [os.path.join(rdir, os.path.basename(u)) for u in urls]

        def sort_key(p):
            b = os.path.basename(p)
            m = re.match(r'(SRR\d+)(?:_(\d))?\.filt\.fastq\.gz$', b)
            srr = m.group(1) if m else b
            tag = int(m.group(2)) if (m and m.group(2)) else 3
            return (srr, tag, b)
        files = sorted(files, key=sort_key)

        parts = []
        for fpath in files:
            if fpath.endswith(".gz"):
                parts.append(f"gunzip -c {shlex.quote(fpath)}")
            else:
                parts.append(f"cat {shlex.quote(fpath)}")

        # Minimal console feedback
        print(f"[{dataset}] Concatenating and decompressing {len(files)} files...")

        cmd = (
            "set -euo pipefail; "
            + f'echo "[INFO] Concatenating {len(files)} files into {shlex.quote(output.combined)}" >> {shlex.quote(log.concat)}; '
            + "mkdir -p " + shlex.quote(os.path.dirname(output.combined)) + "; "
            + "( " + " && ".join(parts) + " )"
            + " > " + shlex.quote(output.combined)
            + " 2>> " + shlex.quote(log.concat)
            + " || ( echo '[ERROR] Concatenation/decompression failed; see "
            + shlex.quote(log.concat)
            + "' >&2; exit 1 ); "
            + "test -s " + shlex.quote(output.combined)
        )
        shell(cmd)

rule ggcat_lighter_correct:
    message: "ggcat_from_reads_lighter: Lighter correction for {wildcards.dataset}"
    input:
        combined=rules.ggcat_lighter_concat_combined.output.combined
    output:
        corrected=os.path.join(DATA_DIR, "{dataset}", "{dataset}.corrected.fastq")
    threads:
        lambda wc: int((ds(wc.dataset).get("lighter", {}) or {}).get("threads", THREADS_LIGHTER_DEF))
    log:
        lighter=os.path.join(OUT_DIR, "logs", "lighter", "{dataset}.log")
    params:
        exe=lambda wc: (ds(wc.dataset).get("lighter", {}) or {}).get("bin", LIGHTER_BIN_DEFAULT),
        K=lambda wc: int((ds(wc.dataset).get("lighter", {}) or {}).get("K", 31)),
        g=lambda wc: int((ds(wc.dataset).get("lighter", {}) or {}).get("genome_size", 3000000000))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.corrected})" "$(dirname {log.lighter})"
        OUTDIR="$(dirname {output.corrected})"
        B="$(basename {input.combined})"

        echo "[lighter] Command: {params.exe} -K {params.K} {params.g} -r {input.combined} -od \"$OUTDIR\" -t {threads}" >> {log.lighter}
        printf "[lighter] Starting at %s\n" "$(date -Is)" >> {log.lighter}
        {params.exe} -K {params.K} {params.g} -r {input.combined} -od "$OUTDIR" -t {threads} >> {log.lighter} 2>&1 \
          || ( echo "[ERROR] Lighter failed; see {log.lighter}" >&2; exit 1 )
        printf "[lighter] Finished at %s\n" "$(date -Is)" >> {log.lighter}

        cf=""
        base="${{B%.*}}"
        for p in \
            "$OUTDIR/$base.cor.fq" "$OUTDIR/$base.cor.fastq" \
            "$OUTDIR/$base.cor.fq.gz" "$OUTDIR/$base.cor.fastq.gz" \
            "$OUTDIR/${{base}}_1.cor.fq" "$OUTDIR/${{base}}_1.cor.fastq" \
            "$OUTDIR/$base.lighter.cor.fq" "$OUTDIR/$base.lighter.cor.fastq" \
            "$OUTDIR/$base.lighter.cor.fq.gz" "$OUTDIR/$base.lighter.cor.fastq.gz" \
            ; do
            [ -f "$p" ] && cf="$p" && break
        done
        if [ -z "$cf" ]; then
            cf=$(ls -1 "$OUTDIR"/*cor*.f*q* 2>/dev/null | head -n1 || true)
        fi
        if [ -z "$cf" ]; then
            echo "[ERROR] Lighter did not produce a corrected file under $OUTDIR" >> {log.lighter}
            echo "[ERROR] Lighter did not produce a corrected file under $OUTDIR; see {log.lighter}" >&2
            exit 1
        fi

        case "$cf" in
            *.gz) gunzip -c "$cf" > {output.corrected} 2>> {log.lighter} ;;
            *)    cp -f "$cf" {output.corrected} 2>> {log.lighter} ;;
        esac
        test -s {output.corrected}
        """

rule ggcat_lighter_raw_gfa:
    message: "ggcat_from_reads_lighter: ggcat build for {wildcards.dataset}"
    input:
        corrected=rules.ggcat_lighter_correct.output.corrected
    output:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.ggcat.lighter.gfa")
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
        (
          echo "[ggcat] Command: ggcat build --gfa-v1 {params.common} -k {params.k} -j {threads} {input.corrected} -o {output.gfa} -t {output.gfa}.temp"
          ggcat build --gfa-v1 {params.common} -k {params.k} -j {threads} {input.corrected} -o {output.gfa} -t {output.gfa}.temp
        ) >> {log.gg} 2>&1 || ( echo "[ERROR] ggcat build failed; see {log.gg}" >&2; exit 1 )
        """