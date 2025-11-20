# modules/pggb_from_fasta.smk
#
# Builder: pggb_from_fasta
# Steps:
#   1) Download a .fa.gz from urls.fa_gz
#   2) samtools faidx <fa.gz>
#   3) pggb -t <threads> -i <fa.gz> -o <outdir>
#   4) Copy a .gfa produced by pggb to data/<dataset>/<dataset>.pggb.gfa
#
# Requires:
#   - PGGB_ENV_YML, THREADS_PGGB, ds(), root_dir(), DATA_DIR, OUT_DIR defined in the Snakefile


PGGB_DEFAULTS = DEFAULTS.get("pggb", {}) or {}


def pggb_fasta_url(name):
    d = ds(name)
    urls = d.get("urls") or {}
    url = urls.get("fa_gz")
    if not url:
        raise WorkflowError(
            f"Dataset {name} (builder=pggb_from_fasta) has no urls.fa_gz defined"
        )
    return url


def pggb_fasta_gz_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.pggb.fa.gz")


def pggb_fasta_fai_path(name):
    return pggb_fasta_gz_path(name) + ".fai"


# ---------------------------------------------------------------------------
# Download .fa.gz
# ---------------------------------------------------------------------------
rule pggb_download_fasta:
    message: "Download FASTA (.fa.gz) for pggb dataset {wildcards.dataset}"
    output:
        # data/<dataset>/<dataset>.pggb.fa.gz
        fa_gz=os.path.join(DATA_DIR, "{dataset}", "{dataset}.pggb.fa.gz")
    run:
        import urllib.request, shutil, os as _os

        dataset = wildcards.dataset
        url = pggb_fasta_url(dataset)
        out = output.fa_gz
        _os.makedirs(_os.path.dirname(out), exist_ok=True)
        tmp = out + ".tmp"

        try:
            with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f)
        except Exception as e:
            raise WorkflowError(
                f"Download failed for dataset {dataset} from {url}: {e}"
            )

        _os.replace(tmp, out)


rule pggb_faidx:
    message: "Index FASTA with samtools faidx for {wildcards.dataset}"
    input:
        fa_gz=os.path.join(DATA_DIR, "{dataset}", "{dataset}.pggb.fa.gz")
    output:
        fai=os.path.join(DATA_DIR, "{dataset}", "{dataset}.pggb.fa.gz.fai")
    conda: PGGB_ENV_YML
    log:
        os.path.join(OUT_DIR, "logs", "pggb", "{dataset}.faidx.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        # samtools crée {input.fa_gz}.fai (et éventuellement .gzi)
        samtools faidx {input.fa_gz} >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# pggb -> GFA
# ---------------------------------------------------------------------------
rule pggb_build:
    message: "Build pangenome GFA with pggb for {wildcards.dataset}"
    input:
        fa_gz=os.path.join(DATA_DIR, "{dataset}", "{dataset}.pggb.fa.gz"),
        fai=os.path.join(DATA_DIR, "{dataset}", "{dataset}.pggb.fa.gz.fai")
    output:
        raw_gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.pggb.gfa")
    conda: PGGB_ENV_YML
    threads:
        THREADS_PGGB
    log:
        os.path.join(OUT_DIR, "logs", "pggb", "{dataset}.pggb.log")
    params:
        extra=lambda wc: (
            (ds(wc.dataset).get("pggb") or {}).get(
                "extra_opts", PGGB_DEFAULTS.get("extra_opts", "")
            )
        ),
        outdir=lambda wc: os.path.join(DATA_DIR, wc.dataset, "pggb_out")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname {log})"

        echo "[pggb] Running pggb on {input.fa_gz} -> {params.outdir}" >> {log}
        # ICI on ajoute explicitement -t {threads}
        pggb -t {threads} -i {input.fa_gz} -o "{params.outdir}" {params.extra} >> {log} 2>&1

        # Sélectionner un GFA produit par pggb
        out_gfa=""
        for pat in \
          "{params.outdir}"/*.smooth.final.gfa \
          "{params.outdir}"/*/*.smooth.final.gfa \
          "{params.outdir}"/*.smooth.gfa \
          "{params.outdir}"/*/*.smooth.gfa \
          "{params.outdir}"/*.gfa \
          "{params.outdir}"/*/*.gfa
        do
          if [ -f "$pat" ]; then
            out_gfa="$pat"
            break
          fi
        done

        if [ -z "$out_gfa" ]; then
          echo "[ERROR] No .gfa produced by pggb in {params.outdir}" >> {log}
          exit 1
        fi

        cp "$out_gfa" {output.raw_gfa}
        """