import os, shlex

# 1) Download VCF (.vcf.gz). If not BGZF, recompress locally to BGZF
rule vg_from_vcf_download_vcf:
    message: "vg_from_vcf: download VCF for {wildcards.dataset}"
    output:
        vcf=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vcf.gz")
    params:
        vcf_url=lambda wc: ds(wc.dataset).get("urls", {}).get("vcf_gz")
    conda: VG_ENV_YML
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.vcf_url}" ]; then
            echo "[ERROR] Missing urls.vcf_gz for dataset {wildcards.dataset}" >&2
            exit 2
        fi
        mkdir -p "$(dirname {output.vcf})"
        part="{output.vcf}.part"
        trap 'rm -f "$part"' INT TERM EXIT
        # Download with wget, fallback to curl
        if ! wget -c -O "$part" {params.vcf_url}; then
            rm -f "$part"
            curl -fSL -o "$part" {params.vcf_url}
        fi
        mv -f "$part" {output.vcf}
        trap - INT TERM EXIT

        # Recompress to BGZF if not already BGZF
        if ! bgzip -t {output.vcf} >/dev/null 2>&1; then
            echo "[INFO] VCF is gzip but not BGZF; recompressing with bgzip..." >&2
            tmp="{output.vcf}.tmp.bgz"
            if ! gzip -cd {output.vcf} | bgzip -c > "$tmp"; then
                echo "[ERR] Failed to recompress VCF as BGZF." >&2
                rm -f "$tmp"
                exit 1
            fi
            mv -f "$tmp" {output.vcf}
        fi
        """

# 2) Tabix index
rule vg_from_vcf_tabix:
    message: "vg_from_vcf: tabix index VCF for {wildcards.dataset}"
    input:
        vcf=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vcf.gz")
    output:
        tbi=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vcf.gz.tbi")
    conda: VG_ENV_YML
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tbi})"
        # Ensure BGZF (safety)
        if ! bgzip -t {input.vcf} >/dev/null 2>&1; then
            echo "[INFO] Recompressing VCF as BGZF before tabix..." >&2
            tmp="{input.vcf}.tmp.bgz"
            gzip -cd {input.vcf} | bgzip -c > "$tmp"
            mv -f "$tmp" {input.vcf}
        fi
        tabix -f -p vcf {input.vcf}
        """

# 3) Download reference FASTA (.fa.gz)
rule vg_from_vcf_download_fa:
    message: "vg_from_vcf: download reference FASTA for {wildcards.dataset}"
    output:
        fa_gz=os.path.join(DATA_DIR, "{dataset}", "{dataset}.fa.gz")
    params:
        fa_url=lambda wc: ds(wc.dataset).get("urls", {}).get("fa_gz")
    conda: VG_ENV_YML
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.fa_url}" ]; then
            echo "[ERROR] Missing urls.fa_gz for dataset {wildcards.dataset}" >&2
            exit 2
        fi
        mkdir -p "$(dirname {output.fa_gz})"

        part="{output.fa_gz}.part"
        # Toujours repartir d'un .part propre pour éviter les concaténations hasardeuses
        rm -f "$part"
        trap 'rm -f "$part"' INT TERM EXIT

        # Téléchargement robuste (pas de -c avec -O)
        if ! wget -q -O "$part" {params.fa_url}; then
            rm -f "$part"
            curl -fSL --retry 5 --retry-delay 2 -o "$part" {params.fa_url}
        fi
        mv -f "$part" {output.fa_gz}
        trap - INT TERM EXIT

        # Valider le gzip en tolérant rc=2 (trailing garbage), géré par l'étape suivante
        set +e
        gzip -t {output.fa_gz} >/dev/null 2>&1
        rc=$?
        set -e
        if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
            echo "[ERR] Downloaded FASTA is not a valid gzip (rc=$rc): {output.fa_gz}" >&2
            exit 1
        fi
        """


# 4) Decompress FASTA to .fa (tolerate trailing garbage from gzip)
rule vg_from_vcf_decompress_fa:
    message: "vg_from_vcf: decompress FASTA for {wildcards.dataset}"
    input:
        fa_gz=os.path.join(DATA_DIR, "{dataset}", "{dataset}.fa.gz")
    output:
        fa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.fa")
    conda: VG_ENV_YML
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fa})"
        # gzip may return 2 on 'trailing garbage' even if decompression is OK; accept rc=0 or rc=2
        set +e
        gzip -cd {input.fa_gz} > {output.fa}
        rc=$?
        set -e
        if [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
            echo "[ERR] gzip returned rc=$rc while decompressing {input.fa_gz}" >&2
            exit $rc
        fi
        test -s {output.fa}
        """

# 5) vg construct (.vg)
rule vg_from_vcf_construct:
    message: "vg_from_vcf: construct VG for {wildcards.dataset}"
    input:
        fa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.fa"),
        vcf=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vcf.gz"),
        tbi=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vcf.gz.tbi")
    output:
        vg=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vg")
    conda: VG_ENV_YML
    threads: THREADS_VG
    params:
        region=lambda wc: (ds(wc.dataset).get("vg", {}) or {}).get("region", VG_DEFAULTS.get("region", "1")),
        max_node_length=lambda wc: int((ds(wc.dataset).get("vg", {}) or {}).get("max_node_length", VG_DEFAULTS.get("max_node_length", 3000000)))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.vg})"
        vg construct -r {input.fa} -v {input.vcf} -R {params.region} -m {params.max_node_length} -t {threads} > {output.vg}
        """

# 6) Convert to GFA (raw GFA consumed by main Snakefile)
rule vg_from_vcf_raw_gfa:
    message: "vg_from_vcf: convert VG -> GFA for {wildcards.dataset}"
    input:
        vg=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vg")
    output:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.vg.gfa")
    conda: VG_ENV_YML
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gfa})"
        vg convert -f {input.vg} > {output.gfa}
        """