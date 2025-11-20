configfile: "datasets.yaml"

import os, re, shlex, shutil, subprocess, threading, sys, json, hashlib
from urllib.request import urlopen
from urllib.error import URLError, HTTPError
from snakemake.shell import shell
from snakemake.exceptions import WorkflowError
from pathlib import Path

shell.executable("/bin/bash")

BASEDIR = Path(workflow.basedir) if "workflow" in globals() else Path(os.getcwd())
def resolve_env(p):
    p = str(p)
    return p if os.path.isabs(p) else str(BASEDIR / p)

BIN_DIR = BASEDIR / "bin"


SPQR_BIN_DEFAULT      = str(BIN_DIR / "BubbleFinder")   
GET_BLUNTED_DEFAULT   = str(BIN_DIR / "get_blunted")
CLSD_BIN_DEFAULT      = str(BIN_DIR / "clsd")
LIGHTER_BIN_DEFAULT   = str(BIN_DIR / "lighter")

DEFAULTS = config.get("defaults", {}) or {}

DATA_DIR = DEFAULTS.get("data_dir", "data")
OUT_DIR = DEFAULTS.get("out_dir", "results")

THREADS_GGCAT = int(DEFAULTS.get("threads", {}).get("ggcat", 8))
THREADS_VG = int(DEFAULTS.get("threads", {}).get("vg", 16))
THREADS_PGGB = int(DEFAULTS.get("threads", {}).get("pggb", 8))

GG_ENV_YML        = resolve_env(DEFAULTS.get("envs", {}).get("ggcat", "config/ggcat.yml"))
BUBBLEGUN_ENV_YML = resolve_env(DEFAULTS.get("envs", {}).get("bubblegun", "config/bubblegun.yml"))
PLOT_ENV_YML      = resolve_env(DEFAULTS.get("envs", {}).get("plot", "config/plot.yml"))
VG_ENV_YML        = resolve_env(DEFAULTS.get("envs", {}).get("vg", "config/vg.yml"))
PYTOOLS_ENV_YML   = resolve_env(DEFAULTS.get("envs", {}).get("pytools", "config/pytools.yml"))
PGGB_ENV_YML      = resolve_env(DEFAULTS.get("envs", {}).get("pggb", "config/pggb.yaml"))
BUILD_TOOLS_ENV_YML = resolve_env(DEFAULTS.get("envs", {}).get("build_tools", "config/build_tools.yml"))

TOOLS = DEFAULTS.get("tools", {}) or {}
TIME_BIN = TOOLS.get("time_bin", "/usr/bin/time")
TIMEOUT_CONF = TOOLS.get("timeout", {}) or {}
TIMEOUT_ENABLED = bool(TIMEOUT_CONF.get("enabled", True))
TIMEOUT_SECS = int(TIMEOUT_CONF.get("seconds", 100))
TIMEOUT_KILL_AFTER = int(TIMEOUT_CONF.get("kill_after", 10))
TIMEOUT_BIN = TIMEOUT_CONF.get("bin", shutil.which("timeout") or "timeout")

SPQR_BIN    = TOOLS.get("spqr_bin",    SPQR_BIN_DEFAULT)
SPQR_GFA_FLAG = TOOLS.get("spqr_gfa_flag", "--gfa")
SPQR_SGRAPH_FLAG = TOOLS.get("spqr_sgraph_flag", "")
SPQR_THREADS_FLAG = TOOLS.get("spqr_threads_flag", "-j")
SPQR_REPORT_FLAG = TOOLS.get("spqr_report_flag", "--report-json")
GET_BLUNTED = TOOLS.get("get_blunted", GET_BLUNTED_DEFAULT)
CLSD_BIN    = TOOLS.get("clsd_bin",    CLSD_BIN_DEFAULT)
LIGHTER_BIN = TOOLS.get("lighter_bin", LIGHTER_BIN_DEFAULT)

GGCAT_COMMON = DEFAULTS.get("ggcat", {}).get("common_opts", "-e -s 1")

SB_CONF = DEFAULTS.get("sb", {}) or {}

BENCH = DEFAULTS.get("bench", {}) or {}
REPS = int(BENCH.get("reps", 2))
DEFAULT_PROGRAMS = list(BENCH.get("programs", ["BubbleGun_gfa", "sbSPQR_gfa", "vg_snarls_gfa"]))
SINK_OUTPUT = bool(BENCH.get("sink_output", False))
DEBUG = bool(BENCH.get("debug", False))
PROGRAM_CLASSES = (BENCH.get("program_classes", {}) or {})
THREADS_BY_PROGRAM_GLOBAL = (BENCH.get("threads_by_program", {}) or {})

VG_DEFAULTS = DEFAULTS.get("vg", {}) or {}

SG_CONF = DEFAULTS.get("sgraph", {}) or {}
SG_BUILD_IN_ALL = bool(SG_CONF.get("build_in_all", False))
SG_GZIP = bool(SG_CONF.get("gzip", True))

SBFIND_GLOBAL = (DEFAULTS.get("sbfind", {}) or {})

workflow.global_resources["benchmark"] = 200

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def all_datasets():
    return list(config.get("datasets", []))

def parse_enabled_flag(val):
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    if s in {"1", "true", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "no", "off", "n"}:
        return False
    return True

def dataset_enabled_auto(d):
    b = d.get("builder", "")
    urls = d.get("urls", {}) or {}
    if b == "ggcat_from_reads_lighter":
        files = urls.get("files") or []
        if len(files) > 0:
            return True
        listing = urls.get("listing")
        pattern = urls.get("pattern", r".*\.fastq(\.gz)?$")
        timeout = int(urls.get("discovery_timeout", 15))
        if not listing:
            return False
        try:
            with urlopen(listing, timeout=timeout) as r:
                html = r.read().decode("utf-8", errors="replace")
            hrefs = re.findall(r'href="([^"]+)"', html)
            rx = re.compile(pattern)
            matches = [h for h in hrefs if rx.search(h or "")]
            return len(matches) > 0
        except Exception:
            return False
    if b == "ggcat_from_fasta":
        return bool(urls.get("fna_gz"))
    if b == "vg_from_vcf":
        return bool(urls.get("fa_gz")) and bool(urls.get("vcf_gz"))
    return True

def is_dataset_enabled(d):
    flag = parse_enabled_flag(d.get("enabled", True))
    if flag == "auto":
        return dataset_enabled_auto(d)
    return bool(flag)

DS = {d["name"]: d for d in all_datasets() if is_dataset_enabled(d)}

def ds(name):
    if name not in DS:
        raise WorkflowError(f"Unknown or disabled dataset: {name}")
    return DS[name]

def builder(name):
    return ds(name).get("builder", "")

def root_dir(name):
    return os.path.join(DATA_DIR, name)

def raw_gfa_path(name):
    b = builder(name)
    if b == "ggcat_from_fasta":
        ext = "ggcat.fasta"
    elif b == "ggcat_from_reads_lighter":
        ext = "ggcat.lighter"
    elif b == "pggb_from_fasta":
        ext = "pggb"
    elif b == "vg_from_vcf":
        ext = "vg"
    else:
        ext = "raw"
    return os.path.join(root_dir(name), f"{name}.{ext}.gfa")

def blunt_gfa_path(name):
    return re.sub(r"\.gfa$", ".bluntified.gfa", raw_gfa_path(name))

def clean_gfa_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.cleaned.gfa")

def bench_programs_for_dataset(name):
    d = ds(name)
    progs = d.get("bench_programs")
    return list(progs) if progs else DEFAULT_PROGRAMS

def program_class(p):
    if p in PROGRAM_CLASSES:
        return PROGRAM_CLASSES[p]
    if p in {"clsd_sb", "sbSPQR_sb"}:
        return "SB"
    if p in {"BubbleGun_gfa", "sbSPQR_gfa", "vg_snarls_gfa", "sbSPQR_snarls_gfa"}:
        return "BiSB"
    return "NA"

def all_enabled_dataset_names():
    return sorted(list(DS.keys()))

def plot_paths():
    return [
        os.path.join(OUT_DIR, "plots", "time_by_dataset_program.png"),
        os.path.join(OUT_DIR, "plots", "rss_by_dataset_program.png"),
    ]

def sgraph_out_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.sgraph{'.gz' if SG_GZIP else ''}")

def sgraph_enabled_for_dataset(name):
    d = ds(name)
    enabled_global = bool(SG_CONF.get("enabled", True))
    local = d.get("sgraph", {}) or {}
    if "enabled" in local:
        return bool(local.get("enabled"))
    return enabled_global

def sgraph_targets():
    if not SG_BUILD_IN_ALL:
        return []
    return [sgraph_out_path(n) for n in all_enabled_dataset_names() if sgraph_enabled_for_dataset(n)]

def is_program_sb(p):
    return p in {"clsd_sb", "sbSPQR_sb"}

def is_program_sbspqr(p):
    return p in {"sbSPQR_gfa", "sbSPQR_sb", "sbSPQR_snarls_gfa"}

def dataset_requests_sb(name):
    progs = bench_programs_for_dataset(name)
    return any(is_program_sb(p) for p in progs)

def sb_force_f_for_dataset(name):
    if builder(name).startswith("ggcat_") and dataset_requests_sb(name):
        return True
    d = ds(name)
    local = (d.get("sb", {}) or {}).get("ggcat_force_f", None)
    if local is not None:
        return bool(local) and builder(name).startswith("ggcat_")
    return bool(SB_CONF.get("ggcat_force_f", False)) and builder(name).startswith("ggcat_")

def sb_ggcat_f_value(name):
    if builder(name).startswith("ggcat_"):
        return "1" if sb_force_f_for_dataset(name) else "0"
    return "NA"

def sb_raw_gfa_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.ggcat.sb.gfa")

def sb_clean_gfa_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.sb.cleaned.gfa")

def clsd_edgelist_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.clsd.edgelist")

def sbspqr_sgraph_path(name):
    return os.path.join(DATA_DIR, name, f"{name}.sbspqr.sgraph")

def sgraph_flags_for_dataset(name):
    local = (ds(name).get("sgraph", {}) or {})
    def flag(key, default_bool, cli):
        val = local.get(key, SG_CONF.get(key, default_bool))
        return cli if bool(val) else ""
    return " ".join([f for f in [
        flag("edges_only", True, "--edges-only"),
        flag("one_based", False, "--one-based"),
        flag("add_rc_edges", False, "--add-rc-edges"),
        flag("dedup_edges", True, "--dedup-edges"),
    ] if f])

def gfa_for_sb_inputs(name):
    return sb_clean_gfa_path(name) if sb_force_f_for_dataset(name) else clean_gfa_path(name)

def all_selected_programs():
    s = set()
    for n in all_enabled_dataset_names():
        s.update(bench_programs_for_dataset(n))
    return s

def threads_for_prog_dataset(prog, dataset):
    d = ds(dataset)
    dt = (d.get("bench_threads", {}) or {}).get(prog)
    if dt:
        return list(dt)
    if is_program_sbspqr(prog):
        ds_sb = (d.get("sbfind", {}) or {}).get("threads")
        if ds_sb:
            return list(ds_sb)
    gp = THREADS_BY_PROGRAM_GLOBAL.get(prog)
    if gp:
        return list(gp)
    if is_program_sbspqr(prog):
        gsb = (SBFIND_GLOBAL.get("threads") or [])
        if gsb:
            return list(gsb)
    return [1]

def all_threads_values():
    vals = set()
    for dname in all_enabled_dataset_names():
        for p in bench_programs_for_dataset(dname):
            for t in threads_for_prog_dataset(p, dname):
                try:
                    vals.add(int(t))
                except Exception:
                    pass
    return sorted(vals or [1])

# -----------------------------------------------------------------------------
# Includes
# -----------------------------------------------------------------------------
include: "modules/ggcat_from_fasta.smk"
include: "modules/ggcat_from_reads_lighter.smk"
include: "modules/vg_from_vcf.smk"
include: "modules/gfa_to_sgraph.smk"
include: "modules/gfa_from_url.smk"
include: "modules/pggb_from_fasta.smk"


if SPQR_BIN == SPQR_BIN_DEFAULT:
    rule build_bubblefinder:
        message: "Clone & build BubbleFinder from Git"
        output:
            SPQR_BIN_DEFAULT
        conda:
            BUILD_TOOLS_ENV_YML
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output})"

            # S'assurer que le dossier build/ existe
            mkdir -p build

            # Cloner explicitement dans build/BubbleFinder si absent
            if [ ! -d build/BubbleFinder ]; then
              git clone https://github.com/algbio/BubbleFinder.git build/BubbleFinder
            fi

            cd build/BubbleFinder
 
            git checkout 244f5ad1a9b258da454eeb5796c1d2e7985cb9aa

            mkdir -p build
            cd build
            cmake .. -DCMAKE_BUILD_TYPE=Release
            make -j"$(nproc)"

            cp BubbleFinder "{output}" 
            """
            
if GET_BLUNTED == GET_BLUNTED_DEFAULT:
    rule build_get_blunted:
        message: "Download precompiled GetBlunted binary"
        output:
            GET_BLUNTED_DEFAULT
        conda:
            BUILD_TOOLS_ENV_YML
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output})"

            curl -L \
              "https://github.com/vgteam/GetBlunted/releases/download/v1.0.0/get_blunted" \
              -o "{output}"

            chmod +x "{output}"
            """


if CLSD_BIN == CLSD_BIN_DEFAULT:
    rule build_clsd:
        message: "Clone & build clsd from Git (pinned commit)"
        output:
            CLSD_BIN_DEFAULT
        conda:
            BUILD_TOOLS_ENV_YML
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output})"

            # Clone the repository if not already present
            if [ ! -d build/clsd ]; then
              git clone https://github.com/Fabianexe/clsd.git build/clsd
            fi

            cd build/clsd

            # Checkout specific commit for reproducibility
            git fetch --all
            git checkout c49598fcb149b2c224a4625e0bf4b870f27ec166

            # Autotools preparation
            echo "==> Running autoreconf -f -i"
            if ! autoreconf -f -i; then
              echo "ERROR: autoreconf failed."
              echo ""
              echo "Please install autoconf and run again:"
              echo "  sudo apt-get install autoconf"
              echo "or (with conda):"
              echo "  conda install -c conda-forge autoconf"
              exit 1
            fi

            # Configure & build
            ./configure
            make -j"$(nproc)"

            # Copy the binary (adapt if needed)
            cp clsd "{output}"
            """


if LIGHTER_BIN == LIGHTER_BIN_DEFAULT:
    rule build_lighter:
        message: "Clone & build Lighter from Git (fixed commit)"
        output:
            LIGHTER_BIN_DEFAULT
        conda:
            BUILD_TOOLS_ENV_YML
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output})"

            # Clone Lighter repository if missing
            if [ ! -d build/Lighter ]; then
                git clone https://github.com/mourisl/Lighter.git build/Lighter
            fi

            cd build/Lighter

            # Check out fixed commit for reproducibility
            git fetch --all
            git checkout d8621db12352895662404b38a4b61862eaf60f6a

            # Build Lighter
            make -j"$(nproc)"

            # Copy the built binary
            cp lighter "{output}"
            chmod +x "{output}"
            """


# -----------------------------------------------------------------------------
# Helpers bench / timing / bluntification etc.
# -----------------------------------------------------------------------------
def time_to_tsv(cmd, bench_tsv, log_file, ensure_path=None, ensure_content=None):
    import os, shlex
    def as_path(x):
        try:
            return os.fspath(x)
        except TypeError:
            try:
                return x[0]
            except Exception:
                return str(x)

    bench_tsv = as_path(bench_tsv)
    log_file = as_path(log_file)

    if ensure_path:
        ep = as_path(ensure_path)
        os.makedirs(os.path.dirname(ep), exist_ok=True)
        if not os.path.exists(ep):
            with open(ep, "w", encoding="utf-8") as f:
                f.write(ensure_content if ensure_content is not None else "")

    dbg = 'export PS4="+ [bench] "; set -x; ' if DEBUG else ''
    inner = f'( {dbg}{cmd} ) > >(tee -a {shlex.quote(log_file)} >/dev/null) 2> >(tee -a {shlex.quote(log_file)} >&2)'

    lines = [
        "set -euo pipefail",
        "export LC_ALL=C LANG=C",
        f"TO_SECS={int(TIMEOUT_SECS)}",
        f"KILL_AFTER={int(TIMEOUT_KILL_AFTER)}",
        f"TB_BIN={shlex.quote(TIMEOUT_BIN)}",
        f"TIME_BIN={shlex.quote(TIME_BIN)}",
        f"TSV_TMP={shlex.quote(bench_tsv)}.tmp",
        f"TSV_OUT={shlex.quote(bench_tsv)}",
        f"LOG_OUT={shlex.quote(log_file)}",
        'mkdir -p "$(dirname "$LOG_OUT")"',
        'mkdir -p "$(dirname "$TSV_OUT")"',
        'mkdir -p "$(dirname "$TSV_TMP")"',
        ': > "$LOG_OUT"',
        ': > "$TSV_OUT" || true',
    ]

    if TIMEOUT_ENABLED:
        lines += [
            'TB=""',
            'if command -v timeout >/dev/null 2>&1; then TB=$(command -v timeout);'
            'elif [ -x "$TB_BIN" ]; then TB="$TB_BIN"; fi',
            'use_timeout=0; [ -n "$TB" ] && use_timeout=1',
            'set +e',
            'if [ "$use_timeout" -eq 1 ]; then ',
            '  KILL_AFTER_S="$KILL_AFTER"s; TO_SECS_S="$TO_SECS"s;',
            f'  "$TIME_BIN" -v -o "$TSV_TMP" "$TB" -k "$KILL_AFTER_S" "$TO_SECS_S" bash -c {shlex.quote(inner)};',
            'else ',
            f'  "$TIME_BIN" -v -o "$TSV_TMP" bash -c {shlex.quote(inner)};',
            'fi',
            'rc=$?',
            'set -e',
            'timed_out=0',
            'if [ "$use_timeout" -eq 1 ] && [ $rc -eq 124 ]; then timed_out=1; fi',
            'if [ $rc -eq 137 ]; then timed_out=1; fi',
        ]
    else:
        lines += [
            'use_timeout=0',
            'set +e',
            f'"$TIME_BIN" -v -o "$TSV_TMP" bash -c {shlex.quote(inner)};',
            'rc=$?',
            'set -e',
            'timed_out=0',
        ]

    lines += [
        'grep -E \'Elapsed \\(wall clock\\) time|Maximum resident set size\' "$TSV_TMP" | '
        "sed 's/ (.*)//' | sed 's/^[[:space:]]*//' | tr -s ' ' '\\t' > \"$TSV_OUT\" || :",
        'printf "Timed out\\t%s\\n" "$timed_out" >> "$TSV_OUT"',
        'printf "Exit status\\t%s\\n" "$rc" >> "$TSV_OUT"',
        'if [ "$use_timeout" -eq 1 ]; then printf "Timeout (s)\\t%s\\n" "$TO_SECS" >> "$TSV_OUT"; fi',
        '[ -s "$TSV_OUT" ] || : > "$TSV_OUT"',
        'if [ "$rc" -ne 0 ]; then printf "[ERROR] Command exited with rc=%s (timed_out=%s).\\n" "$rc" "$timed_out" >> "$LOG_OUT"; fi',
        'if [ "$use_timeout" -eq 1 ] && [ "$timed_out" -eq 1 ]; then printf "[TIMEOUT] Command exceeded %ss (rc=%s).\\n" "$TO_SECS" "$rc" >> "$LOG_OUT"; fi',
        'rm -f "$TSV_TMP" || true',
        'exit 0',
    ]

    script = "\n".join(lines)
    return f"""
    set -euo pipefail
    mkdir -p {shlex.quote(os.path.dirname(bench_tsv))}
    mkdir -p {shlex.quote(os.path.dirname(log_file))}
    bash -c {shlex.quote(script)}
    """

def vg_snarls_from_gfa_cmd(gfa_in, out_path, threads):
    g = shlex.quote(gfa_in)
    o = shlex.quote(out_path)
    return f'( set -euo pipefail; vg snarls -a -T -m 1000000000 -t {int(threads)} {g} > {o} )'

def bubblegun_cmd(gfa_in, out_json):
    g = shlex.quote(gfa_in)
    o = shlex.quote(out_json)
    return (
        'BG=""; '
        'if command -v BubbleGun >/dev/null 2>&1; then BG="$(command -v BubbleGun)"; '
        'elif command -v bubblegun >/dev/null 2>&1; then BG="$(command -v bubblegun)"; '
        'fi; '
        'if [ -z "$BG" ]; then '
        '  echo "[ERROR] BubbleGun not found in the conda environment." >&2; exit 127; '
        'fi; '
        f'"$BG" -g {g} bchains --bubble_json {o}'
    )

# Blunt verification helpers
def dataset_needs_vg_snarls(name):
    try:
        return "vg_snarls_gfa" in bench_programs_for_dataset(name)
    except Exception:
        return False

def any_dataset_needs_vg_snarls():
    return any(dataset_needs_vg_snarls(n) for n in all_enabled_dataset_names())

def need_vg_snarls():
    return any_dataset_needs_vg_snarls()

def gfa_is_blunt(path, max_check=200000):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f, start=1):
                if i > max_check:
                    break
                if not line:
                    continue
                if line.startswith("E\t"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 9:
                        cigar = parts[8]
                        if cigar not in ("*", "0M", "0", "0m", "0S"):
                            return (False, i, line.strip())
                elif line.startswith("L\t"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 6:
                        overlap = parts[5]
                        if overlap not in ("*", "0M", "0", "0m", "0S"):
                            return (False, i, line.strip())
        return (True, None, None)
    except FileNotFoundError:
        return (False, -1, "file not found")

def naive_blunt_gfa(src, dst):
    with open(src, "r", encoding="utf-8", errors="replace") as fi, open(dst, "w", encoding="utf-8") as fo:
        for line in fi:
            if not line:
                continue
            if line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6:
                    parts[5] = "*"
                    line = "\t".join(parts) + "\n"
            elif line.startswith("E\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 9:
                    parts[8] = "*"
                    line = "\t".join(parts) + "\n"
            fo.write(line)

def copy_without_H_lines(src, dst):
    rx = re.compile(r'^H(\t|$)')
    with open(src, "r", encoding="utf-8", errors="replace") as fi, open(dst, "w", encoding="utf-8") as fo:
        for line in fi:
            if rx.match(line or ""):
                continue
            fo.write(line)

def bluntify_clean_impl(src_gfa, blunt_path, clean_path, log_path, require_blunt=False):
    os.makedirs(os.path.dirname(clean_path), exist_ok=True)
    rc = 1
    with open(blunt_path, "wb") as fo, open(log_path, "a", encoding="utf-8") as logf:
        try:
            p = subprocess.Popen([GET_BLUNTED, "-i", src_gfa], stdout=fo, stderr=subprocess.PIPE)
            def reader():
                for line in p.stderr:
                    try:
                        sys.stdout.buffer.write(line); sys.stdout.flush()
                    except Exception:
                        pass
                    try:
                        logf.write(line.decode(errors="replace")); logf.flush()
                    except Exception:
                        pass
            t = threading.Thread(target=reader, daemon=True); t.start()
            rc = p.wait(); t.join(timeout=0.2)
            try:
                logf.write(f"[get_blunted] exit={rc}\n")
            except Exception:
                pass
        except Exception as e:
            rc = 1
            try:
                with open(log_path, "a", encoding="utf-8") as logf2:
                    logf2.write(f"[get_blunted] exception: {e}\n")
            except Exception:
                pass

    use_src = blunt_path if (rc == 0 and os.path.isfile(blunt_path) and os.path.getsize(blunt_path) > 0) else src_gfa

    tmp_noH = clean_path + ".noH.tmp"
    copy_without_H_lines(use_src, tmp_noH)
    shutil.move(tmp_noH, clean_path)

    blunt, bad_line_no, bad_line = gfa_is_blunt(clean_path)
    if not blunt:
        try:
            with open(log_path, "a", encoding="utf-8") as logf:
                logf.write(f"[INFO] GFA not blunt, naive bluntify fallback (first bad at line {bad_line_no}): {bad_line}\n")
        except Exception:
            pass
        try:
            naive_blunt_gfa(clean_path, clean_path + ".naive.tmp")
            shutil.move(clean_path + ".naive.tmp", clean_path)
            blunt2, _, _ = gfa_is_blunt(clean_path)
            if not blunt2:
                with open(log_path, "a", encoding="utf-8") as logf:
                    logf.write("[WARN] Naive bluntify applied but GFA still not blunt.\n")
            else:
                with open(log_path, "a", encoding="utf-8") as logf:
                    logf.write("[INFO] Applied naive bluntify fallback (overlaps set to '*').\n")
        except Exception as e:
            try:
                with open(log_path, "a", encoding="utf-8") as logf:
                    logf.write(f"[WARN] Naive bluntify fallback failed: {e}\n")
            except Exception:
                pass

    if not os.path.isfile(clean_path) or os.path.getsize(clean_path) == 0:
        shutil.copyfile(src_gfa, clean_path)

# -----------------------------------------------------------------------------
# Signature helper for benchmark runs
# -----------------------------------------------------------------------------
def bench_signature(prog, dataset, rep, threads=None):
    conf = {
        "program": prog,
        "dataset": dataset,
        "rep": int(rep),
        "sink_output": SINK_OUTPUT,
        "timeout_enabled": TIMEOUT_ENABLED,
        "timeout_secs": int(TIMEOUT_SECS),
        "timeout_kill_after": int(TIMEOUT_KILL_AFTER),
    }
    if prog in ("vg_snarls_gfa", "sbSPQR_gfa", "sbSPQR_sb", "sbSPQR_snarls_gfa"):
        if threads is not None:
            conf["threads"] = int(threads)
    if prog in ("sbSPQR_gfa", "sbSPQR_snarls_gfa"):
        conf["spqr_bin"] = os.path.abspath(SPQR_BIN)
        conf["spqr_flag"] = SPQR_GFA_FLAG
    if prog == "sbSPQR_sb":
        conf["spqr_bin"] = os.path.abspath(SPQR_BIN)
        conf["spqr_flag"] = SPQR_SGRAPH_FLAG
    if prog == "clsd_sb":
        conf["clsd_bin"] = os.path.abspath(CLSD_BIN)
    blob = json.dumps(conf, sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()

SB_IO_LABELS = ["io/read_graph", "io/write_output"]              
SB_BLOCKS_LABELS = ["snarls/blocks", "sn/phase/SolveBlocks", "sb/phase/SolveBlocks"]    
SB_BUILD_PREFIXES = ["snarls/build", "sn/build", "ogdf"]         
SB_LOGIC_PREFIXES = ["snarls/logic", "sn/logic", "sb/solveSPQR", "sb/checkCutVertices", "sb/findMini"]  
SB_PHASE_PREFIXES = ["sn/phase", "sb/phase"]                    

def _get_num(m, key, default=0.0):
    try:
        v = m.get(key, default)
        if v is None:
            return float(default)
        return float(v)
    except Exception:
        return float(default)

def _bytes_to_gib(x):
    return float(x) / (1024.0**3)

def _label_is_under_prefix(label, pref):
    if not label or not pref:
        return False
    return (label == pref) or label.startswith(pref + "/") or (pref.endswith("/") and label.startswith(pref))

def _labels_under_prefixes(labels, prefixes):
    out = []
    for lab in labels:
        for pref in prefixes:
            if _label_is_under_prefix(lab, pref):
                out.append(lab)
                break
    return out

def _leaf_labels(labels):
    s = set(labels)
    leafs = []
    for lab in labels:
        has_child = any((other != lab and other.startswith(lab + "/")) for other in s)
        if has_child:
            continue
        leafs.append(lab)
    return leafs

def _sum_field_by_labels(marks_by_label, labels, field):
    s = 0.0
    for lab in labels:
        m = marks_by_label.get(lab)
        if m is None:
            continue
        s += _get_num(m, field, 0.0)
    return s

def _max_field_by_labels(marks_by_label, labels, field):
    mx = 0.0
    for lab in labels:
        m = marks_by_label.get(lab)
        if m is None:
            continue
        v = _get_num(m, field, 0.0)
        if v > mx:
            mx = v
    return mx

def _partition_categories(marks):
    labels = [str(m.get("label") or "") for m in marks if isinstance(m.get("label"), str)]
    marks_by_label = {m["label"]: m for m in marks if isinstance(m.get("label"), str)}

    io_present = [lab for lab in SB_IO_LABELS if lab in marks_by_label]
    io_sec = _sum_field_by_labels(marks_by_label, io_present, "seconds")
    io_rss = _sum_field_by_labels(marks_by_label, io_present, "rss_delta_sum_bytes")
    io_peak = _max_field_by_labels(marks_by_label, io_present, "hwm_max_bytes")

    build_all = _labels_under_prefixes(labels, SB_BUILD_PREFIXES)
    build_leaf = _leaf_labels(build_all)
    build_sec = _sum_field_by_labels(marks_by_label, build_leaf, "seconds")
    build_rss = _sum_field_by_labels(marks_by_label, build_leaf, "rss_delta_sum_bytes")
    build_peak = _max_field_by_labels(marks_by_label, build_all, "hwm_max_bytes")

    logic_all = _labels_under_prefixes(labels, SB_LOGIC_PREFIXES)
    logic_leaf = _leaf_labels(logic_all)
    logic_sec = _sum_field_by_labels(marks_by_label, logic_leaf, "seconds")
    logic_rss = _sum_field_by_labels(marks_by_label, logic_leaf, "rss_delta_sum_bytes")
    logic_peak = _max_field_by_labels(marks_by_label, logic_all, "hwm_max_bytes")

    blocks_lab = next((lab for lab in SB_BLOCKS_LABELS if lab in marks_by_label), None)
    blocks_sec = _get_num(marks_by_label.get(blocks_lab, {}), "seconds", 0.0) if blocks_lab else 0.0
    blocks_rss = _get_num(marks_by_label.get(blocks_lab, {}), "rss_delta_sum_bytes", 0.0) if blocks_lab else 0.0

    phase_labels = _labels_under_prefixes(labels, SB_PHASE_PREFIXES)
    phase_sec = _sum_field_by_labels(marks_by_label, phase_labels, "seconds") if phase_labels else 0.0
    phase_rss = _sum_field_by_labels(marks_by_label, phase_labels, "rss_delta_sum_bytes") if phase_labels else 0.0

    if blocks_lab:
        tot_sec = io_sec + blocks_sec
        tot_rss = io_rss + blocks_rss
    elif phase_labels:
        tot_sec = io_sec + phase_sec
        tot_rss = io_rss + phase_rss
    else:
        tot_sec = io_sec + build_sec + logic_sec
        tot_rss = io_rss + build_rss + logic_rss

    tot_peak = 0.0
    for m in marks:
        v = _get_num(m, "hwm_max_bytes", 0.0)
        if v > tot_peak:
            tot_peak = v

    return {
        "io":    {"s": io_sec,    "rss": io_rss,    "peak": io_peak},
        "build": {"s": build_sec, "rss": build_rss, "peak": build_peak},
        "logic": {"s": logic_sec, "rss": logic_rss, "peak": logic_peak},
        "blocks": {"s": blocks_sec, "rss": blocks_rss},
        "total": {"s": tot_sec,   "rss": tot_rss,   "peak": tot_peak},
        "has_blocks": bool(blocks_lab),
        "has_phase": bool(phase_labels),
    }

def sb_summarize_report(json_path):
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        marks = data.get("marks") or []
        if not isinstance(marks, list):
            marks = []
        ok = 1
    except Exception:
        marks = []
        ok = 0

    cats = _partition_categories(marks)

    labels = [str(m.get("label") or "") for m in marks if isinstance(m.get("label"), str)]
    marks_by_label = {m["label"]: m for m in marks if isinstance(m.get("label"), str)}

    def sum_sec(labs):
        return _sum_field_by_labels(marks_by_label, labs, "seconds")

    ogdf_bctree_sec = sum_sec(["ogdf/BCTree::ctor"])
    ogdf_spqr_sec   = sum_sec(["ogdf/StaticSPQRTree::ctor"])
    ogdf_cc_sec     = sum_sec(["ogdf/connectedComponents"])

    ogdf_bctree_rss = _sum_field_by_labels(marks_by_label, ["ogdf/BCTree::ctor"], "rss_delta_sum_bytes")
    ogdf_spqr_rss   = _sum_field_by_labels(marks_by_label, ["ogdf/StaticSPQRTree::ctor"], "rss_delta_sum_bytes")
    ogdf_cc_rss     = _sum_field_by_labels(marks_by_label, ["ogdf/connectedComponents"], "rss_delta_sum_bytes")

    out = {
        "SB_report_ok": ok,
        "SB_report_path": json_path,
        "SB_marks": f"{len(marks)}",

        "SB_tot_phases_s": f"{cats['total']['s']:.6f}",
        "SB_io_s": f"{cats['io']['s']:.6f}",
        "SB_build_s": f"{cats['build']['s']:.6f}",
        "SB_logic_s": f"{cats['logic']['s']:.6f}",

        "SB_tot_phases_rss_delta_sum_bytes": f"{int(cats['total']['rss'])}",
        "SB_io_rss_delta_sum_bytes": f"{int(cats['io']['rss'])}",
        "SB_build_rss_delta_sum_bytes": f"{int(cats['build']['rss'])}",
        "SB_logic_rss_delta_sum_bytes": f"{int(cats['logic']['rss'])}",

        "SB_tot_phases_rss_delta_sum_GiB": f"{_bytes_to_gib(cats['total']['rss']):.6f}",
        "SB_io_rss_delta_sum_GiB": f"{_bytes_to_gib(cats['io']['rss']):.6f}",
        "SB_build_rss_delta_sum_GiB": f"{_bytes_to_gib(cats['build']['rss']):.6f}",
        "SB_logic_rss_delta_sum_GiB": f"{_bytes_to_gib(cats['logic']['rss']):.6f}",

        "SB_io_peak_hwm_bytes": f"{int(cats['io']['peak'])}",
        "SB_build_peak_hwm_bytes": f"{int(cats['build']['peak'])}",
        "SB_logic_peak_hwm_bytes": f"{int(cats['logic']['peak'])}",
        "SB_tot_peak_hwm_bytes": f"{int(cats['total']['peak'])}",

        "SB_io_peak_hwm_GiB": f"{_bytes_to_gib(cats['io']['peak']):.6f}",
        "SB_build_peak_hwm_GiB": f"{_bytes_to_gib(cats['build']['peak']):.6f}",
        "SB_logic_peak_hwm_GiB": f"{_bytes_to_gib(cats['logic']['peak']):.6f}",
        "SB_tot_peak_hwm_GiB": f"{_bytes_to_gib(cats['total']['peak']):.6f}",

        "SB_ogdf_BCTree_ctor_s": f"{ogdf_bctree_sec:.6f}",
        "SB_ogdf_StaticSPQRTree_ctor_s": f"{ogdf_spqr_sec:.6f}",
        "SB_ogdf_connectedComponents_s": f"{ogdf_cc_sec:.6f}",

        "SB_ogdf_BCTree_ctor_rss_delta_sum_bytes": f"{int(ogdf_bctree_rss)}",
        "SB_ogdf_StaticSPQRTree_ctor_rss_delta_sum_bytes": f"{int(ogdf_spqr_rss)}",
        "SB_ogdf_connectedComponents_rss_delta_sum_bytes": f"{int(ogdf_cc_rss)}",

        "SB_ogdf_BCTree_ctor_rss_delta_sum_GiB": f"{_bytes_to_gib(ogdf_bctree_rss):.6f}",
        "SB_ogdf_StaticSPQRTree_ctor_rss_delta_sum_GiB": f"{_bytes_to_gib(ogdf_spqr_rss):.6f}",
        "SB_ogdf_connectedComponents_rss_delta_sum_GiB": f"{_bytes_to_gib(ogdf_cc_rss):.6f}",
    }
    return out

def sb_augment_report_json(json_path):
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        data = {}

    marks = data.get("marks") or []
    if not isinstance(marks, list):
        marks = []

    cats = _partition_categories(marks)
    summary = {
        "version": 2,
        "units": {"time": "seconds", "memory": "bytes"},
        "categories": {
            "io":    {"seconds": cats["io"]["s"],    "rss_delta_sum_bytes": int(cats["io"]["rss"]),    "peak_hwm_bytes": int(cats["io"]["peak"])},
            "build": {"seconds": cats["build"]["s"], "rss_delta_sum_bytes": int(cats["build"]["rss"]), "peak_hwm_bytes": int(cats["build"]["peak"])},
            "logic": {"seconds": cats["logic"]["s"], "rss_delta_sum_bytes": int(cats["logic"]["rss"]), "peak_hwm_bytes": int(cats["logic"]["peak"])},
        },
        "total": {
            "seconds": cats["total"]["s"],
            "rss_delta_sum_bytes": int(cats["total"]["rss"]),
            "peak_hwm_bytes": int(cats["total"]["peak"]),
            "computed_from": ("io+blocks" if cats["has_blocks"] else ("io+sn/phase/*" if cats["has_phase"] else "io+build(leaves)+logic(leaves)")),
        },
        "notes": "BUILD = OGDF build/preproc; LOGIC = solving on the built graph; sum of leaf nodes per category to avoid double counting.",
    }
    data["summary"] = summary
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)

# -----------------------------------------------------------------------------
# Quels binaires doit-on construire automatiquement ?
# -----------------------------------------------------------------------------
def tool_binaries_to_build(wildcards):
    files = []

    # sbSPQR_* (GFA/snarls/sgraph)
    need_sbspqr = any(is_program_sbspqr(p) for p in all_selected_programs())
    if need_sbspqr and SPQR_BIN == SPQR_BIN_DEFAULT:
        files.append(SPQR_BIN_DEFAULT)

    # clsd pour clsd_sb
    need_clsd = "clsd_sb" in all_selected_programs()
    if need_clsd and CLSD_BIN == CLSD_BIN_DEFAULT:
        files.append(CLSD_BIN_DEFAULT)

    # Lighter si un dataset utilise ggcat_from_reads_lighter
    need_lighter = any(d.get("builder") == "ggcat_from_reads_lighter" for d in DS.values())
    if need_lighter and LIGHTER_BIN == LIGHTER_BIN_DEFAULT:
        files.append(LIGHTER_BIN_DEFAULT)

    # get_blunted pour vg_snarls_gfa (et autres vérif blunt)
    need_get_blunted = need_vg_snarls()
    if need_get_blunted and GET_BLUNTED == GET_BLUNTED_DEFAULT:
        files.append(GET_BLUNTED_DEFAULT)

    return files

# -----------------------------------------------------------------------------
# Règle finale
# -----------------------------------------------------------------------------
rule all:
    message: "Final rule: prechecks + aggregation + plots"
    input:
        os.path.join(OUT_DIR, ".prechecks.ok"),
        os.path.join(OUT_DIR, "benchmarks.tsv"),
        plot_paths()

rule bubblegun_smoketest:
    message: "Check BubbleGun availability in its conda environment"
    output:
        os.path.join(OUT_DIR, ".bubblegun.ok")
    conda: BUBBLEGUN_ENV_YML
    shell:
        r"""
        set -euo pipefail
        BG=$(command -v BubbleGun || command -v bubblegun || true)
        if [ -z "$BG" ]; then
            echo "BubbleGun not found in conda env {BUBBLEGUN_ENV_YML}" >&2
            exit 1
        fi
        "$BG" -h >/dev/null 2>&1 || "$BG" --help >/dev/null 2>&1 || true
        touch {output}
        """

rule pggb_smoketest:
    message: "Check pggb availability in its conda environment"
    output:
        os.path.join(OUT_DIR, ".pggb.ok")
    conda: PGGB_ENV_YML
    shell:
        r"""
        set -euo pipefail
        if ! command -v pggb >/dev/null 2>&1; then
            echo "[ERROR] pggb not found in conda env {PGGB_ENV_YML}" >&2
            conda list >&2 || true
            exit 127
        fi
        if ! command -v samtools >/dev/null 2>&1; then
            echo "[ERROR] samtools not found in conda env {PGGB_ENV_YML}" >&2
            conda list >&2 || true
            exit 127
        fi
        pggb -h >/dev/null 2>&1 || pggb --help >/dev/null 2>&1 || true
        samtools --version >/dev/null 2>&1 || samtools --help >/dev/null 2>&1 || true
        touch {output}
        """

rule vg_smoketest:
    message: "Check vg availability in its conda environment"
    output:
        os.path.join(OUT_DIR, ".vg.ok")
    conda: VG_ENV_YML
    shell:
        r"""
        set -euo pipefail
        if ! command -v vg >/dev/null 2>&1; then
            echo "[ERROR] vg not found in conda env {VG_ENV_YML}" >&2
            conda list >&2 || true
            exit 127
        fi
        if ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
            echo "[ERROR] bgzip/tabix not found in conda env {VG_ENV_YML}. Add 'htslib' to config/vg.yml." >&2
            conda list >&2 || true
            exit 127
        fi
        vg version >/dev/null 2>&1 || vg --version >/dev/null 2>&1 || vg -h >/dev/null 2>&1 || true
        touch {output}
        """

rule prechecks:
    message: "Prechecks (binaries, conda/mamba, directories)"
    input:
        os.path.join(OUT_DIR, ".bubblegun.ok"),
        os.path.join(OUT_DIR, ".vg.ok"),
        os.path.join(OUT_DIR, ".pggb.ok"),
        tool_binaries_to_build  # déclenche la construction des outils nécessaires
    output:
        os.path.join(OUT_DIR, ".prechecks.ok")
    run:
        def is_executable(path):
            return os.path.isfile(path) and os.access(path, os.X_OK)

        def is_elf(path):
            try:
                with open(path, "rb") as f:
                    return f.read(4) == b"\x7fELF"
            except Exception:
                return False

        for p in [
            os.path.join(OUT_DIR, "logs", "ggcat"),
            os.path.join(OUT_DIR, "logs", "vg"),
            os.path.join(OUT_DIR, "logs", "lighter"),
            os.path.join(OUT_DIR, "logs", "bench"),
            os.path.join(OUT_DIR, "logs", "sgraph"),
            os.path.join(OUT_DIR, "logs", "pggb"),
            os.path.join(OUT_DIR, "bench"),
            os.path.join(OUT_DIR, "prog_out"),
            os.path.join(OUT_DIR, "plots"),
            os.path.join(OUT_DIR, "summary"),
        ]:
            os.makedirs(p, exist_ok=True)

        if not is_executable(TIME_BIN):
            raise WorkflowError(f"time not found or not executable: {TIME_BIN}")

        if TIMEOUT_ENABLED and not (os.path.isabs(TIMEOUT_BIN) and is_executable(TIMEOUT_BIN)) and shutil.which(TIMEOUT_BIN) is None:
            print(f"[WARN] timeout not found ({TIMEOUT_BIN}); benchmark commands will run without timeout.", flush=True)

        if not is_executable(SPQR_BIN):
            raise WorkflowError(f"sbSPQR (sbfind/BubbleFinder) not found or not executable: {SPQR_BIN}")

        if shutil.which("conda") is None and shutil.which("mamba") is None:
            raise WorkflowError("conda/mamba not found in PATH")

        if not is_executable(GET_BLUNTED) or not is_elf(GET_BLUNTED):
            print(f"[WARN] get_blunted missing or non-ELF ({GET_BLUNTED}); will fall back to input if it fails.", flush=True)

        need_lighter = any(d.get("builder") == "ggcat_from_reads_lighter" for d in DS.values())
        if need_lighter:
            bins = set()
            global_bin = TOOLS.get("lighter_bin", LIGHTER_BIN_DEFAULT)
            if global_bin:
                bins.add(global_bin)
            for d in DS.values():
                if d.get("builder") == "ggcat_from_reads_lighter":
                    b = (d.get("lighter", {}) or {}).get("bin") or global_bin
                    if b:
                        bins.add(b)
            for b in bins:
                if not is_executable(b):
                    raise WorkflowError(f"Lighter not found or not executable: {b}")

        if "clsd_sb" in all_selected_programs():
            if not is_executable(CLSD_BIN):
                raise WorkflowError(f"clsd not found or not executable: {CLSD_BIN}")

        if need_vg_snarls():
            if not (os.path.isfile(GET_BLUNTED) and os.access(GET_BLUNTED, os.X_OK)):
                print(
                    f"[WARN] vg_snarls_gfa: get_blunted not found/executable at {GET_BLUNTED}. "
                    "Naive fallback will be used if necessary. "
                    "For better guarantees, install https://github.com/vgteam/GetBlunted and set defaults.tools.get_blunted.",
                    flush=True
                )

        open(output[0], "w").close()

rule gfa_bluntify_clean:
    message: "Bluntify + clean GFA for {wildcards.dataset}"
    input:
        lambda wc: raw_gfa_path(wc.dataset)
    output:
        os.path.join(DATA_DIR, "{dataset}", "{dataset}.cleaned.gfa")
    run:
        dataset = wildcards.dataset
        src = input[0]
        blunt = os.path.join(DATA_DIR, dataset, f"{dataset}.bluntified.gfa")
        clean = output[0]
        logp = blunt + ".log"
        os.makedirs(os.path.dirname(clean), exist_ok=True)
        bluntify_clean_impl(src, blunt, clean, logp, require_blunt=True)

# -----------------------------------------------------------------------------
# Optional SB reconstruction via ggcat -f
# -----------------------------------------------------------------------------
def _ggcat_sb_input_src(wc):
    b = builder(wc.dataset)
    if b == "ggcat_from_fasta":
        return rules.ggcat_from_fasta_decompress.output.fa
    elif b == "ggcat_from_reads_lighter":
        return rules.ggcat_lighter_correct.output.corrected
    raise WorkflowError(f"ggcat -f requested for non-ggcat dataset: {wc.dataset} (builder={b})")

rule ggcat_sb_raw_gfa:
    message: "ggcat -f (SB) build for {wildcards.dataset}"
    input:
        src=_ggcat_sb_input_src
    output:
        os.path.join(DATA_DIR, "{dataset}", "{dataset}.ggcat.sb.gfa")
    conda: GG_ENV_YML
    threads: THREADS_GGCAT
    log:
        sb=os.path.join(OUT_DIR, "logs", "ggcat", "{dataset}.sb.log")
    params:
        k=lambda wc: int(ds(wc.dataset).get("ggcat", {}).get("k", 31)),
        common=GGCAT_COMMON
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log.sb})"
        ggcat build --gfa-v1 {params.common} -f -k {params.k} -j {threads} {input.src} \
            -o {output} -t {output}.temp 2>&1 | tee -a {log.sb}
        """

rule ggcat_sb_bluntify_clean:
    message: "Bluntify + clean SB GFA for {wildcards.dataset}"
    input:
        os.path.join(DATA_DIR, "{dataset}", "{dataset}.ggcat.sb.gfa")
    output:
        os.path.join(DATA_DIR, "{dataset}", "{dataset}.sb.cleaned.gfa")
    run:
        dataset = wildcards.dataset
        src = input[0]
        blunt = os.path.join(DATA_DIR, dataset, f"{dataset}.sb.bluntified.gfa")
        clean = output[0]
        logp = blunt + ".log"
        os.makedirs(os.path.dirname(clean), exist_ok=True)
        bluntify_clean_impl(src, blunt, clean, logp, require_blunt=True)

# -----------------------------------------------------------------------------
# we prep here clsd edgelist and sbSPQR sgraph
# -----------------------------------------------------------------------------
rule clsd_prepare_edgelist:
    message: "Prepare edgelist for clsd (skip header) for {wildcards.dataset}"
    input:
        gfa=lambda wc: gfa_for_sb_inputs(wc.dataset)
    output:
        os.path.join(DATA_DIR, "{dataset}", "{dataset}.clsd.edgelist")
    conda: PYTOOLS_ENV_YML
    log:
        os.path.join(OUT_DIR, "logs", "sgraph", "{dataset}.clsd_prep.log")
    params:
        script=lambda wc: resolve_env("scripts/gfa_to_sgraph.py"),
        flags=lambda wc: sgraph_flags_for_dataset(wc.dataset)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        tmp="{output}.tmp.sgraph"
        python "{params.script}" "{input.gfa}" "$tmp" --map-out /dev/null {params.flags} \
          2> >(tee -a "{log}" >&2)
        tail -n +2 "$tmp" > "{output}"
        rm -f "$tmp"
        """

rule sbspqr_prepare_sgraph:
    message: "Prepare sgraph for BubbleFinder for {wildcards.dataset}"
    input:
        gfa=lambda wc: gfa_for_sb_inputs(wc.dataset)
    output:
        os.path.join(DATA_DIR, "{dataset}", "{dataset}.sbspqr.sgraph")
    conda: PYTOOLS_ENV_YML
    log:
        os.path.join(OUT_DIR, "logs", "sgraph", "{dataset}.sbspqr_prep.log")
    params:
        script=lambda wc: resolve_env("scripts/gfa_to_sgraph.py"),
        flags=lambda wc: sgraph_flags_for_dataset(wc.dataset)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        python "{params.script}" "{input.gfa}" "{output}" --map-out /dev/null {params.flags} \
          2> >(tee -a "{log}" >&2)
        """

# -----------------------------------------------------------------------------
# Bench rules 
# -----------------------------------------------------------------------------
rule bench_BubbleGun_gfa:
    message: "Bench BubbleGun on {wildcards.dataset} rep={wildcards.rep} t={wildcards.t}"
    input:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.cleaned.gfa"),
        pre=rules.prechecks.output
    output:
        tsv=os.path.join(OUT_DIR, "bench", "BubbleGun_gfa", "{dataset}.t{t}.rep{rep}.tsv"),
        prog=os.path.join(OUT_DIR, "prog_out", "BubbleGun_gfa", "{dataset}.t{t}.rep{rep}.out")
    log:
        os.path.join(OUT_DIR, "logs", "bench", "BubbleGun_gfa", "{dataset}.t{t}.rep{rep}.log")
    conda: BUBBLEGUN_ENV_YML
    resources:
        benchmark = 200
    run:
        os.makedirs(os.path.dirname(output.prog), exist_ok=True)
        cmd = bubblegun_cmd(input.gfa, output.prog if not SINK_OUTPUT else "/dev/null")
        shell(time_to_tsv(cmd, output.tsv, log, ensure_path=(None if SINK_OUTPUT else output.prog), ensure_content=("[]" if not SINK_OUTPUT else None)))
        with open(output.tsv, "a", encoding="utf-8") as f:
            f.write("Program\tBubbleGun_gfa\n")
            f.write(f"Class\t{program_class('BubbleGun_gfa')}\n")
            f.write("Input\tGFA\n")
            f.write(f"Threads\t{wildcards.t}\n")
            sig = bench_signature("BubbleGun_gfa", wildcards.dataset, wildcards.rep)
            f.write("SigV\t1\n")
            f.write(f"Signature\t{sig}\n")

rule bench_sbSPQR_gfa:
    message: "Bench BubbleFinder (bidirectional) on {wildcards.dataset} rep={wildcards.rep} t={wildcards.t}"
    input:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.cleaned.gfa"),
        pre=rules.prechecks.output
    output:
        tsv=os.path.join(OUT_DIR, "bench", "sbSPQR_gfa", "{dataset}.t{t}.rep{rep}.tsv"),
        prog=os.path.join(OUT_DIR, "prog_out", "sbSPQR_gfa", "{dataset}.t{t}.rep{rep}.out"),
        report=os.path.join(OUT_DIR, "prog_out", "sbSPQR_gfa", "{dataset}.t{t}.rep{rep}.report.json")
    log:
        os.path.join(OUT_DIR, "logs", "bench", "sbSPQR_gfa", "{dataset}.t{t}.rep{rep}.log")
    threads:
        lambda wc: int(wc.t)
    resources:
        benchmark = 200
    run:
        os.makedirs(os.path.dirname(output.prog), exist_ok=True)
        os.makedirs(os.path.dirname(output.report), exist_ok=True)
        target = "/dev/null" if SINK_OUTPUT else output.prog
        t = int(wildcards.t)
        thr_flag = (f" {SPQR_THREADS_FLAG} {t}" if SPQR_THREADS_FLAG else "")
        report_clause = (f"{SPQR_REPORT_FLAG}{shlex.quote(output.report)}"
                         if SPQR_REPORT_FLAG.endswith("=")
                         else f"{SPQR_REPORT_FLAG} {shlex.quote(output.report)}")
        cmd = f"OMP_NUM_THREADS={t} {shlex.quote(SPQR_BIN)} -g {shlex.quote(input.gfa)} {SPQR_GFA_FLAG}{thr_flag} {report_clause} -o {shlex.quote(target)}"
        shell(time_to_tsv(cmd, output.tsv, log, ensure_path=(None if SINK_OUTPUT else output.prog)))
        if not os.path.isfile(output.report) or os.path.getsize(output.report) == 0:
            with open(output.report, "w", encoding="utf-8") as jf:
                jf.write("{}")
        with open(output.tsv, "a", encoding="utf-8") as f:
            f.write("Program\tsbSPQR_gfa\n")
            f.write(f"Class\t{program_class('sbSPQR_gfa')}\n")
            f.write("Input\tGFA\n")
            f.write(f"Threads\t{t}\n")
            sig = bench_signature("sbSPQR_gfa", wildcards.dataset, wildcards.rep, threads=t)
            f.write("SigV\t1\n")
            f.write(f"Signature\t{sig}\n")
        vals = sb_summarize_report(output.report)
        with open(output.tsv, "a", encoding="utf-8") as f:
            for k, v in vals.items():
                f.write(f"{k}\t{v}\n")
        try:
            sb_augment_report_json(output.report)
        except Exception:
            pass

rule bench_sbSPQR_snarls_gfa:
    message: "Bench BubbleFinder snarls (GFA) on {wildcards.dataset} rep={wildcards.rep} t={wildcards.t}"
    input:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.cleaned.gfa"),
        pre=rules.prechecks.output
    output:
        tsv=os.path.join(OUT_DIR, "bench", "sbSPQR_snarls_gfa", "{dataset}.t{t}.rep{rep}.tsv"),
        prog=os.path.join(OUT_DIR, "prog_out", "sbSPQR_snarls_gfa", "{dataset}.t{t}.rep{rep}.out"),
        report=os.path.join(OUT_DIR, "prog_out", "sbSPQR_snarls_gfa", "{dataset}.t{t}.rep{rep}.report.json")
    log:
        os.path.join(OUT_DIR, "logs", "bench", "sbSPQR_snarls_gfa", "{dataset}.t{t}.rep{rep}.log")
    threads:
        lambda wc: int(wc.t)
    resources:
        benchmark = 200
    run:
        os.makedirs(os.path.dirname(output.prog), exist_ok=True)
        os.makedirs(os.path.dirname(output.report), exist_ok=True)
        target = "/dev/null" if SINK_OUTPUT else output.prog
        t = int(wildcards.t)
        thr_flag = (f" {SPQR_THREADS_FLAG} {t}" if SPQR_THREADS_FLAG else "")
        report_clause = (f"{SPQR_REPORT_FLAG}{shlex.quote(output.report)}"
                         if SPQR_REPORT_FLAG.endswith("=")
                         else f"{SPQR_REPORT_FLAG} {shlex.quote(output.report)}")
        cmd = f"OMP_NUM_THREADS={t} {shlex.quote(SPQR_BIN)} -g {shlex.quote(input.gfa)} {SPQR_GFA_FLAG} --snarls{thr_flag} {report_clause} -o {shlex.quote(target)}"
        shell(time_to_tsv(cmd, output.tsv, log, ensure_path=(None if SINK_OUTPUT else output.prog)))
        if not os.path.isfile(output.report) or os.path.getsize(output.report) == 0:
            with open(output.report, "w", encoding="utf-8") as jf:
                jf.write("{}")
        with open(output.tsv, "a", encoding="utf-8") as f:
            f.write("Program\tsbSPQR_snarls_gfa\n")
            f.write(f"Class\t{program_class('sbSPQR_snarls_gfa')}\n")
            f.write("Input\tGFA\n")
            f.write(f"Threads\t{t}\n")
            sig = bench_signature("sbSPQR_snarls_gfa", wildcards.dataset, wildcards.rep, threads=t)
            f.write("SigV\t1\n")
            f.write(f"Signature\t{sig}\n")
        vals = sb_summarize_report(output.report)
        with open(output.tsv, "a", encoding="utf-8") as f:
            for k, v in vals.items():
                f.write(f"{k}\t{v}\n")
        try:
            sb_augment_report_json(output.report)
        except Exception:
            pass

rule bench_sbSPQR_sb:
    message: "Bench BubbleFinder (unidirectional, sgraph) on {wildcards.dataset} rep={wildcards.rep} t={wildcards.t}"
    input:
        sgraph=os.path.join(DATA_DIR, "{dataset}", "{dataset}.sbspqr.sgraph"),
        pre=rules.prechecks.output
    output:
        tsv=os.path.join(OUT_DIR, "bench", "sbSPQR_sb", "{dataset}.t{t}.rep{rep}.tsv"),
        prog=os.path.join(OUT_DIR, "prog_out", "sbSPQR_sb", "{dataset}.t{t}.rep{rep}.out"),
        report=os.path.join(OUT_DIR, "prog_out", "sbSPQR_sb", "{dataset}.t{t}.rep{rep}.report.json")
    log:
        os.path.join(OUT_DIR, "logs", "bench", "sbSPQR_sb", "{dataset}.t{t}.rep{rep}.log")
    threads:
        lambda wc: int(wc.t)
    resources:
        benchmark = 200
    run:
        os.makedirs(os.path.dirname(output.prog), exist_ok=True)
        os.makedirs(os.path.dirname(output.report), exist_ok=True)
        target = "/dev/null" if SINK_OUTPUT else output.prog
        t = int(wildcards.t)
        thr_flag = (f" {SPQR_THREADS_FLAG} {t}" if SPQR_THREADS_FLAG else "")
        report_clause = (f"{SPQR_REPORT_FLAG}{shlex.quote(output.report)}"
                         if SPQR_REPORT_FLAG.endswith("=")
                         else f"{SPQR_REPORT_FLAG} {shlex.quote(output.report)}")
        cmd = f"OMP_NUM_THREADS={t} {shlex.quote(SPQR_BIN)} -g {shlex.quote(input.sgraph)}{thr_flag} {report_clause} -o {shlex.quote(target)}"
        shell(time_to_tsv(cmd, output.tsv, log, ensure_path=(None if SINK_OUTPUT else output.prog)))
        if not os.path.isfile(output.report) or os.path.getsize(output.report) == 0:
            with open(output.report, "w", encoding="utf-8") as jf:
                jf.write("{}")
        with open(output.tsv, "a", encoding="utf-8") as f:
            f.write("Program\tsbSPQR_sb\n")
            f.write(f"Class\t{program_class('sbSPQR_sb')}\n")
            f.write("Input\tsgraph\n")
            f.write(f"SB_ggcat_f\t{sb_ggcat_f_value(wildcards.dataset)}\n")
            f.write(f"Threads\t{t}\n")
            sig = bench_signature("sbSPQR_sb", wildcards.dataset, wildcards.rep, threads=t)
            f.write("SigV\t1\n")
            f.write(f"Signature\t{sig}\n")
        vals = sb_summarize_report(output.report)
        with open(output.tsv, "a", encoding="utf-8") as f:
            for k, v in vals.items():
                f.write(f"{k}\t{v}\n")
        try:
            sb_augment_report_json(output.report)
        except Exception:
            pass

rule bench_vg_snarls_gfa:
    message: "Bench vg snarls on {wildcards.dataset} rep={wildcards.rep} t={wildcards.t}"
    input:
        gfa=os.path.join(DATA_DIR, "{dataset}", "{dataset}.cleaned.gfa"),
        pre=rules.prechecks.output
    output:
        tsv=os.path.join(OUT_DIR, "bench", "vg_snarls_gfa", "{dataset}.t{t}.rep{rep}.tsv"),
        prog=os.path.join(OUT_DIR, "prog_out", "vg_snarls_gfa", "{dataset}.t{t}.rep{rep}.out")
    log:
        os.path.join(OUT_DIR, "logs", "bench", "vg_snarls_gfa", "{dataset}.t{t}.rep{rep}.log")
    conda: VG_ENV_YML
    threads:
        lambda wc: int(wc.t)
    resources:
        benchmark = 200
    run:
        import os, shlex
        os.makedirs(os.path.dirname(output.prog), exist_ok=True)
        logfile = log[0] if isinstance(log, (list, tuple)) else str(log)
        os.makedirs(os.path.dirname(logfile), exist_ok=True)

        shell(r"""
            set -euo pipefail
            gfa={gfa}
            log={log}
            (
              echo "=== vg debug prelude ==="
              vg version || vg --version || true
              echo "File: $gfa"
              head -n1 "$gfa" || true
              echo "#H-lines:"; grep -cE '^[H]([[:space:]]|$)' "$gfa" || true
              echo "#W-lines:"; grep -cE '^[W]([[:space:]]|$)' "$gfa" || true
              echo "#P-lines:"; grep -cE '^[P]([[:space:]]|$)' "$gfa" || true
              echo "vg stats -F:"
              vg stats -F "$gfa" || true
              echo "=== end prelude ==="
            ) >> "$log" 2>&1
        """.format(gfa=shlex.quote(input.gfa), log=shlex.quote(logfile)))

        target = "/dev/null" if SINK_OUTPUT else output.prog
        cmd = vg_snarls_from_gfa_cmd(input.gfa, target, threads)
        shell(time_to_tsv(cmd, output.tsv, logfile, ensure_path=(None if SINK_OUTPUT else output.prog)))
        with open(output.tsv, "a", encoding="utf-8") as f:
            f.write("Program\tvg_snarls_gfa\n")
            f.write(f"Class\t{program_class('vg_snarls_gfa')}\n")
            f.write("Input\tGFA\n")
            f.write(f"Threads\t{threads}\n")
            sig = bench_signature("vg_snarls_gfa", wildcards.dataset, wildcards.rep, threads=threads)
            f.write("SigV\t1\n")
            f.write(f"Signature\t{sig}\n")

rule bench_clsd_sb:
    message: "Bench clsd (SB) on {wildcards.dataset} rep={wildcards.rep} t={wildcards.t}"
    input:
        edgelist=os.path.join(DATA_DIR, "{dataset}", "{dataset}.clsd.edgelist"),
        pre=rules.prechecks.output
    output:
        tsv=os.path.join(OUT_DIR, "bench", "clsd_sb", "{dataset}.t{t}.rep{rep}.tsv"),
        prog=os.path.join(OUT_DIR, "prog_out", "clsd_sb", "{dataset}.t{t}.rep{rep}.out")
    log:
        os.path.join(OUT_DIR, "logs", "bench", "clsd_sb", "{dataset}.t{t}.rep{rep}.log")
    resources:
        benchmark = 200
    run:
        os.makedirs(os.path.dirname(output.prog), exist_ok=True)
        target = "/dev/null" if SINK_OUTPUT else output.prog
        cmd = f"{shlex.quote(CLSD_BIN)} {shlex.quote(input.edgelist)} > {shlex.quote(target)}"
        shell(time_to_tsv(cmd, output.tsv, log, ensure_path=(None if SINK_OUTPUT else output.prog)))
        with open(output.tsv, "a", encoding="utf-8") as f:
            f.write("Program\tclsd_sb\n")
            f.write(f"Class\t{program_class('clsd_sb')}\n")
            f.write("Input\tedgelist\n")
            f.write(f"SB_ggcat_f\t{sb_ggcat_f_value(wildcards.dataset)}\n")
            f.write(f"Threads\t{wildcards.t}\n")
            sig = bench_signature("clsd_sb", wildcards.dataset, wildcards.rep)
            f.write("SigV\t1\n")
            f.write(f"Signature\t{sig}\n")


def ALL_TSVS():
    tsvs = []
    for dname in all_enabled_dataset_names():
        progs = bench_programs_for_dataset(dname)
        for p in progs:
            tset = threads_for_prog_dataset(p, dname)
            for t in tset:
                for r in range(1, REPS + 1):
                    tsvs.append(os.path.join(OUT_DIR, "bench", p, f"{dname}.t{t}.rep{r}.tsv"))
    return tsvs


rule aggregate_and_plot:
    message: "Aggregate benchmarks and generate plots"
    input:
        ALL_TSVS(),
        "datasets.yaml"
    output:
        tsv=os.path.join(OUT_DIR, "benchmarks.tsv"),
        time_png=os.path.join(OUT_DIR, "plots", "time_by_dataset_program.png"),
        rss_png=os.path.join(OUT_DIR, "plots", "rss_by_dataset_program.png"),
        rerun=os.path.join(OUT_DIR, "summary", "reruns_planned.tsv")
    params:
        envdir=lambda wc: ".snakemake/conda/plot_cli",
        yml=PLOT_ENV_YML
    shell:
        r"""
        set -euo pipefail
        ENV="{params.envdir}"
        YML="{params.yml}"

        mkdir -p "$(dirname "{output.tsv}")" "$(dirname "{output.time_png}")" \
                 "$(dirname "{output.rss_png}")" "$(dirname "{output.rerun}")"

        if command -v micromamba >/dev/null 2>&1; then
          if [ ! -d "$ENV" ]; then
            micromamba create -y -p "$ENV" -f "$YML"
          else
            micromamba install -y -p "$ENV" -f "$YML" || micromamba update -y -p "$ENV" -f "$YML" || true
          fi
          RUN="micromamba run -p $ENV"
        else
          if [ ! -d "$ENV" ]; then
            conda env create -p "$ENV" -f "$YML" || conda env update -p "$ENV" -f "$YML" --prune
          else
            conda env update -p "$ENV" -f "$YML" --prune || true
          fi
          RUN="conda run -p $ENV"
        fi

        $RUN python scripts/aggregate_and_plot.py \
          --inputs {input} \
          --out-tsv {output.tsv} \
          --time-png {output.time_png} \
          --rss-png {output.rss_png} \
          --rerun-tsv {output.rerun}
        """

wildcard_constraints:
    dataset="|".join(all_enabled_dataset_names()),
    rep="|".join(str(i) for i in range(1, REPS + 1)),
    t="|".join(str(x) for x in all_threads_values())