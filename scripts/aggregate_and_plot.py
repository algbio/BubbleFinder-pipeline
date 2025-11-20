
#!/usr/bin/env python3
import os, re, csv, sys

try:
    snakemake  
except NameError:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--time-png", required=True)
    ap.add_argument("--rss-png", required=True)
    ap.add_argument("--rerun-tsv", required=True)
    args = ap.parse_args()
    bench_files = list(args.inputs)
    out_tsv = args.out_tsv
    time_png = args.time_png
    rss_png = args.rss_png
    rerun_tsv = args.rerun_tsv
    cfg = {}
else:
    bench_files = list(snakemake.input)
    out_tsv = snakemake.output.tsv
    time_png = snakemake.output.time_png
    rss_png = snakemake.output.rss_png
    rerun_tsv = snakemake.output.rerun
    cfg = getattr(snakemake, "config", {}) or {}

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
for bf in bench_files:
    if not os.path.exists(bf) or os.path.getsize(bf) == 0:
        continue
    prog = os.path.basename(os.path.dirname(bf))
    base = os.path.basename(bf)
    m = re.match(r'(.+)\.rep(\d+)\.tsv$', base)
    if not m:
        continue
    dataset, rep = m.group(1), int(m.group(2))
    with open(bf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            c = [x for x in line.rstrip("\n").split("\t") if x != ""]
            if len(c) < 2:
                continue
            metric = " ".join(c[:-1]).strip()
            value = c[-1].strip()
            rows.append((dataset, prog, rep, metric, value))

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
with open(out_tsv, "w", newline="") as w:
    cw = csv.writer(w, delimiter="\t")
    cw.writerow(["dataset", "program", "replicate", "metric", "value"])
    for r in rows:
        cw.writerow(r)

def parse_elapsed(s):
    s = s.strip()
    if re.match(r'^\d+:\d+:\d+(?:\.\d+)?$', s):
        h, m, ss = s.split(":"); return int(h) * 3600 + int(m) * 60 + float(ss)
    if re.match(r'^\d+:\d+(?:\.\d+)?$', s):
        m, ss = s.split(":"); return int(m) * 60 + float(ss)
    m = re.search(r'([0-9.]+)', s)
    return float(m.group(1)) if m else float("nan")

def parse_kbytes(s):
    m = re.search(r'([0-9.]+)', s)
    kb = float(m.group(1)) if m else float("nan")
    return kb / 1024.0

elapsed = {}
rss = {}
for d, p, rep, m, v in rows:
    if m.startswith("Elapsed"):
        elapsed.setdefault((d, p), []).append(parse_elapsed(v))
    elif m.startswith("Maximum resident set size"):
        rss.setdefault((d, p), []).append(parse_kbytes(v))

datasets = sorted(set([d for d, _ in elapsed.keys()] + [d for d, _ in rss.keys()]))


programs = ["BubbleGun_gfa", "sbSPQR_gfa", "vg_snarls_gfa"]
colors = {"BubbleGun_gfa": "#d62728", "sbSPQR_gfa": "#2ca02c", "vg_snarls_gfa": "#1f77b4"}

means = np.array([[np.nanmean(elapsed.get((d, p), [np.nan])) for p in programs] for d in datasets]) if datasets else np.zeros((0, len(programs)))
fig, ax = plt.subplots(figsize=(10, 4.5))
x = np.arange(len(datasets)); w = 0.25
for i, p in enumerate(programs):
    if len(datasets) > 0:
        ax.bar(x + (i - 1) * w, means[:, i], width=w, label=p, color=colors.get(p, "#1f77b4"))
ax.set_xticks(x); ax.set_xticklabels(datasets, rotation=20, ha="right")
ax.set_ylabel("Temps (s)"); ax.set_title("Temps d'exécution")
ax.legend(); fig.tight_layout()
fig.savefig(time_png, dpi=160)


rss_means = np.array([[np.nanmean(rss.get((d, p), [np.nan])) for p in programs] for d in datasets]) if datasets else np.zeros((0, len(programs)))
fig, ax = plt.subplots(figsize=(10, 4.5))
for i, p in enumerate(programs):
    if len(datasets) > 0:
        ax.bar(x + (i - 1) * w, rss_means[:, i], width=w, label=p, color=colors.get(p, "#1f77b4"))
ax.set_xticks(x); ax.set_xticklabels(datasets, rotation=20, ha="right")
ax.set_ylabel("RSS (MB)"); ax.set_title("Mémoire max")
ax.legend(); fig.tight_layout()
fig.savefig(rss_png, dpi=160)

defaults = (cfg.get("defaults") or {})
bench_cfg = (defaults.get("bench") or {})
auto_fail = bool(bench_cfg.get("auto_rerun_failed_on_next", True))
auto_timeout = bool(bench_cfg.get("auto_rerun_timedout_on_next", True))

stats = {}
for d, p, rep, m, v in rows:
    key = (d, p, rep)
    s = stats.setdefault(key, {})
    if m.startswith("Timed out"):
        s["timed_out"] = str(v)
    elif m.startswith("Exit status"):
        s["exit_status"] = str(v)
    elif m.startswith("Timeout (s)"):
        s["timeout_s"] = str(v)

plan = []
for (d, p, rep), s in stats.items():
    reasons = []
    to = s.get("timed_out")
    es = s.get("exit_status")
    if auto_timeout and to in ("1", "true", "True"):
        reasons.append("timeout")
    try:
        if auto_fail and es is not None and int(es) != 0:
            reasons.append(f"exit={es}")
    except Exception:
        pass
    if reasons:
        plan.append((d, p, rep, ",".join(reasons)))

os.makedirs(os.path.dirname(rerun_tsv), exist_ok=True)
with open(rerun_tsv, "w", newline="") as w:
    cw = csv.writer(w, delimiter="\t")
    cw.writerow(["dataset", "program", "replicate", "reason"])
    for r in plan:
        cw.writerow(r)
