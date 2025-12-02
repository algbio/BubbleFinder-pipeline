# Snakemake Pipeline for BubbleFinder 

This repository contains a Snakemake pipeline that:

- Builds graphs (GFA) from FASTA, FASTQ (ggcat+Lighter), VCF (vg), or pre‑existing GFAs.
- Cleans and "bluntifies" the graphs.
- Prepares unidirectional graph representations (sgraph / edgelist).
- Runs benchmarks on several programs (BubbleGun, BubbleFinder, vg snarls, clsd).
- Aggregates results in a single table and produces plots (time / memory).
- Computes 2‑connected block and 3‑connected component (SPQR R‑node) size distributions and plots histograms/boxplots.

---

## 1. Requirements

- Unix‑like OS (Linux, macOS).
- *Snakemake* (recommended via conda/mamba).
- *conda* or *mamba* available in `$PATH`.  
  (The pipeline uses `conda:` directives in the rules.)
- Internet access (to download data and binaries).

The following tools are managed automatically by the Snakefile if no binary is specified in `datasets.yaml`:

- [BubbleFinder](https://github.com/algbio/BubbleFinder/tree/6ba16c9527bbde295a3234b717c594d0178e3cd3), cloned and built from GitHub (commit `6ba16c95`).
- [GetBlunted](https://github.com/vgteam/GetBlunted/releases/tag/v1.0.0), precompiled binary downloaded (release `v1.0.0`).
- [clsd](https://github.com/Fabianexe/clsd/tree/c49598fcb149b2c224a4625e0bf4b870f27ec166), cloned and built (commit `c49598fc`).
- [Lighter](https://github.com/mourisl/Lighter/tree/d8621db12352895662404b38a4b61862eaf60f6a), cloned and built (commit `d8621db1`).

---

## 2. Setup

```bash
# 1) Clone the repository
git clone https://github.com/algbio/BubbleFinder-experiments.git
cd BubbleFinder-experiments

# 2) (Optional) Create a conda/mamba env with Snakemake version 9
mamba env create -f environment.yml
conda activate bench
```

Make sure `conda` *or* `mamba` is in `$PATH` when you run Snakemake.

---

## 3. Dataset configuration (`datasets.yaml`)

The `datasets.yaml` file describes:

- datasets (`datasets:`),
- builders (how to produce raw GFAs),
- tools and conda environments,
- benchmark programs and their parameters.

### Available builders

Each dataset has a `builder` field:

- `ggcat_from_fasta`  
  - Input: FASTA/FNA (or `.tar.gz` archive containing FASTA files).  
  - Example: `coli3682` (from Zenodo).
- `ggcat_from_reads_lighter`  
  - Input: FASTQ(.gz) → Lighter correction → ggcat.
- `vg_from_vcf`  
  - Input: `fa_gz` (reference) + `vcf_gz` → `vg construct` → GFA.
- `gfa_from_url`  
  - Input: pre‑built GFA downloaded from a URL.
- `pggb_from_fasta`  
  - Input: FASTA, graph built with [pggb](https://github.com/pangenome/pggb) (pangenome graphs).

### Enable / disable a dataset

Under `datasets:`:

```yaml
- name: coli3682
  enabled: true
  builder: ggcat_from_fasta
  ...
```

- `enabled: true` → used.
- `enabled: false` → ignored.
- `enabled: auto` → enabled only if input files can be detected (useful for HG00733, etc.).

### Choosing benchmark programs

Global (in `defaults.bench.programs`):

```yaml
defaults:
  bench:
    reps: 2
    programs:
      - BubbleGun_gfa // => BubbleGun, with a .GFA file as input
      - sbSPQR_gfa // => BubbleFinder, bidirectional edges
      - vg_snarls_gfa // => vg snarls, with a .GFA file as input
      - clsd_sb // => clsd only handles unidirectional edges
      - sbSPQR_sb // => BubbleFinder, unidirectional edges
      - sbSPQR_snarls_gfa // => BubbleFinder, snarls mode
```

Per dataset (overrides the global value):

```yaml
- name: coli3682
  ...
  bench_programs:
    - clsd_sb
```

See the comments in `datasets.yaml` for the full list of programs and options.

### Stats / size‑distribution configuration

Block and 3‑connected component size distributions are controlled globally under `defaults.stats`:

```yaml
defaults:
  stats:
    enabled: true          # enable/disable BubbleFinder stats + hist plots
    kinds:
      - blocks             # 2‑connected blocks
      - spqr               # 3‑connected components (SPQR R‑node skeletons)
    threads: 4             # threads for BubbleFinder stats runs
```

- If `enabled: true`, the pipeline will:
  - run `BubbleFinder stats` on each enabled dataset,
  - produce per‑dataset text files with one size per line,
  - allow you to generate global histograms/boxplots.
- `kinds` controls which distributions are collected and plotted:
  - `blocks` → 2‑connected block size distributions.
  - `spqr` → 3‑connected component size distributions.

The actual plot files are:

- `results/plots/block_size_hist.png`
- `results/plots/spqr_size_hist.png`

They are **not** included in the `all` rule by default; you request them explicitly as Snakemake targets (see below).

---

## 4. Running the pipeline

### Dry‑run (no execution)

```bash
snakemake -n -p
```

### Full run (all enabled datasets)

```bash
snakemake all --use-conda -j 8
```

- `--use-conda`: required to create/use the environments defined in `config/*.yml`.
- `-j 8`: number of parallel jobs (adapt to your machine / cluster).

### Targeted examples

- Build and clean only the GFA of `coli3682`:

```bash
snakemake --use-conda -j 4 data/coli3682/coli3682.cleaned.gfa
```

- Run only the benchmarks and aggregation:

```bash
snakemake --use-conda -j 8 results/benchmarks.tsv
```

- Run only the block / 3‑connected component size histogram plots (and their prerequisites):

```bash
snakemake --use-conda -j 8 \
  results/plots/block_size_hist.png \
  results/plots/spqr_size_hist.png
```

If the intermediate stats files already exist (`results/stats/...`), only the plotting steps will be re‑run. Otherwise, Snakemake will first build all required graphs and run `BubbleFinder stats` to populate the size files.

---

## 5. Output structure

### `data/` directory

For each dataset `<name>` (e.g. `coli3682`):

- `data/<name>/<name>.*.gfa`  
  - Raw GFA (`...ggcat.fasta.gfa`, `...vg.gfa`, `...pggb.gfa`, etc.).
- `data/<name>/<name>.bluntified.gfa` (temporary).
- `data/<name>/<name>.cleaned.gfa`  
  - Cleaned GFA (no H‑lines, bluntified).
- `data/<name>/<name>.sb.cleaned.gfa`  
  - unidirectional graph representation version (if `ggcat_force_f` is enabled for this dataset).
- `data/<name>/<name>.sbspqr.sgraph`  
  - sgraph used by BubbleFinder/sbSPQR (unidirectional graph representation mode).
- `data/<name>/<name>.clsd.edgelist`  
  - Edgelist for `clsd`.

### `results/` directory

- `results/.prechecks.ok`  
  - Marker indicating that prechecks have run.
- `results/bench/`  
  - Per‑program, per‑dataset, per‑rep benchmark TSV files.  
  - Name: `<program>/<dataset>.t<threads>.rep<rep>.tsv`.
- `results/prog_out/`  
  - Raw program outputs (JSON, sgraphs, program‑specific logs).
- `results/logs/`  
  - Detailed logs per step (ggcat, vg, sgraph, bench, etc.).
- `results/stats/`  
  - `blocks/<dataset>_block_sizes.txt` (one 2‑connected block size per line).  
  - `spqr/<dataset>_spqr_sizes.txt` (one 3‑connected component size per line).
- `results/benchmarks.tsv`  
  - Aggregated benchmark table (time, memory, signature, etc.).
- `results/plots/`  
  - `time_by_dataset_program.png`
  - `rss_by_dataset_program.png`
  - `block_size_hist.png` (block size histogram + boxplots across datasets)
  - `spqr_size_hist.png` (3‑connected component size histogram + boxplots across datasets)
- `results/summary/reruns_planned.tsv`  
  - Information used for automatic reruns (timeouts, failures).

---

## 6. Customizing tool binaries

In `datasets.yaml`, under `defaults.tools`, you can point to pre‑installed binaries to avoid automatic cloning/building:

```yaml
defaults:
  tools:
    spqr_bin: /path/to/BubbleFinder
    get_blunted: /path/to/get_blunted
    lighter_bin: /path/to/lighter
    clsd_bin: /path/to/clsd
```

If these paths are not set, the Snakefile will:

- clone/build `BubbleFinder`, `clsd`, and `Lighter` under `build/`,
- download `get_blunted` into `bin/`.

---

## 7. Quick troubleshooting

- Error in prechecks (`prechecks`):  
  - Check `results/logs/*` (especially `results/logs/bench`, `results/logs/ggcat`, `results/logs/vg`).
  - Ensure `conda` or `mamba` is available.
- Frequent timeouts: 
  - Adjust `defaults.tools.timeout.seconds` in `datasets.yaml`.
- Problem with a specific dataset:
  - Run a single target, e.g.:  
    `snakemake --use-conda -j 4 data/coli3682/coli3682.cleaned.gfa`
- Bluntification issues:
  - A WARN about `get_blunted` means the pipeline falls back to "naive bluntify" (overlap fields forced to `*`).  
  - Install GetBlunted and set `defaults.tools.get_blunted` for more better behavior.