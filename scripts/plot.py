#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script summary:
# - Input: a "results" directory containing:
#   * prog_out/<program>/*.report.json (phase timing/memory reports, optional)
#   * bench/<program>/*.tsv (per-run benchmark key/value files, optional)
#   * benchmarks.tsv (per-run benchmark table, optional, overrides bench/)
# - Processing:
#   * Load and merge per-run metrics from TSV and JSON
#   * Infer run status (OK, TO, OOM, etc.) and aggregate per dataset/program
#   * Build multi-index tables for total wall time, total peak RSS,
#     and for phase metrics (I/O, BUILD, LOGIC)
# - Output (to --outdir, default "plots"):
#   * For each of 4 tables: .tsv, .html, .png (if dataframe_image is available), .tex
#   * Optional PNG border using Pillow or ImageMagick

import os, sys
os.environ.pop("PYTHONPATH", None)
_prefix = os.environ.get("CONDA_PREFIX", "")
if _prefix:
    sys.path = [p for p in sys.path if p.startswith(_prefix)]

import argparse
from pathlib import Path
import re
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import shutil

try:
    import dataframe_image as dfi
    HAS_DFI = True
except Exception:
    HAS_DFI = False

try:
    from PIL import Image, ImageOps
    HAS_PIL = True
except Exception:
    HAS_PIL = False

FRAME_PX = 28
WHITE_FRAME_PX = FRAME_PX

COL_SPEC = [
    ("SB",    "CSLD",      "clsd_sb",                False),
    ("SB",    "SPQR",      "sbSPQR_sb",              True),
    ("BiSB",  "BubbleGun", "BubbleGun_gfa",          False),
    ("BiSB",  "SPQR",      "sbSPQR_gfa",             True),
    ("Snarl", "vg_snarls", "vg_snarls_gfa",          True),
    ("Snarl", "SPQR",      "sbSPQR_snarls_gfa",      True),
]

FAMILY_COLORS = {
    "SB":    (0.85, 0.92, 0.98),
    "BiSB":  (0.98, 0.88, 0.88),
    "Snarl": (0.88, 0.96, 0.90),
}

DEFAULT_DATASET_RENAME = {
    "GCF_000012685.1_ASM1268v1": "BG1",
    "HG00733": "BG2",
    "Ultrabubble_dataset": "UB",
}

PHASES = ["I/O", "BUILD", "LOGIC"]

PHASE_PROGRAMS_DEFAULT = {"sbSPQR_sb", "sbSPQR_gfa", "sbSPQR_snarls_gfa"}


def _add_frame_to_png(png_path: Path, border_px: int, mode: str = "none") -> bool:
    if mode == "none" or border_px <= 0 or not png_path.exists():
        return False
    if HAS_PIL:
        try:
            with Image.open(png_path) as img:
                if mode == "transparent":
                    if img.mode != "RGBA":
                        img = img.convert("RGBA")
                    fill = (255, 255, 255, 0)
                else:
                    fill = (255, 255, 255, 255) if img.mode in ("RGBA", "LA") else (255, 255, 255)
                ImageOps.expand(img, border=border_px, fill=fill).save(png_path)
            print(f"[INFO] Frame {mode} added: {png_path}")
            return True
        except Exception as e:
            print(f"[WARN] Pillow failed on {png_path}: {e}")
    if mode in ("white", "transparent"):
        for cmd in (["magick", "convert"], ["convert"]):
            try:
                color = "none" if mode == "transparent" else "white"
                subprocess.run(
                    cmd
                    + [
                        str(png_path),
                        "-alpha",
                        "set",
                        "-bordercolor",
                        color,
                        "-border",
                        str(border_px),
                        str(png_path),
                    ],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                print(f"[INFO] Frame {mode} added via ImageMagick: {png_path}")
                return True
            except Exception:
                pass
    print(f"[WARN] Frame not added: {png_path}")
    return False


def _add_frame_in_dir(outdir: Path, border_px: int, mode: str):
    n = 0
    for p in outdir.glob("*.png"):
        if _add_frame_to_png(p, border_px, mode):
            n += 1
    if n:
        print(f"[INFO] Frame {mode} added on {n} PNG file(s) in {outdir}.")
    else:
        print(f"[INFO] No PNG processed for frame in {outdir}.")


def _rgb_to_hex(rgb):
    r, g, b = [int(255 * x) for x in rgb]
    return f"#{r:02x}{g:02x}{b:02x}"


def _lighten(rgb, amount=0.18):
    return tuple(min(1.0, x + (1 - x) * amount) for x in rgb)


def _build_header_table_styles_by_family_3level(columns: pd.MultiIndex):
    ncols = len(columns)
    fams = [c[0] for c in columns]
    tools = [c[1] for c in columns]
    subs = [c[2] for c in columns]

    styles = []
    styles.append(
        {
            "selector": "th.col_heading.level0",
            "props": [
                ("height", "0"),
                ("line-height", "0"),
                ("padding", "0"),
                ("border-top", "0"),
                ("border-bottom", "0"),
                ("font-size", "0"),
                ("color", "transparent"),
            ],
        }
    )
    styles.append(
        {
            "selector": "th.col_heading.level1",
            "props": [("border-top", "1.5px solid black")],
        }
    )
    styles.append(
        {
            "selector": "th.col_heading.level2",
            "props": [("border-bottom", "1.5px solid black")],
        }
    )
    styles.append(
        {
            "selector": "th.blank",
            "props": [("background-color", "white"), ("border", "0")],
        }
    )

    for j in range(ncols):
        fam = fams[j]
        base_rgb = FAMILY_COLORS.get(fam, (1, 1, 1))
        base_hex = _rgb_to_hex(base_rgb)
        lite_hex = _rgb_to_hex(_lighten(base_rgb, 0.18))

        styles.append(
            {
                "selector": f"th.col_heading.level1.col{j}",
                "props": [("background-color", base_hex)],
            }
        )
        styles.append(
            {
                "selector": f"th.col_heading.level2.col{j}",
                "props": [("background-color", lite_hex if subs[j] else base_hex)],
            }
        )

        left_family = (j == 0) or (fams[j] != fams[j - 1])
        if left_family:
            styles += [
                {
                    "selector": f"th.col_heading.level1.col{j}",
                    "props": [("border-left", "1.5px solid black")],
                },
                {
                    "selector": f"th.col_heading.level2.col{j}",
                    "props": [("border-left", "1.5px solid black")],
                },
            ]
        if j == ncols - 1:
            styles += [
                {
                    "selector": f"th.col_heading.level1.col{j}",
                    "props": [("border-right", "1.5px solid black")],
                },
                {
                    "selector": f"th.col_heading.level2.col{j}",
                    "props": [("border-right", "1.5px solid black")],
                },
            ]
        if j > 0 and fams[j] == fams[j - 1] and tools[j] != tools[j - 1]:
            styles += [
                {
                    "selector": f"th.col_heading.level1.col{j}",
                    "props": [("border-left", "1.5px solid black")],
                },
                {
                    "selector": f"th.col_heading.level2.col{j}",
                    "props": [("border-left", "1.5px solid black")],
                },
            ]

    return styles


def _shade_to_cells(val: str) -> str:
    try:
        s = str(val)
    except Exception:
        s = ""
    if s.startswith("TO"):
        return "background-color: #f2f2f2; color: black;"
    return ""


def export_table_png_and_tsv(
    df: pd.DataFrame,
    title: str,
    out_png: Path,
    out_tsv: Path,
    out_html: Path | None = None,
):
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t")

    header_styles = _build_header_table_styles_by_family_3level(df.columns)

    df = df.copy()
    df.index.name = None

    sty = df.style.set_caption(title).applymap(_shade_to_cells).set_table_styles(
        [
            {
                "selector": "caption",
                "props": [
                    ("caption-side", "top"),
                    ("text-align", "center"),
                    ("font-size", "14px"),
                    ("font-weight", "bold"),
                    ("margin-bottom", "8px"),
                ],
            },
            {
                "selector": "th.col_heading",
                "props": [("text-align", "center"), ("padding", "6px 8px")],
            },
            {
                "selector": "th.row_heading",
                "props": [
                    ("text-align", "left"),
                    ("padding", "6px 8px"),
                    ("border-top", "1px solid #444"),
                    ("border-bottom", "1px solid #444"),
                    ("border-left", "1px solid #444"),
                    ("border-right", "1.5px solid black"),
                ],
            },
            {
                "selector": "td",
                "props": [
                    ("text-align", "center"),
                    ("padding", "5px 8px"),
                    ("border", "1px solid #444"),
                ],
            },
            {
                "selector": "table",
                "props": [
                    ("border-collapse", "separate"),
                    ("border-spacing", "0px"),
                    ("border", "1.5px solid black"),
                ],
            },
        ]
        + header_styles
    )

    if out_html is not None:
        out_html.parent.mkdir(parents=True, exist_ok=True)
        out_html.write_text(sty.to_html(), encoding="utf-8")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    if HAS_DFI:
        try:
            dfi.export(sty, out_png)
        except Exception as e:
            print("[WARN] PNG export failed:", e)
            print(
                "       Styled HTML written to:",
                str(out_html) if out_html else "(none). Install/upgrade dataframe_image + kaleido.",
            )
    else:
        print(
            "[INFO] dataframe_image not available: PNG not generated (pip install dataframe_image kaleido)."
        )
        if out_html is None:
            html_fallback = out_png.with_suffix(".html")
            html_fallback.write_text(sty.to_html(), encoding="utf-8")
            print("       Styled HTML written to:", str(html_fallback))


def export_table_latex(
    df: pd.DataFrame,
    title: str,
    out_tex: Path,
    label: str | None = None,
    longtable: bool = False,
):
    def _latex_escape(s) -> str:
        s = "" if s is None else str(s)
        mapping = {
            "\\": r"\textbackslash{}",
            "&": r"\&",
            "%": r"\%",
            "$": r"\$",
            "#": r"\#",
            "_": r"\_",
            "{": r"\{",
            "}": r"\}",
            "~": r"\textasciitilde{}",
            "^": r"\textasciicircum{}",
        }
        return "".join(mapping.get(ch, ch) for ch in s)

    if not isinstance(df.columns, pd.MultiIndex) or df.columns.nlevels != 3:
        raise ValueError(
            "export_table_latex expects a DataFrame with columns as a 3-level MultiIndex."
        )

    cols = list(df.columns)
    n = len(cols)
    fams = [c[0] for c in cols]
    tools = [c[1] for c in cols]
    subs = [c[2] for c in cols]

    def _hex_nohash(rgb):
        return _rgb_to_hex(rgb).lstrip("#").upper()

    uniq_fams = []
    for f in fams:
        if f not in uniq_fams:
            uniq_fams.append(f)
    fam_hex = {
        f: _hex_nohash(FAMILY_COLORS.get(f, (1, 1, 1))) for f in uniq_fams
    }
    fam_hex_lite = {
        f: _hex_nohash(_lighten(FAMILY_COLORS.get(f, (1, 1, 1)), 0.18))
        for f in uniq_fams
    }

    colspec = "l " + " ".join(["c"] * n)
    groups = []
    j = 0
    while j < n:
        f, t = fams[j], tools[j]
        span = 1
        while j + span < n and (fams[j + span], tools[j + span]) == (f, t):
            span += 1
        groups.append((f, t, span))
        j += span

    lines = []
    lines.append("% Auto-generated LaTeX table")
    lines.append("% Requires: \\usepackage[table]{xcolor} and \\usepackage{booktabs}")
    for f in uniq_fams:
        lines.append(f"\\definecolor{{fam{f}}}{{HTML}}{{{fam_hex[f]}}}")
        lines.append(
            f"\\definecolor{{fam{f}Light}}{{HTML}}{{{fam_hex_lite[f]}}}"
        )

    begin_env = "\\begin{longtable}" if longtable else "\\begin{tabular}"
    end_env = "\\end{longtable}" if longtable else "\\end{tabular}"

    if not longtable:
        lines.append("\\begin{table}[ht]")
        lines.append("\\centering")
        if title:
            lines.append(f"\\caption{{{_latex_escape(title)}}}")
        if label:
            lines.append(f"\\label{{{_latex_escape(label)}}}")

    lines.append(f"{begin_env}{{{colspec}}}")
    lines.append("\\toprule")

    hdr1 = _latex_escape("Dataset")
    for f, t, span in groups:
        cell = f"\\cellcolor{{fam{f}}}\\textbf{{{_latex_escape(t)}}}"
        hdr1 += f" & \\multicolumn{{{span}}}{{c}}{{{cell}}}"
    hdr1 += r" \\"
    lines.append(hdr1)

    hdr2 = " "
    for j in range(n):
        f = fams[j]
        sub = subs[j]
        colname = f"fam{f}Light" if (sub and str(sub).strip()) else f"fam{f}"
        text = _latex_escape(sub) if (sub and str(sub).strip()) else " "
        hdr2 += f" & \\cellcolor{{{colname}}}{text}"
    hdr2 += r" \\"
    lines.append(hdr2)
    lines.append("\\midrule")

    for idx, row in df.iterrows():
        line = _latex_escape(idx)
        for _, val in enumerate(row.tolist()):
            s = "" if val is None else str(val)
            if s.startswith("TO"):
                line += (
                    f" & \\cellcolor[HTML]{{F2F2F2}} {_latex_escape(s)}"
                )
            else:
                line += f" & {_latex_escape(s)}"
        line += r" \\"
        lines.append(line)

    lines.append("\\bottomrule")
    lines.append(f"{end_env}")

    if not longtable:
        lines.append("\\end{table}")

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text("\n".join(lines), encoding="utf-8")
    print(f"[INFO] LaTeX written to: {out_tex}")


def parse_time_to_seconds(s: str) -> float:
    if s is None:
        return np.nan
    s = str(s).strip()
    if not s or s.lower() in ("na", "nan"):
        return np.nan
    parts = s.split(":")
    try:
        if len(parts) == 1:
            return float(parts[0])
        elif len(parts) == 2:
            m = int(parts[0])
            sec = float(parts[1])
            return m * 60 + sec
        elif len(parts) == 3:
            h = int(parts[0])
            m = int(parts[1])
            sec = float(parts[2])
            return h * 3600 + m * 60 + sec
        else:
            return np.nan
    except Exception:
        return np.nan


def fmt_seconds(s: float) -> str:
    if pd.isna(s):
        return "N/A"
    s = float(s)
    h = int(s // 3600)
    m = int((s % 3600) // 60)
    sec = int(round(s % 60))
    if sec == 60:
        sec = 0
        m += 1
    if m == 60:
        m = 0
        h += 1
    return f"{h}:{m:02d}:{sec:02d}" if h > 0 else f"{m}:{sec:02d}"


def fmt_mem_kb(kb: float) -> str:
    if pd.isna(kb):
        return "N/A"
    mib = kb / 1024.0
    if mib < 1024:
        return f"{mib:.1f} MiB"
    gib = mib / 1024.0
    return f"{gib:.2f} GiB"


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    if "replicate" not in df.columns and "rep" in df.columns:
        df = df.rename(columns={"rep": "replicate"})
    if "replicate" not in df.columns:
        df["replicate"] = 1
    return df


def _canon_metric_name(metric: str) -> str:
    s = "" if metric is None else str(metric)
    s = s.strip().rstrip(":").replace(" ", "").lower()
    return s


def _parse_bench_kv_file_robust(tsv_path: Path) -> dict:
    kv = {}
    try:
        with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
            for raw in f:
                line = raw.rstrip("\n")
                if not line.strip():
                    continue
                if "\t" in line:
                    parts = line.split("\t")
                    if parts[0].endswith(":"):
                        key = parts[0].rstrip(":").strip()
                        val = parts[1] if len(parts) > 1 else ""
                        kv[key] = val.strip()
                        continue
                    key_tokens, val_tokens, found = [], [], False
                    for tok in parts:
                        if not found:
                            key_tokens.append(tok)
                            if tok.endswith(":"):
                                found = True
                        else:
                            val_tokens.append(tok)
                    if found:
                        key = " ".join(t.rstrip(":") for t in key_tokens).strip()
                        val = "\t".join(val_tokens).strip()
                        kv[key] = val
                        continue
                if ":" in line:
                    k, v = line.split(":", 1)
                    kv[k.strip()] = v.strip()
    except FileNotFoundError:
        return {}
    return kv


def _fill_missing_from_paths(df: pd.DataFrame) -> pd.DataFrame:
    if "path" not in df.columns:
        return df
    df = df.copy()
    for col in ("elapsed_sec", "rss_kb", "timed_out", "exit_status", "timeout_s"):
        if col not in df.columns:
            df[col] = np.nan

    for i, row in df.iterrows():
        need = False
        if pd.isna(row.get("elapsed_sec")) or str(row.get("elapsed_sec")).strip() == "":
            need = True
        if pd.isna(row.get("rss_kb")) or str(row.get("rss_kb")).strip() == "":
            need = True
        if not need:
            continue
        p = row.get("path")
        if not isinstance(p, str) or not p:
            continue
        kv = _parse_bench_kv_file_robust(Path(p))
        el_raw = kv.get("Elapsed") or kv.get("Elapsed (wall clock) time")
        if el_raw is not None and (
            pd.isna(df.at[i, "elapsed_sec"])
            or str(df.at[i, "elapsed_sec"]).strip() == ""
        ):
            df.at[i, "elapsed_sec"] = parse_time_to_seconds(el_raw)
        rss_raw = kv.get("Maximum resident set size")
        if rss_raw is not None and (
            pd.isna(df.at[i, "rss_kb"])
            or str(df.at[i, "rss_kb"]).strip() == ""
        ):
            try:
                df.at[i, "rss_kb"] = float(str(rss_raw).strip())
            except Exception:
                pass
        to_raw = kv.get("Timed out")
        if to_raw is not None and (
            pd.isna(df.at[i, "timed_out"])
            or str(df.at[i, "timed_out"]).strip() == ""
        ):
            try:
                df.at[i, "timed_out"] = int(str(to_raw).strip())
            except Exception:
                pass
        rc_raw = kv.get("Exit status")
        if rc_raw is not None and (
            pd.isna(df.at[i, "exit_status"])
            or str(df.at[i, "exit_status"]).strip() == ""
        ):
            try:
                df.at[i, "exit_status"] = int(str(rc_raw).strip())
            except Exception:
                pass
        tmo_raw = kv.get("Timeout (s)")
        if (
            tmo_raw is not None
            and ("timeout_s" in df.columns)
            and (
                pd.isna(df.at[i, "timeout_s"])
                or str(df.at[i, "timeout_s"]).strip() == ""
            )
        ):
            try:
                df.at[i, "timeout_s"] = float(str(tmo_raw).strip())
            except Exception:
                pass
    return df


def load_runs(tsv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df = _normalize_columns(df)
    cols = set(df.columns)

    if {"dataset", "program", "replicate"}.issubset(cols) and (
        {"elapsed_sec", "rss_kb"} & cols
    ):
        df = _fill_missing_from_paths(df)

        def s_num(col):
            return (
                pd.to_numeric(df[col], errors="coerce")
                if col in df.columns
                else pd.Series([np.nan] * len(df), index=df.index)
            )

        def peak_kb(prefix):
            kb = None
            if f"{prefix}_peak_hwm_bytes" in df.columns:
                kb = s_num(f"{prefix}_peak_hwm_bytes") / 1024.0
            if f"{prefix}_peak_hwm_gib" in df.columns:
                tmp = s_num(f"{prefix}_peak_hwm_gib") * 1024.0 * 1024.0
                kb = tmp if kb is None else kb.fillna(tmp)
            return (
                kb
                if kb is not None
                else pd.Series([np.nan] * len(df), index=df.index)
            )

        runs = pd.DataFrame(
            {
                "dataset": df.get("dataset"),
                "program": df.get("program"),
                "replicate": df.get("replicate"),
                "elapsed_s": pd.to_numeric(
                    df.get("elapsed_sec"), errors="coerce"
                ),
                "max_rss_kb": pd.to_numeric(
                    df.get("rss_kb"), errors="coerce"
                ),
                "timed_out": pd.to_numeric(
                    df.get("timed_out"), errors="coerce"
                )
                .fillna(0)
                .astype(int),
                "timeout_s": pd.to_numeric(
                    df.get("timeout_s")
                    if "timeout_s" in df.columns
                    else df.get("prev_timeout"),
                    errors="coerce",
                ),
                "exit_status": pd.to_numeric(
                    df.get("exit_status"), errors="coerce"
                ),
                "io_s": s_num("sb_io_s"),
                "build_s": s_num("sb_build_s"),
                "logic_s": s_num("sb_logic_s"),
                "io_peak_kb": peak_kb("sb_io"),
                "build_peak_kb": peak_kb("sb_build"),
                "logic_peak_kb": peak_kb("sb_logic"),
            }
        )
        return runs

    if {"dataset", "program", "metric", "value"}.issubset(cols):
        if "replicate" not in df.columns and "rep" in df.columns:
            df = df.rename(columns={"rep": "replicate"})
        if "replicate" not in df.columns:
            df["replicate"] = 1

        rows = []
        for (dset, prog, rep), g in df.groupby(
            ["dataset", "program", "replicate"], dropna=False
        ):
            rec = {
                "dataset": dset,
                "program": prog,
                "replicate": rep,
                "elapsed_s": np.nan,
                "max_rss_kb": np.nan,
                "timed_out": 0,
                "timeout_s": np.nan,
                "exit_status": np.nan,
                "io_s": np.nan,
                "build_s": np.nan,
                "logic_s": np.nan,
                "io_peak_kb": np.nan,
                "build_peak_kb": np.nan,
                "logic_peak_kb": np.nan,
            }
            for _, r in g.iterrows():
                c = _canon_metric_name(r["metric"])
                val = r["value"]

                if c.startswith("elapsed"):
                    rec["elapsed_s"] = parse_time_to_seconds(val)
                    continue
                if c.startswith("maximumresidentsetsiz"):
                    try:
                        rec["max_rss_kb"] = float(val)
                    except Exception:
                        pass
                    continue
                if c == "timedout":
                    try:
                        rec["timed_out"] = int(val)
                    except Exception:
                        rec["timed_out"] = 0
                    continue
                if c == "timeout(s)":
                    try:
                        rec["timeout_s"] = float(val)
                    except Exception:
                        pass
                    continue
                if c == "exitstatus":
                    try:
                        rec["exit_status"] = int(val)
                    except Exception:
                        pass
                    continue

                if c in ("sb_io_s", "sb_build_s", "sb_logic_s"):
                    key = c.replace("sb_", "").replace("_s", "")
                    try:
                        rec[f"{key}_s"] = float(val)
                    except Exception:
                        pass
                    continue

                if c in (
                    "sb_io_peak_hwm_bytes",
                    "sb_build_peak_hwm_bytes",
                    "sb_logic_peak_hwm_bytes",
                ):
                    key = c.replace("sb_", "").replace("_peak_hwm_bytes", "")
                    try:
                        rec[f"{key}_peak_kb"] = float(val) / 1024.0
                    except Exception:
                        pass
                    continue

                if c in (
                    "sb_io_peak_hwm_gib",
                    "sb_build_peak_hwm_gib",
                    "sb_logic_peak_hwm_gib",
                ):
                    key = c.replace("sb_", "").replace("_peak_hwm_gib", "")
                    try:
                        rec[f"{key}_peak_kb"] = float(val) * 1024.0 * 1024.0
                    except Exception:
                        pass
                    continue

            rows.append(rec)
        return pd.DataFrame(rows)

    raise ValueError(
        f"File {tsv_path} not recognized. Columns seen: {sorted(cols)}. "
        "Expected either (dataset, program, replicate, metric, value) or "
        "(dataset, program, rep/replicate, elapsed_sec, rss_kb, ...)."
    )


def _safe_float(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


def _parse_report_name_to_dataset_replicate(
    stem: str,
) -> tuple[str, str] | None:
    if stem.endswith(".report"):
        stem = stem[:-7]
    m = re.match(r"^(?P<dataset>.+)\.rep(?P<rep>\d+)$", stem)
    if not m:
        return None
    return m.group("dataset"), m.group("rep")


def _parse_bench_name_to_dataset_replicate(
    stem: str,
) -> tuple[str, str] | None:
    m = re.match(r"^(?P<dataset>.+)\.rep(?P<rep>\d+)$", stem)
    if not m:
        return None
    return m.group("dataset"), m.group("rep")


def _marks_index_by_label(marks):
    d = {}
    for m in marks or []:
        lbl = m.get("label")
        if isinstance(lbl, str):
            d[lbl] = m
    return d


def _sum_seconds(d, labels):
    vals = []
    for l in labels:
        m = d.get(l)
        if m is not None:
            vals.append(_safe_float(m.get("seconds"), 0.0))
    return sum(vals) if vals else np.nan


def _max_rss_kb(d, labels):
    maxi = -1
    for l in labels:
        m = d.get(l)
        if m is not None:
            v = m.get("rss_max_bytes")
            if isinstance(v, (int, float)):
                maxi = max(maxi, int(v))
    return (maxi / 1024.0) if maxi >= 0 else np.nan


def _incremental_peaks(peaks):
    prev = 0.0
    inc = []
    for p in peaks:
        if p is None:
            inc.append(np.nan)
            continue
        try:
            pf = float(p)
        except Exception:
            inc.append(np.nan)
            continue
        if np.isnan(pf):
            inc.append(np.nan)
            continue
        new_peak = max(prev, pf)
        inc.append(new_peak - prev)
        prev = new_peak
    return inc


def _extract_phases_from_sb_json(
    obj: dict,
) -> tuple[float, float, float, float, float, float]:
    marks = obj.get("marks") or []
    md = _marks_index_by_label(marks)

    io_labels = ["io/read_graph"]
    io_s = _sum_seconds(md, io_labels)
    io_peak_abs = _max_rss_kb(md, io_labels)

    build_children = [
        "sb/phase/GccBuildParallel",
        "sb/phase/BCtrees",
        "sb/phase/BlockDataBuildAll",
    ]
    build_labels = ["sb/io/cc_buckets"] + build_children

    has_any_build = any(l in md for l in build_labels)
    if has_any_build:
        build_s = _sum_seconds(md, build_labels)
        build_peak_abs = _max_rss_kb(md, build_labels)
    elif "sb/build/all" in md:
        build_s = _safe_float(md["sb/build/all"].get("seconds"))
        v = md["sb/build/all"].get("rss_max_bytes")
        build_peak_abs = (
            float(v) / 1024.0 if isinstance(v, (int, float)) else np.nan
        )
    else:
        build_s = np.nan
        build_peak_abs = np.nan

    if "sb/logic/all" in md:
        logic_core = _safe_float(md["sb/logic/all"].get("seconds"))
        extra_io = _sum_seconds(md, ["io/write_output"])
        if pd.isna(logic_core):
            logic_s = extra_io
        elif pd.isna(extra_io):
            logic_s = logic_core
        else:
            logic_s = logic_core + extra_io

        logic_peak_abs = _max_rss_kb(md, ["sb/logic/all", "io/write_output"])
    else:
        logic_labels = [
            "sb/findMini",
            "sb/phase/SolveBlocks",
            "io/write_output",
        ]
        logic_s = _sum_seconds(md, logic_labels)
        logic_peak_abs = _max_rss_kb(md, logic_labels)

    io_peak_kb, build_peak_kb, logic_peak_kb = _incremental_peaks(
        [io_peak_abs, build_peak_abs, logic_peak_abs]
    )

    return io_s, build_s, logic_s, io_peak_kb, build_peak_kb, logic_peak_kb


def _extract_phases_from_snarl_json(
    obj: dict,
) -> tuple[float, float, float, float, float, float]:
    marks = obj.get("marks") or []
    md = _marks_index_by_label(marks)

    io_labels = ["io/read_graph"]
    io_s = _sum_seconds(md, io_labels)
    io_peak_abs = _max_rss_kb(md, io_labels)

    build_labels = [
        "sn/phase/ComputeCC",
        "sn/phase/BucketNodes",
        "sn/phase/BucketEdges",
        "sn/worker_component/gcc_rebuild",
        "sn/worker_bcTree/build",
        "sn/phase/block_SPQR_build",
    ]
    build_s = _sum_seconds(md, build_labels)
    build_peak_abs = _max_rss_kb(md, build_labels)

    logic_labels = [
        "sn/phase/tips_cuts",
        "sn/phase/block_SPQR_solve",
        "io/write_output",
    ]
    logic_s = _sum_seconds(md, logic_labels)
    logic_peak_abs = _max_rss_kb(md, logic_labels)

    io_peak_kb, build_peak_kb, logic_peak_kb = _incremental_peaks(
        [io_peak_abs, build_peak_abs, logic_peak_abs]
    )

    return io_s, build_s, logic_s, io_peak_kb, build_peak_kb, logic_peak_kb


def _extract_total_from_json(obj: dict) -> tuple[float, float]:
    summary = obj.get("summary") or {}
    tot = summary.get("total") or {}
    elapsed_s = _safe_float(tot.get("seconds"))

    marks = obj.get("marks") or []
    maxi = -1
    for m in marks:
        v = m.get("rss_max_bytes")
        if isinstance(v, (int, float)):
            maxi = max(maxi, int(v))
    if maxi >= 0:
        max_rss_kb = maxi / 1024.0
    else:
        peak_hwm = tot.get("peak_hwm_bytes")
        max_rss_kb = (
            float(peak_hwm) / 1024.0 if isinstance(peak_hwm, (int, float)) else np.nan
        )

    return elapsed_s, max_rss_kb


def load_phase_metrics_from_results(results_root: Path) -> pd.DataFrame:
    prog_out = results_root / "prog_out"
    if not prog_out.is_dir():
        return pd.DataFrame(
            columns=[
                "dataset",
                "program",
                "replicate",
                "elapsed_s",
                "max_rss_kb",
                "timed_out",
                "timeout_s",
                "exit_status",
                "io_s",
                "build_s",
                "logic_s",
                "io_peak_kb",
                "build_peak_kb",
                "logic_peak_kb",
            ]
        )

    recs = []
    for prog_dir in prog_out.iterdir():
        if not prog_dir.is_dir():
            continue
        program = prog_dir.name
        if program not in {
            "sbSPQR_gfa",
            "sbSPQR_sb",
            "sbSPQR_snarls_gfa",
        }:
            continue

        for p in prog_dir.glob("*.report.json"):
            parsed = _parse_report_name_to_dataset_replicate(p.stem)
            if not parsed:
                continue
            dataset, rep = parsed

            try:
                with open(p, "r", encoding="utf-8") as fh:
                    obj = json.load(fh)
            except Exception:
                continue

            if program == "sbSPQR_snarls_gfa":
                (
                    io_s,
                    build_s,
                    logic_s,
                    io_peak_kb,
                    build_peak_kb,
                    logic_peak_kb,
                ) = _extract_phases_from_snarl_json(obj)
            else:
                (
                    io_s,
                    build_s,
                    logic_s,
                    io_peak_kb,
                    build_peak_kb,
                    logic_peak_kb,
                ) = _extract_phases_from_sb_json(obj)

            elapsed_s, max_rss_kb = _extract_total_from_json(obj)

            recs.append(
                {
                    "dataset": dataset,
                    "program": program,
                    "replicate": str(rep),
                    "elapsed_s": elapsed_s,
                    "max_rss_kb": max_rss_kb,
                    "timed_out": 0,
                    "timeout_s": np.nan,
                    "exit_status": 0,
                    "io_s": io_s,
                    "build_s": build_s,
                    "logic_s": logic_s,
                    "io_peak_kb": io_peak_kb,
                    "build_peak_kb": build_peak_kb,
                    "logic_peak_kb": logic_peak_kb,
                }
            )

    if not recs:
        return pd.DataFrame(
            columns=[
                "dataset",
                "program",
                "replicate",
                "elapsed_s",
                "max_rss_kb",
                "timed_out",
                "timeout_s",
                "exit_status",
                "io_s",
                "build_s",
                "logic_s",
                "io_peak_kb",
                "build_peak_kb",
                "logic_peak_kb",
            ]
        )
    df = pd.DataFrame(recs)
    df["replicate"] = df["replicate"].astype(str)
    return df


def load_runs_from_bench_dir(results_root: Path) -> pd.DataFrame:
    bench_root = results_root / "bench"
    cols = [
        "dataset",
        "program",
        "replicate",
        "elapsed_s",
        "max_rss_kb",
        "timed_out",
        "timeout_s",
        "exit_status",
        "io_s",
        "build_s",
        "logic_s",
        "io_peak_kb",
        "build_peak_kb",
        "logic_peak_kb",
    ]
    if not bench_root.is_dir():
        return pd.DataFrame(columns=cols)

    recs = []
    for prog_dir in bench_root.iterdir():
        if not prog_dir.is_dir():
            continue
        program = prog_dir.name
        for tsv in prog_dir.glob("*.tsv"):
            parsed = _parse_bench_name_to_dataset_replicate(tsv.stem)
            if not parsed:
                continue
            dataset, rep = parsed
            kv = _parse_bench_kv_file_robust(tsv)

            el_raw = kv.get("Elapsed") or kv.get("Elapsed (wall clock) time")
            elapsed_s = parse_time_to_seconds(el_raw)

            rss_raw = kv.get("Maximum resident set size")
            try:
                max_rss_kb = (
                    float(str(rss_raw))
                    if rss_raw is not None and str(rss_raw).strip() != ""
                    else np.nan
                )
            except Exception:
                max_rss_kb = np.nan

            to_raw = kv.get("Timed out")
            try:
                timed_out = (
                    int(str(to_raw))
                    if to_raw is not None and str(to_raw).strip() != ""
                    else 0
                )
            except Exception:
                timed_out = 0

            tmo_raw = kv.get("Timeout (s)")
            try:
                timeout_s = (
                    float(str(tmo_raw))
                    if tmo_raw is not None and str(tmo_raw).strip() != ""
                    else np.nan
                )
            except Exception:
                timeout_s = np.nan

            rc_raw = kv.get("Exit status")
            try:
                exit_status = (
                    int(str(rc_raw))
                    if rc_raw is not None and str(rc_raw).strip() != ""
                    else (0 if timed_out == 0 and pd.notna(elapsed_s) else np.nan)
                )
            except Exception:
                exit_status = (
                    0 if timed_out == 0 and pd.notna(elapsed_s) else np.nan
                )

            recs.append(
                {
                    "dataset": dataset,
                    "program": program,
                    "replicate": str(rep),
                    "elapsed_s": elapsed_s,
                    "max_rss_kb": max_rss_kb,
                    "timed_out": timed_out,
                    "timeout_s": timeout_s,
                    "exit_status": exit_status,
                    "io_s": np.nan,
                    "build_s": np.nan,
                    "logic_s": np.nan,
                    "io_peak_kb": np.nan,
                    "build_peak_kb": np.nan,
                    "logic_peak_kb": np.nan,
                }
            )

    if not recs:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(recs)


def merge_runs_with_phase_json(
    runs: pd.DataFrame, phases: pd.DataFrame
) -> pd.DataFrame:
    if runs is None or runs.empty:
        return phases.copy() if phases is not None else pd.DataFrame()
    if phases is None or phases.empty:
        return runs.copy()

    left = runs.copy()
    left["replicate"] = left["replicate"].astype(str)
    right = phases.copy()
    right["replicate"] = right["replicate"].astype(str)

    merged = left.merge(
        right,
        on=["dataset", "program", "replicate"],
        how="outer",
        suffixes=("", "_json"),
    )

    for col in [
        "io_s",
        "build_s",
        "logic_s",
        "io_peak_kb",
        "build_peak_kb",
        "logic_peak_kb",
    ]:
        cj = col + "_json"
        if cj in merged.columns:
            if col in merged.columns:
                merged[col] = merged[cj].combine_first(merged[col])
                merged.drop(columns=[cj], inplace=True)
            else:
                merged.rename(columns={cj: col}, inplace=True)

    for col in [
        "elapsed_s",
        "max_rss_kb",
        "exit_status",
        "timed_out",
        "timeout_s",
    ]:
        cj = col + "_json"
        if cj in merged.columns:
            if col in merged.columns:
                merged[col] = merged[col].combine_first(merged[cj])
                merged.drop(columns=[cj], inplace=True)
            else:
                merged.rename(columns={cj: col}, inplace=True)

    return merged


EPS = 0.5


def classify_run(
    exit_status: float | int | None,
    timed_out: int | float | None,
    elapsed_s: float | None,
    timeout_s: float | None,
) -> str:
    rc = np.nan if exit_status is None else float(exit_status)
    to = 0 if timed_out is None or str(timed_out) == "nan" else int(timed_out)
    if not pd.isna(rc):
        irc = int(rc)
        if irc == 0:
            return "OK"
        if irc == 124:
            return "TO"
        if irc == 137:
            return "OOM"
        if irc == 139:
            return "SEGV"
        if irc == 134:
            return "ABRT"
        if irc == 143:
            return "TERM"
        if irc >= 128:
            return "KILL"
        return "ERR"
    if to == 1:
        return "TO"
    if (
        timeout_s is not None
        and not pd.isna(timeout_s)
        and elapsed_s is not None
        and not pd.isna(elapsed_s)
        and float(elapsed_s) >= float(timeout_s) - EPS
    ):
        return "TO"
    return "UNKNOWN"


def aggregate_runs(runs: pd.DataFrame) -> pd.DataFrame:
    if runs.empty:
        return pd.DataFrame(
            columns=[
                "dataset",
                "program",
                "elapsed_ok_mean_s",
                "rss_ok_mean_kb",
                "io_ok_mean_s",
                "build_ok_mean_s",
                "logic_ok_mean_s",
                "io_peak_ok_mean_kb",
                "build_peak_ok_mean_kb",
                "logic_peak_ok_mean_kb",
                "n_rep",
                "n_ok",
                "n_to",
                "n_oom",
                "n_segv",
                "n_abrt",
                "n_term",
                "n_kill",
                "n_err_other",
                "n_unknown",
                "timeout_common_s",
            ]
        )

    runs = runs.copy()
    runs["timed_out"] = (
        pd.to_numeric(runs.get("timed_out", 0), errors="coerce")
        .fillna(0)
        .astype(int)
    )
    runs["elapsed_s"] = pd.to_numeric(
        runs.get("elapsed_s", np.nan), errors="coerce"
    )
    runs["max_rss_kb"] = pd.to_numeric(
        runs.get("max_rss_kb", np.nan), errors="coerce"
    )
    runs["timeout_s"] = pd.to_numeric(
        runs.get("timeout_s", np.nan), errors="coerce"
    )
    runs["exit_status"] = pd.to_numeric(
        runs.get("exit_status", np.nan), errors="coerce"
    )
    for c in [
        "io_s",
        "build_s",
        "logic_s",
        "io_peak_kb",
        "build_peak_kb",
        "logic_peak_kb",
    ]:
        runs[c] = pd.to_numeric(runs.get(c, np.nan), errors="coerce")

    runs["category"] = [
        classify_run(
            row["exit_status"],
            row["timed_out"],
            row["elapsed_s"],
            row["timeout_s"],
        )
        for _, row in runs.iterrows()
    ]

    agg = []
    for (dset, prog), g in runs.groupby(["dataset", "program"], dropna=False):
        counts = g["category"].value_counts()
        n_ok, n_to = int(counts.get("OK", 0)), int(counts.get("TO", 0))
        n_oom = int(counts.get("OOM", 0))
        n_segv = int(counts.get("SEGV", 0))
        n_abrt = int(counts.get("ABRT", 0))
        n_term = int(counts.get("TERM", 0))
        n_kill = int(counts.get("KILL", 0))
        n_err = int(counts.get("ERR", 0))
        n_unk = int(counts.get("UNKNOWN", 0))

        ok_mask = g["category"] == "OK"

        if ok_mask.any():
            elapsed_ok_mean_s = g.loc[ok_mask, "elapsed_s"].mean()
            io_ok_mean_s = g.loc[ok_mask, "io_s"].mean()
            build_ok_mean_s = g.loc[ok_mask, "build_s"].mean()
            logic_ok_mean_s = g.loc[ok_mask, "logic_s"].mean()
        else:
            elapsed_ok_mean_s = io_ok_mean_s = build_ok_mean_s = (
                logic_ok_mean_s
            ) = np.nan

        if ok_mask.any() and g.loc[ok_mask, "max_rss_kb"].notna().any():
            rss_ok_mean_kb = g.loc[ok_mask, "max_rss_kb"].mean()
        elif g["max_rss_kb"].notna().any():
            rss_ok_mean_kb = g["max_rss_kb"].mean()
        else:
            rss_ok_mean_kb = np.nan

        def mean_phase_peak(col):
            if ok_mask.any() and g.loc[ok_mask, col].notna().any():
                return g.loc[ok_mask, col].mean()
            if g[col].notna().any():
                return g[col].mean()
            return np.nan

        io_peak_ok_mean_kb = mean_phase_peak("io_peak_kb")
        build_peak_ok_mean_kb = mean_phase_peak("build_peak_kb")
        logic_peak_ok_mean_kb = mean_phase_peak("logic_peak_kb")

        tmo_candidates = (
            g.loc[g["category"] == "TO", "timeout_s"]
            .dropna()
            .astype(float)
        )
        timeout_common_s = (
            tmo_candidates.median()
            if len(tmo_candidates) > 0
            else g["timeout_s"].dropna().astype(float).median()
            if g["timeout_s"].notna().any()
            else np.nan
        )

        agg.append(
            {
                "dataset": dset,
                "program": prog,
                "elapsed_ok_mean_s": elapsed_ok_mean_s,
                "rss_ok_mean_kb": rss_ok_mean_kb,
                "io_ok_mean_s": io_ok_mean_s,
                "build_ok_mean_s": build_ok_mean_s,
                "logic_ok_mean_s": logic_ok_mean_s,
                "io_peak_ok_mean_kb": io_peak_ok_mean_kb,
                "build_peak_ok_mean_kb": build_peak_ok_mean_kb,
                "logic_peak_ok_mean_kb": logic_peak_ok_mean_kb,
                "n_rep": len(g),
                "n_ok": n_ok,
                "n_to": n_to,
                "n_oom": n_oom,
                "n_segv": n_segv,
                "n_abrt": n_abrt,
                "n_term": n_term,
                "n_kill": n_kill,
                "n_err_other": n_err,
                "n_unknown": n_unk,
                "timeout_common_s": timeout_common_s,
            }
        )

    return pd.DataFrame(agg)


_thr_patterns = [
    re.compile(r"(?i)[._-]thr(\d+)$"),
    re.compile(r"(?i)[._-]t(\d+)$"),
    re.compile(r"(?i)[._-]threads?(\d+)$"),
]


def split_dataset_base_threads(name: str) -> tuple[str, int | None]:
    if not isinstance(name, str):
        return str(name), None
    s = name.strip()
    for pat in _thr_patterns:
        m = pat.search(s)
        if m:
            base = s[: m.start()].rstrip("._-")
            try:
                t = int(m.group(1))
            except Exception:
                t = None
            return base, t
    return s, None


def _build_fail_annotation(row) -> str:
    parts = []
    if int(row.get("n_to", 0)) > 0:
        parts.append(f"{int(row.get('n_to', 0))} TO")
    if int(row.get("n_oom", 0)) > 0:
        parts.append(f"{int(row.get('n_oom', 0))} OOM")
    if int(row.get("n_segv", 0)) > 0:
        parts.append(f"{int(row.get('n_segv', 0))} SEGV")
    if int(row.get("n_abrt", 0)) > 0:
        parts.append(f"{int(row.get('n_abrt', 0))} ABRT")
    if int(row.get("n_term", 0)) > 0:
        parts.append(f"{int(row.get('n_term', 0))} TERM")
    if int(row.get("n_kill", 0)) > 0:
        parts.append(f"{int(row.get('n_kill', 0))} KILL")
    if int(row.get("n_err_other", 0)) > 0:
        parts.append(f"{int(row.get('n_err_other', 0))} ERR")
    if int(row.get("n_unknown", 0)) > 0:
        parts.append(f"{int(row.get('n_unknown', 0))} ?")
    return ", ".join(parts)


def _format_cells_from_row(row, for_time=True):
    n_rep = int(row.get("n_rep", 0) or 0)
    n_ok = int(row.get("n_ok", 0) or 0)
    n_to = int(row.get("n_to", 0) or 0)
    if n_ok > 0 or (not for_time and not pd.isna(row.get("rss_ok_mean_kb"))):
        base = (
            fmt_seconds(row["elapsed_ok_mean_s"])
            if for_time
            else fmt_mem_kb(row["rss_ok_mean_kb"])
        )
        n_fail = max(0, n_rep - n_ok)
        if n_fail > 0 and n_rep > 0:
            counts = {
                "TO": int(row.get("n_to", 0) or 0),
                "OOM": int(row.get("n_oom", 0) or 0),
                "SEGV": int(row.get("n_segv", 0) or 0),
                "ABRT": int(row.get("n_abrt", 0) or 0),
                "TERM": int(row.get("n_term", 0) or 0),
                "KILL": int(row.get("n_kill", 0) or 0),
                "ERR": int(row.get("n_err_other", 0) or 0),
                "?": int(row.get("n_unknown", 0) or 0),
            }
            cats = [k for k, v in counts.items() if v > 0]
            label = cats[0] if len(cats) == 1 else "FAIL"
            return f"{base} ({n_fail}/{n_rep} {label})"
        return base
    if n_rep > 0 and n_to == n_rep:
        to_str = (
            fmt_seconds(row.get("timeout_common_s"))
            if pd.notna(row.get("timeout_common_s"))
            else "N/A"
        )
        return f"TO ({to_str})" if for_time else "TO"
    annot = _build_fail_annotation(row)
    return f"FAIL ({annot})" if annot else "FAIL"


def _format_phase_cell_value(row, val, for_time=True):
    n_rep = int(row.get("n_rep", 0) or 0)
    n_ok = int(row.get("n_ok", 0) or 0)
    n_to = int(row.get("n_to", 0) or 0)
    if n_ok > 0 or (not for_time and pd.notna(val)):
        base = (
            fmt_seconds(val)
            if (for_time and pd.notna(val))
            else (fmt_mem_kb(val) if (not for_time and pd.notna(val)) else "N/A")
        )
        n_fail = max(0, n_rep - n_ok)
        if n_fail > 0 and n_rep > 0:
            counts = {
                "TO": int(row.get("n_to", 0) or 0),
                "OOM": int(row.get("n_oom", 0) or 0),
                "SEGV": int(row.get("n_segv", 0) or 0),
                "ABRT": int(row.get("n_abrt", 0) or 0),
                "TERM": int(row.get("n_term", 0) or 0),
                "KILL": int(row.get("n_kill", 0) or 0),
                "ERR": int(row.get("n_err_other", 0) or 0),
                "?": int(row.get("n_unknown", 0) or 0),
            }
            cats = [k for k, v in counts.items() if v > 0]
            label = cats[0] if len(cats) == 1 else "FAIL"
            return f"{base} ({n_fail}/{n_rep} {label})"
        return base
    if n_rep > 0 and n_to == n_rep:
        return (
            f"TO ({fmt_seconds(row.get('timeout_common_s'))})"
            if for_time
            else "TO"
        )
    annot = _build_fail_annotation(row)
    return f"FAIL ({annot})" if annot else "FAIL"


def build_multiheader_tables(
    agg_df: pd.DataFrame,
    threads: list[int],
    dataset_rename: dict | None = None,
    dataset_order: list[str] | None = None,
):
    if agg_df.empty:
        cols = []
        for fam, tool, _, has_thr in COL_SPEC:
            if has_thr:
                for t in threads:
                    cols.append((fam, tool, f"t={t}"))
            else:
                cols.append((fam, tool, ""))
        mi = pd.MultiIndex.from_tuples(cols, names=[None, None, None])
        empty = pd.DataFrame(columns=mi)
        empty.index.name = None
        return empty, empty.copy()

    idx = {
        (str(r["dataset"]), str(r["program"])): r.to_dict()
        for _, r in agg_df.iterrows()
    }

    base_to_threads = {}
    for d in agg_df["dataset"].unique().tolist():
        base, t = split_dataset_base_threads(str(d))
        base_to_threads.setdefault(base, set()).add(t)

    bases = sorted(base_to_threads.keys())
    pretty = {
        b: (dataset_rename.get(b, b) if dataset_rename else b) for b in bases
    }
    if dataset_order:
        inv = {pretty[b]: b for b in bases}
        ordered = []
        for lbl in dataset_order:
            if lbl in inv:
                ordered.append(inv[lbl])
            elif lbl in bases:
                ordered.append(lbl)
        for b in bases:
            if b not in ordered:
                ordered.append(b)
        bases = ordered

    cols = []
    for fam, tool, _, has_thr in COL_SPEC:
        if has_thr:
            for t in threads:
                cols.append((fam, tool, f"t={t}"))
        else:
            cols.append((fam, tool, ""))

    def get_cell(dataset_name, program, for_time=True):
        r = idx.get((dataset_name, program))
        if r is None:
            return "N/A"
        return _format_cells_from_row(r, for_time)

    def find_candidate(base, t):
        for d in agg_df["dataset"].unique().tolist():
            b2, t2 = split_dataset_base_threads(str(d))
            if b2 == base and t2 == t:
                return d
        return None

    def pick_nonthreaded_candidate(base):
        ts = base_to_threads.get(base, set())
        if None in ts:
            return find_candidate(base, None)
        if 1 in ts:
            return find_candidate(base, 1)
        numeric_ts = sorted([t for t in ts if t is not None])
        return find_candidate(base, numeric_ts[0]) if numeric_ts else None

    time_rows, mem_rows = [], []
    for base in bases:
        t_row, m_row = [], []
        for fam, tool, prog, has_thr in COL_SPEC:
            if has_thr:
                ts = base_to_threads.get(base, set())
                for t in threads:
                    candidate = None
                    if t in ts:
                        candidate = find_candidate(base, t)
                    elif (None in ts) and (t == threads[0]):
                        candidate = find_candidate(base, None)
                    if candidate is None:
                        t_row.append("N/A")
                        m_row.append("N/A")
                    else:
                        t_row.append(get_cell(candidate, prog, True))
                        m_row.append(get_cell(candidate, prog, False))
            else:
                candidate = pick_nonthreaded_candidate(base)
                if candidate is None:
                    t_row.append("N/A")
                    m_row.append("N/A")
                else:
                    t_row.append(get_cell(candidate, prog, True))
                    m_row.append(get_cell(candidate, prog, False))
        time_rows.append(t_row)
        mem_rows.append(m_row)

    mi = pd.MultiIndex.from_tuples(cols, names=[None, None, None])
    return (
        pd.DataFrame(time_rows, index=[pretty[b] for b in bases], columns=mi),
        pd.DataFrame(mem_rows, index=[pretty[b] for b in bases], columns=mi),
    )


def build_multiheader_phase_tables(
    agg_df: pd.DataFrame,
    threads: list[int],
    dataset_rename: dict | None = None,
    dataset_order: list[str] | None = None,
    phase_programs: set[str] | None = None,
):
    if agg_df.empty:
        cols = []
        for fam, tool, prog, has_thr in COL_SPEC:
            if (phase_programs is not None) and (prog not in phase_programs):
                continue
            if has_thr:
                for t in threads:
                    for ph in PHASES:
                        cols.append((fam, tool, f"t={t} {ph}"))
            else:
                for ph in PHASES:
                    cols.append((fam, tool, ph))
        mi = pd.MultiIndex.from_tuples(cols, names=[None, None, None])
        empty = pd.DataFrame(columns=mi)
        empty.index.name = None
        return empty, empty.copy()

    idx = {
        (str(r["dataset"]), str(r["program"])): r.to_dict()
        for _, r in agg_df.iterrows()
    }

    base_to_threads = {}
    for d in agg_df["dataset"].unique().tolist():
        base, t = split_dataset_base_threads(str(d))
        base_to_threads.setdefault(base, set()).add(t)

    bases = sorted(base_to_threads.keys())
    pretty = {
        b: (dataset_rename.get(b, b) if dataset_rename else b) for b in bases
    }
    if dataset_order:
        inv = {pretty[b]: b for b in bases}
        ordered = []
        for lbl in dataset_order:
            if lbl in inv:
                ordered.append(inv[lbl])
            elif lbl in bases:
                ordered.append(lbl)
        for b in bases:
            if b not in ordered:
                ordered.append(b)
        bases = ordered

    cols = []
    for fam, tool, prog, has_thr in COL_SPEC:
        if (phase_programs is not None) and (prog not in phase_programs):
            continue
        if has_thr:
            for t in threads:
                for ph in PHASES:
                    cols.append((fam, tool, f"t={t} {ph}"))
        else:
            for ph in PHASES:
                cols.append((fam, tool, ph))

    def find_candidate(base, t):
        for d in agg_df["dataset"].unique().tolist():
            b2, t2 = split_dataset_base_threads(str(d))
            if b2 == base and t2 == t:
                return d
        return None

    phase_key = {"I/O": "io", "BUILD": "build", "LOGIC": "logic"}

    def get_phase_cell(dataset_name, program, phase, for_time=True):
        r = idx.get((dataset_name, program))
        if r is None:
            return "N/A"
        key = phase_key[phase]
        val = (
            r.get(f"{key}_ok_mean_s")
            if for_time
            else r.get(f"{key}_peak_ok_mean_kb")
        )
        return _format_phase_cell_value(r, val, for_time)

    time_rows, mem_rows = [], []
    for base in bases:
        t_row, m_row = [], []
        for fam, tool, prog, has_thr in COL_SPEC:
            if (phase_programs is not None) and (prog not in phase_programs):
                continue
            if has_thr:
                ts = base_to_threads.get(base, set())
                for t in threads:
                    candidate = None
                    if t in ts:
                        candidate = find_candidate(base, t)
                    elif (None in ts) and (t == threads[0]):
                        candidate = find_candidate(base, None)
                    for ph in PHASES:
                        if candidate is None:
                            t_row.append("N/A")
                            m_row.append("N/A")
                        else:
                            t_row.append(
                                get_phase_cell(candidate, prog, ph, True)
                            )
                            m_row.append(
                                get_phase_cell(candidate, prog, ph, False)
                            )
            else:
                has_base = None in base_to_threads.get(base, set())
                for ph in PHASES:
                    if not has_base:
                        t_row.append("N/A")
                        m_row.append("N/A")
                    else:
                        t_row.append(
                            get_phase_cell(base, prog, ph, True)
                        )
                        m_row.append(
                            get_phase_cell(base, prog, ph, False)
                        )
        time_rows.append(t_row)
        mem_rows.append(m_row)

    mi = pd.MultiIndex.from_tuples([c for c in cols], names=[None, None, None])
    return (
        pd.DataFrame(time_rows, index=[pretty[b] for b in bases], columns=mi),
        pd.DataFrame(mem_rows, index=[pretty[b] for b in bases], columns=mi),
    )


def infer_phase_programs_from_data(agg_df: pd.DataFrame) -> set[str]:
    progs = set()
    if agg_df.empty:
        return progs
    phase_cols = [
        "io_ok_mean_s",
        "build_ok_mean_s",
        "logic_ok_mean_s",
        "io_peak_ok_mean_kb",
        "build_peak_ok_mean_kb",
        "logic_peak_ok_mean_kb",
    ]
    for prog, g in agg_df.groupby("program"):
        if g[phase_cols].notna().any().any():
            progs.add(prog)
    return progs


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Builds time and memory tables from a 'results' directory "
            "containing prog_out/ (JSON) and optionally bench/ (TSV). "
            "JSON phase metrics have priority over TSV."
        )
    )
    ap.add_argument(
        "results",
        type=Path,
        help=(
            "Results directory (containing 'prog_out/'). "
            "Optionally 'bench/' and 'benchmarks.tsv'."
        ),
    )
    ap.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("plots"),
        help="Output directory",
    )
    ap.add_argument(
        "--short-names",
        action="store_true",
        help="Use short dataset names (BG1, BG2, UB)",
    )
    ap.add_argument(
        "--order",
        nargs="*",
        default=None,
        help="Row order after grouping by thread suffix (.thrN).",
    )
    ap.add_argument(
        "--threads",
        nargs="*",
        type=int,
        default=[1, 8, 16],
        help="Sub-columns t=... for multi-threaded tools.",
    )
    ap.add_argument(
        "--frame",
        choices=["none", "white", "transparent"],
        default="none",
        help="Border around the final PNG.",
    )
    ap.add_argument(
        "--frame-px",
        type=int,
        default=FRAME_PX,
        help="Border thickness in pixels.",
    )
    ap.add_argument(
        "--phase-mode",
        choices=["bubblefinder-only", "any-with-data", "all"],
        default="bubblefinder-only",
        help=(
            "Phase columns only for: bubblefinder-only (default), any program "
            "with phase metrics (any-with-data), or all tools (all)."
        ),
    )
    args = ap.parse_args()

    results_root = args.results
    if not results_root.is_dir():
        raise SystemExit(
            f"[ERR] Results directory not found: {results_root}"
        )

    runs_bench = load_runs_from_bench_dir(results_root)

    runs_from_benchmarks = pd.DataFrame()
    tsv_f = results_root / "benchmarks.tsv"
    if tsv_f.is_file():
        try:
            runs_from_benchmarks = load_runs(tsv_f)
            print(
                f"[INFO] Loaded {len(runs_from_benchmarks)} runs from {tsv_f}"
            )
        except Exception as e:
            print(f"[WARN] Could not read {tsv_f}: {e}")

    if runs_from_benchmarks is not None and not runs_from_benchmarks.empty:
        if runs_bench is not None and not runs_bench.empty:
            runs_bench_all = pd.concat(
                [runs_from_benchmarks, runs_bench], ignore_index=True
            )
            runs_bench_all = runs_bench_all.drop_duplicates(
                subset=["dataset", "program", "replicate"], keep="first"
            )
        else:
            runs_bench_all = runs_from_benchmarks
    else:
        runs_bench_all = runs_bench

    json_phase_runs = load_phase_metrics_from_results(results_root)

    runs = merge_runs_with_phase_json(runs_bench_all, json_phase_runs)

    if runs is None or runs.empty:
        print("[WARN] No runs detected.")
        runs = pd.DataFrame(
            columns=[
                "dataset",
                "program",
                "replicate",
                "elapsed_s",
                "max_rss_kb",
                "timed_out",
                "timeout_s",
                "exit_status",
                "io_s",
                "build_s",
                "logic_s",
                "io_peak_kb",
                "build_peak_kb",
                "logic_peak_kb",
            ]
        )

    agg = aggregate_runs(runs)

    dataset_rename = DEFAULT_DATASET_RENAME if args.short_names else None
    time_df, mem_df = build_multiheader_tables(
        agg,
        threads=args.threads,
        dataset_rename=dataset_rename,
        dataset_order=args.order,
    )

    if args.phase_mode == "bubblefinder-only":
        phase_progs = PHASE_PROGRAMS_DEFAULT
    elif args.phase_mode == "any-with-data":
        phase_progs = infer_phase_programs_from_data(agg)
    else:
        phase_progs = None

    time_ph_df, mem_ph_df = build_multiheader_phase_tables(
        agg,
        threads=args.threads,
        dataset_rename=dataset_rename,
        dataset_order=args.order,
        phase_programs=phase_progs,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    export_table_png_and_tsv(
        time_df,
        "Wall clock",
        out_png=args.outdir / "table_time.png",
        out_tsv=args.outdir / "table_time.tsv",
        out_html=args.outdir / "table_time.html",
    )
    export_table_png_and_tsv(
        mem_df,
        "Max RSS",
        out_png=args.outdir / "table_memory.png",
        out_tsv=args.outdir / "table_memory.tsv",
        out_html=args.outdir / "table_memory.html",
    )
    export_table_png_and_tsv(
        time_ph_df,
        "Wall clock (I/O, BUILD, LOGIC)",
        out_png=args.outdir / "table_time_phases.png",
        out_tsv=args.outdir / "table_time_phases.tsv",
        out_html=args.outdir / "table_time_phases.html",
    )
    export_table_png_and_tsv(
        mem_ph_df,
        "Peak RSS (I/O, BUILD, LOGIC)",
        out_png=args.outdir / "table_memory_phases.png",
        out_tsv=args.outdir / "table_memory_phases.tsv",
        out_html=args.outdir / "table_memory_phases.html",
    )

    export_table_latex(
        time_df,
        "Wall clock",
        args.outdir / "table_time.tex",
        label="tab:time",
    )
    export_table_latex(
        mem_df,
        "Max RSS",
        args.outdir / "table_memory.tex",
        label="tab:memory",
    )
    export_table_latex(
        time_ph_df,
        "Wall clock (I/O, BUILD, LOGIC)",
        args.outdir / "table_time_phases.tex",
        label="tab:time_phases",
    )
    export_table_latex(
        mem_ph_df,
        "Peak RSS (I/O, BUILD, LOGIC)",
        args.outdir / "table_memory_phases.tex",
        label="tab:memory_phases",
    )

    _add_frame_in_dir(args.outdir, args.frame_px, args.frame)

    print("Files written to:", args.outdir.resolve())
    print("- table_time(.png/.tsv/.html/.tex)")
    print("- table_memory(.png/.tsv/.html/.tex)")
    print("- table_time_phases(.png/.tsv/.html/.tex)")
    print("- table_memory_phases(.png/.tsv/.html/.tex)")
    if not HAS_DFI:
        print(
            "Note: for PNG output, install dataframe_image and kaleido "
            "(pip install dataframe_image kaleido)."
        )
    if args.frame != "none" and not (
        HAS_PIL or shutil.which("convert") or shutil.which("magick")
    ):
        print(
            "Note: for frame support, install Pillow or ImageMagick "
            "(package 'convert' or 'magick')."
        )
    print(
        "LaTeX note: add \\usepackage[table]{xcolor} and "
        "\\usepackage{booktabs}."
    )


if __name__ == "__main__":
    main()