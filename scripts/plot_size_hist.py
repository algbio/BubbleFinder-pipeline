#!/usr/bin/env python3
import argparse
import pathlib

import numpy as np
import matplotlib.pyplot as plt

# BIN_EDGES_BASE = np.array([1, 31, 61, 101, 201, np.inf])
BIN_EDGES_BASE = np.array([3, 21, 101, 501, 1001, np.inf])
# BIN_LABELS_BASE = ["1–30", "31–60", "61–100", "101–200", "201+"]
BIN_LABELS_BASE = ["3–20", "21–100", "101–500", "501–1000", "1001+"]

KIND_LABELS = {
    "blocks": {
        "xlabel": "Block size (#nodes + #edges)",
        "title": "Block size distribution",
    },
    "spqr": {
        "xlabel": "3-connected component size (#nodes + #edges)",
        "title": "3-connected component size distribution",
    },
}


def get_labels_for_kind(kind: str):
    info = KIND_LABELS.get(kind, KIND_LABELS["blocks"])
    return info["xlabel"], info["title"]


def list_block_files(directory, pattern="*_block_sizes.txt"):
    d = pathlib.Path(directory)
    files = sorted(d.glob(pattern))
    if not files:
        raise SystemExit(f"No files match '{pattern}' in directory '{directory}'")
    return files


def load_dataset_sizes(files, labels=None):
    datasets = {}
    if labels is not None and len(labels) != len(files):
        raise SystemExit("Number of --labels must match number of files")

    for i, path in enumerate(files):
        try:
            sizes = np.loadtxt(path, dtype=int)
        except OSError as e:
            raise SystemExit(f"Error reading {path}: {e}")

        if np.ndim(sizes) == 0:
            sizes = np.array([int(sizes)])

        label = labels[i] if labels is not None else path.stem
        datasets[label] = sizes

    return datasets


def apply_rename_mapping(datasets, rename_pairs):
    if not rename_pairs:
        return datasets

    mapping = {}
    for item in rename_pairs:
        if "=" not in item:
            raise SystemExit(f"Invalid --rename argument '{item}', expected OLD=NEW")
        old, new = item.split("=", 1)
        old = old.strip()
        new = new.strip()
        if not old or not new:
            raise SystemExit(f"Invalid --rename argument '{item}', empty OLD or NEW")
        mapping[old] = new

    renamed = {}
    for label, data in datasets.items():
        new_label = mapping.get(label, label)
        if new_label in renamed:
            raise SystemExit(
                f"Duplicate label after renaming: '{new_label}'. "
                "Please ensure each final label is unique"
            )
        renamed[new_label] = data

    return renamed


def compute_fixed_bins(datasets):
    all_sizes = np.concatenate([np.asarray(v) for v in datasets.values()])
    all_sizes = all_sizes[all_sizes > 0]

    if all_sizes.size == 0:
        return np.array([1, 2]), ["1"]

    edges = BIN_EDGES_BASE.copy()
    if np.isinf(edges[-1]):
        max_size = int(all_sizes.max())
        if max_size >= edges[-2]:
            edges[-1] = max_size + 1
        else:
            edges[-1] = edges[-2] + 1

    return edges, BIN_LABELS_BASE


def compute_histograms(datasets, bin_edges):
    histos = {}
    for label, sizes in datasets.items():
        sizes = np.asarray(sizes)
        sizes = sizes[sizes > 0]

        if sizes.size == 0:
            counts = np.zeros(len(bin_edges) - 1, dtype=int)
        else:
            counts, _ = np.histogram(sizes, bins=bin_edges)

        total = counts.sum()
        probs = counts / total if total > 0 else np.zeros_like(counts, dtype=float)
        histos[label] = (counts, probs)

    return histos


def plot_histograms(bin_labels, histos, kind="blocks", output=None):
    labels = list(histos.keys())
    nbins = len(bin_labels)
    ndatasets = len(labels)

    x_base = np.arange(nbins)
    width = 0.8 / max(ndatasets, 1)

    fig, ax = plt.subplots(figsize=(8, 5))

    hatches = ['/', '\\', 'x', '-', '+', 'o', 'O', '.', '*']

    for i, label in enumerate(labels):
        _, probs = histos[label]
        x = x_base + i * width
        ax.bar(
            x, probs, width=width,
            label=label,
            alpha=0.8,
            edgecolor="black",
            hatch=hatches[i % len(hatches)],
        )

    ax.set_xticks(x_base + width * ndatasets / 2 - width / 2)
    ax.set_xticklabels(bin_labels)

    xlabel, _ = get_labels_for_kind(kind)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability P(size ∈ bin)")
    ax.set_ylim(0, 1)
    ax.legend(title="Dataset")
    ax.grid(axis="y", linestyle=":", alpha=0.5)

    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1)

    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=300)
        print(f"[plot] Figure written to {output}")
    else:
        plt.show()

def plot_histograms_and_boxplots(bin_labels, histos, datasets, kind="blocks", output=None):
    labels = list(histos.keys())
    nbins = len(bin_labels)
    ndatasets = len(labels)

    x_base = np.arange(nbins)
    width = 0.8 / max(ndatasets, 1)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]}
    )

    hatches = ['/', '\\', 'x', '-', '+', 'o', 'O', '.', '*']

    for i, label in enumerate(labels):
        _, probs = histos[label]
        x = x_base + i * width

        ax1.bar(
            x, probs, width=width,
            label=label,
            alpha=0.8,
            edgecolor="black",
            hatch=hatches[i % len(hatches)],
            rasterized=True,
        )

    ax1.set_xticks(x_base + width * ndatasets / 2 - width / 2)
    ax1.set_xticklabels(bin_labels)

    xlabel, title = get_labels_for_kind(kind)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("Probability P(size ∈ bin)")
    ax1.set_ylim(0, 1)
    ax1.legend(title="Dataset", fontsize=7, loc="best")
    ax1.grid(axis="y", linestyle=":", alpha=0.5)
    ax1.set_yscale("log")
    ax1.set_ylim(1e-4, 1)

    # --- Boxplot ---
    data_for_box = [np.asarray(datasets[label]) for label in labels]
    bp = ax2.boxplot(
        data_for_box,
        vert=True,
        labels=labels,
        showfliers=True
    )

    for element in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
        for item in bp[element]:
            item.set_rasterized(True)

    ax2.set_title(title)
    ax2.set_ylabel(xlabel)

    plt.setp(
        ax2.get_xticklabels(),
        rotation=45,
        ha='right',
        rotation_mode='anchor'
    )

    ax2.set_yscale("log")

    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=300)
        print(f"[plot] Figure written to {output}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Histogram of block / 3-connected component sizes"
    )
    parser.add_argument(
        "directory",
        help="Directory containing *_block_sizes.txt or *_spqr_sizes.txt files",
    )
    parser.add_argument(
        "--pattern",
        default="*_block_sizes.txt",
        help="Glob pattern for files inside the directory",
    )
    parser.add_argument(
        "--kind",
        choices=["blocks", "spqr"],
        default="blocks",
        help="What is being plotted; controls axis labels/titles",
    )
    parser.add_argument(
        "--labels",
        nargs="*",
        help=(
            "Override dataset labels in the same order as matched files. "
            "If not provided, labels are taken from file stems."
        ),
    )
    parser.add_argument(
        "--rename",
        nargs="*",
        metavar="OLD=NEW",
        help=(
            "Rename datasets after loading. OLD is the original label (stem or "
            "label from --labels), NEW is the desired name. Example: "
            "--rename sample1=A sample2=B."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file for the figure (png/pdf)",
    )

    args = parser.parse_args()

    files = list_block_files(args.directory, args.pattern)
    datasets = load_dataset_sizes(files, args.labels)
    datasets = apply_rename_mapping(datasets, args.rename)

    bin_edges, bin_labels = compute_fixed_bins(datasets)
    histos = compute_histograms(datasets, bin_edges)

    print("Distribution summary:")
    print("Bins:", bin_labels)
    for label, (counts, probs) in histos.items():
        print(f"  {label}:")
        print("    counts:", counts.tolist())
        print("    probs :", [f"{p:.3f}" for p in probs])

    plot_histograms_and_boxplots(
        bin_labels, histos, datasets, kind=args.kind, output=args.output
    )


if __name__ == "__main__":
    main()