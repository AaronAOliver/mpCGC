#!/usr/bin/env python3
r"""Asymmetric enzyme colocalisation matrix over CGC fingerprints.

The matrix is row-conditioned rather than symmetric:

    M[anchor, colocalized] = P(colocalized family in a CGC | anchor family in it)

so a rare family that almost always travels with a common partner shows a high
value in its own row without dragging down the common family's row. Families
seen in fewer than --min-cgc distinct CGCs are dropped, since a handful of
observations cannot support a conditional probability.

Rows and columns share one ordering, taken from average-linkage hierarchical
clustering of the Jaccard distance between family occurrence vectors, which puts
families that co-occur into contiguous blocks along the diagonal.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fingerprints", required=True, type=Path)
    ap.add_argument("--min-cgc", type=int, default=100)
    ap.add_argument("--out-matrix", required=True, type=Path)
    ap.add_argument("--out-order", required=True, type=Path)
    ap.add_argument("--figure", type=Path, default=None)
    args = ap.parse_args()

    # read each CGC's fingerprint once
    fingerprints = []
    seen = set()
    with args.fingerprints.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr)
        i_cgc = head.index("CGC")
        i_fp = head.index("fingerprint")
        for row in rdr:
            if len(row) <= max(i_cgc, i_fp) or row[i_cgc] in seen:
                continue
            seen.add(row[i_cgc])
            toks = [t for t in row[i_fp].split(",") if t]
            if toks:
                fingerprints.append(frozenset(toks))

    if not fingerprints:
        sys.exit("no fingerprints found in {}".format(args.fingerprints))

    counts = Counter()
    for fp in fingerprints:
        counts.update(fp)
    families = sorted(f for f, n in counts.items() if n >= args.min_cgc)
    if len(families) < 2:
        print("[cooccurrence] only {} family passes --min-cgc {}; nothing to "
              "cluster".format(len(families), args.min_cgc), file=sys.stderr)

    index = {f: i for i, f in enumerate(families)}
    n = len(families)
    pair = np.zeros((n, n), dtype=np.int64)
    solo = np.zeros(n, dtype=np.int64)

    for fp in fingerprints:
        present = sorted(index[f] for f in fp if f in index)
        if not present:
            continue
        for i in present:
            solo[i] += 1
        for a_pos, i in enumerate(present):
            for j in present[a_pos + 1:]:
                pair[i, j] += 1
                pair[j, i] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        matrix = np.where(solo[:, None] > 0, pair / solo[:, None], 0.0)
    np.fill_diagonal(matrix, 1.0)

    order = list(range(n))
    if n >= 3:
        order = cluster_order(pair, solo, n)

    ordered = [families[i] for i in order]
    with args.out_order.open("w", encoding="utf-8") as out:
        out.write("\n".join(ordered) + "\n")

    with args.out_matrix.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["anchor", "n_cgcs"] + ordered)
        for i in order:
            w.writerow([families[i], int(solo[i])]
                       + ["{:.6f}".format(matrix[i, j]) for j in order])

    print("[cooccurrence] {} families x {} CGCs".format(n, len(fingerprints)),
          file=sys.stderr)

    if args.figure and n >= 3:
        try:
            plot(matrix, ordered, order, args.figure)
        except Exception as exc:
            print("[cooccurrence] figure skipped: {}".format(exc), file=sys.stderr)
    return 0


def cluster_order(pair, solo, n):
    """Average-linkage order on Jaccard distance between families."""
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform

    union = solo[:, None] + solo[None, :] - pair
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = np.where(union > 0, pair / union, 0.0)
    dist = 1.0 - jaccard
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0                 # enforce exact symmetry
    link = linkage(squareform(dist, checks=False), method="average")
    return dendrogram(link, no_plot=True)["leaves"]


def plot(matrix, labels, order, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    m = matrix[np.ix_(order, order)]
    size = max(8.0, len(labels) * 0.16)
    fig, ax = plt.subplots(figsize=(size, size))
    # transpose so the anchor family runs along the x-axis
    im = ax.imshow(m.T, cmap="inferno", vmin=0, vmax=1, origin="upper",
                   interpolation="nearest")
    step = 1 if len(labels) <= 80 else max(1, len(labels) // 80)
    ticks = range(0, len(labels), step)
    ax.set_xticks(list(ticks))
    ax.set_xticklabels([labels[i] for i in ticks], rotation=90, fontsize=5.5)
    ax.set_yticks(list(ticks))
    ax.set_yticklabels([labels[i] for i in ticks], fontsize=5.5)
    ax.set_xlabel("Anchor family", fontsize=15)
    ax.set_ylabel("Colocalized family", fontsize=15)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    cbar.set_label("CGCs with colocalized family (%)", fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig("{}.{}".format(prefix, ext),
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
