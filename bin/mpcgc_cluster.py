#!/usr/bin/env python3
r"""Hierarchical clustering of CGCs by their binary enzyme fingerprints.

Distinct fingerprints are clustered on pairwise Hamming distance over the binary
presence/absence vector of S1 / GH / PL / CE families. Two clusters are treated
as overlapping in function when one contains a superset of all annotations of the
other, so the subset relation is reported alongside the distance-based cut: a
fingerprint nested inside another is functionally redundant with it even when the
Hamming distance between them is large.

Outputs the linkage matrix, a flat clustering at --cut, and, for each
fingerprint, the number of fingerprints it is a superset of.
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
    ap.add_argument("--out-linkage", required=True, type=Path)
    ap.add_argument("--out-clusters", required=True, type=Path)
    ap.add_argument("--cut", type=float, default=0.5,
                    help="distance at which to cut the tree for flat clusters")
    ap.add_argument("--max-fingerprints", type=int, default=20000,
                    help="cap on distinct fingerprints clustered; the rarest are "
                         "dropped beyond this to keep the distance matrix feasible")
    ap.add_argument("--figure", type=Path, default=None)
    args = ap.parse_args()

    freq = Counter()
    seen = set()
    with args.fingerprints.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr)
        i_cgc, i_fp = head.index("CGC"), head.index("counted_as")
        for row in rdr:
            if len(row) <= max(i_cgc, i_fp) or row[i_cgc] in seen:
                continue
            seen.add(row[i_cgc])
            if row[i_fp]:
                freq[row[i_fp]] += 1

    if not freq:
        sys.exit("no fingerprints found in {}".format(args.fingerprints))

    kept = [fp for fp, _ in freq.most_common(args.max_fingerprints)]
    dropped = len(freq) - len(kept)
    sets = [frozenset(fp.split(",")) for fp in kept]

    families = sorted({f for s in sets for f in s})
    col = {f: i for i, f in enumerate(families)}
    X = np.zeros((len(sets), len(families)), dtype=bool)
    for r, s in enumerate(sets):
        for f in s:
            X[r, col[f]] = True

    # how many other fingerprints each one fully contains
    token_index = defaultdict(list)
    for i, s in enumerate(sets):
        for f in s:
            token_index[f].append(i)
    n_subsets = [0] * len(sets)
    for i, s in enumerate(sets):
        rarest = min(s, key=lambda f: len(token_index[f]))
        n_subsets[i] = sum(1 for j in token_index[rarest]
                           if j != i and sets[j] < s)

    labels = np.ones(len(sets), dtype=int)
    if len(sets) >= 3:
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import pdist
        d = pdist(X, metric="hamming")
        link = linkage(d, method="average")
        np.save(args.out_linkage, link)
        labels = fcluster(link, t=args.cut, criterion="distance")
        if args.figure:
            try:
                plot(link, args.figure)
            except Exception as exc:
                print("[cluster] figure skipped: {}".format(exc), file=sys.stderr)
    else:
        np.save(args.out_linkage, np.zeros((0, 4)))

    with args.out_clusters.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["fingerprint", "n_cgcs", "n_tokens", "cluster",
                    "n_fingerprints_contained"])
        order = sorted(range(len(sets)),
                       key=lambda i: (labels[i], -freq[kept[i]]))
        for i in order:
            w.writerow([kept[i], freq[kept[i]], len(sets[i]),
                        int(labels[i]), n_subsets[i]])

    print("[cluster] {:,} distinct fingerprints -> {:,} clusters at d={}".format(
        len(sets), len(set(labels.tolist())), args.cut), file=sys.stderr)
    if dropped:
        print("[cluster] {:,} rare fingerprints beyond --max-fingerprints were "
              "not clustered".format(dropped), file=sys.stderr)
    return 0


def plot(link, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram

    fig, ax = plt.subplots(figsize=(11, 6))
    dendrogram(link, no_labels=True, color_threshold=None, ax=ax,
               link_color_func=lambda _k: "#444444")
    ax.set_ylabel("Hamming distance", fontsize=12)
    ax.set_xlabel("CGC fingerprints", fontsize=12)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig("{}.{}".format(prefix, ext),
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
