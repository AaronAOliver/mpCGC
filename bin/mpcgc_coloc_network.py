#!/usr/bin/env python3
r"""Enzyme colocalisation network from the CGC catalog.

Nodes are enzyme families; an edge joins two families that appear together in the
same CGC. Edge weight is the number of distinct CGCs supporting the pair, and
n_genomes records how many different genomes those CGCs come from, so a pair
driven by one over-represented lineage can be told apart from a pair that recurs
across the ocean microbiome.

Writes a plain edge list plus a Cytoscape-compatible JSON graph.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

_FAM = re.compile(r"^(GH|PL|CE|CBM|AA|GT)\d+(_\d+)?$")
_S1 = re.compile(r"^S1(_\d+)?$")
_SPLIT = re.compile(r"[+|;,]")


def families_in(annotation, include_accessory):
    out = set()
    for part in _SPLIT.split(annotation):
        part = part.strip()
        if not part:
            continue
        if _S1.match(part):
            out.add(part)
        elif _FAM.match(part):
            if include_accessory or part[:2] not in ("CB", "AA", "GT"):
                out.add(part)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--edges", required=True, type=Path)
    ap.add_argument("--json", required=True, type=Path)
    ap.add_argument("--min-support", type=int, default=2,
                    help="minimum number of CGCs supporting an edge")
    ap.add_argument("--include-accessory", action="store_true",
                    help="also place CBM, AA and GT modules in the network")
    args = ap.parse_args()

    cgc_fams = defaultdict(set)
    cgc_mag = {}
    with args.catalog.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        for row in rdr:
            if len(row) < 9:
                continue
            cgc, mag, ann = row[0], row[1], row[8]
            cgc_mag[cgc] = mag
            fams = families_in(ann, args.include_accessory)
            if fams:
                cgc_fams[cgc] |= fams

    node_cgcs = Counter()
    node_mags = defaultdict(set)
    edge_cgcs = Counter()
    edge_mags = defaultdict(set)

    for cgc, fams in cgc_fams.items():
        mag = cgc_mag[cgc]
        members = sorted(fams)
        for f in members:
            node_cgcs[f] += 1
            node_mags[f].add(mag)
        for i, a in enumerate(members):
            for b in members[i + 1:]:
                edge_cgcs[(a, b)] += 1
                edge_mags[(a, b)].add(mag)

    edges = [(a, b, n) for (a, b), n in edge_cgcs.items()
             if n >= args.min_support]
    edges.sort(key=lambda e: (-e[2], e[0], e[1]))

    with args.edges.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["family_1", "family_2", "n_cgcs", "n_genomes",
                    "jaccard"])
        for a, b, n in edges:
            union = node_cgcs[a] + node_cgcs[b] - n
            w.writerow([a, b, n, len(edge_mags[(a, b)]),
                        "{:.6f}".format(n / union) if union else "0"])

    in_graph = {f for a, b, _n in edges for f in (a, b)}
    graph = {
        "format_version": "1.0",
        "generated_by": "mpCGC",
        "data": {"name": "mpCGC enzyme colocalization network",
                 "min_support": args.min_support},
        "elements": {
            "nodes": [
                {"data": {"id": f, "family": f, "n_cgcs": node_cgcs[f],
                          "n_genomes": len(node_mags[f]),
                          "class": f[:2] if not f.startswith("S1") else "S1"}}
                for f in sorted(in_graph)
            ],
            "edges": [
                {"data": {"source": a, "target": b, "weight": n,
                          "n_genomes": len(edge_mags[(a, b)])}}
                for a, b, n in edges
            ],
        },
    }
    with args.json.open("w", encoding="utf-8") as out:
        json.dump(graph, out, indent=1)

    print("[coloc] {:,} families, {:,} edges from {:,} CGCs "
          "(min support {})".format(len(in_graph), len(edges), len(cgc_fams),
                                    args.min_support), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
