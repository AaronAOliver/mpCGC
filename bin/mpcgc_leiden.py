#!/usr/bin/env python3
r"""Leiden communities on a per-family sequence similarity network.

The DIAMOND all-vs-all output for one enzyme family becomes an undirected graph:
nodes are proteins, an edge joins two proteins whose alignment passes the E-value
cutoff, and the edge weight is the bitscore. Self-hits are dropped and reciprocal
hits are collapsed onto one canonical edge keeping the better bitscore.

Communities are found with the Leiden algorithm optimising modularity at
--resolution. Each community is a group of sequences more similar to one another
than to the rest of the family, which is what separates, for example, the
characterised members of a family from the uncharacterised ones.

Writes a per-protein community assignment and a Cytoscape-compatible JSON graph.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ssn", required=True, type=Path)
    ap.add_argument("--family", required=True)
    ap.add_argument("--evalue", type=float, default=1e-30)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--out-communities", required=True, type=Path)
    ap.add_argument("--out-json", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    best = {}
    n_lines = 0
    with open_maybe_gzip(args.ssn) as fh:
        for raw in fh:
            n_lines += 1
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            q, s = parts[0], parts[1]
            if q == s:
                continue
            try:
                ev = float(parts[10])
                bs = float(parts[11])
            except ValueError:
                continue
            if ev > args.evalue:
                continue
            key = (q, s) if q < s else (s, q)
            if bs > best.get(key, -1.0):
                best[key] = bs

    nodes = sorted({p for pair in best for p in pair})
    if not nodes:
        # a family with no edge above the cutoff is still reported, as singletons
        write_outputs(args, [], {}, {}, 0)
        print("[leiden] {}: no edges passed E<={:g}".format(
            args.family, args.evalue), file=sys.stderr)
        return 0

    import igraph as ig

    idx = {p: i for i, p in enumerate(nodes)}
    edges = [(idx[a], idx[b]) for a, b in best]
    weights = [best[k] for k in best]

    g = ig.Graph(n=len(nodes), edges=edges, directed=False)
    g.vs["name"] = nodes
    g.es["weight"] = weights

    part = g.community_leiden(objective_function="modularity",
                              weights=weights, resolution=args.resolution,
                              n_iterations=-1)
    membership = part.membership
    sizes = Counter(membership)
    # rank communities largest first so community_rank 1 is the biggest
    rank = {c: r for r, (c, _n) in enumerate(sizes.most_common(), start=1)}
    assign = {nodes[i]: membership[i] for i in range(len(nodes))}

    write_outputs(args, nodes, assign, {"sizes": sizes, "rank": rank,
                                        "edges": edges, "weights": weights},
                  len(set(membership)))
    print("[leiden] {}: {:,} proteins, {:,} edges, {} communities "
          "(modularity {:.3f})".format(args.family, len(nodes), len(best),
                                       len(set(membership)),
                                       g.modularity(membership, weights=weights)),
          file=sys.stderr)
    return 0


def write_outputs(args, nodes, assign, extra, n_comm):
    sizes = extra.get("sizes", Counter())
    rank = extra.get("rank", {})

    with args.out_communities.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["protein", "family", "leiden_community", "community_rank",
                    "community_size"])
        for p in nodes:
            c = assign[p]
            w.writerow([p, args.family, c, rank.get(c, ""), sizes.get(c, 1)])

    graph = {
        "format_version": "1.0",
        "generated_by": "mpCGC",
        "data": {"name": "{} sequence similarity network".format(args.family),
                 "family": args.family, "evalue": args.evalue,
                 "resolution": args.resolution, "n_communities": n_comm},
        "elements": {
            "nodes": [
                {"data": {"id": p, "community": assign[p],
                          "community_rank": rank.get(assign[p], ""),
                          "community_size": sizes.get(assign[p], 1)}}
                for p in nodes
            ],
            "edges": [
                {"data": {"source": nodes[a], "target": nodes[b],
                          "bitscore": wt}}
                for (a, b), wt in zip(extra.get("edges", []),
                                      extra.get("weights", []))
            ],
        },
    }
    with args.out_json.open("w", encoding="utf-8") as out:
        json.dump(graph, out, indent=1)


if __name__ == "__main__":
    sys.exit(main())
