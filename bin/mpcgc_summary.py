#!/usr/bin/env python3
r"""Summary metrics for an mpCGC run: the csv output of the summary dataflow.

Per genome: CGC count, genes inside CGCs, the counts of each signature gene type,
distinct enzyme fingerprints and how many of those fingerprints are unique to that
genome within the run.

Per lineage: genome and CGC totals, observed and Chao2-estimated fingerprint
richness (read from the richness table), and the share of the estimated richness
already recovered.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

GENE_TYPES = ["CAZyme", "TC", "TF", "STP", "Sulfatase", "Peptidase",
              "prodoric", "null"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--fingerprints", required=True, type=Path)
    ap.add_argument("--asymptotes", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-lineage", required=True, type=Path)
    args = ap.parse_args()

    # ---- catalog ------------------------------------------------------
    cgcs = defaultdict(set)
    gene_counts = defaultdict(Counter)
    with args.catalog.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        for row in rdr:
            if len(row) < 9:
                continue
            cgc, mag, gtype = row[0], row[1], row[2]
            cgcs[mag].add(cgc)
            gene_counts[mag][gtype] += 1

    # ---- fingerprints -------------------------------------------------
    mag_fps = defaultdict(set)
    lineage_of = {}
    fp_mags = defaultdict(set)
    lineage_fps = defaultdict(set)
    lineage_mags = defaultdict(set)
    with args.fingerprints.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr)
        i_mag, i_grp, i_fp = (head.index("MAG"), head.index("lineage"),
                              head.index("counted_as"))
        for row in rdr:
            if len(row) <= max(i_mag, i_grp, i_fp) or not row[i_fp]:
                continue
            mag, grp, fp = row[i_mag], row[i_grp], row[i_fp]
            mag_fps[mag].add(fp)
            lineage_of[mag] = grp
            fp_mags[fp].add(mag)
            lineage_fps[grp].add(fp)
            lineage_mags[grp].add(mag)

    with args.out.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(["MAG", "lineage", "n_CGCs", "n_genes_in_CGCs",
                    "n_fingerprints", "n_fingerprints_unique_to_genome"]
                   + ["n_" + t for t in GENE_TYPES])
        for mag in sorted(cgcs):
            fps = mag_fps.get(mag, set())
            uniq = sum(1 for fp in fps if len(fp_mags[fp]) == 1)
            w.writerow([mag, lineage_of.get(mag, "other"), len(cgcs[mag]),
                        sum(gene_counts[mag].values()), len(fps), uniq]
                       + [gene_counts[mag].get(t, 0) for t in GENE_TYPES])

    # ---- lineage roll-up ---------------------------------------------
    asym = {}
    if args.asymptotes.exists():
        with args.asymptotes.open(newline="", encoding="utf-8") as fh:
            rdr = csv.reader(fh, delimiter="\t")
            head = next(rdr, None)
            if head:
                col = {h: i for i, h in enumerate(head)}
                for row in rdr:
                    if row:
                        asym[row[col["lineage"]]] = row

    with args.out_lineage.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(["lineage", "n_genomes", "n_CGCs", "S_obs", "S_chao2",
                    "ci_lower", "ci_upper", "pct_recovered",
                    "mean_CGCs_per_genome"])
        for grp in sorted(lineage_mags, key=lambda g: -len(lineage_mags[g])):
            mags = lineage_mags[grp]
            n_cgc = sum(len(cgcs.get(m, ())) for m in mags)
            row = asym.get(grp)
            if row:
                col = {h: i for i, h in enumerate(
                    ["lineage", "genomes", "S_obs", "Q1", "Q2", "Q0_chao2",
                     "S_chao2", "SE", "ci_lower", "ci_upper", "pct_recovered"])}
                s_obs = row[col["S_obs"]]
                s_est = row[col["S_chao2"]]
                lo, hi = row[col["ci_lower"]], row[col["ci_upper"]]
                pct = row[col["pct_recovered"]]
            else:
                s_obs, s_est, lo, hi, pct = len(lineage_fps[grp]), "", "", "", ""
            w.writerow([grp, len(mags), n_cgc, s_obs, s_est, lo, hi, pct,
                        round(n_cgc / len(mags), 2) if mags else 0])

    print("[summary] {:,} genomes, {:,} CGCs, {:,} lineages".format(
        len(cgcs), sum(len(v) for v in cgcs.values()), len(lineage_mags)),
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
