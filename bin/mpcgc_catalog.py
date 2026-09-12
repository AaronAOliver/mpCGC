#!/usr/bin/env python3
r"""Merge per-genome dbCAN cgc_standard_out.tsv files into one mpCGC catalog.

The genome name is taken from the filename (<genome>.cgc_standard_out.tsv) and
CGC identifiers are made globally unique as <genome>_CGC<n>, matching the
catalog schema published at https://mpcgcdb.com.

Output columns:
    CGC#  MAG  Gene Type  Contig ID  Protein ID
    Gene Start  Gene Stop  Gene Strand  Gene Annotation
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

# dbCAN-sub appends its internal eCAMI cluster index to family names, e.g.
# "GH140_e33" or "AA1_e". That index is not a CAZy subfamily, and leaving it in
# would scatter one family across dozens of pseudo-families in every downstream
# step, so it is removed. Real CAZy subfamilies (GH13_3) and sulfatase
# subfamilies (S1_16) are untouched because they never carry the _e marker.
_ECAMI = re.compile(r"\b((?:GH|PL|CE|CBM|AA|GT)\d+(?:_\d+)?)_e\d*")


def strip_ecami(annotation):
    return _ECAMI.sub(r"\1", annotation)

HEADER = ["CGC#", "MAG", "Gene Type", "Contig ID", "Protein ID",
          "Gene Start", "Gene Stop", "Gene Strand", "Gene Annotation"]

SUFFIX = ".cgc_standard_out.tsv"


def genome_name(path: Path) -> str:
    name = path.name
    return name[: -len(SUFFIX)] if name.endswith(SUFFIX) else path.stem


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--per-genome", required=True, type=Path)
    ap.add_argument("--keep-ecami", action="store_true",
                    help="keep dbCAN-sub eCAMI cluster indices such as GH140_e33")
    args = ap.parse_args()

    n_rows = 0
    per_genome: dict[str, Counter] = {}
    cgc_counts: dict[str, int] = {}

    with args.catalog.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)

        for path in sorted(args.inputs):
            mag = genome_name(path)
            tally = per_genome.setdefault(mag, Counter())
            seen: set[str] = set()

            with path.open(newline="", encoding="utf-8") as fh:
                rdr = csv.reader(fh, delimiter="\t")
                head = next(rdr, None)
                if head is None:
                    continue
                # tolerate dbCAN writing the header more than once
                for row in rdr:
                    if len(row) < 8 or row[0] == "CGC#":
                        continue
                    cgc, gtype, contig, pid, start, stop, strand, ann = row[:8]
                    gid = f"{mag}_{cgc}"
                    seen.add(gid)
                    tally[gtype] += 1
                    if not args.keep_ecami:
                        ann = strip_ecami(ann)
                    w.writerow([gid, mag, gtype, contig, pid,
                                start, stop, strand, ann])
                    n_rows += 1
            cgc_counts[mag] = len(seen)

    with args.per_genome.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        types = ["CAZyme", "TC", "TF", "STP", "Sulfatase", "Peptidase",
                 "prodoric", "null"]
        w.writerow(["MAG", "n_CGCs", "n_genes_in_CGCs"] + [f"n_{t}" for t in types])
        for mag in sorted(cgc_counts):
            tally = per_genome[mag]
            w.writerow([mag, cgc_counts[mag], sum(tally.values())]
                       + [tally.get(t, 0) for t in types])

    print(f"[catalog] {len(cgc_counts):,} genomes  "
          f"{sum(cgc_counts.values()):,} CGCs  {n_rows:,} gene rows",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
