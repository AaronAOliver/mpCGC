#!/usr/bin/env python3
r"""Write one protein FASTA per enzyme family from the mpCGC catalog.

Every CGC member gene carrying a GH / PL / CE / CBM / AA / GT or S1 annotation
contributes its protein sequence to that family's fasta. Sequences are pulled
from the per-genome uniInput.faa files produced during identification. These
fastas are the input to the mining dataflow (all-vs-all DIAMOND, Leiden
communities, MUSCLE alignment, FastTree phylogeny).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# one annotation token may hold several families, e.g. "CAZyme|GH13_3+CBM48"
_FAM = re.compile(r"\b((?:GH|PL|CE|CBM|AA|GT)\d+(?:_\d+)?|S1(?:_\d+)?)\b")


def read_fasta(path: Path):
    """Yield (header_id, sequence) pairs without loading the whole file."""
    pid, chunks = None, []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(chunks)
                pid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if pid is not None:
        yield pid, "".join(chunks)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--faa", nargs="+", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--counts", required=True, type=Path)
    ap.add_argument("--protein-map", type=Path, default=None,
                    help="TSV of protein/genome/CGC/family, used to colour "
                         "family phylogenies by host taxonomy")
    ap.add_argument("--min-seqs", type=int, default=1)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # protein id -> families, plus protein id -> CGC id for the fasta header
    want: dict[str, set[str]] = defaultdict(set)
    origin: dict[str, tuple[str, str]] = {}
    with args.catalog.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        for row in rdr:
            if len(row) < 9:
                continue
            cgc, mag, gtype, _contig, pid, _s, _e, _str, ann = row[:9]
            fams = set(_FAM.findall(ann))
            if not fams:
                continue
            want[pid] |= fams
            origin[pid] = (mag, cgc)

    # stream each genome's proteins once, appending to the family files
    handles: dict[str, object] = {}
    counts: dict[str, int] = defaultdict(int)
    n_found = 0
    emitted: list[tuple[str, str, str, str]] = []
    try:
        for faa in sorted(args.faa):
            for pid, seq in read_fasta(faa):
                fams = want.get(pid)
                if not fams:
                    continue
                n_found += 1
                mag, cgc = origin[pid]
                emitted.append((pid, mag, cgc, ",".join(sorted(fams))))
                for fam in fams:
                    fh = handles.get(fam)
                    if fh is None:
                        fh = handles[fam] = (args.outdir / f"{fam}.faa").open(
                            "w", encoding="utf-8")
                    fh.write(f">{pid} mag={mag} cgc={cgc} family={fam}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                    counts[fam] += 1
    finally:
        for fh in handles.values():
            fh.close()

    # drop families below the threshold so downstream steps are not started
    for fam, n in list(counts.items()):
        if n < args.min_seqs:
            (args.outdir / f"{fam}.faa").unlink(missing_ok=True)
            del counts[fam]

    if args.protein_map:
        kept = set(counts)
        with args.protein_map.open("w", newline="", encoding="utf-8") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["protein", "mag", "cgc", "families"])
            for pid, mag, cgc, fams in emitted:
                # only proteins whose family survived --min-seqs are in a fasta
                if any(f in kept for f in fams.split(",")):
                    w.writerow([pid, mag, cgc, fams])

    with args.counts.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "n_proteins"])
        for fam in sorted(counts, key=lambda f: (-counts[f], f)):
            w.writerow([fam, counts[fam]])

    missing = len(want) - n_found
    print(f"[proteins] {len(counts):,} families  {n_found:,} sequences written"
          + (f"  ({missing:,} catalog proteins not found in the fastas)" if missing else ""),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
