#!/usr/bin/env python3
"""Give dbCAN short sequence IDs when the input's IDs are too long for it.

dbCAN 5.2.9 cuts protein IDs longer than 80 characters (and renumbers them to
keep them unique) when it reads a protein FASTA, without warning and without
changing the GFF. The FASTA and GFF then disagree, no gene is annotated and
the genome yields zero CGCs. In nucleotide mode, Pyrodigal names proteins
<contig>_<n>, so a long contig name leads to the same cut.

This script writes the copy of the input that dbCAN works on:

  * IDs short enough: an unchanged copy, and an empty map (so dbCAN, which
    rewrites its input in place, never touches the user's file);
  * otherwise: every protein ID (protein mode) or contig name (nucleotide
    mode) replaced by a fixed-width placeholder such as mpcgcP000000001,
    with the FASTA and GFF changed together, plus a two-column map
    placeholder -> original.

mpcgc_restore_ids.py puts the original IDs back in every output file.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

DBCAN_MAX_ID = 80          # dbCAN 5.2.9 ID_MAX_LENGTH (dbcan/IO/fasta.py)
PYRODIGAL_SUFFIX = 7       # room for "_<gene number>" after a contig name


def opener(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def fasta_ids(path: Path) -> list[str]:
    with opener(path) as fh:
        return [ln[1:].split()[0] if ln[1:].split() else "" for ln in fh if ln.startswith(">")]


def placeholder(kind: str, i: int) -> str:
    return f"mpcgc{kind}{i:09d}"


def write_fasta(src: Path, dst: Path, mapping: dict[str, str]) -> None:
    """Copy a FASTA, replacing header IDs found in mapping (descriptions dropped
    for renamed records, so nothing long reaches the gene caller)."""
    with opener(src) as fin, open(dst, "w") as fout:
        for ln in fin:
            if ln.startswith(">"):
                parts = ln[1:].split(maxsplit=1)
                sid = parts[0] if parts else ""
                if sid in mapping:
                    ln = f">{mapping[sid]}\n"
            fout.write(ln)


def write_gff(src: Path, dst: Path, proteins: dict[str, str], contigs: dict[str, str]) -> None:
    """Copy a GFF, renaming contigs in column 1 and any attribute value that is
    exactly a renamed protein ID (ID=, protein_id=, Name=, Parent=, ...)."""
    with opener(src) as fin, open(dst, "w") as fout:
        for ln in fin:
            if ln.startswith("#") or not ln.strip():
                fout.write(ln)
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                fout.write(ln)
                continue
            f[0] = contigs.get(f[0], f[0])
            if proteins:
                attrs = []
                for kv in f[8].split(";"):
                    k, eq, v = kv.partition("=")
                    attrs.append(f"{k}={proteins.get(v, v)}" if eq else kv)
                f[8] = ";".join(attrs)
            fout.write("\t".join(f) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["prok", "meta", "protein"])
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--gff", type=Path, help="protein mode: gene coordinates matching --fasta")
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-gff", type=Path)
    ap.add_argument("--map", required=True, type=Path, help="output: placeholder -> original ID")
    ap.add_argument("--sample", default="", help="sample name, for the log message")
    args = ap.parse_args()

    ids = fasta_ids(args.fasta)
    if len(set(ids)) != len(ids):
        sys.exit(f"[rename] {args.fasta.name}: duplicate sequence IDs; each FASTA ID must be unique")

    proteins: dict[str, str] = {}
    contigs: dict[str, str] = {}
    if args.mode == "protein":
        if ids and max(map(len, ids)) > DBCAN_MAX_ID:
            proteins = {sid: placeholder("P", i) for i, sid in enumerate(ids, 1)}
    else:
        if ids and max(map(len, ids)) > DBCAN_MAX_ID - PYRODIGAL_SUFFIX:
            contigs = {sid: placeholder("C", i) for i, sid in enumerate(ids, 1)}

    write_fasta(args.fasta, args.out_fasta, proteins or contigs)
    if args.gff and args.out_gff:
        write_gff(args.gff, args.out_gff, proteins, contigs)

    with open(args.map, "w") as fh:
        fh.write("placeholder\toriginal\n")
        for orig, ph in {**proteins, **contigs}.items():
            fh.write(f"{ph}\t{orig}\n")

    renamed = proteins or contigs
    if renamed:
        what = "protein IDs" if proteins else "contig names"
        longest = max(renamed, key=len)
        print(f"WARNING [{args.sample or args.fasta.name}] {len(renamed):,} {what} replaced with short "
              f"placeholders for dbCAN (longest original: {len(longest)} characters, dbCAN's limit is "
              f"{DBCAN_MAX_ID}). Original IDs are restored in every output.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
