#!/usr/bin/env python3
r"""Cluster encoding: turn every CGC into a binary enzyme fingerprint.

One token per gene is emitted, using the priority GH > PL > CE > S1, so a gene
annotated "CAZyme|GH13_3+CBM48" contributes GH13_3 and a gene annotated
"Sulfatase|S1_7" contributes S1_7. AA, GT and CBM modules are not encoded: they
are accessory to, rather than diagnostic of, the polysaccharide being targeted.
The fingerprint is the sorted, deduplicated set of those tokens.

With --collapse, fingerprints are reduced to the maximal sets of the inclusion
order (an antichain): a fingerprint that is a perfect subset of another is
folded into its superset and does not count as a distinct cluster. When a subset
has more than one maximal superset the --ambiguous rule decides:
    most_frequent  assign to the most abundant superset (ties broken by name)
    all            assign to every superset
    drop           leave the subset uncollapsed
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

_FAM = re.compile(r"^(GH|PL|CE)\d+(_\d+)?$")
_S1 = re.compile(r"^S1(_\d+)?$")
_PRIORITY = ["GH", "PL", "CE", "S1"]
_RANK = {pre: i for i, pre in enumerate(_PRIORITY)}
_SPLIT = re.compile(r"[+|;,]")

_PHYLA = {"Actinobacteriota", "Bacteroidota", "Chloroflexota", "Cyanobacteria",
          "Desulfobacterota", "Firmicutes", "Planctomycetota", "Verrucomicrobiota"}
_CLASSES = {"Alphaproteobacteria", "Gammaproteobacteria"}


def _rank(token):
    return _RANK.get(token[:2], len(_PRIORITY))


def gene_token(annotation):
    """Highest-priority family token carried by one gene, or None."""
    best = None
    best_key = None
    for part in _SPLIT.split(annotation):
        part = part.strip()
        if not part:
            continue
        if _FAM.match(part) or _S1.match(part):
            key = (_rank(part), part)
            if best_key is None or key < best_key:
                best, best_key = part, key
    return best


def _short(name):
    return name.split("_", 1)[0] if "_" in name else name


def lineage_of(taxonomy):
    """Collapse a GTDB taxonomy string onto the mpCGC lineage groups."""
    names = [t.strip()[3:] for t in taxonomy.split(";")
             if len(t.strip()) >= 4 and t.strip()[1:3] == "__"]
    if not names:
        return "other"
    if names[0].startswith("Archaea"):
        return "Archaea"
    if len(names) >= 2 and _short(names[1]) in _PHYLA:
        return _short(names[1])
    if len(names) >= 3 and names[2] in _CLASSES:
        return names[2]
    return "other"


def load_groups(path, group_by):
    """Return genome -> lineage label, read from a metadata TSV."""
    groups = {}
    if path is None or not path.exists():
        return groups
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr, None)
        if head is None:
            return groups
        idx = {h.strip(): i for i, h in enumerate(head)}
        id_col = next((idx[k] for k in ("Bin_id", "sample", "MAG", "genome")
                       if k in idx), 0)
        tax_col = next((i for h, i in idx.items() if "taxonomy" in h.lower()), None)
        grp_col = idx.get(group_by)
        for row in rdr:
            if len(row) <= id_col or not row[id_col]:
                continue
            if grp_col is not None and len(row) > grp_col and row[grp_col]:
                groups[row[id_col]] = row[grp_col]
            elif tax_col is not None and len(row) > tax_col:
                groups[row[id_col]] = lineage_of(row[tax_col])
    return groups


def collapse_maximal(fingerprints, rule):
    """Map each fingerprint onto its maximal superset(s).

    Returns (assignment, maximal) where assignment maps a fingerprint to the
    list of fingerprints it should be counted as.
    """
    freq = Counter(fingerprints)
    uniq = sorted(set(fingerprints), key=lambda s: (-len(s), sorted(s)))

    by_size = defaultdict(list)
    for s in uniq:
        by_size[len(s)].append(s)

    # token -> maximal sets containing it, so a subset only has to be compared
    # against the sets that share its rarest token
    token_index = defaultdict(list)
    maximal = set()
    assignment = {}

    for size in sorted(by_size, reverse=True):
        for s in by_size[size]:
            rarest = min(s, key=lambda t: len(token_index[t]))
            supers = [c for c in token_index[rarest] if len(c) > len(s) and s < c]
            if not supers or rule == "drop":
                maximal.add(s)
                assignment[s] = [s]
                for tok in s:
                    token_index[tok].append(s)
            elif rule == "all":
                assignment[s] = sorted(supers, key=sorted)
            else:  # most_frequent
                best = max(supers, key=lambda c: (freq[c], -len(c),
                                                  tuple(sorted(c, reverse=True))))
                assignment[s] = [best]
    return assignment, maximal


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--metadata", type=Path, default=None)
    ap.add_argument("--group-by", default="phylum")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--collapse-map", required=True, type=Path)
    ap.add_argument("--collapse", action="store_true")
    ap.add_argument("--ambiguous", default="most_frequent",
                    choices=["most_frequent", "all", "drop"])
    args = ap.parse_args()

    groups = load_groups(args.metadata, args.group_by)

    cgc_tokens = defaultdict(set)
    cgc_mag = {}
    with args.catalog.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)
        for row in rdr:
            if len(row) < 9:
                continue
            cgc, mag, ann = row[0], row[1], row[8]
            cgc_mag[cgc] = mag
            token = gene_token(ann)
            if token:
                cgc_tokens[cgc].add(token)

    records = [(cgc, cgc_mag[cgc], frozenset(toks))
               for cgc, toks in cgc_tokens.items() if toks]
    records.sort()

    n_distinct = len({r[2] for r in records})
    if args.collapse:
        assignment, maximal = collapse_maximal([r[2] for r in records],
                                               args.ambiguous)
        print("[fingerprint] {:,} distinct -> {:,} maximal".format(
            n_distinct, len(maximal)), file=sys.stderr)
    else:
        assignment = {r[2]: [r[2]] for r in records}
        maximal = set(assignment)

    with args.out.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["CGC", "MAG", "lineage", "n_tokens", "fingerprint",
                    "counted_as"])
        for cgc, mag, fp in records:
            joined = ",".join(sorted(fp))
            for target in assignment[fp]:
                w.writerow([cgc, mag, groups.get(mag, "other"), len(fp),
                            joined, ",".join(sorted(target))])

    with args.collapse_map.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["fingerprint", "is_maximal", "counted_as"])
        for fp in sorted(assignment, key=lambda s: (-len(s), sorted(s))):
            w.writerow([",".join(sorted(fp)), int(fp in maximal),
                        ";".join(",".join(sorted(t)) for t in assignment[fp])])

    print("[fingerprint] {:,} CGCs encoded across {:,} genomes".format(
        len(records), len({r[1] for r in records})), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
