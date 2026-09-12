#!/usr/bin/env python3
r"""Compare an mpCGC catalog against a reference catalog, gene by gene.

Written to check a pipeline run against the published mpCGCdb record for the same
genome, but it will compare any two catalogs sharing the mpCGC schema.

Genes are matched on (contig, start, stop) rather than on protein identifier,
because identifiers are only unique within the run that produced them (mpCGCdb
prefixes them with the genome name, a fresh run does not). Coordinates come
straight from the gene caller, so identical calls match exactly.

Reported:
  * cluster counts and gene-row counts on each side
  * per-gene agreement on Gene Type and on the enzyme families in Gene Annotation
  * cluster-level agreement, pairing each reference cluster with the run cluster
    it shares the most genes with and reporting the Jaccard overlap
  * genes present on only one side, split by whether they are signature or null

Usage:
    validate_against_mpcgcdb.py \
        --run results/catalog/all_cgc_catalog.tsv \
        --reference test/reference/mpcgcdb_GCA_000325705.1_ASM32570v1_genomic.tsv.gz \
        --genome GCA_000325705.1_ASM32570v1_genomic

Either file may be gzipped.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

_FAM = re.compile(r"\b((?:GH|PL|CE|CBM|AA|GT)\d+(?:_\d+)?|S1(?:_\d+)?)\b")


def _open(path):
    """Open a catalog, transparently handling gzip."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="", encoding="utf-8")
    return path.open(newline="", encoding="utf-8")


def load(path, genome=None):
    """(contig, start, stop) -> row, plus cluster -> set of those keys."""
    genes = {}
    clusters = defaultdict(set)
    with _open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        for row in rdr:
            if len(row) < 9 or row[0] in ("CGC#", ""):
                continue
            cgc, mag, gtype, contig, pid, start, stop, strand, ann = row[:9]
            if genome and mag != genome:
                continue
            try:
                key = (contig, int(start), int(stop))
            except ValueError:
                continue
            genes[key] = {"cgc": cgc, "type": gtype, "protein": pid,
                          "strand": strand, "annotation": ann}
            clusters[cgc].add(key)
    return genes, clusters


def families(annotation):
    return frozenset(_FAM.findall(annotation or ""))


def pct(n, d):
    return "{:.1f}%".format(100.0 * n / d) if d else "n/a"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--genome", default=None)
    ap.add_argument("--report", type=Path, default=None)
    ap.add_argument("--max-examples", type=int, default=10)
    args = ap.parse_args()

    run_genes, run_clusters = load(args.run, args.genome)
    ref_genes, ref_clusters = load(args.reference, args.genome)
    if not run_genes:
        sys.exit("no genes loaded from {}".format(args.run))
    if not ref_genes:
        sys.exit("no genes loaded from {}".format(args.reference))

    out = []
    def say(line=""):
        out.append(line)
        print(line)

    say("=" * 72)
    say("mpCGC catalog comparison" + ("  [{}]".format(args.genome) if args.genome else ""))
    say("=" * 72)
    say("{:<34} {:>14} {:>14}".format("", "run", "reference"))
    say("{:<34} {:>14,} {:>14,}".format("clusters", len(run_clusters), len(ref_clusters)))
    say("{:<34} {:>14,} {:>14,}".format("gene rows", len(run_genes), len(ref_genes)))
    for t in ("CAZyme", "TC", "TF", "STP", "Sulfatase", "Peptidase", "prodoric", "null"):
        a = sum(1 for g in run_genes.values() if g["type"] == t)
        b = sum(1 for g in ref_genes.values() if g["type"] == t)
        if a or b:
            say("{:<34} {:>14,} {:>14,}".format("  " + t, a, b))

    shared = set(run_genes) & set(ref_genes)
    only_run = set(run_genes) - set(ref_genes)
    only_ref = set(ref_genes) - set(run_genes)

    say("")
    say("-" * 72)
    say("GENE-LEVEL AGREEMENT (matched on contig + coordinates)")
    say("-" * 72)
    say("genes in both            : {:,}  ({} of reference)".format(
        len(shared), pct(len(shared), len(ref_genes))))
    say("only in run              : {:,}".format(len(only_run)))
    say("only in reference        : {:,}".format(len(only_ref)))

    # signature vs null among the non-shared genes
    for label, keys, src in (("only in run", only_run, run_genes),
                             ("only in reference", only_ref, ref_genes)):
        if keys:
            sig = sum(1 for k in keys if src[k]["type"] != "null")
            say("  {:<22}: {:,} signature, {:,} null".format(
                label, sig, len(keys) - sig))

    type_same = sum(1 for k in shared
                    if run_genes[k]["type"] == ref_genes[k]["type"])
    fam_same = sum(1 for k in shared
                   if families(run_genes[k]["annotation"])
                   == families(ref_genes[k]["annotation"]))
    ann_same = sum(1 for k in shared
                   if run_genes[k]["annotation"] == ref_genes[k]["annotation"])
    say("")
    say("of the {:,} shared genes:".format(len(shared)))
    say("  identical Gene Type          : {:,}  ({})".format(
        type_same, pct(type_same, len(shared))))
    say("  identical enzyme families    : {:,}  ({})".format(
        fam_same, pct(fam_same, len(shared))))
    say("  byte-identical annotation    : {:,}  ({})".format(
        ann_same, pct(ann_same, len(shared))))

    mism_type = Counter()
    for k in shared:
        a, b = run_genes[k]["type"], ref_genes[k]["type"]
        if a != b:
            mism_type[(b, a)] += 1
    if mism_type:
        say("")
        say("Gene Type differences (reference -> run):")
        for (b, a), n in mism_type.most_common(args.max_examples):
            say("  {:>12} -> {:<12} {:,}".format(b, a, n))

    mism_fam = []
    for k in sorted(shared):
        fa, fb = (families(run_genes[k]["annotation"]),
                  families(ref_genes[k]["annotation"]))
        if fa != fb:
            mism_fam.append((k, sorted(fb), sorted(fa)))
    if mism_fam:
        say("")
        say("enzyme-family differences: {:,}".format(len(mism_fam)))
        for k, fb, fa in mism_fam[:args.max_examples]:
            say("  {}:{}-{}   reference {} | run {}".format(
                k[0], k[1], k[2], ",".join(fb) or "-", ",".join(fa) or "-"))

    # ---- cluster level -------------------------------------------------
    say("")
    say("-" * 72)
    say("CLUSTER-LEVEL AGREEMENT")
    say("-" * 72)
    gene_to_run = {}
    for cgc, keys in run_clusters.items():
        for k in keys:
            gene_to_run[k] = cgc

    exact = 0
    jaccs = []
    unmatched = 0
    merges = Counter()
    for ref_cgc, ref_keys in ref_clusters.items():
        hits = Counter(gene_to_run[k] for k in ref_keys if k in gene_to_run)
        if not hits:
            unmatched += 1
            continue
        best, n = hits.most_common(1)[0]
        union = len(ref_keys | run_clusters[best])
        j = n / union if union else 0.0
        jaccs.append(j)
        if ref_keys == run_clusters[best]:
            exact += 1
        merges[len(hits)] += 1

    say("reference clusters              : {:,}".format(len(ref_clusters)))
    say("  gene sets identical           : {:,}  ({})".format(
        exact, pct(exact, len(ref_clusters))))
    if jaccs:
        jaccs.sort()
        mean = sum(jaccs) / len(jaccs)
        say("  mean Jaccard with best match  : {:.3f}".format(mean))
        say("  median Jaccard                : {:.3f}".format(
            jaccs[len(jaccs) // 2]))
        say("  >=0.90 overlap                : {:,}  ({})".format(
            sum(1 for j in jaccs if j >= 0.90),
            pct(sum(1 for j in jaccs if j >= 0.90), len(jaccs))))
        say("  >=0.50 overlap                : {:,}  ({})".format(
            sum(1 for j in jaccs if j >= 0.50),
            pct(sum(1 for j in jaccs if j >= 0.50), len(jaccs))))
    if unmatched:
        say("  no overlapping gene at all    : {:,}".format(unmatched))
    if merges:
        split = {k: v for k, v in merges.items() if k > 1}
        if split:
            say("  reference clusters split across several run clusters:")
            for k in sorted(split):
                say("    into {} run clusters          : {:,}".format(k, split[k]))

    # families recovered overall, which is what most downstream analysis uses
    run_fams = Counter()
    ref_fams = Counter()
    for g in run_genes.values():
        run_fams.update(families(g["annotation"]))
    for g in ref_genes.values():
        ref_fams.update(families(g["annotation"]))
    both = set(run_fams) & set(ref_fams)
    say("")
    say("-" * 72)
    say("ENZYME FAMILY INVENTORY")
    say("-" * 72)
    say("families in run                 : {:,}".format(len(run_fams)))
    say("families in reference           : {:,}".format(len(ref_fams)))
    say("families in both                : {:,}  ({} of reference)".format(
        len(both), pct(len(both), len(ref_fams))))
    miss = sorted(set(ref_fams) - set(run_fams))
    extra = sorted(set(run_fams) - set(ref_fams))
    if miss:
        say("only in reference               : " + ", ".join(miss[:20])
            + (" ..." if len(miss) > 20 else ""))
    if extra:
        say("only in run                     : " + ", ".join(extra[:20])
            + (" ..." if len(extra) > 20 else ""))
    same_count = sum(1 for f in both if run_fams[f] == ref_fams[f])
    say("identical gene counts           : {:,} of {:,} shared families".format(
        same_count, len(both)))

    if args.report:
        args.report.write_text("\n".join(out) + "\n", encoding="utf-8")
        print("\nwrote {}".format(args.report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
