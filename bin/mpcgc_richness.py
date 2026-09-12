#!/usr/bin/env python3
r"""CGC diversity: rarefaction, extrapolation and asymptotic richness per lineage.

Genomes are the sampling units and distinct CGC fingerprints are the species, so
a fingerprint is "detected" in a genome when that genome encodes at least one CGC
with that fingerprint. Repeated CGCs within a genome are not counted twice, which
is what makes incidence (rather than abundance) data the right framework.

For every lineage the script writes
  * the interpolated curve from 1 to T genomes and the extrapolated curve out to
    --extrapolate x T, with 95% confidence bands
  * the Chao2 asymptote with its analytic log-transformed interval

Curves and asymptotes use the estimators in mpcgc_inext.py.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mpcgc_inext import Z, chao2_ci, inext_ci, inext_point  # noqa: E402

# 11-group mpCGC palette plus the grey catch-all
GROUP_COLORS = {
    "Actinobacteriota": "#1f77b4", "Alphaproteobacteria": "#ff7f0e",
    "Archaea": "#2ca02c", "Bacteroidota": "#d62728", "Chloroflexota": "#9467bd",
    "Cyanobacteria": "#8c564b", "Desulfobacterota": "#e377c2",
    "Firmicutes": "#7f7f7f", "Gammaproteobacteria": "#bcbd22",
    "Planctomycetota": "#17becf", "Verrucomicrobiota": "#1a3f70",
    "other": "#c9c9c9",
}


def build_grid(T, extrapolate, knots):
    """Knots from 1 to T, then on to extrapolate*T, always including T."""
    end = max(int(round(T * extrapolate)), T)
    n_in = max(int(knots) // 2, 2)
    interp = np.unique(np.linspace(1, T, n_in, dtype=int))
    if end > T:
        n_ex = max(int(knots) - len(interp), 2)
        extrap = np.unique(np.linspace(T, end, n_ex, dtype=int))
    else:
        extrap = np.array([T], dtype=int)
    return np.unique(np.concatenate([interp, extrap, [T]]))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fingerprints", required=True, type=Path)
    ap.add_argument("--knots", type=int, default=40)
    ap.add_argument("--extrapolate", type=float, default=2.0)
    ap.add_argument("--out-curves", required=True, type=Path)
    ap.add_argument("--out-asymptotes", required=True, type=Path)
    ap.add_argument("--figure", type=Path, default=None)
    ap.add_argument("--ci", action="store_true")
    ap.add_argument("--nboot", type=int, default=0,
                    help="0 uses the analytic Chao2 interval for the asymptote "
                         "and skips bootstrap bands")
    ap.add_argument("--min-genomes", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    # lineage -> fingerprint -> set of genomes encoding it
    incidence = defaultdict(lambda: defaultdict(set))
    genomes = defaultdict(set)
    with args.fingerprints.open(newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr)
        i_mag = head.index("MAG")
        i_grp = head.index("lineage")
        i_fp = head.index("counted_as")
        for row in rdr:
            if len(row) <= max(i_mag, i_grp, i_fp):
                continue
            lineage, mag, fp = row[i_grp], row[i_mag], row[i_fp]
            if not fp:
                continue
            incidence[lineage][fp].add(mag)
            genomes[lineage].add(mag)

    curves = []
    asymptotes = []
    for lineage in sorted(incidence, key=lambda g: -len(genomes[g])):
        T = len(genomes[lineage])
        if T < args.min_genomes:
            print("[richness] skipping {} ({} genomes < --min-genomes)".format(
                lineage, T), file=sys.stderr)
            continue

        y = np.array([len(m) for m in incidence[lineage].values()], float)
        grid = build_grid(T, args.extrapolate, args.knots)
        curve, S_obs, Q1, Q2, Q0 = inext_point(T, y, grid)

        se = None
        if args.ci and args.nboot > 0:
            se = inext_ci(T, y, grid, nboot=args.nboot, seed=args.seed)

        for t, val, idx in zip(grid, curve, range(len(grid))):
            lo = hi = ""
            if se is not None:
                lo = max(val - Z * se[idx], 0.0)
                hi = val + Z * se[idx]
            curves.append([lineage, int(t), "observed" if t <= T else "extrapolated",
                           round(float(val), 4),
                           round(float(lo), 4) if lo != "" else "",
                           round(float(hi), 4) if hi != "" else ""])

        S_est, se_a, lo_a, hi_a = chao2_ci(T, Q1, Q2, S_obs)
        asymptotes.append([lineage, T, S_obs, Q1, Q2,
                           round(Q0, 2), round(S_est, 1), round(se_a, 1),
                           round(lo_a, 1), round(hi_a, 1),
                           round(S_obs / S_est * 100, 1) if S_est > 0 else ""])

    with args.out_curves.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["lineage", "genomes", "method", "richness",
                    "ci_lower", "ci_upper"])
        w.writerows(curves)

    with args.out_asymptotes.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["lineage", "genomes", "S_obs", "Q1", "Q2", "Q0_chao2",
                    "S_chao2", "SE", "ci_lower", "ci_upper", "pct_recovered"])
        w.writerows(asymptotes)

    print("[richness] {} lineages; asymptotes in {}".format(
        len(asymptotes), args.out_asymptotes.name), file=sys.stderr)

    if args.figure and curves:
        try:
            plot(curves, asymptotes, args.figure)
        except Exception as exc:                      # figures are optional
            print("[richness] figure skipped: {}".format(exc), file=sys.stderr)
    return 0


def plot(curves, asymptotes, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_lineage = defaultdict(list)
    for lineage, t, method, val, lo, hi in curves:
        by_lineage[lineage].append((t, method, val, lo, hi))
    obs_T = {a[0]: a[1] for a in asymptotes}

    fig, ax = plt.subplots(figsize=(9, 6))
    for lineage, rows in by_lineage.items():
        rows.sort()
        color = GROUP_COLORS.get(lineage, "#888888")
        T = obs_T.get(lineage)
        solid = [(t, v) for t, m, v, _l, _h in rows if m == "observed"]
        dashed = [(t, v) for t, m, v, _l, _h in rows if m == "extrapolated"]
        if solid:
            ax.plot(*zip(*solid), color=color, lw=2, label=lineage)
        if dashed:
            if solid:
                dashed = [solid[-1]] + dashed
            ax.plot(*zip(*dashed), color=color, lw=2, ls="--")
        band = [(t, l, h) for t, _m, _v, l, h in rows if l != "" and h != ""]
        if band:
            t_b, lo_b, hi_b = zip(*band)
            ax.fill_between(t_b, lo_b, hi_b, color=color, alpha=0.15, lw=0)
        if T is not None:
            here = [v for t, _m, v, _l, _h in rows if t == T]
            if here:
                ax.plot([T], [here[0]], "o", color=color, ms=5, zorder=5)

    ax.set_xlabel("Number of genomes", fontsize=13)
    ax.set_ylabel("CGC richness", fontsize=13)
    ax.set_xscale("log")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(color="0.92", lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=10,
              frameon=False, title="Lineage")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig("{}.{}".format(prefix, ext),
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
