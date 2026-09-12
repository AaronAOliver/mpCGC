#!/usr/bin/env python3
r"""Check mpcgc_inext against iNEXT in R.

Expected values were produced by iNEXT 3.0.2 / ChaoRichness(datatype =
"incidence_freq") under R 4.5.2, so this runs without R installed. Regenerate
them with test/chao_ref.R if the estimators are ever changed.

    python test/test_inext.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from mpcgc_inext import chao2_ci, chao2_Q0, inext_point  # noqa: E402

# incidence counts per species, and T = number of sampling units
DATASETS = {
    "low": (40, [1] * 30 + [2] * 12 + [3] * 6 + [5] * 3 + [9, 14, 22]),
    "mid": (120, [1] * 45 + [2] * 30 + [3] * 20 + [6] * 14 + [11] * 8
            + [25] * 5 + [60, 88]),
    "high": (500, [1] * 20 + [2] * 25 + [4] * 40 + [9] * 50 + [30] * 35
             + [120] * 18 + [300, 410, 455]),
}

# iNEXT 3.0.2 ChaoRichness: S_obs, Q1, Q2, S_chao2, SE, ci_lower, ci_upper
R_REFERENCE = {
    "low": (54, 30, 12, 90.5620, 18.0610, 68.6320, 145.3650),
    "mid": (124, 45, 30, 157.4690, 13.0530, 140.0080, 193.9730),
    "high": (191, 20, 25, 198.9840, 4.8250, 193.6750, 214.8330),
}

TOL = 0.02   # absolute tolerance against R's 4-decimal output


def main() -> int:
    failures = []
    print("{:<6} {:<12} {:>12} {:>12} {:>10}".format(
        "set", "quantity", "mpCGC", "R iNEXT", "diff"))
    print("-" * 56)

    for name, (T, y) in DATASETS.items():
        y = np.array(y, float)
        exp_sobs, exp_q1, exp_q2, exp_s, exp_se, exp_lo, exp_hi = R_REFERENCE[name]

        _curve, S_obs, Q1, Q2, _Q0 = inext_point(T, y, [T])
        S_est, se, lo, hi = chao2_ci(T, Q1, Q2, S_obs)

        checks = [
            ("S_obs", S_obs, exp_sobs),
            ("Q1", Q1, exp_q1),
            ("Q2", Q2, exp_q2),
            ("S_chao2", S_est, exp_s),
            ("SE", se, exp_se),
            ("ci_lower", lo, exp_lo),
            ("ci_upper", hi, exp_hi),
        ]
        for label, got, want in checks:
            diff = abs(got - want)
            ok = diff <= TOL
            if not ok:
                failures.append((name, label, got, want))
            print("{:<6} {:<12} {:>12.4f} {:>12.4f} {:>10.4f} {}".format(
                name, label, float(got), float(want), diff,
                "" if ok else "  <-- MISMATCH"))
        print("-" * 56)

    # structural properties that must hold for any input
    for name, (T, y) in DATASETS.items():
        y = np.array(y, float)
        _c, S_obs, Q1, Q2, Q0 = inext_point(T, y, [T])
        S_est, _se, lo, hi = chao2_ci(T, Q1, Q2, S_obs)
        if not lo >= S_obs:
            failures.append((name, "ci_lower >= S_obs", lo, S_obs))
        if not lo <= S_est <= hi:
            failures.append((name, "estimate inside interval", S_est, (lo, hi)))
        if abs(S_est - (S_obs + Q0)) > 1e-9:
            failures.append((name, "S_chao2 == S_obs + Q0", S_est, S_obs + Q0))
        # the interval is asymmetric by construction, not estimate +/- 1.96 SE
        if abs((S_est - lo) - (hi - S_est)) < 1e-6:
            failures.append((name, "interval should be asymmetric", lo, hi))

    # rarefaction must pass through S_obs at t = T and be monotonic
    for name, (T, y) in DATASETS.items():
        y = np.array(y, float)
        grid = np.array([1, T // 4, T // 2, T, 2 * T, 4 * T])
        curve, S_obs, Q1, Q2, Q0 = inext_point(T, y, grid)
        if abs(curve[3] - S_obs) > 1e-6:
            failures.append((name, "curve(T) == S_obs", curve[3], S_obs))
        if not np.all(np.diff(curve) >= -1e-9):
            failures.append((name, "curve is monotonic", curve, None))
        if curve[-1] > S_obs + Q0 + 1e-6:
            failures.append((name, "extrapolation <= Chao2 asymptote",
                             curve[-1], S_obs + Q0))

    if failures:
        print("\n{} FAILURE(S):".format(len(failures)))
        for f in failures:
            print("  {}: {} got {!r}, expected {!r}".format(*f))
        return 1

    print("\nall checks passed: point estimates, analytic Chao2 interval, "
          "monotonicity and asymptote bound agree with iNEXT 3.0.2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
