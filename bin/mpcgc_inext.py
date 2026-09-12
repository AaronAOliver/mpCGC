#!/usr/bin/env python3
r"""Sample-based (incidence) species richness interpolation and extrapolation.

Vectorised NumPy implementation of the q=0 estimators of the iNEXT framework
(Colwell et al. 2012; Chao et al. 2014), applied here with genomes as the
sampling units and distinct CGC fingerprints as the species.

Interpolation (t <= T):
    S(t) = sum_j Q_j * [1 - C(T-j, t) / C(T, t)]

Extrapolation (t = T + t*):
    S(T + t*) = S_obs + Q0 * [1 - (1 - Q1 / (T*Q0 + Q1))**t*]

Chao2 asymptotic richness:
    Q0 = (T-1)/T * Q1**2 / (2*Q2)          when Q2 > 0
    Q0 = (T-1)/T * Q1*(Q1-1) / 2           when Q2 == 0

Two interval methods are provided:
    chao2_ci   analytic variance of Chao2 with a log-transformed interval, which
               is what iNEXT::ChaoRichness reports and is asymmetric about the
               point estimate
    inext_ci   parametric bootstrap SE on an estimated detection-probability
               community, for the interpolated/extrapolated curve

Point estimates and analytic asymptotes from this module were checked against
iNEXT 3.0.2 / ChaoRichness in R.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammaln
from scipy.stats import norm

Z = norm.ppf(0.975)   # 1.959964


def _lchoose(n, k):
    n = np.asarray(n, float)
    k = np.asarray(k, float)
    out = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    return np.where((k < 0) | (k > n), -np.inf, out)


def chao2_Q0(T, Q1, Q2):
    """Estimated number of undetected fingerprints (the Chao2 correction)."""
    if Q2 > 0:
        return (T - 1) / T * Q1 * Q1 / (2 * Q2)
    return (T - 1) / T * Q1 * (Q1 - 1) / 2


def inext_point(T, y, t_grid):
    """Rarefaction / extrapolation curve.

    T       number of sampling units (genomes)
    y       incidence counts per fingerprint (Y_i >= 1)
    t_grid  sampling-unit counts at which to evaluate

    Returns (curve, S_obs, Q1, Q2, Q0).
    """
    y = np.asarray(y, float)
    S_obs = int((y > 0).sum())
    Q1 = int((y == 1).sum())
    Q2 = int((y == 2).sum())
    j = np.arange(1, T + 1)
    Qj = np.bincount(y.astype(int), minlength=T + 1)[1:T + 1].astype(float)
    Q0 = chao2_Q0(T, Q1, Q2)

    t_grid = np.asarray(t_grid)
    out = np.empty(len(t_grid), float)
    lchoose_T = _lchoose(T, t_grid)
    for i, t in enumerate(t_grid):
        if t <= T:
            ratio = np.exp(_lchoose(T - j, t) - lchoose_T[i])
            ratio = np.where(np.isfinite(ratio), ratio, 0.0)
            out[i] = float(np.sum(Qj * (1.0 - ratio)))
        else:
            tstar = t - T
            if Q1 == 0 or Q0 <= 0:
                out[i] = S_obs
            else:
                out[i] = S_obs + Q0 * (1.0 - (1.0 - Q1 / (T * Q0 + Q1)) ** tstar)
    return out, S_obs, Q1, Q2, Q0


def chao2_ci(T, Q1, Q2, S_obs):
    """Chao2 point estimate, SE and log-transformed 95% interval.

    Variance follows Chao (1987) for incidence data; the interval is the
    log-transformed form used by iNEXT, so it is asymmetric and never falls
    below S_obs.
    """
    k = (T - 1) / T
    Q0 = chao2_Q0(T, Q1, Q2)
    S_est = S_obs + Q0

    if Q2 > 0:
        r = Q1 / Q2
        var = Q2 * (0.5 * k * r ** 2 + k ** 2 * r ** 3 + 0.25 * k ** 2 * r ** 4)
    elif Q1 > 1:
        var = (0.25 * k ** 2 * Q1 * (2 * Q1 - 1) ** 2
               + 0.5 * k * Q1 * (Q1 - 1) - 0.25 * k ** 2 * Q1 ** 4 / max(S_est, 1e-12))
    else:
        var = 0.0
    var = max(var, 0.0)
    se = float(np.sqrt(var))

    if Q0 > 0 and var > 0:
        ratio = float(np.exp(Z * np.sqrt(np.log1p(var / (Q0 ** 2)))))
        lower = S_obs + Q0 / ratio
        upper = S_obs + Q0 * ratio
    else:
        lower = upper = S_est
    return S_est, se, lower, upper


def _boot_probs(T, y):
    """Detection probabilities of the iNEXT incidence bootstrap community."""
    y = np.asarray(y, float)
    U = y.sum()
    Q1 = (y == 1).sum()
    Q2 = (y == 2).sum()
    Q0 = chao2_Q0(T, Q1, Q2)

    if U > 0 and ((T - 1) * Q1 + 2 * Q2) > 0:
        C_hat = 1 - (Q1 / U) * ((T - 1) * Q1 / ((T - 1) * Q1 + 2 * Q2))
    else:
        C_hat = 1.0

    p = y / T
    denom = np.sum(p * (1 - p) ** T)
    w = (1 - C_hat) / denom if denom > 0 else 0.0
    prob_obs = p * (1 - w * (1 - p) ** T)

    n0 = int(np.ceil(Q0))
    a = (U / T) * (1 - C_hat)
    prob_uns = np.full(n0, a / n0) if n0 > 0 else np.array([])
    return np.concatenate([prob_obs, prob_uns])


def inext_ci(T, y, t_grid, nboot=50, seed=0):
    """Bootstrap standard error of the rarefaction / extrapolation curve."""
    rng = np.random.default_rng(seed)
    probs = _boot_probs(T, y)
    boot = np.empty((nboot, len(t_grid)), float)
    for b in range(nboot):
        Yb = rng.binomial(T, probs)
        Yb = Yb[Yb > 0]
        boot[b], *_ = inext_point(T, Yb, t_grid)
    return boot.std(axis=0, ddof=1)
