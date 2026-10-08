"""Gaussian constraint-bit estimate for symmetric slices over an abelian group (same model as the repo's
scripts/b1_multiplier_enum.py::cost, generalised to any G and symmetry group S).

For each row sum a used by the problem, S-invariant +-1 sequences with that row sum are sampled
exactly uniformly (DP over orbit sizes), and the mean and covariance of their PAF vector on the lag-class
representatives are measured.  The sum over blocks is treated as Gaussian on a lattice of spacing 4
(PAF = |G| mod 4), with one dimension removed (implied by the row sums):
   n     = -log2 P(sum PAF = target)                     [constraint bits]
   log2E = sum_k log2 N(a_k) - n                         [expected ordered solutions]
kinds: 'b1' (|G| even: row sums (2,0,0,0), target -4), 'lp' (row sums (1,1), target -2).
"""
from math import log, log2, pi

import numpy as np


def rowsum_dp(sizes):
    """counts[i][s+off] = number of sign choices on orbits i.. with sum s (float, exact for moderate r)"""
    r = len(sizes)
    tot = int(sum(sizes))
    off = tot
    cnt = np.zeros((r + 1, 2 * tot + 1))
    cnt[r, off] = 1.0
    for i in range(r - 1, -1, -1):
        s = sizes[i]
        cnt[i, s:] += cnt[i + 1, :2 * tot + 1 - s]
        cnt[i, :2 * tot + 1 - s] += cnt[i + 1, s:]
    return cnt, off


def sample_rowsum(sizes, a, ns, rng, cnt=None, off=None):
    if cnt is None:
        cnt, off = rowsum_dp(sizes)
    r = len(sizes)
    if cnt[0, a + off] == 0:
        return None
    out = np.empty((ns, r), dtype=np.int8)
    for t in range(ns):
        rem = a
        for i in range(r):
            s = sizes[i]
            w_plus = cnt[i + 1, rem - s + off] if 0 <= rem - s + off < cnt.shape[1] else 0.0
            w_minus = cnt[i + 1, rem + s + off] if 0 <= rem + s + off < cnt.shape[1] else 0.0
            if rng.random() * (w_plus + w_minus) < w_plus:
                out[t, i] = 1
                rem -= s
            else:
                out[t, i] = -1
                rem += s
    return out


def paf_stats(P, Y):
    G = P['G']
    X = Y[:, P['orb']].astype(float).reshape((len(Y),) + G.mods)
    ax = tuple(range(1, len(G.mods) + 1))
    F = np.fft.fftn(X, axes=ax)
    R = np.fft.ifftn(np.abs(F) ** 2, axes=ax).real.reshape(len(Y), -1)
    Pf = np.rint(R[:, P['lrep']])
    return Pf.mean(0), np.cov(Pf.T).reshape(len(P['lrep']), len(P['lrep']))


def gauss_cost(P, kind, ns=3000, seed=0):
    rng = np.random.default_rng(seed)
    sizes = [int(s) for s in P['sizes']]
    cnt, off = rowsum_dp(sizes)
    if kind == 'b1':
        assert P['n'] % 2 == 0
        profile = [2, 0, 0, 0]
        target = -4.0
    elif kind == 'lp':
        profile = [1, 1]
        target = -2.0
    else:
        raise ValueError(kind)
    stats = {}
    logN = 0.0
    for a in profile:
        if a not in stats:
            if cnt[0, a + off] == 0:
                return dict(feasible=False)
            Y = sample_rowsum(sizes, a, ns, rng, cnt, off)
            stats[a] = paf_stats(P, Y)
        logN += log2(cnt[0, a + off])
    mu = sum(stats[a][0] for a in profile)
    C = sum(stats[a][1] for a in profile)
    w, V = np.linalg.eigh(np.atleast_2d(C))
    keep = w > 1e-6 * max(w.max(), 1e-12)
    proj = V.T @ (np.full(len(mu), target) - mu)
    null_mis = float(np.abs(proj[~keep]).max()) if (~keep).any() else 0.0
    D = int(keep.sum())  # nondegenerate directions (the row-sum direction is already degenerate)
    idx = np.argsort(-w)[:max(D, 0)]
    lp = D * log(4) - 0.5 * np.sum(np.log(2 * pi * w[idx])) - 0.5 * np.sum(proj[idx] ** 2 / w[idx])
    n = -lp / log(2)
    return dict(feasible=True, D=D, n=round(float(n), 1), log2E=round(float(logN - n), 1),
                null_mismatch=round(null_mis, 3))
