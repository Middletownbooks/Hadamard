"""Gaussian cost of a b=0 (Goethals-Seidel / Williamson) slice with a fixed row-sum profile.
usage: cost_b0.py p unit a1,a2,a3,a4      (G = Z_p, S = <unit>)"""
import sys
from math import log, log2, pi
import numpy as np
from abgroup import AbGroup, AutGroup
from slicesolve import prepare
from cost import rowsum_dp, sample_rowsum, paf_stats
p = int(sys.argv[1]); u = int(sys.argv[2]); prof = [int(x) for x in sys.argv[3].split(',')]
G = AbGroup([('C', p)]); AG = AutGroup(G)
auts = [[u]]
P = prepare(G, auts)
sizes = [int(s) for s in P['sizes']]
cnt, off = rowsum_dp(sizes)
rng = np.random.default_rng(0)
stats = {}; logN = 0.0
for a in prof:
    if a not in stats:
        assert cnt[0, a + off] > 0, f'row sum {a} infeasible'
        stats[a] = paf_stats(P, sample_rowsum(sizes, a, 3000, rng, cnt, off))
    logN += log2(cnt[0, a + off])
mu = sum(stats[a][0] for a in prof); C = sum(stats[a][1] for a in prof)
w, V = np.linalg.eigh(np.atleast_2d(C)); keep = w > 1e-6 * w.max()
proj = V.T @ (0.0 - mu)
D = int(keep.sum()); idx = np.argsort(-w)[:D]
# lattice spacing: estimate from the data (gcd of PAF-sum differences)
lp = D * log(4) - 0.5 * np.sum(np.log(2 * pi * w[idx])) - 0.5 * np.sum(proj[idx] ** 2 / w[idx])
n = -lp / log(2)
print(f'p={p} S=<{u}> |S|={len(set(pow(u,k,p) for k in range(p)))} r={P["r"]} L={len(P["lrep"])} profile={prof}: n={n:.1f} bits, sum log2 N(a)={logN:.1f}, log2E(ordered, this profile)={logN-n:.1f}, null_mismatch={float(np.abs(proj[~keep]).max()) if (~keep).any() else 0:.3f}')
