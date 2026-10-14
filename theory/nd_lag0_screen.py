r"""Lag-0 screen of a compressed stage 1 [M; the argument is P].

For the compression of a cyclic slice (Z_m, M = <mu>) along a Hall factor Z_K, every entry y(b) of every block is a
letter of its Z_F-orbit (signed sums of the M_b-orbit sizes of the fibre) and sum_k sum_b y_k(b)^2 = 4m - 4(K - 1).
This script computes the letter sets exactly (theory/nd_s1_sd.structure) and decides by a DP over the orbits (weights =
orbit sizes, 4 blocks) whether that energy is attainable at all; it also reports the minimum attainable energy.
usage: python3 theory/nd_lag0_screen.py m mu K [m mu K ...]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import structure
import numpy as np

def screen(m, mu, K):
    F, Mg, Fo, reps, letters = structure(m, mu, K)
    w = np.bincount(Fo)
    T0 = 4 * m - 4 * (K - 1)
    reach = {0}
    for _ in range(4):
        for j, L in enumerate(letters):
            sq = sorted({int(y) * int(y) * int(w[j]) for y in L})
            reach = {r + s for r in reach for s in sq if r + s <= T0}
    minE = 4 * sum(int(w[j]) * min(int(y) * int(y) for y in L) for j, L in enumerate(letters))
    return dict(m=m, v=m + 1, mu=mu, K=K, F=F, orbits=len(reps), orbit_sizes=w.tolist(),
                letters=[sorted({abs(int(y)) for y in L}) for L in letters], target=T0, min_energy=minE,
                lag0_attainable=T0 in reach)

if __name__ == "__main__":
    a = list(map(int, sys.argv[1:]))
    for i in range(0, len(a), 3):
        print(json.dumps(screen(*a[i:i + 3])), flush=True)
