"""Solver fingerprints of the b = 0 cores (published and ours), at the level of single blocks.

For each block of each b = 0 core (cyclic Z_v, v prime, U = common multiplier stabiliser, fixed point moved to 0):
  * flatness: mean PAF^2 over the nonzero lags, divided by v;
  * the fraction of lag classes whose PAF has the least possible magnitude (PAF = v mod 4, so {-1, 3} for v = 3 mod 4,
    {1, -3} for v = 1 mod 4);
  * local minimality: the number of row-sum-preserving swaps (flip one +1 U-orbit and one -1 U-orbit) that lower the
    block's own PAF energy sum_{t != 0} PAF(t)^2.  A block taken from a pool built by per-block descent on that energy
    is a local minimum (0 improving swaps); a block of a solution found without such pools need not be.
  * pair sums: for each pairing (AB|CD, AC|BD, AD|BC), the number of lag classes with |PAF_A + PAF_B| = 2 (the least
    possible value, since every PAF is = v mod 4) and the rms of the pair sum.
usage: python3 forensics/b0_fingerprints.py"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import near_structure as N
from census_originals import stabiliser


def paf(x):
    F = np.fft.fft(x.astype(float)); return np.rint(np.fft.ifft(np.abs(F) ** 2).real).astype(np.int64)


def main():
    for c in [c for c in N.load_all() if c['b'] == 0]:
        X = c['X']; v = c['g'].m
        d, fp, U = stabiliser(X, c['g']); X = X[:, (np.arange(v) + fp) % v]
        orb = -np.ones(v, np.int64); reps = []
        for y in range(v):
            if orb[y] < 0:
                o = sorted(set((u * y) % v for u in U)); orb[o] = len(reps); reps.append(o)
        seen = np.zeros(v, bool); lreps = []
        for t in range(1, v):
            if not seen[t]:
                o = sorted(set((u * t) % v for u in U) | set((-u * t) % v for u in U)); seen[o] = True; lreps.append(o[0])
        lreps = np.array(lreps); P = np.array([paf(x) for x in X]); Pc = P[:, lreps]
        minvals = [-1, 3] if v % 4 == 3 else [1, -3]
        flat = [round(float((P[k, 1:] ** 2).mean()) / v, 3) for k in range(4)]
        fmin = [round(float(np.isin(Pc[k], minvals).mean()), 2) for k in range(4)]
        loc = []
        sz = np.array([len(o) for o in reps])
        for k in range(4):
            x = X[k]; E0 = float((P[k, 1:] ** 2).sum()); lab = np.array([x[o[0]] for o in reps])
            plus = [i for i in range(len(reps)) if lab[i] > 0 and sz[i] > 1]; minus = [i for i in range(len(reps)) if lab[i] < 0 and sz[i] > 1]
            better = 0
            for i in plus:
                for j in minus:
                    y = x.copy(); y[reps[i]] *= -1; y[reps[j]] *= -1
                    if float((paf(y)[1:] ** 2).sum()) < E0 - 1e-6: better += 1
            loc.append(f'{better}/{len(plus) * len(minus)}')
        pairs = []
        for name, (i, j) in [('AB', (0, 1)), ('AC', (0, 2)), ('AD', (0, 3))]:
            s = Pc[i] + Pc[j]; pairs.append(f'{name} {int((np.abs(s) == 2).sum())}/{len(lreps)} rms {float(np.sqrt((s ** 2).mean())):.1f}')
        print(f"{c['label']:>16} v={v} |U|={d} rows={X.sum(1).tolist()} | mean PAF^2/v {flat} | min-|PAF| fraction {fmin} | "
              f"improving swaps {loc} | pair sums: {'; '.join(pairs)}", flush=True)


if __name__ == '__main__':
    main()
