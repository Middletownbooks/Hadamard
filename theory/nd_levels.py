r"""Constraint budget of stage 2 by Z_F-character level for a cyclic multiplier slice, given a stage-1 y [E].

Generalises theory/q523_levels.py (523: Z_18 x quartic Z_29) to any slice (m, mu, K) with Z_F = Z_m / Z_K.
Stage-2 configurations are drawn uniformly given y (as in theory/nd_stage2_est.py); the full PAF on Z_m is written as
P(h, z) (h in Z_F, z in Z_K, CRT) and projected onto the characters of Z_F of order d by Ramanujan sums
    P_d(z) = sum_h c_d(h) P(h, z),   c_d(h) = sum_{chi of order d} chi(h)  (integers),
(the projection keeps every h: P_d(h, z) = sum_h' c_d(h - h') P(h', z)), integer combinations of the lag sums.  A level is a set of orders d; its constraint bits are the
lattice-coordinate Gaussian (theory/q523_gauss.py) of the projected vectors (all z) at the projected target.
Conditional bits C(B | A) = bits(A u B) - bits(A).  The per-level freedom h_D = H(D-data | S-data) is computed
exactly per coset pair from the pattern lists (S-level data of a coset pair = the sum of the two lifts' x-values).
usage: python3 theory/nd_levels.py m mu K Yfile [--idx 0] [--N 20000] [--levels "1,5;2,10"]"""
import sys, os, json, math, itertools
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits
from nd_stage2_est import Slice

def ramanujan(d, n):
    """c_d(h) for h in Z_n: sum over the characters of Z_n of order d"""
    out = np.zeros(n)
    for j in range(n):
        if n // math.gcd(j, n) == d: out += np.cos(2 * np.pi * j * np.arange(n) / n)
    return np.rint(out).astype(np.int64)

def main(m, mu, K, Y, N=20000, levels=None, seed=1):
    S = Slice(m, mu, K); F = S.F; rng = np.random.default_rng(seed); Y = np.array(Y)
    pats = [[S.patterns(j, int(Y[k][b])) for j, b in enumerate(S.freps)] for k in range(4)]
    H2 = sum(math.log2(len(p)) for pk in pats for p in pk)
    # CRT index: element with (h mod F, z mod K)
    idx = np.zeros((F, K), np.int64)
    for x in range(m): idx[x % F, x % K] = x
    ds = sorted({F // math.gcd(j, F) for j in range(F)})
    R = {d: ramanujan(d, F) for d in ds}
    hh = np.arange(F)
    C = {d: R[d][(hh[:, None] - hh[None, :]) % F] for d in ds}           # projection onto order-d characters
    proj = {d: np.zeros((N, F * K), np.int64) for d in ds}
    chunk = 1000
    for st in range(0, N, chunk):
        c = min(chunk, N - st); acc = np.zeros((c, m))
        for k in range(4):
            lab = np.zeros((c, S.nor), np.int8)
            for j in range(len(S.freps)):
                P = pats[k][j]; pick = P[rng.integers(len(P), size=c)]
                for i, (o, s) in enumerate(S.parts[j]): lab[:, o] = pick[:, i]
            x = lab[:, S.orb].astype(np.float64); fx = np.fft.rfft(x, axis=1)
            acc += np.fft.irfft(fx * fx.conj(), n=m, axis=1)
        Pm = np.rint(acc).astype(np.int64)[:, idx]                 # (c, F, K)
        for d in ds: proj[d][st:st + c] = np.einsum('hj,njk->nhk', C[d], Pm).reshape(c, F * K)
    Tm = np.full(m, -4, np.int64); Tm[0] = 4 * m
    TP = Tm[idx]
    Tproj = {d: (C[d] @ TP).ravel() for d in ds}
    def bits(level):
        V = np.concatenate([proj[d] for d in level], 1); T = np.concatenate([Tproj[d] for d in level])
        # drop exactly duplicated columns (M-symmetry) to keep the lattice small
        _, keep = np.unique(V[:200].T, axis=0, return_index=True); keep = np.sort(keep)
        return gauss_bits(V[:, keep], T[keep])
    out = dict(m=m, mu=mu, K=K, F=F, H2=round(H2, 2), orders=ds)
    allb = bits(ds); out['all'] = allb
    if levels:
        for lv in levels:
            out['level_' + ','.join(map(str, lv))] = bits(lv)
    print(json.dumps(out), flush=True)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--idx": "0", "--N": "20000", "--levels": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    m, mu, K = map(int, pos[:3]); d = json.load(open(pos[3]))
    sols = d['solutions'] if 'solutions' in d else [d]
    lv = [list(map(int, t.split(','))) for t in opt["--levels"].split(';')] if opt["--levels"] else None
    main(m, mu, K, sols[int(opt["--idx"])]['y'], int(opt["--N"]), lv)
