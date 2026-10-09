"""Z2 x Zh structure of the published b=1 cores (m = v-1 = 2h, h odd): P = (X0+X1)/2, R = (X0-X1)/2 where
X0(y) = x(i), X1(y) = x(i+h) for the CRT index i <-> (i mod 2, i mod h).  Measures, per block:
  * |S| (P-support), symmetry of S under y -> -y (broken pairs), symmetry of S under the multipliers of Zh
  * on symmetric pairs inside S: fraction with P(-y) = P(y) ; inside S^c: fraction with R(-y) = R(y)
Null: random +-1 blocks with the same row sum (symmetry of S), 2000 draws."""
import sys, json
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/scripts')
from hadlib import Core
rng = np.random.default_rng(1)
def to_cyclic(c):
    n = c.m
    if len(c.mods) == 1: return c.X
    idx = np.array([c.idx(np.array([i % mm for mm in c.mods])) for i in range(n)])
    return c.X[:, idx]
def split(x, h):
    i = np.arange(2 * h)
    X0 = np.zeros(h, int); X1 = np.zeros(h, int)
    for t in range(2 * h):
        (X0 if t % 2 == 0 else X1)[t % h] = x[t]
    P = (X0 + X1) // 2; R = (X0 - X1) // 2
    return P, R
def sym_stats(P, R, h):
    S = P != 0
    neg = (-np.arange(h)) % h
    broken = int(np.sum(S[1:] != S[neg][1:]) // 2)
    pairs_S = [y for y in range(1, h) if S[y] and S[neg[y]] and y < neg[y]]
    pairs_C = [y for y in range(1, h) if (not S[y]) and (not S[neg[y]]) and y < neg[y]]
    Psym = np.mean([P[y] == P[neg[y]] for y in pairs_S]) if pairs_S else np.nan
    Rsym = np.mean([R[y] == R[neg[y]] for y in pairs_C]) if pairs_C else np.nan
    return int(S.sum()), broken, float(Psym), float(Rsym)
for f, v in [('/home/user/Hadamard/data/core_v167_b1.txt', 167), ('/home/user/Hadamard/data/core_v179_b1.txt', 179),
             ('/home/user/Hadamard/data/core_v419_b1.txt', 419), ('/home/user/Hadamard/data/core_v443_b1.txt', 443)]:
    c = Core(f); X = to_cyclic(c); h = c.m // 2
    for k in range(4):
        P, R = split(X[k], h)
        sz, br, ps, rs = sym_stats(P, R, h)
        nb = []
        for _ in range(2000):
            y = rng.permutation(X[k]); Pn, Rn = split(y, h); nb.append(sym_stats(Pn, Rn, h)[1])
        nb = np.array(nb)
        print(json.dumps(dict(v=v, block='ABCD'[k], h=h, S_size=sz, broken_pairs=br, null_broken_mean=round(float(nb.mean()), 1),
                              null_broken_min=int(nb.min()), P_symmetric_frac=round(ps, 3), R_symmetric_frac=round(rs, 3))), flush=True)
