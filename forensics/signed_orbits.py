"""Signed multiplier relations beyond the stabiliser, for every available original (natural coordinates).
U = largest scalar-multiplier group fixing all four blocks about a common fixed point c (census_originals.py).
For each multiplier a with a not in U but a^2 in U (an involution mod U), the U-orbits pair up as (O, aO); we count
pairs with x(aO) = x(O) ('same') vs x(aO) = -x(O) ('anti').  Null: each pair is a fair coin (two-sided binomial tail).
The one strong hit is v=359, a = 17 (identity on Z2 x Z2 x Z8, a non-square on Z11): blocks are almost exactly of the
form x(g, z) = eps(g) * Legendre_11(z) for z != 0, i.e. anti-invariant under the non-squares of Z11."""
import os, sys, glob, re
from math import gcd, comb, log10
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census_originals import natural_core, stabiliser
from gbuilder import read_seq
from load_originals import load_any
def two_sided(n, k):
    k = max(k, n - k); return min(1.0, 2 * sum(comb(n, j) for j in range(k, n + 1)) / 2 ** n)
def run(v, b, X, g):
    d, c, U = stabiliser(X, g)
    sh = g.idx(g.coords + g.coords[c]); X = X[:, sh]                    # move the fixed point to 0
    e = int(np.lcm.reduce(g.mods)); m = g.m
    orb = -np.ones(m, int); reps = []
    for y in range(m):
        if orb[y] < 0:
            o = sorted(set(int(g.idx(a * g.coords[y])) for a in U)); orb[o] = len(reps); reps.append(o)
    seen, out = set(), []
    for a in range(2, e):
        if gcd(a, e) != 1 or a in U or (a * a) % e not in U: continue
        key = tuple(sorted(a * u % e for u in U))
        if key in seen: continue
        seen.add(key)
        img = orb[g.idx(a * g.coords)]
        pairs = sorted({(min(orb[y], img[y]), max(orb[y], img[y])) for y in range(m) if orb[y] != img[y]})
        per = []
        for k in range(4):
            xo = np.array([X[k][o[0]] for o in reps])
            per.append(sum(1 for p, q in pairs if xo[p] == xo[q]))
        n, s = 4 * len(pairs), sum(per)
        out.append((log10(max(two_sided(n, s), 1e-300)), a, s, n, per, len(pairs)))
    out.sort()
    tests = len(out)
    if out:
        lp, a, s, n, per, npairs = out[0]
        print(f'v={v}: |U|={d}, {len(reps)} orbits/block, {tests} involution class(es) tested; strongest: a={a}: '
              f'{s}/{n} pairs "same" (per block {per} of {npairs}); log10 two-sided p = {lp:.1f}'
              f' (x{tests} for multiple testing: {lp + log10(tests):.1f})')
    else:
        print(f'v={v}: |U|={d}: no involution classes mod U')
if __name__ == '__main__':
    for f in sorted(glob.glob(os.path.join(HERE, 'originals', 'matrix_H*.txt*')) + glob.glob(os.path.join(HERE, 'originals', 'seq_v*_g.txt'))):
        base = os.path.basename(f).replace('.gz', '').replace('.txt', '')
        if base.startswith('seq'):
            hdr, seqs = read_seq(f); H = None
        else:
            hdr = (open(f).readline() if f.endswith('.txt') else __import__('gzip').open(f, 'rt').readline()); seqs = None; H = load_any(base)
        meta = dict(re.findall(r'(\w+)=(\w+)', hdr)); v, b = int(meta['v']), int(meta['b'])
        X, g = natural_core(v, b, meta['func'], H, seqs)
        run(v, b, X, g)
