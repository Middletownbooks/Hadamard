"""Symmetry census of the authors' own records (all available originals), in natural Goethals-Seidel coordinates:
the back-identity / inversion R is undone for blocks B, C, D (cyclic: f(y) = g(-1-y); group-developed: f(y) = g(-y)).
For each core: largest scalar-multiplier group U fixing all four blocks about one common fixed point c (searched over
all c), the number of U-orbits per block, and the row sums."""
import os, sys, re, glob
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from halfshift_b1 import core_seqs
from gbuilder import read_seq
from load_originals import originals, load_any, header
class G:
    def __init__(s, mods):
        s.mods = np.array(mods); s.m = int(np.prod(mods))
        s.coords = np.array(np.unravel_index(np.arange(s.m), mods)).T
    def idx(s, c):
        c = np.asarray(c) % s.mods
        return np.ravel_multi_index(tuple(np.moveaxis(c, -1, 0)), s.mods)
GROUPS = {347: (2, 3, 57), 359: (2, 2, 8, 11), 479: (2, 2, 119)}
def natural_core(v, b, func, H=None, seqs=None):
    m = v - b
    if b >= 2:
        g = G(GROUPS[v]); off = 4 * b; neg = g.idx(-g.coords)
        X = np.array([H[off, off + q * m: off + (q + 1) * m] for q in range(4)])
        return np.array([X[0]] + [X[q][neg] for q in (1, 2, 3)]), g
    g = G((m,)); I = np.arange(m)
    if seqs is not None:
        X = np.array(seqs)                                   # records: row 0 of the first block row
    elif b == 1:
        X = np.array([x for _, x in core_seqs(H, m)])
    else:
        X = np.array([H[0, q * m:(q + 1) * m] for q in range(4)])
    return np.array([X[0]] + [X[q][(-1 - I) % m] for q in (1, 2, 3)]), g
def stabiliser(X, g):
    e = int(np.lcm.reduce(g.mods)); scal = [a for a in range(2, e) if gcd(a, e) == 1]
    best = (1, 0, [1])
    for c in range(g.m):
        sh = g.idx(g.coords + g.coords[c]); Y = X[:, sh]            # Y(y) = X(y + c)
        U = [1] + [a for a in scal if all(np.array_equal(Y[k][g.idx(a * g.coords)], Y[k]) for k in range(4))]
        if len(U) > best[0]: best = (len(U), c, U)
    return best
def orbits(g, U, c):
    seen = np.zeros(g.m, bool); r = 0
    for y in range(g.m):
        if not seen[y]:
            r += 1
            for a in U:
                seen[int(g.idx(a * (g.coords[y] - g.coords[c]) + g.coords[c]))] = True
    return r
if __name__ == '__main__':
    rows = []
    srcs = [('m', n) for n in originals()] + [('s', f) for f in sorted(glob.glob(os.path.join(HERE, 'originals', 'seq_v*.txt')))]
    for kind, f in srcs:
        if kind == 's':
            hdr, seqs = read_seq(f); H = None; meta = dict(re.findall(r'(\w+)=(\w+)', hdr)); f = os.path.basename(f)
        else:
            H = load_any(f); seqs = None; meta = dict(re.findall(r'(\w+)=(\w+)', header(f)))
        v, b, func = int(meta['v']), int(meta['b']), meta['func']
        X, g = natural_core(v, b, func, H, seqs)
        d, c, U = stabiliser(X, g)
        r = orbits(g, U, c)
        rows.append((v, b, func, 'x'.join(map(str, g.mods)), d, r, tuple(int(t) for t in g.coords[c]), X.sum(1).tolist(), f))
    print('| v | N | b | builder | group | d = |U| | orbits/block | common fixed point | row sums (natural) | source |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for v, b, func, grp, d, r, c, rs, src in sorted(rows):
        cs = c[0] if len(c) == 1 else c
        print(f'| {v} | {4*v} | {b} | {func} | Z{grp.replace("x", "×Z")} | {d} | {r} | {cs} | {tuple(rs)} | `{src}` |')
