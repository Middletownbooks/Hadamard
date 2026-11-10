#!/usr/bin/env python3
"""Monomial complex orthogonal designs of order 4 for typed quadruples.

Rows j = 1..4:  z_j = sum_i  e_{ji} * x_i * u_{ji},  u_{ji} a letter y_k or conj(y_k), each k once per row.
Requirement: sum_j conj(M_{ji}) M_{jl} = delta_il |y|^2 formally, on the type.
If such M exists, T-matrices (any, order t) x typed GS quadruple -> GS quadruple (x_i = transform of T_i).
Exhaustive [M].
"""
import sys, itertools, collections
import numpy as np
sys.path.insert(0, 'theory')
from fable_tmult import symbols_for_type

def run(typ):
    sym, conj, letters = symbols_for_type(typ)
    nsym = len(conj)
    mono = {}
    for a in range(nsym):
        for b in range(a, nsym):
            mono[(a, b)] = len(mono)
    nm = len(mono)
    pairs = [(i, l) for i in range(4) for l in range(i + 1, 4)]
    # a row: assignment slot i -> (sign, k, s)
    rows = []
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((1, -1), repeat=4):
            if signs[0] != 1:
                continue      # row sign normalisation
            for cj in itertools.product((0, 1), repeat=4):
                row = tuple((signs[i], perm[i], cj[i]) for i in range(4))
                # skip duplicates caused by type identifications
                key = tuple((e, sym[(k, s)]) for e, k, s in row)
                rows.append((key, row))
    uniq = {}
    for key, row in rows:
        uniq.setdefault(key, row)
    rows = list(uniq.values())
    def vec(row):
        v = np.zeros((len(pairs), 2, nm), dtype=np.int16)   # conj(M_ji) M_jl : keep ordered (need both P and its formal conj -> store as is)
        for pi, (i, l) in enumerate(pairs):
            ei, ki, si = row[i]
            el, kl, sl = row[l]
            a = conj[sym[(ki, si)]]
            b = sym[(kl, sl)]
            if a > b: a, b = b, a
            v[pi, 0, mono[(a, b)]] += ei * el
        c = np.zeros((4, 4), dtype=np.int16)                 # column i contains index k
        for i in range(4):
            c[i, row[i][1]] += 1
        return np.concatenate([v.reshape(-1), c.reshape(-1)])
    V = np.array([vec(r) for r in rows])
    tgt = np.concatenate([np.zeros(len(pairs) * 2 * nm, dtype=np.int16), np.ones(16, dtype=np.int16)])
    n = len(rows)
    d = collections.defaultdict(list)
    for i in range(n):
        S = V[i] + V[i:]
        for jj in range(S.shape[0]):
            d[S[jj].tobytes()].append((i, i + jj))
    nsol = 0; ex = None
    for key, prs in d.items():
        s = np.frombuffer(key, dtype=np.int16)
        comp = (tgt - s).astype(np.int16).tobytes()
        if comp in d and key <= comp:
            for p in prs:
                for q in d[comp]:
                    if key != comp or p <= q:
                        nsol += 1
                        if ex is None: ex = (p, q)
    print(f"type={typ}: rows {n}, designs (unordered row-pair pairs) {nsol}")
    if ex:
        names = 'ABCD'
        for i in ex[0] + ex[1]:
            print('   ', ' '.join(('+' if e > 0 else '-') + names[k] + ('*' if s and sym[(k,1)] != sym[(k,0)] else ' ') for e, k, s in rows[i]))
    return nsol

if __name__ == '__main__':
    for typ in sys.argv[1:]:
        run(typ)
