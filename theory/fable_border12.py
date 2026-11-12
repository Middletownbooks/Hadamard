#!/usr/bin/env python3
"""General (any block function with Walsh support in T) sigma-invariant borders for the 2396 slice: the first exact
gate.  Every sigma-invariant border has NFIX >= 1 fixed rows (28 = 1 mod 3) and orbit rows.
 (1) <z, sigma z> = sum_c <f_c, f_c o sigma>: computed per block function (0 or 16 for each of the 30 functions); an
     orbit row needs <z, sigma z> = 0.
 (2) A fixed E-row f and an orbit Z-row z: K(f, sigma^i z) = K(f, z) for all i (f is fixed), so the condition is |K(f, z)| = 1.
 (3) Symmetrically a fixed Z-row and an orbit E-row.
Counts, exhaustive over all 30^4 rows.  If no orbit row is compatible with any fixed row of the other side, there is no
sigma-invariant border at all.  [M]
"""
import sys, itertools, time
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2, fable_border6 as B6, fable_tcensus2 as C2

def main():
    idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    d, sol, M, F, t = B6.load('logs/fable_fresh/q599_solutions.json', idx)
    k = 16; b = 28
    sig = [x for x in C2.gl_odd_classes() if (x[0], x[1]) == (3, 1)][0]; img = sig[2]
    Ps = np.zeros((k, k), dtype=np.int64)
    for q in range(k): Ps[q, img[q]] = 1
    no = len(F); key = {tuple(f): i for i, f in enumerate(F)}
    pi = np.array([key[tuple(F[o] @ Ps.T)] for o in range(no)])
    w = np.array([int(F[o] @ F[pi[o]]) for o in range(no)])          # <f, f o sigma>
    fixed = np.nonzero(pi == np.arange(no))[0]
    print(f'case {idx}: {no} block functions; <f, f o sigma> values {sorted(set(w.tolist()))}; sigma-fixed: {fixed.tolist()} '
          f'with <f, f o sigma> = {w[fixed].tolist()}', flush=True)
    # tables: tab[c][c'][o][o'] = F_o (block c) M^T F_o' (block c')
    MT = M.T
    tab = np.array([[F @ MT[c * k:(c + 1) * k, c2 * k:(c2 + 1) * k] @ F.T for c2 in range(4)] for c in range(4)])   # 4,4,no,no
    rows = np.array(list(itertools.product(range(no), repeat=4)), dtype=np.int16)        # all 30^4
    self_sig = w[rows].sum(1)                                                         # <z, sigma z>
    orbit_ok = self_sig == 0
    is_fixed = np.all(pi[rows] == rows, axis=1)
    print(f'{len(rows)} rows; sigma-fixed rows {int(is_fixed.sum())}; non-fixed rows with <z, sigma z> = 0: '
          f'{int((orbit_ok & ~is_fixed).sum())}', flush=True)
    frows = rows[is_fixed]; orows = rows[orbit_ok & ~is_fixed]
    # K(f, z) * 64 = sum_{c,c'} tab[c][c'][f_c][z_c']; need |.| = 64
    best = 0; ncomp_E = 0; ncomp_Z = 0
    t0 = time.time()
    compE = np.zeros(len(orows), dtype=bool); compZ = np.zeros(len(orows), dtype=bool)
    for fi, f in enumerate(frows):
        # f as E-row against orbit Z-rows
        v = np.zeros(len(orows), dtype=np.int64)
        for c in range(4):
            for c2 in range(4):
                v += tab[c, c2][f[c]][orows[:, c2]]
        compZ |= (np.abs(v) == 64)
        # f as Z-row against orbit E-rows
        u = np.zeros(len(orows), dtype=np.int64)
        for c in range(4):
            for c2 in range(4):
                u += tab[c, c2][orows[:, c], f[c2]]
        compE |= (np.abs(u) == 64)
    print(f'orbit Z-rows compatible with some fixed E-row: {int(compZ.sum())}; orbit E-rows compatible with some fixed '
          f'Z-row: {int(compE.sum())}  ({time.time()-t0:.0f}s)', flush=True)
    if compZ.sum() == 0 or compE.sum() == 0:
        print(f'=> NO sigma-invariant border (general block functions) exists for case {idx}', flush=True)
    else:
        print('=> the first gate does not close the general sigma-invariant case', flush=True)

if __name__ == '__main__':
    main()
