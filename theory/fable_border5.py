#!/usr/bin/env python3
"""Border for a Q-level solution (no core needed): Q = Z_2^a, character set T, compressions y (4 x 2^a integers).

Rows of E and Z: block by block a +-1 function on Q with Walsh support inside T (all such functions are listed by
brute force).  CP-SAT: choose the functions; K = E M^T Z^T / 4k must be +-1; rows of E pairwise orthogonal, rows of
Z pairwise orthogonal.  Check afterwards: K Hadamard, and the compressed matrix [[K, sqrt(s) E], [sqrt(s) F, M]]
with F = -Z^T has Gram n I (this is the whole border condition; a core with these compressions then completes to a
Hadamard matrix).  [M]
Usage: fable_border5.py solutions.json index s time_limit   (file written by fable_q599.py)
"""
import sys, json, time, itertools
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2

def main():
    d = json.load(open(sys.argv[1])); idx = int(sys.argv[2]); s_ = int(sys.argv[3]); tl = float(sys.argv[4])
    T = d['T']; t = len(T); K_ = 16; a = 4; b = 4 * t
    sol = d['solutions'][idx]
    cellorb = d['cell_orbits']
    ys = []
    for c in sol['cells']:
        y = [0] * K_
        for j, ob in enumerate(cellorb):
            for q in ob: y[q] = c[j]
        ys.append(np.array(y, dtype=np.int64))
    G = B2.Grp2(a, 1)                       # Q itself
    ch = G.chars()
    # element order of G.Qel is the binary expansion; index q of the solution = integer with the same bits
    perm = [G.Qel.index(tuple((q >> (a - 1 - i)) & 1 for i in range(a))) for q in range(K_)]
    ys_g = [np.array([y[[perm.index(i) for i in range(K_)][j]] for j in range(K_)]) for y in ys]
    M = FB.gs_array(G, np.array(ys_g))
    n = 4 * (K_ * s_ + t)
    # Walsh check in G's own character table
    Tm = [m for m in range(K_) if m in set(T)]
    E_ = [sum(int(ch[m] @ y) ** 2 for y in ys_g) for m in range(K_)]
    lowset = [m for m in range(K_) if E_[m] == 4 * t]
    assert sorted(E_) == sorted([4 * t] * t + [n] * (K_ - t)), E_
    print('characters with energy 4t (in the table of fable_border2):', lowset, flush=True)
    # +-1 functions with Walsh support inside the low set
    funcs = []
    for bits in range(1 << K_):
        f = np.array([1 - 2 * ((bits >> q) & 1) for q in range(K_)])
        if all(int(ch[m] @ f) == 0 for m in range(K_) if m not in lowset): funcs.append(f)
    F = np.array(funcs); no = len(F)
    print(len(funcs), 'admissible block functions', flush=True)
    # tables: tab[a][a2][o][o2] = f_o (block a) M^T f_o2 (block a2)
    tabs = {}
    for a1 in range(4):
        for a2 in range(4):
            blk = M.T[a1 * K_:(a1 + 1) * K_, a2 * K_:(a2 + 1) * K_]
            tb = (F @ blk @ F.T)
            tabs[a1, a2] = [int(x) for x in tb.reshape(-1)]
    orth = [int(x) for x in (F @ F.T).reshape(-1)]
    mdl = cp_model.CpModel()
    E = [[mdl.NewIntVar(0, no - 1, '') for _ in range(4)] for r in range(b)]
    Z = [[mdl.NewIntVar(0, no - 1, '') for _ in range(4)] for r in range(b)]
    lo = min(min(tb) for tb in tabs.values()); hi = max(max(tb) for tb in tabs.values())
    for r in range(b):
        for j in range(b):
            terms = []
            for a1 in range(4):
                for a2 in range(4):
                    ix = mdl.NewIntVar(0, no * no - 1, ''); mdl.Add(ix == E[r][a1] * no + Z[j][a2])
                    v = mdl.NewIntVar(lo, hi, ''); mdl.AddElement(ix, tabs[a1, a2], v); terms.append(v)
            sg = mdl.NewBoolVar('')
            mdl.Add(sum(terms) == (8 * K_) * sg - 4 * K_)          # +-4k
    for Xv in (E, Z):
        for r in range(b):
            for r2 in range(r + 1, b):
                terms = []
                for a1 in range(4):
                    ix = mdl.NewIntVar(0, no * no - 1, ''); mdl.Add(ix == Xv[r][a1] * no + Xv[r2][a1])
                    v = mdl.NewIntVar(-K_, K_, ''); mdl.AddElement(ix, orth, v); terms.append(v)
                mdl.Add(sum(terms) == 0)
    slv = cp_model.CpSolver(); slv.parameters.max_time_in_seconds = tl; slv.parameters.num_workers = 4
    t0 = time.time(); r_ = slv.Solve(mdl); st = slv.StatusName(r_)
    print(f'border for solution {idx}: {st} ({time.time()-t0:.0f}s)', flush=True)
    if r_ not in (cp_model.OPTIMAL, cp_model.FEASIBLE): return
    Em = np.array([np.concatenate([F[slv.Value(E[r][a1])] for a1 in range(4)]) for r in range(b)])
    Zm = np.array([np.concatenate([F[slv.Value(Z[r][a1])] for a1 in range(4)]) for r in range(b)])
    Km = (Em @ M.T @ Zm.T) // (4 * K_)
    okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=int))
    Fm = -Zm.T
    # Gram of the compressed matrix, in integers: rows [K, E] weighted (1, s) and rows [F, M] weighted (s, 1)
    top = Km @ Km.T + s_ * (Em @ Em.T)
    cross = Km @ Fm.T + Em @ M.T
    bot = s_ * (Fm @ Fm.T) + M @ M.T
    good = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * K_, dtype=int))
    print(f'  K Hadamard: {okK}; compressed matrix of order {b + 4 * K_} has Gram {n} I: {good}', flush=True)
    json.dump(dict(T=T, s=s_, cells=sol['cells'], K=Km.tolist(), E=Em.tolist(), Z=Zm.tolist()),
              open(f'logs/fable_fresh/q599_border_{idx}.json', 'w'))
if __name__ == '__main__':
    main()
