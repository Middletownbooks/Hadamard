#!/usr/bin/env python3
"""CP-SAT search for single-character borders (see fable_border3.py for the setting and fable_border2.py for the
family).  Variables: for every row of E and of Z and every block, an option o in [0, 2t) (character o % t of T,
sign + if o < t).  Constraints: K[r, j] = (1/4) sum_{a, a'} tab[a][a'][e_ra][z_ja'] in {+1, -1};
rows of E pairwise orthogonal, rows of Z pairwise orthogonal (then K is Hadamard).  A found border is assembled
and H H^T = n I checked.  [M]     Usage: fable_border4.py a s T ncores time_limit [seed]
"""
import sys, time, itertools
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2

def solve(G, T, M, tl, workers=4, seed=0):
    k = G.k; t = len(T); b = 4 * t; no = 2 * t
    ch = G.chars()
    V = np.zeros((4, t, 4 * k), dtype=np.int64)
    for a in range(4):
        for i, m in enumerate(T):
            V[a, i, a * k:(a + 1) * k] = ch[m]
    D = np.einsum('aic,dc,bjd->aibj', V, M, V) // k
    sg = lambda o: 1 if o < t else -1
    mdl = cp_model.CpModel()
    E = [[mdl.NewIntVar(0, no - 1, '') for a in range(4)] for r in range(b)]
    Z = [[mdl.NewIntVar(0, no - 1, '') for a in range(4)] for r in range(b)]
    # tables
    tabs = {}
    for a in range(4):
        for a2 in range(4):
            tab = [int(sg(o) * sg(o2) * D[a, o % t, a2, o2 % t]) if (o % t) == (o2 % t) else 0 for o in range(no) for o2 in range(no)]
            tabs[a, a2] = tab
    orth_tab = [sg(o) * sg(o2) if (o % t) == (o2 % t) else 0 for o in range(no) for o2 in range(no)]
    lo = min(min(tb) for tb in tabs.values()); hi = max(max(tb) for tb in tabs.values())
    for r in range(b):
        for j in range(b):
            terms = []
            for a in range(4):
                for a2 in range(4):
                    tb = tabs[a, a2]
                    if not any(tb): continue
                    idx = mdl.NewIntVar(0, no * no - 1, '')
                    mdl.Add(idx == E[r][a] * no + Z[j][a2])
                    v = mdl.NewIntVar(lo, hi, '')
                    mdl.AddElement(idx, tb, v)
                    terms.append(v)
            sgn = mdl.NewBoolVar('')
            mdl.Add(sum(terms) == 8 * sgn - 4)
    for Xv in (E, Z):
        for r in range(b):
            for r2 in range(r + 1, b):
                terms = []
                for a in range(4):
                    idx = mdl.NewIntVar(0, no * no - 1, '')
                    mdl.Add(idx == Xv[r][a] * no + Xv[r2][a])
                    v = mdl.NewIntVar(-1, 1, '')
                    mdl.AddElement(idx, orth_tab, v)
                    terms.append(v)
                mdl.Add(sum(terms) == 0)
    for a in range(4): mdl.Add(E[0][a] < t)        # first row of E: all signs +
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    s.parameters.random_seed = seed
    r = s.Solve(mdl); st = s.StatusName(r)
    if r not in (cp_model.OPTIMAL, cp_model.FEASIBLE): return st, None
    def rows(Xv):
        R = np.zeros((b, 4 * k), dtype=np.int64)
        for rr in range(b):
            for a in range(4):
                o = s.Value(Xv[rr][a]); R[rr, a * k:(a + 1) * k] = sg(o) * ch[T[o % t]]
        return R
    Em = rows(E); Zm = rows(Z)
    K = (Em @ M.T @ Zm.T) // (4 * k)
    return st, (K, Em, Zm)

def main():
    a = int(sys.argv[1]); s_ = int(sys.argv[2]); T = [int(x) for x in sys.argv[3].split(',')]
    ncores = int(sys.argv[4]); tl = float(sys.argv[5]); seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    G = B2.Grp2(a, s_); rng = np.random.default_rng(seed)
    k = G.k; b = 4 * len(T); n = 4 * G.v + b
    print(f'Q = Z_2^{a}, S = Z_{s_}, T = {T}, |G| = {G.v}, b = {b}, n = {n}', flush=True)
    for ci_ in range(ncores):
        t0 = time.time(); X = None; tries = 0
        while X is None and tries < 40:
            X = B2.find_core(G, T, rng); tries += 1
        if X is None:
            print(f'core {ci_}: none found'); continue
        Cm = FB.gs_array(G, X); M = FB.compress(G, Cm); ch = G.chars()
        prof = [[int(sum(ch[m][q] * M[0, j * k + q] for q in range(k))) for j in range(4)] for m in T]
        print(f'core {ci_}: found in {time.time()-t0:.0f}s; Walsh values on T: {prof}; sequences: ' + ' '.join(''.join('+' if v > 0 else '-' for v in x) for x in X), flush=True)
        t1 = time.time()
        st, res = solve(G, T, M, tl, seed=seed)
        if res is None:
            print(f'   single-character border: {st} ({time.time()-t1:.0f}s)', flush=True); continue
        K, E, Z = res
        okK = np.array_equal(np.abs(K), np.ones_like(K)) and np.array_equal(K @ K.T, b * np.eye(b, dtype=int))
        H = FB.assemble(G, Cm, K, E, Z)
        good = np.array_equal(H @ H.T, n * np.eye(n, dtype=int))
        print(f'   single-character border: {st} ({time.time()-t1:.0f}s); K Hadamard {okK}; assembled H({n}): H H^T = n I is {good}', flush=True)

if __name__ == '__main__':
    main()
