#!/usr/bin/env python3
"""sigma-invariant borders for the 2396 slice (report/fable_fresh_answer.md section 7.1), CP-SAT.

sigma: the element of order 3 of GL(4,2) without nonzero fixed points that defines the slice; the compressions are
constant on its orbits, so the Goethals-Seidel array M of the compressions commutes with P = diag(P_s, P_s, P_s, P_s),
P_s f = f o sigma (checked below).  If the row sets of E and Z are unions of sigma-orbits (row r -> r o sigma), then
  K(sigma e, sigma z) = K(e, z)  and  <sigma a, sigma b> = <a, b>,
so it suffices to impose the K-conditions for E-representatives against ALL rows of Z, and orthogonality of each
representative against all rows of its own side.  Rows: NORB orbit representatives r_i (their images sigma r_i,
sigma^2 r_i are derived through the permutation that sigma induces on the block functions) and NFIX sigma-fixed rows
(every block a sigma-fixed function); 3 NORB + NFIX = 28.
Modes: 'single' (block = one signed character of T) or 'general' (block = any +-1 function with Walsh support in T).
Every border found is checked with integers (K Hadamard, Gram of the compressed matrix) and saved.  [M]
Usage: fable_border10.py mode nfix time_limit workers [case_index] [seed]
"""
import sys, json, time
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2, fable_border6 as B6, fable_tcensus2 as C2

def main():
    mode = sys.argv[1]; nfix = int(sys.argv[2]); tl = float(sys.argv[3]); workers = int(sys.argv[4])
    idx = int(sys.argv[5]) if len(sys.argv) > 5 else 0   # -1 = control
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    if idx >= 0:
        d, sol, M, Fgen, t = B6.load('logs/fable_fresh/q599_solutions.json', idx)
        k = 16; s_ = 37; G = B2.Grp2(4, 1); ch = G.chars(); T = d['T']
        sig = [x for x in C2.gl_odd_classes() if (x[0], x[1]) == (3, 1)][0]
        img = sig[2]
    else:
        # control (case -1): Q = Z_2^2, sigma = Singer cycle 1 -> 2 -> 3 -> 1, spike compressions, s = 7, n = 124
        k = 4; s_ = 7; G = B2.Grp2(2, 1); ch = G.chars(); T = [1, 2, 3]; t = 3
        ys = [np.array(y, dtype=np.int64) for y in ([1, 1, 1, 1], [3, 1, 1, 1], [3, 1, 1, 1], [3, 1, 1, 1])]
        assert all(sum(int(ch[mm] @ y) ** 2 for y in ys) == (12 if mm in T else 124) for mm in range(4))
        M = FB.gs_array(G, np.array(ys)); img = [0, 2, 3, 1]
        Fgen = np.array([f for f in (np.array([1 - 2 * ((bits >> q) & 1) for q in range(k)]) for bits in range(1 << k))
                         if all(int(ch[mm] @ f) == 0 for mm in range(k) if mm not in T)], dtype=np.int64)
        sol = dict(cells=[y.tolist() for y in ys]); d = dict(T=T)
    b = 4 * t; n = 4 * (k * s_ + t)
    assert (b - nfix) % 3 == 0
    norb = (b - nfix) // 3
    Ps = np.zeros((k, k), dtype=np.int64)
    for q in range(k): Ps[q, img[q]] = 1          # (f Ps^T)(q) = f(img(q))  i.e. row vector f -> f @ Ps.T
    P = np.kron(np.eye(4, dtype=np.int64), Ps)
    assert np.array_equal(P @ M, M @ P), 'M does not commute with sigma'
    if mode == 'single':
        F = np.array([e * ch[m] for e in (1, -1) for m in T], dtype=np.int64)
    else:
        F = Fgen
    no = len(F)
    key = {tuple(f): i for i, f in enumerate(F)}
    pi = [key[tuple(F[o] @ Ps.T)] for o in range(no)]            # sigma on block functions
    assert all(pi[pi[pi[o]]] == o for o in range(no))
    fixed = [o for o in range(no) if pi[o] == o]
    print(f'mode {mode}: {no} block functions, {len(fixed)} sigma-fixed; NORB = {norb}, NFIX = {nfix}', flush=True)
    if nfix and not fixed: print('no sigma-fixed block functions: NFIX must be 0'); return
    MT = M.T
    tabs = {}
    for c in range(4):
        for c2 in range(4):
            tabs[c, c2] = [int(v) for v in (F @ MT[c * k:(c + 1) * k, c2 * k:(c2 + 1) * k] @ F.T).reshape(-1)]
    orth = [int(v) for v in (F @ F.T).reshape(-1)]
    lo = min(min(v) for v in tabs.values()); hi = max(max(v) for v in tabs.values())
    m = cp_model.CpModel()
    def side(name):
        # [P] <z, sigma z> = sum_c <f_c, f_c o sigma>, each term 16 for a sigma-fixed block function and 0 otherwise,
        # so the rows of a 3-orbit use only non-fixed block functions (domain reduction, exact)
        nonfix = cp_model.Domain.FromValues([o for o in range(no) if pi[o] != o])
        reps = [[m.NewIntVarFromDomain(nonfix, f'{name}r{i}_{c}') for c in range(4)] for i in range(norb)]
        rows = []
        for r in reps:
            r1 = [m.NewIntVarFromDomain(nonfix, '') for c in range(4)]
            r2 = [m.NewIntVarFromDomain(nonfix, '') for c in range(4)]
            for c in range(4):
                m.AddElement(r[c], pi, r1[c]); m.AddElement(r1[c], pi, r2[c])
            rows += [r, r1, r2]
        fx = []
        for i in range(nfix):
            f = [m.NewIntVarFromDomain(cp_model.Domain.FromValues(fixed), f'{name}f{i}_{c}') for c in range(4)]
            fx.append(f); rows.append(f)
        return reps + fx, rows
    Erep, Erows = side('E'); Zrep, Zrows = side('Z')
    nodes = 0
    for e in Erep:                                          # K conditions: representatives x all Z rows
        for z in Zrows:
            terms = []
            for c in range(4):
                for c2 in range(4):
                    tb = tabs[c, c2]
                    if not any(tb): continue
                    ix = m.NewIntVar(0, no * no - 1, ''); m.Add(ix == e[c] * no + z[c2])
                    v = m.NewIntVar(lo, hi, ''); m.AddElement(ix, tb, v); terms.append(v)
            sg = m.NewBoolVar(''); m.Add(sum(terms) == 8 * k * sg - 4 * k); nodes += 1
    for reps, rows in ((Erep, Erows), (Zrep, Zrows)):      # orthogonality: representatives x all other rows
        for a in reps:
            for bb in rows:
                if bb is a: continue
                terms = []
                for c in range(4):
                    ix = m.NewIntVar(0, no * no - 1, ''); m.Add(ix == a[c] * no + bb[c])
                    v = m.NewIntVar(-k, k, ''); m.AddElement(ix, orth, v); terms.append(v)
                m.Add(sum(terms) == 0)
    # symmetry breaking: orbit representatives of E ordered by their first block
    for i in range(norb - 1): m.Add(Erep[i][0] <= Erep[i + 1][0])
    print(f'{nodes} K-conditions; building done', flush=True)
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = tl; sv.parameters.num_workers = workers
    sv.parameters.random_seed = seed
    t0 = time.time(); st = sv.Solve(m); name = sv.StatusName(st)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        print(f'sigma-invariant border ({mode}, NFIX {nfix}, case {idx}): {name} ({time.time()-t0:.0f}s)', flush=True); return
    Em = np.array([np.concatenate([F[sv.Value(r[c])] for c in range(4)]) for r in Erows])
    Zm = np.array([np.concatenate([F[sv.Value(r[c])] for c in range(4)]) for r in Zrows])
    Km = (Em @ M.T @ Zm.T) // (4 * k)
    okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=int))
    Fm = -Zm.T
    top = Km @ Km.T + s_ * (Em @ Em.T); cross = Km @ Fm.T + Em @ M.T; bot = s_ * (Fm @ Fm.T) + M @ M.T
    good = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * k, dtype=int))
    print(f'sigma-invariant border ({mode}, NFIX {nfix}, case {idx}): {name} ({time.time()-t0:.0f}s); K Hadamard {okK}; '
          f'compressed matrix of order {b + 4*k}: Gram = {n} I is {good}', flush=True)
    json.dump(dict(case=idx, mode=mode, nfix=nfix, T=T, s=s_, cells=sol['cells'], M=M.tolist(), K=Km.tolist(),
                   E=Em.tolist(), Z=Zm.tolist(), checks=dict(K_hadamard=bool(okK), compressed_gram=bool(good))),
              open(f'logs/fable_fresh/q599_sigma_border_{mode}_f{nfix}_case{idx}.json', 'w'))

if __name__ == '__main__':
    main()
