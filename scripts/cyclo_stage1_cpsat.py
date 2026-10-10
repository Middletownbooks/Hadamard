"""EXACT stage 1 by CP-SAT for a generalised-Paley slice (b = 1): Z_F fixed pointwise, index-s classes on Z_q.
Find y_k in L^F (L = {a + f * (sum of s signs)}, f = (q - 1)/s), k = 0..3, with
    sum_k PAF_{Z_F}(y_k)(h) = T(h)   (T(0) = 4m - 4(q - 1), T(h != 0) = -4q),  rows (2, 0, 0, 0).
(The row sums follow from the lag equations up to order and sign; fixing them breaks that symmetry.)
Products y_{k,g} y_{k,g+h} via AddMultiplicationEquality.  Lags h = 1..F/2 suffice (PAF symmetric); lag 0 included.
Status FEASIBLE/OPTIMAL gives a verified solution; INFEASIBLE is a proof (CP-SAT is exact); UNKNOWN = undecided.
usage: cyclo_stage1_cpsat.py F q s [--time 600] [--workers 4] [--seed 0] [--fix-comp n13,n15,...]"""
import sys, json, time, itertools as it
import numpy as np
from ortools.sat.python import cp_model

def solve(F, q, s, time_limit=600, workers=4, seed=0, extra_lag0=True):
    f = (q - 1) // s; m = F * q
    L = sorted({a + f * sum(c) for a in (1, -1) for c in it.product((1, -1), repeat=s)})
    T = [4 * m - 4 * (q - 1)] + [-4 * q] * (F // 2)
    M = cp_model.CpModel()
    dom = cp_model.Domain.FromValues(L)
    y = [[M.NewIntVarFromDomain(dom, f"y{k}_{g}") for g in range(F)] for k in range(4)]
    big = max(abs(x) for x in L) ** 2
    sq_dom = cp_model.Domain.FromValues(sorted({a * b for a in L for b in L}))
    for h in range(F // 2 + 1):
        terms = []
        for k in range(4):
            for g in range(F):
                if h == 0:
                    p = M.NewIntVarFromDomain(cp_model.Domain.FromValues(sorted({a * a for a in L})), f"p{k}_{g}_{h}")
                    M.AddMultiplicationEquality(p, [y[k][g], y[k][g]])
                else:
                    p = M.NewIntVarFromDomain(sq_dom, f"p{k}_{g}_{h}")
                    M.AddMultiplicationEquality(p, [y[k][g], y[k][(g + h) % F]])
                terms.append(p)
        M.Add(sum(terms) == T[h])
    rows = [2, 0, 0, 0]
    for k in range(4): M.Add(sum(y[k]) == rows[k])
    # symmetry: translate block 0 so that its largest-|y| entry ... (kept simple: none)
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = time_limit; S.parameters.num_workers = workers
    S.parameters.random_seed = seed
    t0 = time.time(); st = S.Solve(M); secs = round(time.time() - t0, 1)
    name = S.StatusName(st)
    out = dict(F=F, q=q, s=s, letters=L, status=name, secs=secs)
    if name in ("FEASIBLE", "OPTIMAL"):
        Y = [[S.Value(y[k][g]) for g in range(F)] for k in range(4)]
        P = [sum(sum(Y[k][g] * Y[k][(g + h) % F] for g in range(F)) for k in range(4)) for h in range(F)]
        out.update(y=Y, verified=(P[0] == T[0] and all(P[h] == -4 * q for h in range(1, F))))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--time": "600", "--workers": "4", "--seed": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos
    print(json.dumps(solve(F, q, s, float(opt["--time"]), int(opt["--workers"]), int(opt["--seed"]))), flush=True)
