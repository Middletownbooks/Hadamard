r"""Stage 1 by CP-SAT with a fixed lag-0 composition (report/paley_523_answer.md) [M when FEASIBLE].

Same model as scripts/cyclo_stage1_cpsat.py (y_k in the letter set L on Z_F, sum_k PAF(y_k) = T_J, rows (2, 0, 0, 0)),
plus: the number of entries with |y| = a is fixed for every letter magnitude a (the lag-0 equation only fixes
sum n_a a^2; fixing the n_a removes most of the combinatorial freedom), optional per-block counts, and an optional
symmetry break (block 0 starts with its largest |y|).  Letters are encoded by one boolean per (position, letter).
usage: python3 theory/q523_stage1.py F q s --comp a1:n1,a2:n2,... [--time 600] [--workers 4] [--seed 0]
       python3 theory/q523_stage1.py F q s --list-comps        (all lag-0 compositions, by multinomial weight)"""
import sys, json, time, math, itertools as it
import numpy as np
from ortools.sat.python import cp_model

def letters(q, s):
    f = (q - 1) // s
    return sorted({a + f * sum(c) for a in (1, -1) for c in it.product((1, -1), repeat=s)}), f

def compositions(F, q, s):
    L, f = letters(q, s); m = F * q; T0 = 4 * m - 4 * (q - 1)
    mags = sorted({abs(x) for x in L}); big = [a for a in mags if a > 1]
    out = []
    def rec(i, left, rem, cur):
        if i == len(big):
            if rem == left: out.append({**cur, 1: left})
            return
        a = big[i]; c = 0
        while c <= left and c * a * a <= rem:
            rec(i + 1, left - c, rem - c * a * a, {**cur, a: c}); c += 1
    rec(0, 4 * F, T0, {})
    # weight: multinomial count of placements times sign choices (2^{4F}) -- log2
    def lw(c):
        n = 4 * F; r = math.lgamma(n + 1) - sum(math.lgamma(v + 1) for v in c.values())
        return r / math.log(2) + n
    return sorted(out, key=lambda c: -lw(c)), lw

def solve(F, q, s, comp, time_limit=600, workers=4, seed=0):
    L, f = letters(q, s); m = F * q
    T = [4 * m - 4 * (q - 1)] + [-4 * q] * (F // 2)
    M = cp_model.CpModel()
    b = [[[M.NewBoolVar(f"b{k}_{g}_{i}") for i in range(len(L))] for g in range(F)] for k in range(4)]
    y = [[M.NewIntVar(-max(L), max(L), f"y{k}_{g}") for g in range(F)] for k in range(4)]
    for k in range(4):
        for g in range(F):
            M.AddExactlyOne(b[k][g]); M.Add(y[k][g] == sum(L[i] * b[k][g][i] for i in range(len(L))))
    for a, n in comp.items():
        M.Add(sum(b[k][g][i] for k in range(4) for g in range(F) for i in range(len(L)) if abs(L[i]) == a) == n)
    sq = cp_model.Domain.FromValues(sorted({x * z for x in L for z in L}))
    for h in range(1, F // 2 + 1):
        terms = []
        for k in range(4):
            for g in range(F):
                p = M.NewIntVarFromDomain(sq, ""); M.AddMultiplicationEquality(p, [y[k][g], y[k][(g + h) % F]]); terms.append(p)
        M.Add(sum(terms) == T[h])
    for k, r in enumerate((2, 0, 0, 0)): M.Add(sum(y[k]) == r)
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = time_limit; S.parameters.num_workers = workers
    S.parameters.random_seed = seed
    t0 = time.time(); st = S.Solve(M); name = S.StatusName(st)
    out = dict(F=F, q=q, s=s, comp={str(k): v for k, v in comp.items()}, status=name, secs=round(time.time() - t0, 1))
    if name in ("FEASIBLE", "OPTIMAL"):
        Y = [[S.Value(y[k][g]) for g in range(F)] for k in range(4)]
        P = [sum(sum(Y[k][g] * Y[k][(g + h) % F] for g in range(F)) for k in range(4)) for h in range(F)]
        out.update(y=Y, verified=(P[0] == T[0] and all(P[h] == -4 * q for h in range(1, F)) and [sum(r) for r in Y] == [2, 0, 0, 0]))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--time": "600", "--workers": "4", "--seed": "0", "--comp": ""}; pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith('--'): flags.add(a[i]); i += 1
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos
    if '--list-comps' in flags:
        cs, lw = compositions(F, q, s)
        for c in cs[:40]: print(json.dumps(dict(comp=c, log2_weight=round(lw(c), 1))))
        sys.exit()
    comp = {int(t.split(':')[0]): int(t.split(':')[1]) for t in opt["--comp"].split(',')}
    print(json.dumps(solve(F, q, s, comp, float(opt["--time"]), int(opt["--workers"]), int(opt["--seed"]))), flush=True)
