#!/usr/bin/env python3
"""Exact Q-level gate for partial character borders over Z_2^a x Z_s (report section 7), all T at once.

Unknowns: the 4 * 2^a cell sums y_k(q), each in the set Y of sums of a U-invariant +-1 function over a coset of Z_s.
Constraints: for every character rho of Z_2^a,  sum_k yhat_k(rho)^2 = n  or  = 4t, and exactly t characters take
the value 4t  (t = p mod 2^a, n = 4p).  CP-SAT with integer squares; the set T is read off a solution.
INFEASIBLE = no core with this multiplier group for ANY T [M].  FEASIBLE says only that the Q-level lives.
Usage: fable_qlevel.py p a g time_limit [forbid_classical]     (g generates the multiplier group on Z_s)
"""
import sys, time, itertools
from ortools.sat.python import cp_model
sys.path.insert(0, 'theory')
from fable_tcensus import orbits
from fable_tclasses import kind

def main():
    p = int(sys.argv[1]); a = int(sys.argv[2]); g = int(sys.argv[3]); tl = float(sys.argv[4])
    forbid = len(sys.argv) > 5 and sys.argv[5] == '1'
    k = 2 ** a; t = p % k; s = (p - t) // k; n = 4 * p
    sizes = orbits(s, g)
    Y = {0}
    for z in sizes: Y = {y + z for y in Y} | {y - z for y in Y}
    Y = sorted(Y)
    E = ((k - t) * n + 4 * t * t) // k
    print(f'p={p} a={a} t={t} s={s} g={g}: {len(sizes)} orbits per coset, r={k*len(sizes)}, |Y|={len(Y)}, sum y^2 = {E}', flush=True)
    m = cp_model.CpModel()
    y = [[m.NewIntVarFromDomain(cp_model.Domain.FromValues(Y), '') for q in range(k)] for b in range(4)]
    ymax = max(abs(v) for v in Y)
    chars = [[(-1) ** bin(r & q).count('1') for q in range(k)] for r in range(k)]
    low = [m.NewBoolVar('') for r in range(k)]
    W = k * ymax
    for r in range(k):
        sqs = []
        for b in range(4):
            w = m.NewIntVar(-W, W, '')
            m.Add(w == sum(chars[r][q] * y[b][q] for q in range(k)))
            w2 = m.NewIntVar(0, n, '')
            m.AddMultiplicationEquality(w2, [w, w])
            sqs.append(w2)
        m.Add(sum(sqs) == n - (n - 4 * t) * low[r])
    m.Add(sum(low) == t)
    # Parseval (redundant, helps propagation)
    y2 = []
    for b in range(4):
        for q in range(k):
            v = m.NewIntVar(0, ymax * ymax, '')
            m.AddMultiplicationEquality(v, [y[b][q], y[b][q]])
            y2.append(v)
    m.Add(sum(y2) == E)
    # symmetry: block signs
    for b in range(4): m.Add(y[b][0] >= 0) if 0 not in Y else None
    if forbid:
        # forbid T of classical kind: exclude by no-good cuts found iteratively
        pass
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = tl; sol.parameters.num_workers = 4
    t0 = time.time(); found = []
    while True:
        r_ = sol.Solve(m); st = sol.StatusName(r_)
        if r_ not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            print(f'  {st} ({time.time()-t0:.0f}s)' + (f' after {len(found)} character sets' if found else ''), flush=True); break
        T = [r for r in range(k) if sol.Value(low[r])]
        kd = kind(T, k)
        yy = [[sol.Value(y[b][q]) for q in range(k)] for b in range(4)]
        print(f'  FEASIBLE ({time.time()-t0:.0f}s): T = {T} [{kd}]; compressions {yy}', flush=True)
        found.append(T)
        if not forbid or len(found) >= 6: break
        # forbid this T and continue
        m.AddBoolOr([low[r].Not() for r in T])
        remaining = tl - (time.time() - t0)
        if remaining <= 5: break
        sol.parameters.max_time_in_seconds = remaining

if __name__ == '__main__':
    main()
