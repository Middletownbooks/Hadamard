#!/usr/bin/env python3
"""CP-SAT version of theory/fable_tmult.py: formal multiplication of typed GS quadruples by Z_t.

Variables: one letter (sign x symbol) per slot (j, tau), j = 0..3, tau in Z_t.
Constraints: for d = 1..(t-1)/2 and every formal monomial, the signed number of slot pairs
(tau+d, tau) producing it is zero; every index k occupies exactly t slots (with type multiplicity).
Symmetry breaking: slot (0,0) carries +A (plain) when the type allows.
Usage: fable_tmult_sat.py t type [time_limit_s] [workers]
"""
import sys, itertools, time
from ortools.sat.python import cp_model
sys.path.insert(0, 'theory')
from fable_tmult import symbols_for_type

def solve(t, typ, tl=600, workers=4, log=False):
    sym, conj, letters = symbols_for_type(typ)
    nsym = len(conj)
    mono = {}
    for a in range(nsym):
        for b in range(a, nsym):
            mono[(a, b)] = len(mono)
    L = [(e, l) for e in (1, -1) for l in range(len(letters))]
    nL = len(L)
    m = cp_model.CpModel()
    x = {(j, tau, i): m.NewBoolVar(f"x{j}_{tau}_{i}") for j in range(4) for tau in range(t) for i in range(nL)}
    for j in range(4):
        for tau in range(t):
            m.AddExactlyOne(x[j, tau, i] for i in range(nL))
    # counts
    cls = {}
    for s in range(nsym):
        cls.setdefault(min(s, conj[s]), len(cls))
    target = [0] * len(cls)
    for k in range(4):
        s = sym[(k, 0)]
        target[cls[min(s, conj[s])]] += t
    for c in range(len(cls)):
        m.Add(sum(x[j, tau, i] for j in range(4) for tau in range(t) for i in range(nL)
                  if cls[min(sym[letters[L[i][1]]], conj[sym[letters[L[i][1]]]])] == c) == target[c])
    half = (t - 1) // 2
    for d in range(1, half + 1):
        terms = {mm: [] for mm in mono.values()}
        for j in range(4):
            for tau in range(t):
                for i1 in range(nL):
                    for i0 in range(nL):
                        e1, l1 = L[i1]; e0, l0 = L[i0]
                        a = sym[letters[l1]]; b = conj[sym[letters[l0]]]
                        if a > b: a, b = b, a
                        z = m.NewBoolVar('')
                        m.AddBoolAnd([x[j, (tau + d) % t, i1], x[j, tau, i0]]).OnlyEnforceIf(z)
                        m.AddBoolOr([x[j, (tau + d) % t, i1].Not(), x[j, tau, i0].Not(), z])
                        terms[mono[(a, b)]].append((e1 * e0, z))
        for mm, lst in terms.items():
            if lst:
                m.Add(sum(c * z for c, z in lst) == 0)
    m.Add(x[0, 0, 0] == 1)       # +A plain at slot (0,0): sign and (if A is not constrained) nothing else
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl
    s.parameters.num_workers = workers
    s.parameters.log_search_progress = log
    t0 = time.time()
    r = s.Solve(m)
    st = s.StatusName(r)
    print(f"t={t} type={typ}: {st} in {time.time()-t0:.1f}s", flush=True)
    if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        names = 'ABCD'
        for j in range(4):
            row = []
            for tau in range(t):
                for i in range(nL):
                    if s.Value(x[j, tau, i]):
                        e, l = L[i]
                        row.append(('+' if e > 0 else '-') + names[letters[l][0]] + ('*' if letters[l][1] else ' '))
            print('   ', ' '.join(row))
    return st

if __name__ == '__main__':
    t = int(sys.argv[1]); typ = sys.argv[2]
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    w = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    solve(t, typ, tl, w)
