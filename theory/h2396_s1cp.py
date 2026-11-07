r"""Stage 1 of the Paley plug-in Legendre-quadruple slice on Z_n x F_q by CP-SAT (report/prompt_599_plugin.md, Part A).

Independent of theory/open_screen_stage1_exact.py (method C there uses one integer letter per cell; here the plug-in
variables are used).  Variables per block i and k in Z_n: x0 in {+-1} (Boolean b: x0 = 2b - 1) and the defect sign
e in {-1, 0, 1} (Booleans P, M with P + M <= 1, e = P - M).  u = x0 + (q - 1) e.  Constraints [P, all exact]:
  - for every lag d = 0..n/2:  sum_i sum_k u_i(k) u_i(k + d) = T(d),  T(0) = 4m + 4 - 4q, T(d != 0) = -4q  (m = qn);
    each product u(a) u(b) is a table constraint over the 36 letter pairs;
  - the allocation: sum of defects = e_plus + e_minus, and #{x0 = e at a defect} = e_plus (the defect identity
    (q + 1) e_plus + (q - 3) e_minus = 4(n - 1) is checked first);
  - redundant: row sums r_i = sum_k u_i(k) with sum r_i^2 = 4 (n even: one block +-2), and the order-2 character
    (n even) A_i = sum_k (-1)^k u_i(k) with sum A_i^2 = 4m + 4;
  - symmetry breaking (solutions map to solutions under block negation and block permutation): r_0 = 2 (n even) or
    r_i = 1 (n odd) and the blocks 1..3 ordered by their number of defects.
INFEASIBLE decides emptiness [M]; FEASIBLE returns a quadruple, re-checked from scratch; UNKNOWN decides nothing.
usage: python3 theory/h2396_s1cp.py n q e_plus e_minus [--secs S] [--workers W] [--seed s] [--nsol K] [--out FILE]"""
import sys, json, time, itertools


def target(n, q):
    m = q * n
    return [4 * m + 4 - 4 * q] + [-4 * q] * (n - 1)


def check(us, n, q):
    T = target(n, q)
    for d in range(n):
        if sum(u[k] * u[(k + d) % n] for u in us for k in range(n)) != T[d]: return False
    return True


def solve(n, q, ep, em, secs=600, workers=4, seed=0, exclude=()):
    from ortools.sat.python import cp_model
    assert (q + 1) * ep + (q - 3) * em == 4 * (n - 1), "allocation violates the defect identity"
    E = ep + em; L = sorted({a + (q - 1) * e for a in (1, -1) for e in (-1, 0, 1)})
    M = cp_model.CpModel()
    B = [[M.NewBoolVar(f"b{i}_{k}") for k in range(n)] for i in range(4)]
    P = [[M.NewBoolVar(f"P{i}_{k}") for k in range(n)] for i in range(4)]
    Mi = [[M.NewBoolVar(f"M{i}_{k}") for k in range(n)] for i in range(4)]
    U = [[M.NewIntVarFromDomain(cp_model.Domain.FromValues(L), f"u{i}_{k}") for k in range(n)] for i in range(4)]
    same = []                                   # defect with x0 = e  (type e_plus)
    for i in range(4):
        for k in range(n):
            M.AddBoolOr([P[i][k].Not(), Mi[i][k].Not()])
            M.Add(U[i][k] == 2 * B[i][k] - 1 + (q - 1) * (P[i][k] - Mi[i][k]))
            s = M.NewBoolVar(f"s{i}_{k}")       # s <=> (P and b) or (M and not b)
            a1 = M.NewBoolVar(""); a2 = M.NewBoolVar("")
            M.AddBoolAnd([P[i][k], B[i][k]]).OnlyEnforceIf(a1); M.AddBoolOr([P[i][k].Not(), B[i][k].Not()]).OnlyEnforceIf(a1.Not())
            M.AddBoolAnd([Mi[i][k], B[i][k].Not()]).OnlyEnforceIf(a2); M.AddBoolOr([Mi[i][k].Not(), B[i][k]]).OnlyEnforceIf(a2.Not())
            M.AddMaxEquality(s, [a1, a2])
            same.append(s)
    M.Add(sum(P[i][k] + Mi[i][k] for i in range(4) for k in range(n)) == E)
    M.Add(sum(same) == ep)
    pairs = [(a, b, a * b) for a in L for b in L]
    sq = sorted({a * a for a in L})
    T = target(n, q)
    for d in range(0, n // 2 + 1):
        terms = []
        for i in range(4):
            for k in range(n):
                if d == 0:
                    z = M.NewIntVarFromDomain(cp_model.Domain.FromValues(sq), "")
                    M.AddAllowedAssignments([U[i][k], z], [(a, a * a) for a in L])
                else:
                    z = M.NewIntVarFromDomain(cp_model.Domain.FromValues(sorted({p for _, _, p in pairs})), "")
                    M.AddAllowedAssignments([U[i][k], U[i][(k + d) % n], z], pairs)
                terms.append(z)
        M.Add(sum(terms) == T[d])
    # redundant rows / order-2 character
    R = [M.NewIntVar(-4 * n * q, 4 * n * q, f"r{i}") for i in range(4)]
    for i in range(4): M.Add(R[i] == sum(U[i]))
    if n % 2 == 0:
        M.Add(R[0] == 2)
        for i in range(1, 4): M.Add(R[i] == 0)
        A = [M.NewIntVar(-4 * n * q, 4 * n * q, f"A{i}") for i in range(4)]
        A2 = [M.NewIntVar(0, (4 * n * q) ** 2, f"A2{i}") for i in range(4)]
        for i in range(4):
            M.Add(A[i] == sum((-1) ** k * U[i][k] for k in range(n)))
            M.AddMultiplicationEquality(A2[i], [A[i], A[i]])
        M.Add(sum(A2) == 4 * q * n + 4)
    else:
        for i in range(4): M.Add(R[i] == 1)
    # blocks 1..3 have equal rows: order them by defect count
    Dc = [sum(P[i][k] + Mi[i][k] for k in range(n)) for i in range(4)]
    for i in (1, 2): M.Add(Dc[i] >= Dc[i + 1])
    for ex in exclude:                         # forbid earlier solutions (for --nsol)
        lits = []
        for i in range(4):
            for k in range(n):
                v = ex[i][k]; bv = M.NewBoolVar("")
                M.Add(U[i][k] != v).OnlyEnforceIf(bv); M.Add(U[i][k] == v).OnlyEnforceIf(bv.Not())
                lits.append(bv)
        M.AddBoolOr(lits)
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = secs
    S.parameters.num_search_workers = workers; S.parameters.random_seed = seed
    r = S.Solve(M)
    if r == cp_model.INFEASIBLE: return 'INFEASIBLE', None, S.WallTime()
    if r in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        us = [[S.Value(U[i][k]) for k in range(n)] for i in range(4)]
        assert check(us, n, q)
        return 'FEASIBLE', us, S.WallTime()
    return 'UNKNOWN', None, S.WallTime()


if __name__ == "__main__":
    a = sys.argv[1:]
    n, q, ep, em = map(int, a[:4])
    secs = float(a[a.index('--secs') + 1]) if '--secs' in a else 600
    workers = int(a[a.index('--workers') + 1]) if '--workers' in a else 4
    seed = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    nsol = int(a[a.index('--nsol') + 1]) if '--nsol' in a else 1
    outp = a[a.index('--out') + 1] if '--out' in a else None
    found = []
    t0 = time.time()
    while len(found) < nsol:
        st, us, wt = solve(n, q, ep, em, secs, workers, seed, exclude=found)
        rec = dict(n=n, q=q, alloc=[ep, em], status=st, wall=round(wt, 1), index=len(found))
        if us: rec['u'] = us; found.append(us)
        print(json.dumps(rec), flush=True)
        if outp: open(outp, 'a').write(json.dumps(rec) + "\n")
        if st != 'FEASIBLE': break
