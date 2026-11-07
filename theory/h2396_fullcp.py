r"""The whole Paley plug-in Legendre-quadruple system on Z_n x F_q by CP-SAT (report/prompt_599_plugin.md, gate A2).

A solution is four blocks (x0, x1, x2): Z_n -> {+-1}^3 with sum_i PAF_{X_i} = (4m + 4) delta - 4 on G = Z_n x F_q
(m = qn).  By theory/h2396_plugin_system.py [P, checked 66/66] this is equivalent to
   stage 1   sum_i PAF_{u_i} = (4m + 4) delta - 4q,        u = s + q e = x0 + (q - 1) e,
   (I)       sum_i (PAF_{s_i} + q PAF_{t_i}) = (4m + 4) delta,
   (Q)       sum_i (CCF(s_i, t_i) + eps CCF(t_i, s_i)) = 0,  eps = chi(-1),
with e = (x1 + x2)/2, t = (x1 - x2)/2, s = x0 - e.  Model: one integer state z in 0..7 per (block, k) encoding
(x0, x1, x2); for every block, lag d = 0..n/2 (all d for (Q) when eps = -1 is not needed: (Q) is antisymmetric then)
and k a table constraint over (z(k), z(k+d)) gives the products u u', s s', t t', s t', t s'.  Redundant exact
constraints: the defect allocation ((q+1) e_plus + (q-3) e_minus = 4(n-1)), the row sums (sum r_i^2 = 4), the trivial
character of K for t (t^_i = sum t_i: sum (s^_i^2 + q t^_i^2) = 4m + 4, sum s^_i t^_i = 0: the Z_q level), and for even n
the order-2 character.  Symmetry breaking: block negation and permutation (block 0 carries the row sum 2 for even n).
INFEASIBLE decides that the analogue has no plug-in core [M]; FEASIBLE returns a core, checked at the G level by
direct PAF (code shared with nothing else); UNKNOWN decides nothing.
(The real-character equation sum s^ t^ = 0 holds only for eps = +1; for eps = -1, (Q) is antisymmetric and gives
nothing at a real character.)
usage: python3 theory/h2396_fullcp.py n q [--secs S] [--workers W] [--seed s] [--alloc ep,em] [--out FILE]"""
import sys, json, time, itertools


def legendre(q):
    sq = {(x * x) % q for x in range(1, q)}
    return [0] + [1 if y in sq else -1 for y in range(1, q)]


def g_check(n, q, X):
    """X[i] = (x0, x1, x2) lists on Z_n; direct G-level PAF on Z_n x F_q against (4m+4) delta - 4"""
    chi = legendre(q); m = n * q
    f = [[(X[i][0][k] if z == 0 else (X[i][1][k] if chi[z] == 1 else X[i][2][k])) for k in range(n) for z in range(q)]
         for i in range(4)]
    idx = lambda k, z: k * q + z
    for h in range(n):
        for w in range(q):
            tot = 0
            for i in range(4):
                for k in range(n):
                    for z in range(q):
                        tot += f[i][idx(k, z)] * f[i][idx((k + h) % n, (z + w) % q)]
            if tot != (4 * m if (h == 0 and w == 0) else -4): return False
    return True


def allocations(q, n):
    R = 4 * (n - 1)
    return [(a, (R - (q + 1) * a) // (q - 3)) for a in range(R // (q + 1) + 1) if (R - (q + 1) * a) % (q - 3) == 0]


def solve(n, q, secs=600, workers=4, seed=0, alloc=None):
    from ortools.sat.python import cp_model
    eps = legendre(q)[q - 1]; m = n * q
    states = list(itertools.product((1, -1), repeat=3))           # (x0, x1, x2)
    def dec(st):
        x0, x1, x2 = st; e = (x1 + x2) // 2; t = (x1 - x2) // 2; s = x0 - e
        return dict(u=s + q * e, s=s, t=t, e=e)
    D = [dec(st) for st in states]
    als = allocations(q, n) if alloc is None else [alloc]
    if not als: return 'INFEASIBLE', None, 0.0, 'no allocation'
    M = cp_model.CpModel()
    Z = [[M.NewIntVar(0, 7, f"z{i}_{k}") for k in range(n)] for i in range(4)]
    def lin(i, k, fld):
        v = M.NewIntVar(-2 * q, 2 * q, ""); M.AddElement(Z[i][k], [D[j][fld] for j in range(8)], v); return v
    U = [[lin(i, k, 'u') for k in range(n)] for i in range(4)]
    S = [[lin(i, k, 's') for k in range(n)] for i in range(4)]
    T = [[lin(i, k, 't') for k in range(n)] for i in range(4)]
    E = [[lin(i, k, 'e') for k in range(n)] for i in range(4)]
    # defect allocation
    isdef = [[M.NewBoolVar("") for k in range(n)] for i in range(4)]
    isplus = [[M.NewBoolVar("") for k in range(n)] for i in range(4)]
    for i in range(4):
        for k in range(n):
            M.AddAllowedAssignments([Z[i][k], isdef[i][k], isplus[i][k]],
                                    [(j, int(D[j]['e'] != 0), int(D[j]['e'] != 0 and D[j]['s'] == 0)) for j in range(8)])
    nd = sum(isdef[i][k] for i in range(4) for k in range(n)); npl = sum(isplus[i][k] for i in range(4) for k in range(n))
    if len(als) == 1:
        M.Add(nd == als[0][0] + als[0][1]); M.Add(npl == als[0][0])
    else:
        ch = [M.NewBoolVar("") for _ in als]; M.AddExactlyOne(ch)
        for c, (a, b) in zip(ch, als):
            M.Add(nd == a + b).OnlyEnforceIf(c); M.Add(npl == a).OnlyEnforceIf(c)
    # products by table over (z(k), z(k+d))
    tab = [(a, b, D[a]['u'] * D[b]['u'], D[a]['s'] * D[b]['s'], D[a]['t'] * D[b]['t'], D[a]['s'] * D[b]['t'],
            D[a]['t'] * D[b]['s']) for a in range(8) for b in range(8)]
    lagset = range(0, n // 2 + 1)
    for d in lagset:
        su = []; ss = []; tt = []; qq = []
        for i in range(4):
            for k in range(n):
                vs = [M.NewIntVar(-4 * q * q, 4 * q * q, "") for _ in range(5)]
                if d == 0: M.AddAllowedAssignments([Z[i][k]] + vs, [(a,) + r[2:] for r in tab for a in [r[0]] if r[0] == r[1]])
                else: M.AddAllowedAssignments([Z[i][k], Z[i][(k + d) % n]] + vs, tab)
                su.append(vs[0]); ss.append(vs[1]); tt.append(vs[2]); qq.append(vs[3] + eps * vs[4])
        M.Add(sum(su) == (4 * m + 4 - 4 * q if d == 0 else -4 * q))
        M.Add(sum(ss) + q * sum(tt) == (4 * m + 4 if d == 0 else 0))
        M.Add(sum(qq) == 0)
    if eps == -1:                       # (Q) antisymmetric: lags n/2+1..n-1 are the negatives, nothing to add
        pass
    # redundant: rows, Z_q level, order-2 character
    big = 4 * n * q
    R = [M.NewIntVar(-big, big, "") for i in range(4)]
    Sh = [M.NewIntVar(-big, big, "") for i in range(4)]
    Th = [M.NewIntVar(-big, big, "") for i in range(4)]
    for i in range(4):
        M.Add(R[i] == sum(U[i])); M.Add(Sh[i] == sum(S[i])); M.Add(Th[i] == sum(T[i]))
    sq = lambda v: (lambda z: (M.AddMultiplicationEquality(z, [v, v]), z)[1])(M.NewIntVar(0, big * big, ""))
    pr = lambda a, b: (lambda z: (M.AddMultiplicationEquality(z, [a, b]), z)[1])(M.NewIntVar(-big * big, big * big, ""))
    M.Add(sum(sq(R[i]) for i in range(4)) == 4)
    M.Add(sum(sq(Sh[i]) for i in range(4)) + q * sum(sq(Th[i]) for i in range(4)) == 4 * m + 4)
    if eps == 1: M.Add(sum(pr(Sh[i], Th[i]) for i in range(4)) == 0)   # (Q) at a real character: (1 + eps) sum s^ t^ = 0
    if n % 2 == 0:
        M.Add(R[0] == 2)
        for i in range(1, 4): M.Add(R[i] == 0)
        Sp = [M.NewIntVar(-big, big, "") for i in range(4)]; Tp = [M.NewIntVar(-big, big, "") for i in range(4)]
        Up = [M.NewIntVar(-big, big, "") for i in range(4)]
        for i in range(4):
            M.Add(Sp[i] == sum((-1) ** k * S[i][k] for k in range(n)))
            M.Add(Tp[i] == sum((-1) ** k * T[i][k] for k in range(n)))
            M.Add(Up[i] == sum((-1) ** k * U[i][k] for k in range(n)))
        M.Add(sum(sq(Up[i]) for i in range(4)) == 4 * m + 4)
        M.Add(sum(sq(Sp[i]) for i in range(4)) + q * sum(sq(Tp[i]) for i in range(4)) == 4 * m + 4)
        if eps == 1: M.Add(sum(pr(Sp[i], Tp[i]) for i in range(4)) == 0)
        Dc = [sum(isdef[i][k] for k in range(n)) for i in range(4)]
        for i in (1, 2): M.Add(Dc[i] >= Dc[i + 1])
    else:
        for i in range(4): M.Add(R[i] == 1)
        Dc = [sum(isdef[i][k] for k in range(n)) for i in range(4)]
        for i in range(3): M.Add(Dc[i] >= Dc[i + 1])
    Sv = cp_model.CpSolver(); Sv.parameters.max_time_in_seconds = secs
    Sv.parameters.num_search_workers = workers; Sv.parameters.random_seed = seed
    r = Sv.Solve(M)
    if r == cp_model.INFEASIBLE: return 'INFEASIBLE', None, Sv.WallTime(), None
    if r in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        X = [[[states[Sv.Value(Z[i][k])][c] for k in range(n)] for c in range(3)] for i in range(4)]
        assert g_check(n, q, X), "CP-SAT core fails the G-level check"
        return 'FEASIBLE', X, Sv.WallTime(), None
    return 'UNKNOWN', None, Sv.WallTime(), None


if __name__ == "__main__":
    a = sys.argv[1:]
    n, q = int(a[0]), int(a[1])
    secs = float(a[a.index('--secs') + 1]) if '--secs' in a else 600
    workers = int(a[a.index('--workers') + 1]) if '--workers' in a else 4
    seed = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    alloc = tuple(map(int, a[a.index('--alloc') + 1].split(','))) if '--alloc' in a else None
    outp = a[a.index('--out') + 1] if '--out' in a else None
    st, X, wt, why = solve(n, q, secs, workers, seed, alloc)
    rec = dict(n=n, q=q, v=n * q + 1, alloc=alloc, allocations=allocations(q, n), status=st, wall=round(wt, 1))
    if why: rec['why'] = why
    if X: rec['core_x0x1x2'] = X
    print(json.dumps(rec), flush=True)
    if outp: open(outp, 'a').write(json.dumps(rec) + "\n")
