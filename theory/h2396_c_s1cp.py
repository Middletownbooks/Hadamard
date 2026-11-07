r"""Part C, gate C2: stage 1 (compression along F_q to K = Z_nK) and the Z_q level (compression along K to F_q) of
Kumar's index-4 Szolloesi plug-in, jointly, by CP-SAT [M when INFEASIBLE; FEASIBLE gives instances].

Credit: the index-r plug-in class is K. Kumar's (search/NOTES.md section 82; Apache-2.0).  Conventions as his
verify.py (_szollosi_zero_diagonal_core): A, B symmetric, A(0) = B(0) = 0, +-1 elsewhere; C, D +-1;
sum PAF = (4N - 2) delta on G = Z_nK x F_q, N = nK q.  New code.

Variables per block i, k in Z_nK: x0(k) in {+-1} (0 at k = 0 for A, B), class values x_j(k) in {+-1}, j = 0..3
(C_j = g^j S, S the index-4 subgroup of F_q^*; -1 in S, so A, B symmetric <=> x0, x_j even in k).  Then
  stage 1:   y = x0 + L sigma, sigma = sum_j x_j, L = (q - 1)/4;  sum_i PAF_{y_i}(d) = (4N - 2) delta(d), d = 0..nK/2;
  Z_q level: a = sum_k x0(k), W_j = sum_k x_j(k);  sum_i [a_i^2 + L sum_j W_ij^2] = 4N - 2 and, for l = 0..3,
             sum_i [2 a_i W_il + sum_{j,j'} W_ij W_ij' (j - l, j' - l)] = 0   ((i, j) cyclotomic numbers of order 4).
Products y(k) y(k + d) by table constraints over the letters; aggregate products by AddMultiplicationEquality.
Symmetry breaking: none beyond fixing nothing (sound).  Each FEASIBLE answer is re-checked: stage 1 by direct PAF on
Z_nK and the Z_q level by direct PAF on F_q of the compressed functions.
usage: python3 theory/h2396_c_s1cp.py nK q [--secs S] [--workers W] [--seed s] [--nozq] [--out FILE]"""
import sys, json, time


def classes(q, r=4):
    g = next(x for x in range(2, q) if all(pow(x, (q - 1) // p, q) != 1 for p in {2, 3, 5, 7, 11, 13} if (q - 1) % p == 0))
    S = sorted({pow(x, r, q) for x in range(1, q)})
    cls = {}
    for j in range(r):
        for s in S: cls[(pow(g, j, q) * s) % q] = j
    cyc = [[sum(1 for y in range(1, q) if cls[y] == i and (y + 1) % q and cls[(y + 1) % q] == j) for j in range(r)] for i in range(r)]
    return S, cls, cyc


def solve(nK, q, secs=600, workers=4, seed=0, zq=True):
    from ortools.sat.python import cp_model
    S, cls, cyc = classes(q); L = (q - 1) // 4; N = nK * q; T0 = 4 * N - 2
    assert (q - 1) % 4 == 0 and (q - 1) // 4 in [len(S)] and q - 1 in S
    M = cp_model.CpModel()
    X0 = [[None] * nK for _ in range(4)]; XJ = [[[None] * 4 for _ in range(nK)] for _ in range(4)]
    for i in range(4):
        for k in range(nK):
            sym = i < 2 and k > nK // 2
            if sym:
                X0[i][k] = X0[i][nK - k]; XJ[i][k] = XJ[i][nK - k]; continue
            if i < 2 and k == 0: X0[i][k] = 0
            else:
                b = M.NewBoolVar(""); v = M.NewIntVar(-1, 1, ""); M.Add(v == 2 * b - 1); X0[i][k] = v
            for j in range(4):
                b = M.NewBoolVar(""); v = M.NewIntVar(-1, 1, ""); M.Add(v == 2 * b - 1); XJ[i][k][j] = v
    letters = sorted({a + L * s for a in (1, -1, 0) for s in range(-4, 5, 2)})
    Y = [[None] * nK for _ in range(4)]
    for i in range(4):
        for k in range(nK):
            if i < 2 and k > nK // 2: Y[i][k] = Y[i][nK - k]; continue
            yv = M.NewIntVarFromDomain(cp_model.Domain.FromValues(letters), "")
            M.Add(yv == X0[i][k] + L * sum(XJ[i][k])); Y[i][k] = yv
    pairs = [(a, b, a * b) for a in letters for b in letters]
    pv = sorted({p for _, _, p in pairs})
    for d in range(0, nK // 2 + 1):
        terms = []
        for i in range(4):
            for k in range(nK):
                z = M.NewIntVarFromDomain(cp_model.Domain.FromValues(pv), "")
                if d == 0 or Y[i][k] is Y[i][(k + d) % nK]:      # y(k) y(-k) = y(k)^2 for a symmetric block
                    M.AddAllowedAssignments([Y[i][k], z], [(a, a * a) for a in letters])
                else: M.AddAllowedAssignments([Y[i][k], Y[i][(k + d) % nK], z], pairs)
                terms.append(z)
        M.Add(sum(terms) == (T0 if d == 0 else 0))
    big = nK * 4
    R = []
    for i in range(4):
        r = M.NewIntVar(-big * (L * 4 + 1), big * (L * 4 + 1), ""); M.Add(r == sum(Y[i])); R.append(r)
    R2 = []
    for r in R:
        z = M.NewIntVar(0, (big * (4 * L + 1)) ** 2, ""); M.AddMultiplicationEquality(z, [r, r]); R2.append(z)
    M.Add(sum(R2) == T0)
    if zq:
        Aa = []; Ww = []
        for i in range(4):
            a = M.NewIntVar(-nK, nK, ""); M.Add(a == sum(X0[i][k] for k in range(nK) if not isinstance(X0[i][k], int)) +
                                                sum(0 for k in range(nK) if isinstance(X0[i][k], int))); Aa.append(a)
            Ws = []
            for j in range(4):
                w = M.NewIntVar(-nK, nK, ""); M.Add(w == sum(XJ[i][k][j] for k in range(nK))); Ws.append(w)
            Ww.append(Ws)
        def prod(u, v):
            z = M.NewIntVar(-nK * nK, nK * nK, ""); M.AddMultiplicationEquality(z, [u, v]); return z
        e0 = []
        for i in range(4):
            e0.append(prod(Aa[i], Aa[i])); e0 += [L * prod(Ww[i][j], Ww[i][j]) for j in range(4)]
        M.Add(sum(e0) == T0)
        WW = [[[prod(Ww[i][j], Ww[i][jj]) for jj in range(4)] for j in range(4)] for i in range(4)]
        for l in range(4):
            terms = []
            for i in range(4):
                terms.append(2 * prod(Aa[i], Ww[i][l]))
                for j in range(4):
                    for jj in range(4):
                        c = cyc[(j - l) % 4][(jj - l) % 4]
                        if c: terms.append(c * WW[i][j][jj])
            M.Add(sum(terms) == 0)
    Sv = cp_model.CpSolver(); Sv.parameters.max_time_in_seconds = secs
    Sv.parameters.num_search_workers = workers; Sv.parameters.random_seed = seed
    st = Sv.Solve(M)
    if st == cp_model.INFEASIBLE: return 'INFEASIBLE', None, Sv.WallTime()
    if st not in (cp_model.FEASIBLE, cp_model.OPTIMAL): return 'UNKNOWN', None, Sv.WallTime()
    val = lambda v: v if isinstance(v, int) else Sv.Value(v)
    sol = dict(x0=[[val(X0[i][k]) for k in range(nK)] for i in range(4)],
               xj=[[[val(XJ[i][k][j]) for j in range(4)] for k in range(nK)] for i in range(4)])
    assert check(nK, q, sol, zq)
    return 'FEASIBLE', sol, Sv.WallTime()


def check(nK, q, sol, zq=True):
    S, cls, cyc = classes(q); L = (q - 1) // 4; N = nK * q
    y = [[sol['x0'][i][k] + L * sum(sol['xj'][i][k]) for k in range(nK)] for i in range(4)]
    for d in range(nK):
        if sum(y[i][k] * y[i][(k + d) % nK] for i in range(4) for k in range(nK)) != (4 * N - 2 if d == 0 else 0): return False
    for i in range(2):
        if sol['x0'][i][0] != 0: return False
        for k in range(nK):
            if sol['x0'][i][k] != sol['x0'][i][(-k) % nK] or sol['xj'][i][k] != sol['xj'][i][(-k) % nK]: return False
    if zq:
        for z in range(q):
            tot = 0
            for i in range(4):
                w = [sum(sol['x0'][i])] + [sum(sol['xj'][i][k][cls[t]] for k in range(nK)) for t in range(1, q)]
                tot += sum(w[t] * w[(t + z) % q] for t in range(q))
            if tot != (4 * N - 2 if z == 0 else 0): return False
    return True


if __name__ == "__main__":
    a = sys.argv[1:]
    nK, q = int(a[0]), int(a[1])
    secs = float(a[a.index('--secs') + 1]) if '--secs' in a else 600
    workers = int(a[a.index('--workers') + 1]) if '--workers' in a else 4
    seed = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    zq = '--nozq' not in a
    outp = a[a.index('--out') + 1] if '--out' in a else None
    st, sol, wt = solve(nK, q, secs, workers, seed, zq)
    rec = dict(nK=nK, q=q, N=nK * q, Zq_level=zq, status=st, wall=round(wt, 1))
    if sol: rec['solution'] = sol
    print(json.dumps(rec)); sys.stdout.flush()
    if outp: open(outp, 'a').write(json.dumps(rec) + "\n")
