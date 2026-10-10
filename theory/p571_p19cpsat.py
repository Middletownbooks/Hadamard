r"""Exact stage 1 (and stage 2) of the <61> = Paley(19) slice of v = 571 per gadget class, with OR-tools CP-SAT
(report section 7).  Input: the gadget classes of theory/p571_p19gadgets.py.

Stage 1 model for a class (W, eps): Booleans b_k(g) on the twisted cosets (A = 2b - 1), A = eps W on gadgets;
  row sums  sum_g A_k(g) = r_k - 18 t_k                      (r = (2, 0, 0, 0));
  lags      sum_k [PAF(A_k)(h) + 18 X_k(h)] = -76 - 324 s(h),   h = 1..15,
with A(g) A(g') = 1 - 2 XOR(b_g, b_g') for two twisted cosets (one XOR literal per pair), and
  Lemma 6 at psi = (-1)^g (necessary for stage 2):  sum_k What_k (2 Ahat_k + 17 What_k) = sum_k shat_k^2, shat_k = 30 - K_k
  (mod 2), |shat_k| <= 30 - K_k  (hats at psi; a domain constraint on a linear expression).
INFEASIBLE for every class proves the slice empty (S1 and Lemma 6 are necessary for a core).
Stage 2 (for any stage-1 solution): u_k(g) = +-1 on the twisted cosets with
  sum_k PAF(u_k)(h) = 4 + x_h + 17 s(h)  (h = 1..15),   sum_k sum_g u_k(g) [beta_k(g - h) - beta_k(g + h)] = 0 (h = 1..14),
beta = A on twisted cosets, A - W on gadgets (real and imaginary parts of the S2 condition); then the core is built
and checked by theory/paley_b1pipe.finish (bordered GS, H H^T = N I).
usage: python3 theory/p571_p19cpsat.py classes.pkl [--start i] [--count c] [--time 60] [--workers 1] [--log file]"""
import sys, os, json, time, pickle, math
import numpy as np
from ortools.sat.python import cp_model

n = 30; V = 571; P19 = 19
R = (2, 0, 0, 0)

def psi_allowed(Ks):
    """{sum_k shat_k^2 : shat_k = 30 - K_k (mod 2), |shat_k| <= 30 - K_k}"""
    sets = [{q * q for q in range(-(n - K), n - K + 1) if (q - (n - K)) % 2 == 0} for K in Ks]
    tot = {0}
    for s in sets: tot = {a + b for a in tot for b in s}
    return sorted(tot)

def stage1_model(W, eps, s, lag_targets=None, row_targets=None, lemma6=True, lags=range(1, 16), fix=None, comp=(),
                 split=False, dp_domains=False, upsd=False):
    """W, eps: (4, n) int arrays; s: (16,) gadget PAF sum. Returns (model, b vars, A-expression helper).
    lag_targets / row_targets override -76 - 324 s(h) / r_k - 18 t_k (planted tests); lemma6=False drops the psi test"""
    m = cp_model.CpModel()
    b = [[None] * n for _ in range(4)]
    for k in range(4):
        for g in range(n):
            if W[k][g] == 0: b[k][g] = m.NewBoolVar(f'b{k}_{g}')
    def A(k, g):                                            # linear expression or int
        g %= n
        return eps[k][g] * W[k][g] if W[k][g] != 0 else 2 * b[k][g] - 1
    t = [int(sum(W[k])) for k in range(4)]
    for k in range(4):
        m.Add(sum(A(k, g) for g in range(n)) == (R[k] - 18 * t[k] if row_targets is None else int(row_targets[k])))
    xor_cache = {}
    def prod(k, g1, g2):
        g1 %= n; g2 %= n
        f1, f2 = W[k][g1] != 0, W[k][g2] != 0
        if f1 and f2: return int(eps[k][g1] * W[k][g1] * eps[k][g2] * W[k][g2])
        if f1: return int(eps[k][g1] * W[k][g1]) * (2 * b[k][g2] - 1)
        if f2: return int(eps[k][g2] * W[k][g2]) * (2 * b[k][g1] - 1)
        key = (k, min(g1, g2), max(g1, g2))
        if key not in xor_cache:
            z = m.NewBoolVar(f'z{k}_{key[1]}_{key[2]}')
            m.AddBoolXOr([b[k][g1], b[k][g2], z.Not()])
            xor_cache[key] = z
        return 1 - 2 * xor_cache[key]
    if fix is not None:                                     # planted tests: fix some A-values
        for (k, g), val in fix.items(): m.Add(b[k][g] == (1 if val > 0 else 0))
    for h in lags:
        terms = []
        if split or dp_domains:                              # per-block PAF and X parts as integer variables
            Pk, Xk = [], []
            for k in range(4):
                pv = m.NewIntVar(-30, 30, ''); m.Add(pv == sum(prod(k, g, g + h) for g in range(n)))
                xs_ = [int(W[k][g]) * A(k, g + h) for g in np.nonzero(W[k])[0]] + [int(W[k][g]) * A(k, g - h) for g in np.nonzero(W[k])[0]]
                xv = m.NewIntVar(-14, 14, ''); m.Add(xv == sum(xs_) if xs_ else xv == 0)
                Pk.append(pv); Xk.append(xv)
                if dp_domains:
                    import p571_p19lagdp as LD
                    dom = sorted(LD.block_sigma(np.asarray(W), np.asarray(eps), k, h))
                    m.AddLinearExpressionInDomain(pv + 18 * xv, cp_model.Domain.FromValues(dom))
            tgt = (-76 - 324 * int(s[h]) if lag_targets is None else int(lag_targets[h]))
            xt = m.NewIntVar(-14, 14, ''); m.Add(xt == sum(Xk))
            if lag_targets is None:
                if s[h] == -1: m.Add(xt >= 8)
                else: m.Add(xt <= 2); m.Add(xt >= -10)
            m.Add(sum(Pk) + 18 * xt == tgt)
            continue
        for k in range(4):
            for g in range(n): terms.append(prod(k, g, g + h))
            for g in np.nonzero(W[k])[0]:
                terms.append(18 * int(W[k][g]) * A(k, g + h)); terms.append(18 * int(W[k][g]) * A(k, g - h))
        m.Add(sum(terms) == (-76 - 324 * int(s[h]) if lag_targets is None else int(lag_targets[h])))
    if comp:
        add_compressions(m, lambda k, g: A(k, g) + 18 * int(W[k][g % n]), comp)
    if upsd:
        add_upsd(m, W, A, s)
    if not lemma6: return m, b, A
    # Lemma 6 at psi
    sg = [1 - 2 * (g % 2) for g in range(n)]
    What = [int(sum(sg[g] * W[k][g] for g in range(n))) for k in range(4)]
    expr = sum(What[k] * (2 * sum(sg[g] * A(k, g) for g in range(n)) + 17 * What[k]) for k in range(4))
    Ks = [int(np.count_nonzero(W[k])) for k in range(4)]
    m.AddLinearExpressionInDomain(expr, cp_model.Domain.FromValues(psi_allowed(Ks)))
    return m, b, A

def solve_stage1(W, eps, s, tlim=60.0, workers=1, all_solutions=False, max_sols=50, **kw):
    m, b, A = stage1_model(W, eps, s, **kw)
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tlim
    solver.parameters.num_search_workers = workers
    sols = []
    if all_solutions:
        class CB(cp_model.CpSolverSolutionCallback):
            def __init__(s_): super().__init__()
            def on_solution_callback(s_):
                sols.append(np.array([[s_.Value(A(k, g)) if not isinstance(A(k, g), int) else A(k, g) for g in range(n)] for k in range(4)]))
                if len(sols) >= max_sols: s_.StopSearch()
        solver.parameters.enumerate_all_solutions = True; solver.parameters.num_search_workers = 1
        st = solver.Solve(m, CB())
    else:
        st = solver.Solve(m)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            sols.append(np.array([[solver.Value(A(k, g)) if not isinstance(A(k, g), int) else A(k, g) for g in range(n)] for k in range(4)]))
    return solver.StatusName(st), sols, solver.WallTime()

def lag_sums(Aarr, W):
    """sum_k [PAF(A_k)(h) + 18 X_k(h)], h = 0..15 (= sum_k PAF(y_k)(h) - 324 s(h))"""
    W = np.asarray(W); out = np.zeros(16, np.int64)
    for h in range(16):
        for k in range(4):
            out[h] += int((Aarr[k] * np.roll(Aarr[k], -h)).sum())
            out[h] += 18 * int((W[k] * (np.roll(Aarr[k], -h) + np.roll(Aarr[k], h))).sum())
    return out

def planted_partial_test(trials=5, seed=2, nfree=24):
    """planted A with all but nfree twisted A-values fixed: CP-SAT must return FEASIBLE, and every returned solution
    must reproduce the planted lag sums and row sums when recomputed with numpy"""
    rng = np.random.default_rng(seed); res = []
    for t_ in range(trials):
        W = np.zeros((4, n), np.int64); eps = np.zeros((4, n), np.int64)
        for k in range(4):
            pos = rng.choice(n, size=int(rng.integers(1, 3)), replace=False)
            W[k, pos] = rng.choice((-1, 1), size=len(pos)); eps[k, pos] = rng.choice((-1, 1), size=len(pos))
        Aarr = rng.choice((-1, 1), size=(4, n)); Aarr[W != 0] = (eps * W)[W != 0]
        Lt = lag_sums(Aarr, W); rows = Aarr.sum(1)
        tw = [(k, g) for k in range(4) for g in range(n) if W[k, g] == 0]
        keep = set(map(tuple, np.array(tw)[rng.choice(len(tw), size=nfree, replace=False)].tolist()))
        fix = {(k, g): int(Aarr[k, g]) for (k, g) in tw if (k, g) not in keep}
        st, sols, wt = solve_stage1(W, eps, np.zeros(16, np.int64), tlim=60, workers=4, lag_targets=Lt, row_targets=rows,
                                    lemma6=False, fix=fix, all_solutions=True, max_sols=200)
        ok = bool(sols) and all((lag_sums(S, W)[1:] == Lt[1:]).all() and (S.sum(1) == rows).all() for S in sols)
        res.append(dict(trial=t_, status=st, solutions=len(sols), all_consistent=bool(ok), secs=round(wt, 2)))
        print(json.dumps(res[-1]), flush=True)
    return res

def planted_test(trials=5, seed=1):
    """random gadget sets and random A: the model with the planted lag/row targets must be FEASIBLE and return an A
    with the same lag sums (tests the XOR / X-term encoding)"""
    rng = np.random.default_rng(seed); res = []
    for t_ in range(trials):
        W = np.zeros((4, n), np.int64); eps = np.zeros((4, n), np.int64)
        for k in range(4):
            pos = rng.choice(n, size=int(rng.integers(0, 3)), replace=False)
            W[k, pos] = rng.choice((-1, 1), size=len(pos)); eps[k, pos] = rng.choice((-1, 1), size=len(pos))
        Aarr = rng.choice((-1, 1), size=(4, n)); Aarr[W != 0] = (eps * W)[W != 0]
        Lt = lag_sums(Aarr, W); rows = Aarr.sum(1)
        s = np.zeros(16, np.int64)
        st, sols, wt = solve_stage1(W, eps, s, tlim=60, workers=4, lag_targets=Lt, row_targets=rows, lemma6=False)
        ok = bool(sols) and (lag_sums(sols[0], W)[1:] == Lt[1:]).all() and (sols[0].sum(1) == rows).all()
        res.append(dict(trial=t_, status=st, consistent=bool(ok), secs=round(wt, 2)))
        print(json.dumps(res[-1]), flush=True)
    return res

def stage1_model_y(W, eps, lag_targets=None, row_targets=None, fix=None):
    """cross-check formulation (pure S1, no X/s algebra, no Lemma 6): y_k(g) = 2b - 1 on twisted cosets,
    y = (eps + 18) W on gadgets (19 W for J, 17 W for J'); sum_g y_k = r_k; sum_k PAF(y_k)(h) = -76 (h = 1..15)"""
    m = cp_model.CpModel()
    b = [[m.NewBoolVar('') if W[k][g] == 0 else None for g in range(n)] for k in range(4)]
    yv = lambda k, g: int((eps[k][g % n] + 18) * W[k][g % n]) if W[k][g % n] != 0 else 2 * b[k][g % n] - 1
    for k in range(4): m.Add(sum(yv(k, g) for g in range(n)) == (R[k] if row_targets is None else int(row_targets[k])))
    if fix is not None:
        for (k, g), val in fix.items(): m.Add(b[k][g] == (1 if val > 0 else 0))
    xc = {}
    def pr(k, g1, g2):
        g1 %= n; g2 %= n
        f1, f2 = W[k][g1] != 0, W[k][g2] != 0
        if f1 and f2: return yv(k, g1) * yv(k, g2)
        if f1: return yv(k, g1) * (2 * b[k][g2] - 1)
        if f2: return yv(k, g2) * (2 * b[k][g1] - 1)
        key = (k, min(g1, g2), max(g1, g2))
        if key not in xc:
            z = m.NewBoolVar(''); m.AddBoolXOr([b[k][g1], b[k][g2], z.Not()]); xc[key] = z
        return 1 - 2 * xc[key]
    for h in range(1, 16):
        m.Add(sum(pr(k, g, g + h) for k in range(4) for g in range(n)) == (-76 if lag_targets is None else int(lag_targets[h])))
    return m

def planted_partial_test_y(trials=5, seed=3, nfree=24):
    """model y: planted A (gadgets random), all but nfree twisted values fixed; must be feasible and consistent"""
    rng = np.random.default_rng(seed); res = []
    for t_ in range(trials):
        W = np.zeros((4, n), np.int64); eps = np.zeros((4, n), np.int64)
        for k in range(4):
            pos = rng.choice(n, size=int(rng.integers(1, 3)), replace=False)
            W[k, pos] = rng.choice((-1, 1), size=len(pos)); eps[k, pos] = rng.choice((-1, 1), size=len(pos))
        Aarr = rng.choice((-1, 1), size=(4, n)); Aarr[W != 0] = (eps * W)[W != 0]
        y = Aarr + 18 * W
        Lt = np.array([sum(int((y[k] * np.roll(y[k], -h)).sum()) for k in range(4)) for h in range(16)])
        tw = [(k, g) for k in range(4) for g in range(n) if W[k, g] == 0]
        keep = set(map(tuple, np.array(tw)[rng.choice(len(tw), size=nfree, replace=False)].tolist()))
        fix = {(k, g): int(Aarr[k, g]) for (k, g) in tw if (k, g) not in keep}
        m = stage1_model_y(W, eps, lag_targets=Lt, row_targets=y.sum(1), fix=fix)
        solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = 60
        st = solver.Solve(m)
        res.append(dict(trial=t_, status=solver.StatusName(st)))
        print(json.dumps(res[-1]), flush=True)
    return res

def solve_stage1_y(W, eps, tlim=60.0, workers=1):
    m = stage1_model_y(W, eps)
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tlim; solver.parameters.num_search_workers = workers
    st = solver.Solve(m)
    return solver.StatusName(st), solver.WallTime()

def add_compressions(m, yexpr, divisors=(2, 3, 5, 6, 10)):
    """redundant DK-style constraints: compressed sequences c_k(u) = sum_{g = u mod d} y_k(g) as integer variables and
    the exact compressed PAF equations sum_k sum_u c_k(u) c_k(u + h') = T_d(h') (T(0) = 2208, T(h != 0) = -76)"""
    for d in divisors:
        f = n // d
        T0 = 2208 - 76 * (f - 1); Th = -76 * f
        c = [[None] * d for _ in range(4)]
        for k in range(4):
            for u in range(d):
                c[k][u] = m.NewIntVar(-19 * f, 19 * f, '')
                m.Add(c[k][u] == sum(yexpr(k, u + d * j) for j in range(f)))
        for hp in range(0, d // 2 + 1):
            prods = []
            for k in range(4):
                for u in range(d):
                    pv = m.NewIntVar(-(19 * f) ** 2, (19 * f) ** 2, '')
                    m.AddMultiplicationEquality(pv, [c[k][u], c[k][(u + hp) % d]])
                    prods.append(pv)
            m.Add(sum(prods) == (T0 if hp == 0 else Th))

def add_upsd(m, W, A, s, scale=10 ** 6):
    """necessary S2 condition on stage 1: U(h) = sum_k PAF(u_k)(h) = 4 + x_h + 17 s(h) (h != 0), U(0) = #twisted cosets,
    so Uhat(lambda) = sum_k |uhat_k(lambda)|^2 >= 0 for every character lambda of Z30.  Real coefficients are scaled by
    10^6 and rounded; the right-hand side is relaxed by the worst-case rounding error (sum |U(h)| / 2), so the integer
    constraint is implied by the real one (no feasible point is cut)."""
    W = np.asarray(W); ntw = int((W == 0).sum())
    xh = {}
    for h in range(1, 16):
        terms = [int(W[k][g]) * A(k, g + h) for k in range(4) for g in np.nonzero(W[k])[0]] + \
                [int(W[k][g]) * A(k, g - h) for k in range(4) for g in np.nonzero(W[k])[0]]
        xv = m.NewIntVar(-14, 14, ''); m.Add(xv == sum(terms)); xh[h] = xv
    Umax = {h: 4 + 14 + 17 * abs(int(s[h])) for h in range(1, 16)}
    for lam in range(16):
        coef = {h: (2.0 if h < 15 else 1.0) * math.cos(2 * math.pi * lam * h / n) for h in range(1, 16)}
        c_int = {h: int(round(scale * coef[h])) for h in coef}
        const = scale * ntw + sum(c_int[h] * (4 + 17 * int(s[h])) for h in coef)
        slack = sum(abs(c_int[h] - scale * coef[h]) * Umax[h] for h in coef) + 1
        # exact:  ntw + sum_h coef_h U(h) >= 0,  U(h) = 4 + 17 s(h) + x_h
        m.Add(const + sum(c_int[h] * xh[h] for h in coef) >= -int(math.ceil(slack)))

def joint_model(W, eps, s):
    """stage 1 and stage 2 together: A and u (both +-1 on twisted cosets); S1 (15 lags), the real part of S2
    (sum_k PAF(u_k)(h) = 4 + x_h + 17 s(h)), and the imaginary part (sum_k sum_g u_k(g) [beta_k(g-h) - beta_k(g+h)] = 0,
    beta = A - W), with XOR literals for A A', u u' and u A products; plus Lemma 6 at psi and the row sums."""
    m, b, A = stage1_model(W, eps, s)
    W = np.asarray(W)
    u = [[m.NewBoolVar('') if W[k][g] == 0 else None for g in range(n)] for k in range(4)]
    U = lambda k, g: (2 * u[k][g % n] - 1) if W[k][g % n] == 0 else 0
    cache = {}
    def xor_lit(a, b_, key):
        if key not in cache:
            z = m.NewBoolVar(''); m.AddBoolXOr([a, b_, z.Not()]); cache[key] = z
        return cache[key]
    for h in range(1, 16):
        terms = []
        for k in range(4):
            for g in range(n):
                g1, g2 = g, (g + h) % n
                if W[k][g1] or W[k][g2]: continue
                key = ('uu', k, min(g1, g2), max(g1, g2))
                terms.append(1 - 2 * xor_lit(u[k][g1], u[k][g2], key))
        xs = [int(W[k][g]) * A(k, g + h) for k in range(4) for g in np.nonzero(W[k])[0]] + \
             [int(W[k][g]) * A(k, g - h) for k in range(4) for g in np.nonzero(W[k])[0]]
        m.Add(sum(terms) == 4 + sum(xs) + 17 * int(s[h]))
    def beta(k, g):
        g %= n
        return int(eps[k][g] * W[k][g] - W[k][g]) if W[k][g] else ('A', k, g)
    for h in range(1, 15):
        terms = []
        for k in range(4):
            for g in range(n):
                if W[k][g]: continue
                for sign, g2 in ((1, (g - h) % n), (-1, (g + h) % n)):
                    bt = beta(k, g2)
                    if isinstance(bt, tuple):          # u(g) * A(g2), both +-1
                        key = ('uA', k, g, g2)
                        terms.append(sign * (1 - 2 * xor_lit(u[k][g], b[k][g2], key)))
                    elif bt != 0:
                        terms.append(sign * bt * (2 * u[k][g] - 1))
        m.Add(sum(terms) == 0)
    return m

def solve_joint(W, eps, s, tlim=120.0, workers=4):
    m = joint_model(W, eps, s)
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tlim; solver.parameters.num_search_workers = workers
    st = solver.Solve(m)
    return solver.StatusName(st), solver.WallTime()

def check_s1(Aarr, W):
    y = Aarr + 18 * np.asarray(W)
    P = np.stack([[int((y[k] * np.roll(y[k], -h)).sum()) for h in range(n)] for k in range(4)]).sum(0)
    return bool(P[0] == 2208 and (P[1:] == -76).all() and sorted(np.abs(y.sum(1)).tolist()) == [0, 0, 0, 2])

def stage2_model(Aarr, W, eps, s):
    W = np.asarray(W); m = cp_model.CpModel()
    beta = Aarr - W                                          # A on twisted, A - W on gadgets
    u = [[m.NewBoolVar(f'u{k}_{g}') if W[k][g] == 0 else None for g in range(n)] for k in range(4)]
    U = lambda k, g: (2 * u[k][g % n] - 1) if W[k][g % n] == 0 else 0
    yy = Aarr + 18 * W
    x = {}
    for h in range(1, 16):
        x[h] = sum(int(W[k][g]) * (int(Aarr[k][(g + h) % n]) + int(Aarr[k][(g - h) % n])) for k in range(4) for g in range(n) if W[k][g])
    xc = {}
    def uu(k, g1, g2):
        g1 %= n; g2 %= n
        if W[k][g1] or W[k][g2]: return 0
        key = (k, min(g1, g2), max(g1, g2))
        if key not in xc:
            z = m.NewBoolVar(''); m.AddBoolXOr([u[k][g1], u[k][g2], z.Not()]); xc[key] = z
        return 1 - 2 * xc[key]
    for h in range(1, 16):
        m.Add(sum(uu(k, g, g + h) for k in range(4) for g in range(n)) == 4 + x[h] + 17 * int(s[h]))
    for h in range(1, 15):
        m.Add(sum(U(k, g) * int(beta[k][(g - h) % n] - beta[k][(g + h) % n]) for k in range(4) for g in range(n)) == 0)
    return m, U

def solve_stage2(Aarr, W, eps, s, tlim=300.0, workers=4):
    m, U = stage2_model(Aarr, W, eps, s)
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tlim; solver.parameters.num_search_workers = workers
    st = solver.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        Uarr = np.array([[solver.Value(U(k, g)) if not isinstance(U(k, g), int) else 0 for g in range(n)] for k in range(4)])
        return solver.StatusName(st), Uarr
    return solver.StatusName(st), None

def build_and_verify(Aarr, Uarr, W):
    """x_k(g, z) on Z30 x Z19: z = 0 -> A; z != 0 -> u chi(z) (twisted) or W (gadget); bordered GS check"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import paley_b1pipe as PB
    W = np.asarray(W)
    return PB.finish(Aarr, Uarr, W, V, P19, n)

if __name__ == "__main__" and sys.argv[1] == 'test':
    planted_test(int(sys.argv[2]) if len(sys.argv) > 2 else 5); sys.exit()
if __name__ == "__main__" and sys.argv[1] == 'testy':
    planted_partial_test_y(int(sys.argv[2]) if len(sys.argv) > 2 else 5); sys.exit()
if __name__ == "__main__" and sys.argv[1] == 'testpartial':
    planted_partial_test(int(sys.argv[2]) if len(sys.argv) > 2 else 5); sys.exit()

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--start": "0", "--count": "1000000", "--time": "60", "--workers": "1", "--log": "", "--model": "A",
                             "--only": ""}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    classes = pickle.load(open(pos[0], 'rb')); keys = sorted(classes, key=lambda c: (classes[c]['K'], str(classes[c]['dist']), str(c)))
    st0 = int(opt['--start']); cnt = int(opt['--count'])
    logf = open(opt['--log'], 'a') if opt['--log'] else None
    tot = {}
    todo = range(st0, min(len(keys), st0 + cnt))
    if opt['--only']: todo = [int(x) for x in open(opt['--only']).read().split()]
    for idx in todo:
        info = classes[keys[idx]]
        W = np.array(info['W']); eps = np.array(info['eps'])
        if opt['--model'] == 'y':
            status, wt = solve_stage1_y(W, eps, tlim=float(opt['--time']), workers=int(opt['--workers']))
            tot[status] = tot.get(status, 0) + 1
            line = json.dumps(dict(i=idx, st=status, t=round(wt, 2)))
            print(line, flush=True)
            if logf: logf.write(line + '\n'); logf.flush()
            continue
        s = np.zeros(16, np.int64)
        for h in range(1, 16):
            s[h] = sum(int((W[k] * np.roll(W[k], -h)).sum()) for k in range(4))
        status, sols, wt = solve_stage1(W, eps, s, tlim=float(opt['--time']), workers=int(opt['--workers']))
        rec = dict(i=idx, st=status, t=round(wt, 2))
        if status not in ('INFEASIBLE',): rec.update(K=info['K'], dist=info['dist'], spec=info['spec'])
        if sols:
            A1 = sols[0]; rec['s1_exact'] = check_s1(A1, W); rec['A'] = A1.tolist()
            st2, U = solve_stage2(A1, W, eps, s)
            rec['stage2'] = st2
            if U is not None:
                res = build_and_verify(A1, U, W); rec['core_ok'] = res['core_ok']; rec['hadamard'] = res['hadamard']
        tot[status] = tot.get(status, 0) + 1
        line = json.dumps(rec)
        print(line, flush=True)
        if logf: logf.write(line + '\n'); logf.flush()
    print(json.dumps(dict(summary=True, range=[st0, min(len(keys), st0 + cnt)], statuses=tot)), flush=True)
