r"""b >= 3 border (template) search by CP-SAT, with an imposed symmetry group [M].

A template for compressions z (4 x c, any integers) in a given GS-type layout is a pair of frames K, Kc (4b x 4c,
b = c - 1, +-1, rows orthogonal, every row balanced against the block's null character) and a +-1 corner C with
    C Kc + K Cb^T = 0        (Cb = the 4c x 4c array of z in the layout; report section 3.1)
(C C^T = 4b I then follows).  We write N = [[C, K], [Kc^T, Cb]] and solve for the +-1 entries of C, K, Kc:
  * balance   sum_x n_J(x) K[i, (J, x)] = 0, same for Kc (n_J = all-ones unless the layout is character-twisted);
  * frames    K K^T = 4c I, Kc Kc^T = 4c I;
  * linking   sum_l C[i, l] Kc[l, r] + sum_s K[i, s] Cb[r, s] = 0 for every border row i and core row r;
  * corner    C C^T = 4b I (implied; kept as a redundant constraint).
Products of two unknowns are XNOR auxiliaries.  A symmetry group is imposed by generators (pr, sr, pc, sc): signed
permutations of the rows and columns of N that fix Cb; entries in one orbit share one variable (with signs), and
duplicate constraints are dropped.  Every solution is checked exactly (C = -K Cb^T Kc^T/(4c), +-1, frames).
library: model(...) / solve(...); see theory/ab_border_sat_run.py for the drivers."""
import sys, os, time, json
import numpy as np
from ortools.sat.python import cp_model


class UF:
    """union-find over entries with parity (entry = parity * rep)"""
    def __init__(self, n):
        self.p = list(range(n)); self.par = [0] * n
    def find(self, a):
        path = []
        while self.p[a] != a:
            path.append(a); a = self.p[a]
        root = a; acc = 0
        for x in reversed(path):
            acc ^= self.par[x]; self.par[x] = acc; self.p[x] = root
        return root
    def get(self, a):
        r = self.find(a)
        return r, (self.par[a] if a != r else 0)
    def union(self, a, b, odd):
        """impose x_a = (-1)^odd x_b; returns False on contradiction"""
        ra, pa = self.get(a); rb, pb = self.get(b)
        if ra == rb: return (pa ^ pb) == odd
        self.p[ra] = rb; self.par[ra] = pa ^ pb ^ odd
        return True


def model(Cb, b, c, gens=(), nullK=None, nullKc=None, fix=None, redundant=True):
    """Cb: 4c x 4c integer core; gens: list of (pr, sr, pc, sc) on N's 4b + 4c rows/cols (N[pr[r], pc[s]] =
    sr[r] sc[s] N[r, s]); nullK/nullKc: 4 x c null characters (default all-ones); fix: dict {(r, s): +-1}"""
    nb, nc = 4 * b, 4 * c; n = nb + nc
    Cb = np.asarray(Cb, np.int64)
    # entries: (r, s) unknown iff r < nb or s < nb; index them
    idx = -np.ones((n, n), np.int64); ent = []
    for r in range(n):
        for s in range(n):
            if r < nb or s < nb: idx[r, s] = len(ent); ent.append((r, s))
    uf = UF(len(ent))
    for (pr, sr, pc, sc) in gens:
        for e, (r, s) in enumerate(ent):
            r2, s2 = pr[r], pc[s]; sg = sr[r] * sc[s]
            if r2 >= nb and s2 >= nb: raise ValueError("generator maps an unknown entry into the core")
            if not uf.union(idx[r2, s2], e, 0 if sg > 0 else 1): return None
        # the generator must fix Cb
        for r in range(nb, n):
            for s in range(nb, n):
                if Cb[pr[r] - nb, pc[s] - nb] != sr[r] * sc[s] * Cb[r - nb, s - nb]:
                    raise ValueError("generator does not fix Cb")
    m = cp_model.CpModel()
    reps = {}
    def lit(e):
        """(var, parity) for entry e: value = (-1)^parity * (2 var - 1)"""
        r, p = uf.get(e)
        if r not in reps: reps[r] = m.NewBoolVar(f"x{r}")
        return reps[r], p
    if fix:
        for (r, s), v in fix.items():
            x, p = lit(idx[r, s]); m.Add(x == (1 if (v > 0) != bool(p) else 0))
    xn = {}
    def prod(e1, e2):
        """returns (const, var, sign): product N[e1] N[e2] = const if var is None else sign * (2 var - 1)"""
        x1, p1 = lit(e1); x2, p2 = lit(e2)
        sg = -1 if (p1 ^ p2) else 1
        if x1.Index() == x2.Index(): return sg, None, 0
        key = (min(x1.Index(), x2.Index()), max(x1.Index(), x2.Index()))
        if key not in xn:
            y = m.NewBoolVar(f"p{key[0]}_{key[1]}")
            m.AddBoolXOr([x1, x2, y])          # y = XNOR(x1, x2)
            xn[key] = y
        return 0, xn[key], sg
    seen = set(); ncons = [0]
    def add_sum(terms, rhs):
        """terms: list of (coef, kind, a, b): kind 'l' = coef * N[a]; kind 'p' = coef * N[a] N[b]; sum == rhs"""
        lin = {}; const = 0
        for (cf, kind, a, b_) in terms:
            if cf == 0: continue
            if kind == 'l':
                x, p = lit(a); sg = -1 if p else 1
                lin[x.Index()] = lin.get(x.Index(), [x, 0])
                lin[x.Index()][1] += 2 * cf * sg; const -= cf * sg
            else:
                k0, y, sg = prod(a, b_)
                if y is None: const += cf * k0
                else:
                    lin[y.Index()] = lin.get(y.Index(), [y, 0])
                    lin[y.Index()][1] += 2 * cf * sg; const -= cf * sg
        items = tuple(sorted((k, v[1]) for k, v in lin.items() if v[1] != 0))
        key = (items, rhs - const)
        if key in seen: return
        seen.add(key)
        if not items:
            if rhs - const != 0: m.AddBoolOr([])  # infeasible
            return
        m.Add(sum(v[1] * v[0] for k, v in lin.items() if v[1] != 0) == rhs - const); ncons[0] += 1
    ones = np.ones((4, c), np.int64)
    nK = ones if nullK is None else np.asarray(nullK); nKc = ones if nullKc is None else np.asarray(nullKc)
    for i in range(nb):                                   # K balance (rows i, core cols)
        for J in range(4):
            add_sum([(int(nK[J, x]), 'l', idx[i, nb + J * c + x], None) for x in range(c)], 0)
    for l in range(nb):                                   # Kc balance (core rows, border col l)
        for J in range(4):
            add_sum([(int(nKc[J, x]), 'l', idx[nb + J * c + x, l], None) for x in range(c)], 0)
    for i in range(nb):                                   # K frame
        for i2 in range(i + 1, nb):
            add_sum([(1, 'p', idx[i, s], idx[i2, s]) for s in range(nb, n)], 0)
    for l in range(nb):                                   # Kc frame
        for l2 in range(l + 1, nb):
            add_sum([(1, 'p', idx[r, l], idx[r, l2]) for r in range(nb, n)], 0)
    for i in range(nb):                                   # linking: top row i . core row r = 0
        for r in range(nb, n):
            t = [(1, 'p', idx[i, l], idx[r, l]) for l in range(nb)]
            t += [(int(Cb[r - nb, s - nb]), 'l', idx[i, s], None) for s in range(nb, n)]
            add_sum(t, 0)
    if redundant:
        for i in range(nb):                               # corner rows orthogonal
            for i2 in range(i + 1, nb):
                add_sum([(1, 'p', idx[i, l], idx[i2, l]) for l in range(nb)], 0)
    info = dict(free_vars=len(reps), xnor=len(xn), constraints=ncons[0])
    return m, (idx, uf, reps, ent), info


def extract(solver, data, nb, n):
    idx, uf, reps, ent = data
    N = np.zeros((n, n), np.int64)
    for e, (r, s) in enumerate(ent):
        rr, p = uf.get(e)
        v = 1 if solver.Value(reps[rr]) else -1
        N[r, s] = -v if p else v
    return N


def check(N, Cb, b, c, nullK=None, nullKc=None):
    nb = 4 * b; C = N[:nb, :nb]; K = N[:nb, nb:]; Kc = N[nb:, :nb].T
    ok = bool((np.abs(C) == 1).all() and (np.abs(K) == 1).all() and (np.abs(Kc) == 1).all())
    ok &= bool((K @ K.T == 4 * c * np.eye(nb, dtype=np.int64)).all() and (Kc @ Kc.T == 4 * c * np.eye(nb, dtype=np.int64)).all())
    nK = np.ones((4, c), np.int64) if nullK is None else np.asarray(nullK)
    nKc = np.ones((4, c), np.int64) if nullKc is None else np.asarray(nullKc)
    ok &= bool(all((K[:, J * c:(J + 1) * c] @ nK[J] == 0).all() and (Kc[:, J * c:(J + 1) * c] @ nKc[J] == 0).all() for J in range(4)))
    num = -(K @ np.asarray(Cb).T @ Kc.T)
    ok &= bool((num == 4 * c * C).all())
    ok &= bool((C @ C.T == 4 * b * np.eye(nb, dtype=np.int64)).all())
    return ok


def solve(Cb, b, c, gens=(), nullK=None, nullKc=None, fix=None, time_limit=600.0, workers=4, seed=0, log=False,
          redundant=True):
    t0 = time.time()
    res = model(Cb, b, c, gens, nullK, nullKc, fix, redundant)
    if res is None: return dict(status='SYMMETRY_CONTRADICTION', secs=round(time.time() - t0, 1)), None
    m, data, info = res
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    s.parameters.random_seed = seed; s.parameters.log_search_progress = log
    st = s.Solve(m)
    out = dict(info, status=s.StatusName(st), secs=round(time.time() - t0, 1))
    N = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        n = 4 * b + 4 * c
        N = extract(s, data, 4 * b, n)
        N[4 * b:, 4 * b:] = np.asarray(Cb)
        out['verified'] = check(N, Cb, b, c, nullK, nullKc)
    return out, N
