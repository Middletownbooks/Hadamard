r"""Exact stage 1 of a b >= 3 (or b = 0) staged slice by CP-SAT [M when INFEASIBLE / every FEASIBLE answer is
re-verified exactly].

Stage 1 (theory/ab_lib.py): four integer vectors y_k on G' (constant on the G'-orbits of the slice), y_k(b) in the
letter set of b's orbit, with sum_k PAF_{G'}(y_k) = T1 at every lag-class representative.
Model: one integer variable per (block, G'-orbit) with the letters as domain; products p_ij = y_i y_j per block
(AddMultiplicationEquality); one linear equation per lag class, PAF_y(c) = sum_{i,j} N_c(i, j) y_i y_j with
N_c(i, j) = #{x in orbit i : x + c in orbit j}.
Implied constraints added for propagation [P]: block sums s_k with sum_k s_k^2 = sum_c T1 (all lags); for b >= 3 and
Q = the first nQ axes, the Q-coset sums z_k(q) satisfy |zhat_k(psi)|^2 <= 4b for psi != 1 (characters of Q = Z_2^2:
zhat in {0, +-2}), and sum_k zhat_k(psi)^2 = 4b.
Symmetry breaking (solution sets are closed under them): each block's first coordinate is positive (y and -y have the
same PAF), and the block energies are non-decreasing (block order is free).
Modes: --decide (feasibility; INFEASIBLE is a proof that stage 1, hence the configuration, is empty);
       --enumerate N (up to N solutions, each verified exactly with the exact PAF, written to --out as JSON lines).
usage: python3 theory/ab_s1_cpsat.py CFG.json [--decide] [--enumerate N] [--time S] [--workers W] [--out FILE]"""
import sys, os, json, time
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD


CHAR = True


def build_model(cfg, sl=None):
    sl = W.build(cfg) if sl is None else sl
    Gp = sl.Gp; L1 = len(sl.breps); borb = np.asarray(sl.borb); w = np.asarray(sl.bsz, np.int64)
    Up = [Gp.act(u) for u in sl.Sp] if len(Gp.mods) else [np.zeros(1, np.int64)]
    reps, _ = AGs.lag_reps(Gp, Up)
    T1 = np.asarray(sl.T1, np.int64)
    # N_c(i, j)
    N = np.zeros((len(reps), L1, L1), np.int64)
    for ci, c in enumerate(reps):
        xc = Gp.add(np.arange(Gp.n), np.full(Gp.n, c))
        np.add.at(N[ci], (borb, borb[xc]), 1)
    m = cp_model.CpModel()
    Y = [[m.NewIntVarFromDomain(cp_model.Domain.FromValues(sorted(int(v) for v in sl.letters[j])), f"y{k}_{j}")
          for j in range(L1)] for k in range(4)]
    P = []
    for k in range(4):
        Pk = {}
        for i in range(L1):
            for j in range(i, L1):
                if (N[:, i, j] + N[:, j, i]).any() or i == j:
                    li = sl.letters[i]; lj = sl.letters[j]
                    prods = sorted({int(a) * int(b) for a in li for b in lj})
                    p = m.NewIntVarFromDomain(cp_model.Domain.FromValues(prods), f"p{k}_{i}_{j}")
                    m.AddMultiplicationEquality(p, [Y[k][i], Y[k][j]])
                    Pk[(i, j)] = p
        P.append(Pk)
    def paf_expr(k, ci):
        terms = []
        for (i, j), p in P[k].items():
            coef = int(N[ci, i, j]) if i == j else int(N[ci, i, j] + N[ci, j, i])
            if coef: terms.append(coef * p)
        return sum(terms)
    for ci in range(len(reps)):
        m.Add(sum(paf_expr(k, ci) for k in range(4)) == int(T1[reps[ci]]))
    # implied: block sums and sum of squares of block sums
    R4 = int(T1.sum())
    S = []; S2 = []
    for k in range(4):
        smax = int(sum(w[j] * max(abs(int(v)) for v in sl.letters[j]) for j in range(L1)))
        s = m.NewIntVar(-smax, smax, f"s{k}"); m.Add(s == sum(int(w[j]) * Y[k][j] for j in range(L1)))
        s2 = m.NewIntVar(0, smax * smax, f"s2_{k}"); m.AddMultiplicationEquality(s2, [s, s])
        S.append(s); S2.append(s2)
    m.Add(sum(S2) == R4)
    # implied Q-level (b >= 3, Q = Z_2^2 on the first two axes)
    nQ = cfg.get('nQ', 0); b = cfg.get('b', 0)
    if nQ == 2 and tuple(Gp.mods[:2]) == (2, 2) and b >= 3:
        qof = Gp.co[np.asarray(sl.breps)][:, :2]
        qidx = qof[:, 0] * 2 + qof[:, 1]
        chars = [(1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)]          # psi((a, b)) on q index a*2+b
        for ch in chars:
            Zh = []
            for k in range(4):
                expr = sum(int(ch[qidx[j]]) * int(w[j]) * Y[k][j] for j in range(L1))
                zk = m.NewIntVarFromDomain(cp_model.Domain.FromValues([-2, 0, 2]), "")
                m.Add(zk == expr)
                z2 = m.NewIntVarFromDomain(cp_model.Domain.FromValues([0, 4]), "")
                m.AddMultiplicationEquality(z2, [zk, zk])
                Zh.append(z2)
            m.Add(sum(Zh) == 4 * b)
    # redundant character formulation over Q = Z_2^2 [P]: Z_k[psi, o] = sum_q psi(q) y_k(q, o) on the H'-orbits o;
    # for every psi: sum_k PAF_{H'}(Z_k[psi]) = T^psi, T^psi(c) = sum_q' psi(q') T1(q', c)  (Fourier over Q)
    if CHAR and nQ == 2 and tuple(Gp.mods[:2]) == (2, 2):
        n = Gp.n; co = Gp.co
        hkey = lambda x: tuple(co[x, 2:].tolist())
        base = {}                                               # H'-element -> G'-orbit of (0, 0, h)
        for x in range(n):
            if co[x, 0] == 0 and co[x, 1] == 0: base[hkey(x)] = int(borb[x])
        olab = sorted(set(base.values())); oid = {g: i for i, g in enumerate(olab)}; no = len(olab)
        j_of = {}                                               # (q index, o) -> G'-orbit j
        hor = np.array([oid[base[hkey(x)]] for x in range(n)])  # H'-orbit id of every G' element
        qix = co[:, 0] * 2 + co[:, 1]
        for x in range(n): j_of[(int(qix[x]), int(hor[x]))] = int(borb[x])
        wo = np.array([int((hor[(qix == 0)] == o).sum()) for o in range(no)], np.int64)
        H0 = np.nonzero(qix == 0)[0]                            # H' as the coset q = 0
        crep = [c for c in reps if qix[c] == 0]                 # H'-lag classes (G' lag reps with Q-part 0)
        NH = np.zeros((len(crep), no, no), np.int64)
        for ci, c in enumerate(crep):
            xc = Gp.add(H0, np.full(len(H0), c))
            np.add.at(NH[ci], (hor[H0], hor[xc]), 1)
        allchars = [(1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)]
        for ch in allchars:
            # target T^psi on H'-lag c: sum over q' of psi(q') T1(q', c)
            Tpsi = []
            for c in crep:
                tot = 0
                for x in range(n):
                    if hor[x] == hor[c] and co[x, 2:].tolist() == co[c, 2:].tolist():
                        tot += ch[int(qix[x])] * int(T1[x])
                Tpsi.append(tot)
            Zs = []
            for k in range(4):
                zr = []
                for o in range(no):
                    lin = sum(int(ch[qq]) * Y[k][j_of[(qq, o)]] for qq in range(4))
                    zmax = sum(max(abs(int(v)) for v in sl.letters[j_of[(qq, o)]]) for qq in range(4))
                    z = m.NewIntVar(-zmax, zmax, ""); m.Add(z == lin); zr.append(z)
                Zs.append(zr)
            for k in range(4):
                pass
            prods = []
            for k in range(4):
                pk = {}
                for o in range(no):
                    for o2 in range(o, no):
                        if o == o2 or (NH[:, o, o2] + NH[:, o2, o]).any():
                            zm = sum(max(abs(int(v)) for v in sl.letters[j_of[(qq, o)]]) for qq in range(4))
                            zm2 = sum(max(abs(int(v)) for v in sl.letters[j_of[(qq, o2)]]) for qq in range(4))
                            pv = m.NewIntVar(-zm * zm2 if o != o2 else 0, zm * zm2, "")
                            m.AddMultiplicationEquality(pv, [Zs[k][o], Zs[k][o2]]); pk[(o, o2)] = pv
                prods.append(pk)
            for ci in range(len(crep)):
                expr = []
                for k in range(4):
                    for (o, o2), pv in prods[k].items():
                        cf = int(NH[ci, o, o2]) if o == o2 else int(NH[ci, o, o2] + NH[ci, o2, o])
                        if cf: expr.append(cf * pv)
                m.Add(sum(expr) == int(Tpsi[ci]))
    # symmetry breaking
    for k in range(4):
        pos = [int(v) for v in sl.letters[0] if v > 0]
        m.AddAllowedAssignments([Y[k][0]], [[v] for v in pos])
    E = [paf_expr(k, 0) for k in range(4)]                     # rep 0 is lag 0 (energy)
    assert reps[0] == 0
    for k in range(3): m.Add(E[k] <= E[k + 1])
    return m, Y, sl, reps, N


def verify(sl, Yv):
    Gp = sl.Gp
    full = np.array([[Yv[k][int(sl.borb[x])] for x in range(Gp.n)] for k in range(4)], np.int64)
    return bool((Gp.paf(full).sum(0) == sl.T1).all()), full


class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, Y, sl, limit, out):
        super().__init__(); self.Y, self.sl, self.limit, self.out = Y, sl, limit, out; self.n = 0; self.t0 = time.time()
    def on_solution_callback(self):
        Yv = [[self.Value(v) for v in row] for row in self.Y]
        ok, full = verify(self.sl, Yv)
        self.n += 1
        if self.out:
            with open(self.out, 'a') as f: f.write(json.dumps(dict(n=self.n, verified=ok, y=Yv, secs=round(time.time() - self.t0, 1))) + "\n")
        if self.n >= self.limit: self.StopSearch()


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read())
    tlim = float(a[a.index('--time') + 1]) if '--time' in a else 600.0
    workers = int(a[a.index('--workers') + 1]) if '--workers' in a else 4
    out = a[a.index('--out') + 1] if '--out' in a else None
    t0 = time.time()
    global CHAR
    CHAR = '--nochar' not in a
    m, Y, sl, reps, N = build_model(cfg)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tlim; s.parameters.num_workers = workers
    info = dict(cfg=dict(v=cfg['v'], q=cfg['q'], s=cfg['s'], U=[list(u) for u in cfg['U']]), L1=len(sl.breps),
                lag_classes=len(reps), letters=sorted({int(x) for l in sl.letters for x in l}))
    if '--enumerate' in a:
        lim = int(a[a.index('--enumerate') + 1])
        s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
        cb = Collector(Y, sl, lim, out)
        st = s.Solve(m, cb)
        info.update(mode='enumerate', status=s.StatusName(st), solutions=cb.n, complete=(st == cp_model.OPTIMAL),
                    secs=round(time.time() - t0, 1))
    else:
        st = s.Solve(m)
        info.update(mode='decide', status=s.StatusName(st), secs=round(time.time() - t0, 1))
        if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            Yv = [[s.Value(v) for v in row] for row in Y]
            ok, full = verify(sl, Yv); info.update(verified=ok, y=Yv)
            if out:
                with open(out, 'a') as f: f.write(json.dumps(dict(n=1, verified=ok, y=Yv)) + "\n")
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    main()
