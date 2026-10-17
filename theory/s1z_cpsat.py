r"""Stage-1 solutions of a b >= 3 GP slice with a prescribed Q-compression quadruple, by CP-SAT [M for every solution
returned: each is re-verified with the exact integer PAF; the search is a sample, not exhaustive].

Model (as theory/ab_s1_cpsat.py): one integer variable per (block, G'-orbit) with the letters as domain; products
y_i y_j per block; one linear equation per lag class of G'.  Added: the coset sums sum_j prof[j, q] y_kj = z_k(q) of
the prescribed compression quadruple.  Diversity: every solve uses a fresh random seed, a random linear objective over
the letters and random value hints; solutions already found are excluded by no-good cuts on their block patterns
(optionally only the first DIFF blocks, so the next solution differs in a whole block).
usage: python3 theory/s1z_cpsat.py CFG.json '<z json 4 x c>' OUT.jsonl [--want 20] [--secs 120] [--workers 4]
                                   [--seed 1] [--total 3600]"""
import sys, os, json, time
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD


def build(cfg, Z, sl=None, found=(), rng=None, objective=True):
    sl = W.build(cfg) if sl is None else sl
    Gp = sl.Gp; L1 = len(sl.breps); borb = np.asarray(sl.borb); w = np.asarray(sl.bsz, np.int64)
    nQ = cfg['nQ']
    from ab_lib import AG
    Q = AG(Gp.mods[:nQ]); qfull = Q.idx(Gp.co[:, :nQ])
    prof = np.zeros((L1, Q.n), np.int64); np.add.at(prof, (borb, qfull), 1)
    Up = [Gp.act(u) for u in sl.Sp]
    reps, _ = AGs.lag_reps(Gp, Up)
    T1 = np.asarray(sl.T1, np.int64)
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
                    prods = sorted({int(a) * int(b) for a in sl.letters[i] for b in sl.letters[j]})
                    p = m.NewIntVarFromDomain(cp_model.Domain.FromValues(prods), f"p{k}_{i}_{j}")
                    m.AddMultiplicationEquality(p, [Y[k][i], Y[k][j]])
                    Pk[(i, j)] = p
        P.append(Pk)
    for ci in range(len(reps)):
        terms = []
        for k in range(4):
            for (i, j), p in P[k].items():
                coef = int(N[ci, i, j]) if i == j else int(N[ci, i, j] + N[ci, j, i])
                if coef: terms.append(coef * p)
        m.Add(sum(terms) == int(T1[reps[ci]]))
    for k in range(4):
        for q in range(Q.n):
            m.Add(sum(int(prof[j, q]) * Y[k][j] for j in range(L1) if prof[j, q]) == int(Z[k][q]))
    # no-goods: exclude found block patterns (block k of a found solution, as a whole vector)
    for (k, yv) in found:
        bs = []
        for j in range(L1):
            b = m.NewBoolVar("")
            m.Add(Y[k][j] != int(yv[j])).OnlyEnforceIf(b); m.Add(Y[k][j] == int(yv[j])).OnlyEnforceIf(b.Not())
            bs.append(b)
        m.AddBoolOr(bs)
    if rng is not None:
        if objective:
            m.Maximize(sum(int(rng.integers(-3, 4)) * Y[k][j] for k in range(4) for j in range(L1)))
        for k in range(4):
            for j in range(L1):
                L = sorted(int(v) for v in sl.letters[j])
                m.AddHint(Y[k][j], int(L[rng.integers(len(L))]))
    return m, Y, sl


def verify(sl, Yv):
    full = np.asarray(Yv, np.int64)[:, np.asarray(sl.borb)]
    return bool((sl.Gp.paf(full).sum(0) == sl.T1).all())


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read() if a[0].endswith('.json') else a[0])
    Z = json.loads(a[1]); outp = a[2]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    want = int(opt('--want', '20')); secs = float(opt('--secs', '120')); workers = int(opt('--workers', '4'))
    seed = int(opt('--seed', '1')); total = float(opt('--total', '3600')); diff = int(opt('--diff', '1'))
    rng = np.random.default_rng(seed)
    sl = W.build(cfg); found = []; nsol = 0; t0 = time.time(); tries = 0
    seen = set()
    if os.path.exists(outp):
        for l in open(outp):
            r = json.loads(l)
            if r.get('y') and r.get('z') == Z:
                seen.add(json.dumps(r['y'])); nsol += 1
                for k in range(diff): found.append((k, r['y'][k]))
    while nsol < want and time.time() - t0 < total:
        tries += 1
        m, Y, _ = build(cfg, Z, sl=sl, found=found, rng=rng, objective=True)
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = secs; s.parameters.num_workers = workers
        s.parameters.random_seed = int(rng.integers(1 << 30)); s.parameters.stop_after_first_solution = True
        t1 = time.time(); st = s.Solve(m)
        rec = dict(stage='solve', try_=tries, status=s.StatusName(st), secs=round(time.time() - t1, 1))
        if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            Yv = [[s.Value(v) for v in row] for row in Y]
            ok = verify(sl, Yv); key = json.dumps(Yv)
            rec['verified'] = ok
            if ok and key not in seen:
                seen.add(key); nsol += 1
                with open(outp, 'a') as f: f.write(json.dumps(dict(y=Yv, z=Z, verified=True, secs=rec['secs'])) + "\n")
                for k in range(diff): found.append((k, Yv[k]))
        elif st == cp_model.INFEASIBLE:
            rec['note'] = 'infeasible with the current no-goods'
            print(json.dumps(rec), flush=True); break
        print(json.dumps(rec), flush=True)
    print(json.dumps(dict(stage='end', solutions=nsol, tries=tries, secs=round(time.time() - t0, 1))), flush=True)


if __name__ == '__main__':
    main()
