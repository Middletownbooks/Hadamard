r"""Border problem of the Z_q all-character family as one CP-SAT model [M] (Kumar's border3.c reduction,
search/border3.c in https://github.com/vatskishan03/hadamard, Apache-2.0; this code is new).

With W = Mc / 2 (W W^T = q I), the border condition X Ch^T + Rh Mc^T = 0 is  Rh = -X Ch^T W / (2q).
X may be fixed to one matrix of each Hadamard equivalence class (its left signed permutations go into Rh, its right
ones into Ch).  Unknowns: Ch in {+-1}^{4q x 4q} with mutually orthogonal rows; constraints: every entry of
X Ch^T W equals +-2q (linear in Ch).  Then Rh is +-1 and Hadamard automatically.  X: the Paley matrix of order 4q
(theory/zq_border.paley); a model infeasible for it decides only that class.
usage: python3 theory/zq_border_sat.py logs/zq/patterns_qQ.json [--only i] [--tmax 600] [--out FILE]
"""
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB, zq_patterns as ZP, zq_border as ZB


def solve(Qg, z, X=None, tmax=600.0, workers=4, seed=1):
    from ortools.sat.python import cp_model
    q = Qg.n; m = 4 * q
    W = BB.gs_array([np.asarray(v, np.int64) for v in z], Qg.add, Qg.neg)
    assert (W @ W.T == q * np.eye(m)).all()
    X = ZB.paley(m) if X is None else X
    md = cp_model.CpModel()
    c = [[md.NewBoolVar(f'c{k}_{b}') for b in range(m)] for k in range(m)]      # Ch[k,b] = 2c - 1
    s = [[md.NewBoolVar(f's{a}_{j}') for j in range(m)] for a in range(m)]      # Rh[a,j] = 1 - 2 s  (sign)
    # sum_{k,b} X[a,b] W[k,j] Ch[k,b] = -2q Rh[a,j]
    for a in range(m):
        for j in range(m):
            ks = np.flatnonzero(W[:, j])
            terms = []; const = 0
            for k in ks:
                for b in range(m):
                    co = int(X[a, b] * W[k, j])
                    terms.append(2 * co * c[k][b]); const -= co
            md.Add(sum(terms) + const == -2 * q * (1 - 2 * s[a][j]))
    # orthogonal rows of Ch: exactly 2q disagreements
    for k in range(m):
        for l in range(k + 1, m):
            d = []
            for b in range(m):
                x = md.NewBoolVar('')
                md.AddBoolXOr([c[k][b], c[l][b], x.Not()])      # x = c_k xor c_l
                d.append(x)
            md.Add(sum(d) == 2 * q)
    md.Add(c[0][0] == 1)                                        # global sign of Ch (and Rh)
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = tmax
    sol.parameters.num_workers = workers; sol.parameters.random_seed = seed
    t0 = time.time(); st = sol.Solve(md)
    res = dict(status=sol.StatusName(st), secs=round(time.time() - t0, 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        Ch = np.array([[2 * sol.Value(v) - 1 for v in row] for row in c], np.int64)
        num = -(X @ Ch.T @ W)
        assert not (num % (2 * q)).any()
        Rh = num // (2 * q)
        Mc = 2 * W
        ok = (np.abs(Rh) == 1).all() and (Rh @ Rh.T == m * np.eye(m)).all() and \
            (Ch @ Ch.T == m * np.eye(m)).all() and (X @ Ch.T + Rh @ Mc.T == 0).all()
        res.update(ok=bool(ok), X=X.tolist(), Rh=Rh.tolist(), Ch=Ch.tolist())
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('pfile'); ap.add_argument('--only', type=int, default=None)
    ap.add_argument('--tmax', type=float, default=600.0); ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    P = json.load(open(a.pfile)); Qg = ZP.QGroup(tuple(P['mods']))
    done = set()
    if a.out and os.path.exists(a.out):
        done = {json.loads(l)['index'] for l in open(a.out)}
    for i, p in enumerate(P['patterns']):
        if (a.only is not None and i != a.only) or i in done:
            continue
        r = solve(Qg, p['z'], tmax=a.tmax, workers=a.workers)
        rec = dict(q=P['q'], mods=P['mods'], index=i, z=p['z'], **r)
        print(json.dumps({k: rec.get(k) for k in ('q', 'mods', 'index', 'status', 'ok', 'secs')}), flush=True)
        if a.out:
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
