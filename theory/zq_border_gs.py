r"""Border problem of the Z_q all-character family with a Goethals-Seidel ansatz [M] (Kumar's reduction,
search/border3.c in https://github.com/vatskishan03/hadamard, Apache-2.0; this code is new).

Border condition (theory/zq_border_sat.py): Rh = -X Ch^T W / (2q) must be +-1, X and Ch Hadamard of order 4q,
W = Mc / 2 the GS array over Q of the pattern z.  Ansatz: X and Ch are themselves GS arrays over Q of +-1 sequences
x_1..x_4, c_1..c_4 with sum PAF = 4q delta.  Every entry of X Ch^T W is linear in the products x_i(g) c_l(h), so the
model is: bits for x and c, one XOR bit per product that occurs, PAF constraints on x and c, and 16 q^2 linear
equalities.  A solution is a border (then Rh is Hadamard automatically); infeasibility decides only the ansatz.
usage: python3 theory/zq_border_gs.py logs/zq/patterns_qQ.json [--only i] [--tmax 600] [--out FILE]
"""
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB, zq_patterns as ZP


def sym_gs(Qg):
    """GS array over Q of symbolic sequences: (sign, variable index i*q + g)"""
    q = Qg.n
    ids = [np.arange(1, q + 1) + i * q for i in range(4)]
    S = BB.gs_array(ids, Qg.add, Qg.neg)
    return np.sign(S), np.abs(S) - 1


def solve(Qg, z, tmax=600.0, workers=4, seed=1, sym=False, same=False, mid=False):
    from ortools.sat.python import cp_model
    q = Qg.n; m = 4 * q
    W = BB.gs_array([np.asarray(v, np.int64) for v in z], Qg.add, Qg.neg)
    assert (W @ W.T == q * np.eye(m)).all()
    sg, vid = sym_gs(Qg)
    md = cp_model.CpModel()
    xb = [md.NewBoolVar(f'x{u}') for u in range(m)]
    cb = [md.NewBoolVar(f'c{u}') for u in range(m)]
    xor = {}

    def xv(A, B, ka, kb, key):
        if key not in xor:
            d = md.NewBoolVar('')
            md.AddBoolXOr([A[ka], B[kb], d.Not()])
            xor[key] = d
        return xor[key]
    # PAF of x and of c
    for nm, V in (('x', xb), ('c', cb)):
        for g in range(1, q):
            terms = []
            for i in range(4):
                for h in range(q):
                    a, b = i * q + h, i * q + Qg.add[h, g]
                    terms.append(xv(V, V, a, b, (nm, min(a, b), max(a, b), a == b)))
            # sum (1 - 2 d) = 0
            md.Add(len(terms) - 2 * sum(terms) == 0)
    # entries: (X Ch^T W)[a, j] = sum_k W[k, j] sum_b X[a, b] Ch[k, b]            (default ansatz: X, Ch GS)
    # mid: x = the sequences of Rh, c those of Ch, entries of -Rh W^T Ch = 2q X  (ansatz: Rh, Ch GS)
    s = [[md.NewBoolVar('') for j in range(m)] for a in range(m)]
    WT = W.T
    for a in range(m):
        for j in range(m):
            coef = {}
            if not mid:
                for k in np.flatnonzero(W[:, j]):
                    w = int(W[k, j])
                    for b in range(m):
                        key = (int(vid[a, b]), int(vid[k, b]))
                        coef[key] = coef.get(key, 0) + w * int(sg[a, b]) * int(sg[k, b])
            else:
                for k in range(m):
                    for l in np.flatnonzero(WT[k]):
                        key = (int(vid[a, k]), int(vid[l, j]))
                        coef[key] = coef.get(key, 0) + int(WT[k, l]) * int(sg[a, k]) * int(sg[l, j])
            terms = []; const = 0
            for (u, v), co in coef.items():
                if co:
                    d = xv(xb, cb, u, v, ('p', u, v))
                    const += co; terms.append(-2 * co * d)
            # value = const + sum terms  = -2q (1 - 2 s)
            md.Add(const + sum(terms) == -2 * q + 4 * q * s[a][j])
    md.Add(xb[0] == 1); md.Add(cb[0] == 1)
    if sym:                                                     # Williamson-type: symmetric sequences
        for V in (xb, cb):
            for i in range(4):
                for g in range(q):
                    md.Add(V[i * q + g] == V[i * q + int(Qg.neg[g])])
    if same:                                                    # X = Ch
        for u in range(m):
            md.Add(xb[u] == cb[u])
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = tmax
    sol.parameters.num_workers = workers; sol.parameters.random_seed = seed
    t0 = time.time(); st = sol.Solve(md)
    res = dict(status=sol.StatusName(st), secs=round(time.time() - t0, 1), nxor=len(xor))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        x = np.array([2 * sol.Value(v) - 1 for v in xb]).reshape(4, q)
        c = np.array([2 * sol.Value(v) - 1 for v in cb]).reshape(4, q)
        A = BB.gs_array(list(x), Qg.add, Qg.neg); Ch = BB.gs_array(list(c), Qg.add, Qg.neg)
        if not mid:
            X = A; num = -(X @ Ch.T @ W); Rh = num // (2 * q)
        else:
            Rh = A; num = -(Rh @ W.T @ Ch); X = num // (2 * q)
        ok = not (num % (2 * q)).any() and (np.abs(Rh) == 1).all() and (np.abs(X) == 1).all() and \
            (Rh @ Rh.T == m * np.eye(m)).all() and \
            (X @ X.T == m * np.eye(m)).all() and (Ch @ Ch.T == m * np.eye(m)).all() and \
            (X @ Ch.T + Rh @ (2 * W).T == 0).all()
        res.update(ok=bool(ok), x=x.tolist(), c=c.tolist(), X=X.tolist(), Rh=Rh.tolist(), Ch=Ch.tolist())
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('pfile'); ap.add_argument('--only', type=int, default=None)
    ap.add_argument('--tmax', type=float, default=600.0); ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--sym', action='store_true'); ap.add_argument('--same', action='store_true'); ap.add_argument('--mid', action='store_true')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    P = json.load(open(a.pfile)); Qg = ZP.QGroup(tuple(P['mods']))
    done = set()
    if a.out and os.path.exists(a.out):
        done = {json.loads(l)['index'] for l in open(a.out)}
    for i, p in enumerate(P['patterns']):
        if (a.only is not None and i != a.only) or i in done:
            continue
        r = solve(Qg, p['z'], tmax=a.tmax, workers=a.workers, sym=a.sym, same=a.same, seed=a.seed, mid=a.mid)
        r['ansatz'] = 'GS' + ('+sym' if a.sym else '') + ('+same' if a.same else '') + ('+mid' if a.mid else '')
        rec = dict(q=P['q'], mods=P['mods'], index=i, z=p['z'], **r)
        print(json.dumps({k: rec.get(k) for k in ('q', 'mods', 'index', 'ansatz', 'status', 'ok', 'secs', 'nxor')}), flush=True)
        if a.out:
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
