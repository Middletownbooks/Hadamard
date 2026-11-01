r"""Border problem of the Z_q all-character family: GS ansatz with X fixed [M] (Kumar's reduction, search/border3.c in
https://github.com/vatskishan03/hadamard, Apache-2.0; this code is new).

Rh = -X Ch^T W / (2q) must be +-1 (theory/zq_border_sat.py).  Ansatz: X = GS array over Q of a perfect quadruple x
(sum PAF = 4q delta), Ch = GS array of an unknown perfect quadruple c.  For fixed x every entry of X Ch^T W is linear
in the bits of c, so each x is a small CP-SAT model (4q bits, PAF constraints on c).  The x are sampled uniformly
from the perfect quadruples by meet in the middle on the PAF vectors of pairs.  A solution is a border, checked in
full (Rh, X, Ch Hadamard, X Ch^T + Rh Mc^T = 0); 'not found' decides only the sampled x.
usage: python3 theory/zq_border_fixx.py logs/zq/patterns_qQ.json [--only i] [--tries 2000] [--tmax 3600] [--out F]
"""
import sys, os, json, time, argparse, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB, zq_patterns as ZP, zq_border_gs as ZG


class Perfect:
    """sampler of +-1 quadruples on Q with sum PAF = 4q delta"""

    def __init__(self, Qg, seed=1):
        q = Qg.n
        self.Qg = Qg
        seqs = np.array(list(itertools.product((1, -1), repeat=q)), np.int64)
        P = np.stack([(seqs * seqs[:, Qg.add[:, d]]).sum(1) for d in range(1, q)], 1)   # 2^q x (q-1)
        self.seqs, self.P = seqs, P
        base = 8 * q + 1; off = 4 * q
        w = base ** np.arange(q - 1, dtype=np.int64)
        self.enc = lambda V: ((V + off) * w).sum(-1)
        n = len(seqs)
        I = np.repeat(np.arange(n), n); J = np.tile(np.arange(n), n)
        keys = self.enc(P[I] + P[J])
        o = np.argsort(keys, kind='stable')
        self.keys, self.I, self.J = keys[o], I[o].astype(np.int32), J[o].astype(np.int32)
        self.rng = np.random.default_rng(seed)

    def sample(self):
        n = len(self.seqs)
        while True:
            a = int(self.rng.integers(len(self.keys)))
            k = self.enc(-(self.P[self.I[a]] + self.P[self.J[a]]))
            lo, hi = np.searchsorted(self.keys, k, 'left'), np.searchsorted(self.keys, k, 'right')
            if hi > lo:
                b = int(self.rng.integers(lo, hi))
                x = self.seqs[[self.I[a], self.J[a], self.I[b], self.J[b]]]
                tot = sum((x[i] * x[i][self.Qg.add[:, d]]).sum() for i in range(4) for d in [1])
                return x


def solve_c(Qg, W, x, tmax=10.0, X=None):
    from ortools.sat.python import cp_model
    q = Qg.n; m = 4 * q
    if X is None:
        X = BB.gs_array(list(x), Qg.add, Qg.neg)
    sg, vid = ZG.sym_gs(Qg)
    md = cp_model.CpModel()
    cb = [md.NewBoolVar(f'c{u}') for u in range(m)]
    for g in range(1, q):
        terms = []
        for i in range(4):
            for h in range(q):
                a, b = i * q + h, i * q + int(Qg.add[h, g])
                d = md.NewBoolVar(''); md.AddBoolXOr([cb[a], cb[b], d.Not()]); terms.append(d)
        md.Add(len(terms) - 2 * sum(terms) == 0)
    # (X Ch^T W)[a, j] = sum_k W[k, j] sum_b X[a, b] sg[k, b] c[vid[k, b]]
    for a in range(m):
        for j in range(m):
            coef = np.zeros(m, np.int64)
            for k in np.flatnonzero(W[:, j]):
                np.add.at(coef, vid[k], int(W[k, j]) * X[a] * sg[k])
            nz = np.flatnonzero(coef)
            s = md.NewBoolVar('')
            # sum coef_u (2c_u - 1) = -2q (1 - 2 s)
            md.Add(sum(2 * int(coef[u]) * cb[u] for u in nz) - int(coef[nz].sum()) == -2 * q + 4 * q * s)
    md.Add(cb[0] == 1)
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = tmax; sol.parameters.num_workers = 1
    st = sol.Solve(md)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return sol.StatusName(st), None
    c = np.array([2 * sol.Value(v) - 1 for v in cb]).reshape(4, q)
    Ch = BB.gs_array(list(c), Qg.add, Qg.neg)
    num = -(X @ Ch.T @ W)
    Rh = num // (2 * q)
    ok = not (num % (2 * q)).any() and (np.abs(Rh) == 1).all() and (Rh @ Rh.T == m * np.eye(m)).all() and \
        (X @ X.T == m * np.eye(m)).all() and (Ch @ Ch.T == m * np.eye(m)).all() and \
        (X @ Ch.T + Rh @ (2 * W).T == 0).all()
    return ('FOUND' if ok else 'BAD'), dict(x=x.tolist(), c=c.tolist(), X=X.tolist(), Rh=Rh.tolist(), Ch=Ch.tolist())


def run(Qg, z, sampler, tries=2000, tmax_total=3600.0, tmax_one=10.0, xpaley=False, seed=1):
    W = BB.gs_array([np.asarray(v, np.int64) for v in z], Qg.add, Qg.neg)
    t0 = time.time(); stat = {}
    if xpaley:                           # X = Paley matrix with a random signed column permutation (Ch still GS)
        import zq_border as ZB
        H0 = ZB.paley(4 * Qg.n); rng = np.random.default_rng(seed)
    for k in range(tries):
        if time.time() - t0 > tmax_total:
            break
        if xpaley:
            m = 4 * Qg.n
            Xm = H0[:, rng.permutation(m)] * rng.choice([-1, 1], m)[None, :]
            st, res = solve_c(Qg, W, None, tmax=tmax_one, X=Xm)
            if res is not None:
                res['x'] = None
        else:
            x = sampler.sample()
            st, res = solve_c(Qg, W, x, tmax=tmax_one)
        stat[st] = stat.get(st, 0) + 1
        if st == 'FOUND':
            return dict(found=True, tries=k + 1, stat=stat, secs=round(time.time() - t0, 1), **res)
    return dict(found=False, tries=sum(stat.values()), stat=stat, secs=round(time.time() - t0, 1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('pfile'); ap.add_argument('--only', default=None)
    ap.add_argument('--tries', type=int, default=2000); ap.add_argument('--tmax', type=float, default=3600.0)
    ap.add_argument('--tone', type=float, default=10.0); ap.add_argument('--out', default=None)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--xpaley', action='store_true')
    ap.add_argument('--stop-first', action='store_true', help='stop after the first bordered pattern')
    a = ap.parse_args()
    P = json.load(open(a.pfile)); Qg = ZP.QGroup(tuple(P['mods']))
    only = None if a.only is None else {int(i) for i in a.only.split(',')}
    done = set()
    if a.out and os.path.exists(a.out):
        done = {json.loads(l)['index'] for l in open(a.out)}
    S = Perfect(Qg, seed=a.seed)
    for i, p in enumerate(P['patterns']):
        if (only is not None and i not in only) or i in done:
            continue
        r = run(Qg, p['z'], S, tries=a.tries, tmax_total=a.tmax, tmax_one=a.tone, xpaley=a.xpaley, seed=a.seed)
        rec = dict(q=P['q'], mods=P['mods'], index=i, z=p['z'], ansatz='Ch GS, X ' + ('Paley, random signed column permutation' if a.xpaley else 'GS (sampled)'), **r)
        print(json.dumps({k: rec.get(k) for k in ('q', 'mods', 'index', 'found', 'tries', 'stat', 'secs')}),
              flush=True)
        if a.out:
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
        if a.stop_first and r['found']:
            break
