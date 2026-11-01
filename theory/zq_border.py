r"""Border problem for the Z_q all-character bordered family [M] (Kumar, search/NOTES.md section 18, border3.c /
border_probe.py in https://github.com/vatskishan03/hadamard, Apache-2.0; theory/bt_border.py; this code is new).

T = all characters of Q, t = q, so V_T = R^{4q} and a frame Rh is any Hadamard matrix of order 4q.  For a pattern z
(compressions y = 2z), Mc = GS array over Q of y; a border exists for the frame Rh iff
  S_R = {eps in {+-1}^{4q} : Mc Rh^T eps in {+-4q}^{4q}}
contains 4q mutually orthogonal vectors (the columns of X); then Ch = -Mc Rh^T X / (4q)  (theory/bt_border.py, [P]).
The left signed permutations of Rh are absorbed by X, so frames = Paley matrices with random signed column
permutations.  S_R: brute force (numba) for 4q <= 28, CP-SAT enumeration above.
usage: python3 theory/zq_border.py logs/zq/patterns_qQ.json [--frames 200] [--out FILE]
"""
import sys, os, json, time, argparse, itertools, math
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_border as BB, paley_multi as PM, zq_patterns as ZP


def paley(order):
    """Hadamard matrix of the given order: Paley I (order - 1 a prime power = 3 mod 4) or Paley II
    (order / 2 - 1 a prime power = 1 mod 4)"""
    def jac(Q):
        mods, cl = PM.gf(Q)
        co = np.array(list(itertools.product(*[range(m) for m in mods])))
        idx = lambda c: int(np.ravel_multi_index(tuple(np.mod(c, mods)), mods))
        chi = np.where(cl == 0, 0, np.where(cl == 1, 1, -1))
        return np.array([[chi[idx(co[j] - co[i])] for j in range(Q)] for i in range(Q)], np.int64)
    def ppow(Q):
        for r in range(2, Q + 1):
            if Q % r == 0:
                k = round(math.log(Q, r)); return r ** k == Q
        return False
    Q = order - 1
    if ppow(Q) and Q % 4 == 3:
        J = jac(Q)
        S = np.zeros((order, order), np.int64); S[0, 1:] = 1; S[1:, 0] = -1; S[1:, 1:] = J
        H = S + np.eye(order, dtype=np.int64)
    else:
        Q = order // 2 - 1
        assert ppow(Q) and Q % 4 == 1
        J = jac(Q)
        S = np.zeros((Q + 1, Q + 1), np.int64); S[0, 1:] = 1; S[1:, 0] = 1; S[1:, 1:] = J
        A = np.array([[1, -1], [-1, -1]]); B = np.array([[1, 1], [1, -1]])
        H = np.kron(S, A) + np.kron(np.eye(Q + 1, dtype=np.int64), B)
    assert (H @ H.T == order * np.eye(order)).all(), order
    return H


@nb.njit(cache=True)
def _brute(A, target, cap):
    n = A.shape[1]; m = A.shape[0]
    out = np.zeros((cap, n), np.int8); k = 0
    v = np.zeros(m, np.int64)
    for j in range(m):
        s = 0
        for c in range(n):
            s += A[j, c]
        v[j] = s                                    # eps = all +1
    eps = np.ones(n, np.int8)
    # Gray code over coordinates 1..n-1
    N = 1 << (n - 1)
    for it in range(N):
        if it > 0:
            b = 0; x = it
            while (x & 1) == 0:
                x >>= 1; b += 1
            c = b + 1
            eps[c] = -eps[c]
            for j in range(m):
                v[j] += 2 * eps[c] * A[j, c]
        ok = True
        for j in range(m):
            if v[j] != target and v[j] != -target:
                ok = False; break
        if ok:
            if k < cap:
                out[k] = eps
            k += 1
    return out[:min(k, cap)], k


def admissible(A, target, cap=200000, tmax=60.0):
    n = A.shape[1]
    if n <= 28:
        S, k = _brute(A.astype(np.int64), target, cap)
        return S.astype(np.int64), k
    from ortools.sat.python import cp_model
    md = cp_model.CpModel()
    e = [md.NewBoolVar(f'e{i}') for i in range(n)]
    s = [md.NewBoolVar(f's{i}') for i in range(A.shape[0])]
    md.Add(e[0] == 1)
    for j in range(A.shape[0]):
        # sum A_jc (2e_c - 1) = target (2 s_j - 1)
        md.Add(sum(2 * int(A[j, c]) * e[c] for c in range(n) if A[j, c]) - int(A[j].sum())
               == 2 * target * s[j] - target)
    sol = cp_model.CpSolver(); sol.parameters.enumerate_all_solutions = True
    sol.parameters.max_time_in_seconds = tmax; sol.parameters.num_workers = 1
    res = []

    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(s_):
            super().__init__()
        def on_solution_callback(s_):
            res.append([2 * s_.Value(x) - 1 for x in e])
            if len(res) >= cap:
                s_.StopSearch()
    st = sol.Solve(md, CB())
    return np.array(res, np.int64).reshape(-1, n), (len(res) if st == cp_model.OPTIMAL else -len(res))


def clique(S, size, tmax=10.0):
    """size mutually orthogonal rows of S (greedy-ordered branch and bound with a time limit)"""
    if len(S) < size:
        return None
    G = (S @ S.T == 0)
    deg = G.sum(1)
    order = np.argsort(-deg)
    S = S[order]; G = G[np.ix_(order, order)]
    t0 = time.time(); best = [None]
    def rec(cands, chosen):
        if len(chosen) == size:
            best[0] = list(chosen); return True
        if time.time() - t0 > tmax or len(chosen) + len(cands) < size:
            return False
        for k, u in enumerate(cands):
            if len(chosen) + len(cands) - k < size:
                return False
            nc = [x for x in cands[k + 1:] if G[u, x]]
            chosen.append(u)
            if rec(nc, chosen):
                return True
            chosen.pop()
        return False
    rec(list(range(len(S))), [])
    return None if best[0] is None else S[best[0]]


def border(Qg, z, frames=200, seed=1, tclique=10.0, tadm=60.0, verbose=False):
    q = Qg.n
    y = 2 * np.asarray(z, np.int64)
    Mc = BB.gs_array(list(y), Qg.add, Qg.neg)
    assert (Mc @ Mc.T == 4 * q * np.eye(4 * q)).all()
    H = paley(4 * q)
    rng = np.random.default_rng(seed)
    stats = []
    for k in range(frames):
        if k == 0:
            Rh = H.copy()
        else:
            Rh = H[:, rng.permutation(4 * q)] * rng.choice([-1, 1], 4 * q)[None, :]
        A = Mc @ Rh.T
        S, cnt = admissible(A, 4 * q, tmax=tadm)
        stats.append(int(cnt))
        if verbose:
            print(k, cnt, flush=True)
        if len(S) < 4 * q:
            continue
        Xc = clique(S, 4 * q, tmax=tclique)
        if Xc is None:
            continue
        X = Xc.T.copy()
        num = -(A @ X)
        if (num % (4 * q)).any():
            continue
        Ch = num // (4 * q)
        if (np.abs(Ch) == 1).all() and (Ch @ Ch.T == 4 * q * np.eye(4 * q)).all():
            return dict(found=True, frame=k, X=X.tolist(), Rh=Rh.tolist(), Ch=Ch.tolist(), adm=stats)
    return dict(found=False, frames=frames, adm=stats)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('pfile'); ap.add_argument('--frames', type=int, default=200)
    ap.add_argument('--out', default=None); ap.add_argument('--only', type=int, default=None)
    ap.add_argument('--tadm', type=float, default=60.0)
    a = ap.parse_args()
    P = json.load(open(a.pfile))
    Qg = ZP.QGroup(tuple(P['mods']))
    done = set()
    if a.out and os.path.exists(a.out):
        done = {json.loads(l)['index'] for l in open(a.out)}
    for i, p in enumerate(P['patterns']):
        if (a.only is not None and i != a.only) or i in done:
            continue
        t0 = time.time()
        r = border(Qg, p['z'], frames=a.frames, tadm=a.tadm)
        rec = dict(q=P['q'], mods=P['mods'], index=i, z=p['z'], found=r['found'], frame=r.get('frame'),
                   adm_first=r['adm'][:10], n_frames=len(r['adm']), secs=round(time.time() - t0, 1))
        if r['found']:
            rec.update(X=r['X'], Rh=r['Rh'], Ch=r['Ch'])
        print(json.dumps({k: rec[k] for k in ('q', 'index', 'found', 'frame', 'adm_first', 'n_frames', 'secs')}),
              flush=True)
        if a.out:
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
