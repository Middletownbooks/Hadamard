r"""Borders for c >= 10: search for a flat corner (report §3.1) [M].

A b >= 3 template for c is (z, K, Kc, C): z = four Q-compressions satisfying the Q-level equation, K and Kc frames
(4b x 4c, +-1, K K^T = 4c I, rows balanced on each block), and the forced corner C = -K Cb^T Kc^T/(4c) a +-1 matrix
(Cb = GS array of z, layout of the 347 template).  By Lemma 2, C C^T = 4b I always, and
    k in rows(K)  <=>  k block-balanced +-1  and  c(k) := -k Cb^T Kc^T/(4c) in {+-1}^{4b}     (then C's row = c(k)).
So, for fixed z and Kc, the admissible rows S = {k} are enumerated by a meet in the middle over the four blocks
(the condition is linear in k: residues mod 8c, then magnitudes), and a template is a clique of 4b mutually
orthogonal elements of S.  The h = 1 instance (z_k in {+-1}^Q, sum_k PAF_Q(z_k) = 4 off 0) is the natural source of
compressions: every template transfers along the all-ones direction to any h (report §3).
usage: python3 theory/ab_border10.py --c 6 --published        (validation: 347's z and Kc must give S containing K)
       python3 theory/ab_border10.py --c 10 [--trials N] [--seed S]"""
import sys, os, json, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import AG
import ab_assemble as AS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAYOUT = AS.load_template347()['layout']          # 347's GS layout (untwisted), reused on Q = Z_c


def qlevel_h1(c):
    """all (z_0..z_3) in {+-1}^{Z_c} with sum_k PAF = 4 off 0; one representative per (block-wise shift, negation)
    class is not taken: returns all solutions as an array (n x 4 x c) together with the PAF table"""
    Q = AG((c,))
    X = np.array(list(itertools.product((1, -1), repeat=c)), np.int64)
    P = Q.paf(X)[:, 1:c // 2 + 1]                  # lags 1..c/2 determine the rest
    key = defaultdict(list)
    for i, row in enumerate(map(tuple, P)): key[row].append(i)
    uk = list(key)
    U = np.array(uk, np.int64)
    pairs = defaultdict(list)
    for i in range(len(uk)):
        s = U[i] + U[i:]
        for j, r in enumerate(map(tuple, s)): pairs[r].append((i, i + j))
    target = np.full(c // 2, 4, np.int64)
    sols = []
    for k1, lst in pairs.items():
        need = tuple(target - np.array(k1))
        if need not in pairs or need < k1: continue
        for (a, b) in lst:
            for (cc, d) in pairs[need]:
                sols.append((uk[a], uk[b], uk[cc], uk[d]))
    return X, key, sols


def Cb_of(z, c):
    return AS.compressed_core(np.asarray(z), (c,), LAYOUT)


def frame_lemma1(c, seed=0):
    """frame from H2 x H(2c) (Paley I needs 2c - 1 prime = 3 mod 4), columns grouped into the four pattern blocks; the
    order of columns inside each block is randomised by seed (seed 0: natural order)"""
    q = 2 * c - 1
    chi = np.array([0] + [1 if pow(x, (q - 1) // 2, q) == 1 else -1 for x in range(1, q)])
    S = np.zeros((q + 1, q + 1), int); S[0, 1:] = 1; S[1:, 0] = -1
    S[1:, 1:] = np.array([[chi[(j - i) % q] for j in range(q)] for i in range(q)])
    H = S + np.eye(q + 1, dtype=int)
    H = H * H[0][None, :]
    H4c = np.kron(np.array([[1, 1], [1, -1]]), H)
    rows = [0, 1, 2 * c, 2 * c + 1]
    R = H4c[rows]
    K = np.delete(H4c, rows, axis=0)
    # blocks: patterns ordered as H4 = [[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]] (columns)
    pats = [(1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)]
    rng = np.random.default_rng(seed)
    cols = []
    for p in pats:
        idx = [j for j in range(4 * c) if tuple(R[:, j]) == p]
        assert len(idx) == c
        if seed: idx = list(rng.permutation(idx))
        cols += idx
    F = K[:, cols]
    assert (F @ F.T == 4 * c * np.eye(4 * c - 4)).all() and (F.reshape(4 * c - 4, 4, c).sum(2) == 0).all()
    return F


def admissible(Cb, Kc, c):
    """S = {k block-balanced +-1 : -k Cb^T Kc^T/(4c) in {+-1}^{4b}} by a meet in the middle over blocks (0,1) | (2,3)"""
    b4 = Kc.shape[0]; M = Cb.T @ Kc.T                # 4c x 4b; condition on k M
    bal = np.array([v for v in itertools.product((1, -1), repeat=c) if sum(v) == 0], np.int64)
    def side(J0, J1):
        A0 = bal @ M[J0 * c:(J0 + 1) * c]; A1 = bal @ M[J1 * c:(J1 + 1) * c]
        U = (A0[:, None, :] + A1[None, :, :]).reshape(-1, b4)
        return U
    L = side(0, 1); R = side(2, 3)
    mod = 8 * c
    rng = np.random.default_rng(1); w = rng.integers(1, 2 ** 61, size=b4, dtype=np.int64)
    hl = ((L % mod) @ w); hr = (((4 * c - R) % mod) @ w)
    o = np.argsort(hl); hs = hl[o]
    lo = np.searchsorted(hs, hr, 'left'); hi = np.searchsorted(hs, hr, 'right')
    out = []
    nb = len(bal)
    for r in np.nonzero(hi > lo)[0]:
        for p in o[lo[r]:hi[r]]:
            u = L[p] + R[r]
            if (np.abs(u) == 4 * c).all():
                i0, i1 = divmod(int(p), nb); i2, i3 = divmod(int(r), nb)
                out.append(np.concatenate([bal[i0], bal[i1], bal[i2], bal[i3]]))
    return np.array(out, np.int64).reshape(-1, 4 * c)


def clique(S, size, tmax=60.0):
    """a set of `size` mutually orthogonal rows of S (up to sign), by backtracking with a time limit"""
    if len(S) == 0: return None
    # identify k and -k
    keep = {}
    for v in S:
        t = tuple(v if v[0] > 0 else -v)
        keep[t] = 1
    V = np.array(list(keep), np.int64)
    G = (V @ V.T == 0)
    t0 = time.time(); best = [None]
    def rec(cands, chosen):
        if len(chosen) == size: best[0] = list(chosen); return True
        if time.time() - t0 > tmax: return False
        if len(chosen) + len(cands) < size: return False
        for i, u in enumerate(cands):
            nc = [x for x in cands[i + 1:] if G[u, x]]
            chosen.append(u)
            if rec(nc, chosen): return True
            chosen.pop()
            if len(chosen) + len(cands) - i - 1 < size: return False
        return False
    rec(list(range(len(V))), [])
    return V[best[0]] if best[0] is not None else None, len(V), G.sum(1).mean() if len(V) else 0


def main():
    a = sys.argv[1:]
    c = int(a[a.index('--c') + 1]) if '--c' in a else 10
    b = c - 1
    log = lambda r: print(json.dumps(r), flush=True)
    if '--published' in a:
        T = AS.load_template347(); assert c == 6
        # shift each block of 347's compression to +-1 values (same deviations)
        comp = T['comp']; z = np.zeros_like(comp)
        for k in range(4):
            lo_ = comp[k].min(); z[k] = np.where(comp[k] == lo_, -1, 1) if len(set(comp[k])) > 1 else -np.ones(c, int)
        Q = AG((2, 3)); Pz = Q.paf(z).sum(0)
        log(dict(stage="z", z=z.tolist(), paf=Pz.tolist()))
        Cb = AS.compressed_core(z, (2, 3), T['layout'])
        S = admissible(Cb, T['Kc'], c)
        inS = sum(any((s == k).all() or (s == -k).all() for s in S) for k in T['K'])
        log(dict(stage="published", admissible=len(S), K_rows_in_S=int(inS)))
        return
    t0 = time.time()
    X, key, sols = qlevel_h1(c)
    log(dict(stage="qlevel", c=c, solutions=len(sols), secs=round(time.time() - t0, 1)))


if __name__ == "__main__":
    main()
