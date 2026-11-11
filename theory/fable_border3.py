#!/usr/bin/env python3
"""Single-character borders for partial character borders over Q = Z_2^a (see fable_border2.py).

Every row of E and of Z is, block by block, one signed character of T.  In the basis e_a (x) rho the compressed
core acts as D = sum_rho Gamma_rho (4 x 4 integer blocks), and  K = C_E D C_Z^T / 4  must be a +-1 matrix with
C_E, C_Z the 0/+-1 coefficient matrices (one nonzero per row and block), rows of each pairwise orthogonal.
Search: candidates = all (2t)^4 rows; compat[e, z] = (entry is +-4); depth-first choice of 4t orthogonal rows for E
keeping the common compatible set for Z, then 4t orthogonal rows for Z inside it.  A found border is assembled and
H H^T = n I is checked on the full matrix.  [M]   Usage: fable_border3.py a s T ncores [seed] [maxnodes]
"""
import sys, time, itertools
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2

def search(G, T, M, maxnodes=2_000_000, rng=None):
    k = G.k; t = len(T); b = 4 * t
    ch = G.chars()
    # basis vectors v[a][i] = e_a (x) rho_i  (length 4k)
    V = np.zeros((4, t, 4 * k), dtype=np.int64)
    for a in range(4):
        for i, m in enumerate(T):
            V[a, i, a * k:(a + 1) * k] = ch[m]
    D = np.einsum('aic,dc,bjd->aibj', V, M, V) // k        # D[a,i,b,j] = v_ai M^T v_bj^T / k   (M^T[c,d] = M[d,c])
    assert all(D[a, i, b_, j] == 0 for a in range(4) for b_ in range(4) for i in range(t) for j in range(t) if i != j), 'not block diagonal'
    cands = list(itertools.product(range(2 * t), repeat=4))            # option o: character o % t, sign + if o < t
    nc = len(cands)
    ci = np.array([[o % t for o in c] for c in cands]); cs = np.array([[1 if o < t else -1 for o in c] for c in cands])
    # dense coefficient vectors (nc x 4t)
    C = np.zeros((nc, 4 * t), dtype=np.int8)
    for a in range(4):
        C[np.arange(nc), a * t + ci[:, a]] = cs[:, a]
    Dm = D.reshape(4 * t, 4 * t)
    val = (C.astype(np.int32) @ Dm.astype(np.int32)) @ C.T.astype(np.int32)       # e D z^T
    compat = (np.abs(val) == 4)
    orth = (C.astype(np.int32) @ C.T.astype(np.int32)) == 0
    nodes = [0]; best = [0]
    order = np.argsort(-compat.sum(1))                                  # rows with many compatible partners first
    if rng is not None: order = rng.permutation(nc)
    def find_orth(cset, need, chosen):
        """need more mutually orthogonal rows inside cset (boolean mask), orthogonal to chosen."""
        if need == 0: return chosen
        idx = np.nonzero(cset)[0]
        if len(idx) < need: return None
        nodes[0] += 1
        if nodes[0] > maxnodes: return None
        for x in idx:
            r = find_orth(cset & orth[x] & (np.arange(nc) > x), need - 1, chosen + [x])
            if r is not None: return r
            if nodes[0] > maxnodes: return None
        return None
    def dfs(Eset, emask, zmask):
        if len(Eset) == b:
            Zs = find_orth(zmask, b, [])
            return (Eset, Zs) if Zs is not None else None
        nodes[0] += 1
        if nodes[0] > maxnodes: return None
        if zmask.sum() < b: return None
        best[0] = max(best[0], len(Eset))
        idx = np.nonzero(emask)[0]
        if len(idx) < b - len(Eset): return None
        for x in idx:
            r = dfs(Eset + [x], emask & orth[x] & (np.arange(nc) > x), zmask & compat[x])
            if r is not None: return r
            if nodes[0] > maxnodes: return None
        return None
    res = dfs([], np.ones(nc, dtype=bool), np.ones(nc, dtype=bool))
    if res is None: return None, nodes[0], best[0]
    Eset, Zs = res
    # function-space rows
    def rows(sel):
        R = np.zeros((b, 4 * k), dtype=np.int64)
        for r, x in enumerate(sel):
            for a in range(4):
                R[r, a * k:(a + 1) * k] = cs[x, a] * ch[T[ci[x, a]]]
        return R
    E = rows(Eset); Z = rows(Zs)
    K = (E @ M.T @ Z.T) // (4 * k)
    return (K, E, Z), nodes[0], best[0]

def main():
    a = int(sys.argv[1]); s = int(sys.argv[2]); T = [int(x) for x in sys.argv[3].split(',')]
    ncores = int(sys.argv[4]); seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    maxnodes = int(float(sys.argv[6])) if len(sys.argv) > 6 else 2_000_000
    G = B2.Grp2(a, s); rng = np.random.default_rng(seed)
    k = G.k; b = 4 * len(T); n = 4 * G.v + b
    print(f'Q = Z_2^{a}, S = Z_{s}, T = {T}, |G| = {G.v}, b = {b}, n = {n}', flush=True)
    for ci_ in range(ncores):
        t0 = time.time(); X = None; tries = 0
        while X is None and tries < 30:
            X = B2.find_core(G, T, rng); tries += 1
        if X is None:
            print(f'core {ci_}: none found'); continue
        Cm = FB.gs_array(G, X); M = FB.compress(G, Cm)
        ch = G.chars()
        prof = [[int(sum(ch[m][q] * M[0, j * k + q] for q in range(k))) for j in range(4)] for m in T]
        print(f'core {ci_}: found in {time.time()-t0:.0f}s; Walsh values of the compressions on T: {prof}', flush=True)
        t1 = time.time()
        res, nodes, best = search(G, T, M, maxnodes)
        if res is None:
            print(f'   single-character border: none found ({nodes} nodes, deepest E-set {best}, {time.time()-t1:.0f}s){" [node cap: undecided]" if nodes > maxnodes else " [exhaustive]"}', flush=True)
            continue
        K, E, Z = res
        okK = np.array_equal(np.abs(K), np.ones_like(K)) and np.array_equal(K @ K.T, b * np.eye(b, dtype=int))
        H = FB.assemble(G, Cm, K, E, Z)
        good = np.array_equal(H @ H.T, n * np.eye(n, dtype=int))
        print(f'   single-character border found ({nodes} nodes, {time.time()-t1:.0f}s): K Hadamard {okK}; assembled H({n}): H H^T = n I is {good}', flush=True)

if __name__ == '__main__':
    main()
