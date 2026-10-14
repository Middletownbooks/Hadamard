r"""Exact dual-compression test for b = 0 GP slices on Z_F x Z_q (the analogue of nd_fq.dual_row_screen) [P, M].

[P] Compressing x along Z_F gives xbar_k on Z_q: xbar_k(0) = A_k = sum_g alpha_k(g), xbar_k(C_j) = T_{k,j} =
sum_g c_{k,j}(g) (a sum of F signs), M-invariant, with sum_k PAF_{Z_q}(xbar_k) = 4v delta (b = 0).  Given the stage-1
y, the per-coset class signs of block k form an F x s sign matrix with prescribed row sums sigma_k(g); its column sums
are the T_{k,j}.  So every stage-2 completion of y induces a solution of this "dual" system with T_k a column-sum
vector of such a matrix.  This is exact and the Gaussian cannot see it (the T's are discrete).
Method: per block, all reachable (A_k, T_k) (column-sum vectors realisable with the given row sums, counted with
multiplicity = number of sign matrices); PAF signature on Z_q at the <M, -1>-lag classes; 2+2 hashed join.
usage: python3 theory/ab_b0_dual.py CFGJSON [--ys all|N]"""
import sys, os, json, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json
from ab_b0_stage1 import stage1
from ab_lib import AG

def col_vectors(rows_sigma, s):
    """all column-sum vectors (t_1..t_s) of F x s sign matrices with row sums rows_sigma, with multiplicities"""
    F = len(rows_sigma)
    # DP over columns: state = tuple of remaining row sums
    cols = [np.array(c) for c in itertools.product((1, -1), repeat=F)]
    start = tuple(rows_sigma)
    layer = {start: {(): 1}}
    for j in range(s):
        new = defaultdict(lambda: defaultdict(int))
        rem = s - j - 1
        for st, d in layer.items():
            for c in cols:
                st2 = tuple(a - b for a, b in zip(st, c))
                if any(abs(x) > rem for x in st2) or any((x - rem) % 2 for x in st2): continue
                t = int(c.sum())
                for pre, m in d.items(): new[st2][pre + (t,)] += m
        layer = new
    out = layer.get(tuple([0] * F), {})
    return out

def dual_test(sl, Y, cfg):
    q = sl.q; F = sl.Gp.n; s = cfg['s']; f = cfg['f']
    # class structure on Z_q: M = index-s subgroup; classes C_j
    from ab_sweep import mgen
    g = mgen(q, s)
    M = [1]
    while True:
        x = M[-1] * g % q
        if x == 1: break
        M.append(x)
    from ab_lib import primitive_root
    gam = primitive_root(q)
    cls = -np.ones(q, np.int64); cls[0] = s       # label s = the zero class
    for j in range(s):
        rep = pow(gam, j, q)
        for m_ in M: cls[rep * m_ % q] = j
    assert (cls >= 0).all()
    Zq = AG((q,))
    perms = [Zq.act((m_,)) for m_ in M]
    reps, _ = AGs.lag_reps(Zq, perms)
    T = np.zeros(q, np.int64); T[0] = 4 * sl.v
    per_block = []
    for k in range(4):
        alpha = np.array([int(Y[k][b]) - f * 0 for b in range(F)])  # placeholder; alpha from letters below
        # alpha(g) = y(g) mod 2f -> +-1 ; sigma(g) = (y - alpha)/f
        al = np.array([1 if (int(Y[k][b]) - 1) % (2 * f) == 0 else -1 for b in range(F)])
        sig = (np.array(Y[k]) - al) // f
        A = int(al.sum())
        cv = col_vectors([int(x) for x in sig], s)
        vecs = np.array(list(cv.keys()), np.int64).reshape(-1, s)
        mult = np.array(list(cv.values()), np.float64)
        xb = np.zeros((len(vecs), q), np.int64)
        xb[:, 0] = A
        for j in range(s): xb[:, cls == j] = vecs[:, j:j + 1]
        P = Zq.paf(xb)[:, reps]
        per_block.append((P, mult, vecs, A))
    # dedupe per block by signature (weights summed), then an energy-split 2+2 join (coordinate 0 = energy)
    U = []
    for (P, mult, vecs, A) in per_block:
        up, inv = np.unique(P, axis=0, return_inverse=True)
        w = np.bincount(inv.ravel(), weights=mult)
        U.append((up, w))
    import importlib
    from ab_s1sd import ejoin4
    hits = ejoin4(U[0][0], U[1][0], U[2][0], U[3][0], T[reps].astype(np.int64), max_pairs=6 * 10 ** 8, max_out=10 ** 7)
    sols = len(hits); weighted = 0.0
    for (a_, b_, c_, d_) in hits:
        weighted += U[0][1][a_] * U[1][1][b_] * U[2][1][c_] * U[3][1][d_]
    return dict(dual_vectors=[len(pb[0]) for pb in per_block], dual_signatures=[len(u[0]) for u in U], dual_solutions=sols,
                log2_weighted=round(float(np.log2(weighted)), 2) if weighted else None)

if __name__ == "__main__":
    a = sys.argv[1:]
    cfg = cfg_from_json(a[0]); sl = W.build(cfg)
    sols, nb = stage1(sl)
    print(json.dumps(dict(v=cfg['v'], q=cfg['q'], s=cfg['s'], stage1_solutions=len(sols))), flush=True)
    tot = 0
    for i, Y in enumerate(sols):
        t0 = time.time()
        r = dual_test(sl, Y, cfg)
        tot += r['dual_solutions'] > 0
        print(json.dumps(dict(i=i, rows=Y.sum(1).tolist(), **r, secs=round(time.time() - t0, 1))), flush=True)
    print(json.dumps(dict(y_with_dual_solutions=tot, of=len(sols))))
