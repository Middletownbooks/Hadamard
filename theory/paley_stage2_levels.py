r"""Stage 2 (u given y) at the Q-level, v = 359 with the published y: exact count of Q-level solutions.

For Q = G'/H' elementary abelian (all characters real) the Q-level part of S2 reads, for every xi in Q^,
    sum_k uhat_Q,k(xi)^2 = (4v - sum_k betahat_Q,k(xi)^2) / p,
where u_Q,k(q) = sum of u over the twisted positions of coset q (beta_Q is fixed by stage 1).  We count the 4-tuples
(u_Q,A..u_Q,D) solving this, unweighted and weighted by the number of u-configurations realising each u_Q.
usage: python3 theory/paley_stage2_levels.py"""
import sys, os, json, math, itertools as it
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
import paley_ls as LSm
from paley_algebra import load_core, components
from paley_stage1 import coset_index

def main():
    v, b, mods, X = load_core(os.path.join(ROOT, 'data/core_v359_b7.txt')); p = 11
    P = LSm.Prob(359, 7, 11, (2, 2, 8), [(0, 0, 2)])
    comps = [components(x, p) for x in X]
    A0 = [c[0].reshape(-1) for c in comps]; U0 = [c[1].reshape(-1) for c in comps]; W0 = [c[2].reshape(-1) for c in comps]
    rep = coset_index(P); nq = rep.max() + 1
    co = P.G.co
    # Q characters: Q = G'/H' with coset index q = a*4 + b*2 + (e mod 2); characters (-1)^(x.a + y.b + z.(e mod 2))
    qco = np.array([[q // 4, (q // 2) % 2, q % 2] for q in range(nq)])
    chars = np.array([[(-1) ** int((xi * qc).sum()) for qc in qco] for xi in qco])          # (8 chars, 8 cosets)
    beta = [A0[k] - W0[k] for k in range(4)]
    betaQ = np.array([[beta[k][rep == q].sum() for q in range(nq)] for k in range(4)])
    tgt = (4 * v - ((chars @ betaQ.T) ** 2).sum(1))
    assert (tgt % p == 0).all(); tgt //= p
    per = []
    for k in range(4):
        ntw = [int(((rep == q) & (W0[k] == 0)).sum()) for q in range(nq)]
        rng = [list(range(-n, n + 1, 2)) for n in ntw]
        vecs = np.array(list(it.product(*rng)), np.int64)
        mult = np.prod([[math.comb(n, (n - s) // 2) for s in vecs[:, q]] for q, n in enumerate(ntw)], axis=0)
        spec = (vecs @ chars.T) ** 2
        keep = spec[:, 0] == 0                                                   # U_k = 0 (trivial character)
        per.append((vecs[keep], mult[keep], spec[keep]))
    # 4-list join on the spectrum vectors, pruned by S <= target (all spectra are >= 0), vectorised counting
    def pair_counts(i, j):
        Si, Mi = per[i][2], per[i][1]; Sj, Mj = per[j][2], per[j][1]
        keys_all, cnt_all, w_all = [], [], []
        base = np.array([1 << (10 * t) for t in range(8)], dtype=object)
        for a0 in range(0, len(Si), 2000):
            S = Si[a0:a0 + 2000, None, :] + Sj[None, :, :]
            ok = (S <= tgt[None, None, :]).all(2)
            ia, ib = np.nonzero(ok)
            if len(ia) == 0: continue
            Sv = S[ia, ib]                                                      # (M, 8), entries < 1024
            key = (Sv * (1 << (10 * np.arange(8, dtype=np.int64))) [None, :]).sum(1)
            keys_all.append(key); cnt_all.append(np.ones(len(key), np.int64))
            w_all.append(Mi[a0 + ia].astype(np.float64) * Mj[ib].astype(np.float64))
        if not keys_all: return {}, {}
        key = np.concatenate(keys_all); w = np.concatenate(w_all)
        uk, inv = np.unique(key, return_inverse=True)
        c = np.bincount(inv); ww = np.bincount(inv, weights=w)
        return dict(zip(uk.tolist(), c.tolist())), dict(zip(uk.tolist(), ww.tolist()))
    assert (tgt < 1024).all()
    c01, w01 = pair_counts(0, 1); c23, w23 = pair_counts(2, 3)
    tkey = int((tgt * (1 << (10 * np.arange(8, dtype=np.int64)))).sum())
    n_unw = 0; n_w = 0.0
    for k_, c1 in c01.items():
        k2 = tkey - k_
        if k2 in c23: n_unw += c1 * c23[k2]; n_w += w01[k_] * w23[k2]
    pubQ = [tuple(int(U0[k][(rep == q) & (W0[k] == 0)].sum()) for q in range(nq)) for k in range(4)]
    out = dict(per_block_vectors=[len(x[0]) for x in per], Q_level_solutions=n_unw, log2_Q_solutions=round(math.log2(max(1, n_unw)), 1),
               u_configs_passing_Q_level=float(n_w), log2_u_configs_passing=round(math.log2(max(1.0, n_w)), 1),
               log2_u_configs_total=sum(int((W0[k] == 0).sum()) for k in range(4)), target_per_char=tgt.tolist(),
               published_uQ=pubQ)
    print(json.dumps(out))

if __name__ == "__main__":
    main()
