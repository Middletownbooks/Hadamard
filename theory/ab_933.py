r"""b = 0 at v = 933 = 3 x 311 in the GP slice of index 10 on Z_311 (f = 31, G' = Z_3): staged solver [P, M].

Group G = Z_3 x Z_311 (= Z_933 by CRT).  Slice: x_k(g, 0) = alpha_k(g), x_k(g, z) = c_{k,j}(g) on the 10 cosets C_j of the
order-31 subgroup M of F_311^*.  Levels (all exact identities):
  stage 1  (compression along Z_311 onto Z_3): y_k(g) = alpha + 31 sigma, sum_k PAF_{Z_3}(y_k) = 4v delta; 36 solutions
           up to block order and sign (theory/ab_b0_stage1.py);
  dual     (compression along Z_3 onto Z_311): xbar_k = (A_k at 0, T_{k,j} on C_j), T_{k,j} = sum_g c_{k,j}(g) in
           {+-1, +-3}, sum_k PAF_{Z_311}(xbar_k) = 4v delta  (6 integer equations at the <M, -1>-lag classes);
  cross    the remaining characters (phi != 1 on Z_3, psi != 1 on Z_311).
Given y, block k's class signs form a 3 x 10 sign matrix with row sums sigma_k(g) (stage 1) and column sums T_k (dual).
Pipeline: per block all dual vectors T (with the number m(T) of sign matrices) under the spectral bound
|xbar^(psi)|^2 <= 4v; a Wagner merge on the dual signature (first level: energy and lag class 1 summed over a pair)
streams dual solutions; for each, every block's fibre (all sign matrices with those margins) is listed and an exact
2+2 join on the full PAF over G (all <M, -1>-lag classes of G) finds the cores; every hit is re-verified, mapped to
Z_933 and completed by the Goethals-Seidel array (H(3732), H H^T = 3732 I checked).
usage: python3 theory/ab_933.py --y I [--middles N] [--tag T] [--seed S] [--stop-on-core]"""
import sys, os, json, time, itertools, hashlib
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json
from ab_b0_stage1 import stage1
from ab_b0_dual import col_vectors
from ab_lib import AG, primitive_root
from ab_sweep import mgen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(ROOT, 'logs/ab/cfg_933_s10.json')
V, F, Q_, S_, FF = 933, 3, 311, 10, 31


class Slice933:
    def __init__(self):
        self.cfg = cfg_from_json('@' + CFG); self.sl = W.build(self.cfg)
        sl = self.sl
        # classes on Z_q
        g = mgen(Q_, S_); M = [1]
        while True:
            x = M[-1] * g % Q_
            if x == 1: break
            M.append(x)
        self.M = M
        gam = primitive_root(Q_)
        cls = -np.ones(Q_, np.int64); cls[0] = S_
        for j in range(S_):
            rep = pow(gam, j, Q_)
            for m_ in M: cls[rep * m_ % Q_] = j
        self.cls = cls
        self.Zq = AG((Q_,))
        self.qreps, _ = AGs.lag_reps(self.Zq, [self.Zq.act((m_,)) for m_ in M])
        self.Tq = np.zeros(Q_, np.int64); self.Tq[0] = 4 * V
        # Gauss periods for the spectral bound
        zeta = np.exp(2j * np.pi * np.arange(Q_) / Q_)
        reps_a = [pow(gam, i, Q_) for i in range(S_)]
        self.eta = np.array([[zeta[(a * np.nonzero(cls == j)[0]) % Q_].sum() for j in range(S_)] for a in reps_a])
        # G lag classes (for the full PAF join)
        self.G = sl.G
        self.greps, _ = AGs.lag_reps(sl.G, sl.perms)
        self.T = sl.T[self.greps]
        # map orbit labels: for block assembly, element (g, z) -> x value from alpha(g) or c_j(g)
        self.gco = sl.G.co                       # (g, z)

    def alpha_sigma(self, yk):
        al = np.array([1 if (int(yk[b]) - 1) % (2 * FF) == 0 else -1 for b in range(F)])
        sig = (np.array(yk) - al) // FF
        return al, sig

    def dual_lists(self, Y):
        """per block: dual vectors T (n x 10), multiplicities, dual PAF signatures, after the spectral bound"""
        out = []
        for k in range(4):
            al, sig = self.alpha_sigma(Y[k])
            cv = col_vectors([int(x) for x in sig], S_)
            vecs = np.array(list(cv.keys()), np.int64); mult = np.array(list(cv.values()), np.float64)
            A = int(al.sum())
            X = A + vecs.astype(float) @ self.eta.T
            ok = ((np.abs(X) ** 2) <= 4 * V + 1e-6).all(1)
            vecs, mult = vecs[ok], mult[ok]
            xb = np.zeros((len(vecs), Q_), np.int64); xb[:, 0] = A
            for j in range(S_): xb[:, self.cls == j] = vecs[:, j:j + 1]
            P = self.Zq.paf(xb)[:, self.qreps]
            out.append(dict(alpha=al, sigma=sig, A=A, vecs=vecs, mult=mult, sig=P))
        return out

    def fibre(self, al, sig, T):
        """all 3 x 10 sign matrices with row sums sig and column sums T -> list of x on G (each 933 values)"""
        cols_opts = []
        for t in T:
            opts = [np.array(c) for c in itertools.product((1, -1), repeat=F) if sum(c) == t]
            cols_opts.append(opts)
        mats = []
        def rec(j, rows):
            if j == S_:
                if all(r == s_ for r, s_ in zip(rows, sig)): mats.append(list(cur))
                return
            rem = S_ - j - 1
            for c in cols_opts[j]:
                r2 = [r + int(x) for r, x in zip(rows, c)]
                if all(abs(s_ - r) <= rem for s_, r in zip(sig, r2)):
                    cur.append(c); rec(j + 1, r2); cur.pop()
        cur = []
        rec(0, [0] * F)
        if not mats: return np.zeros((0, self.G.n), np.int64)
        X = np.zeros((len(mats), self.G.n), np.int64)
        g_ = self.gco[:, 0]; z_ = self.gco[:, 1]
        zc = self.cls[z_]
        for i, m in enumerate(mats):
            Mx = np.array(m).T                      # F x 10
            vals = np.where(z_ == 0, al[g_], Mx[g_, np.minimum(zc, S_ - 1)])
            X[i] = vals
        return X


def wagner_dual(D, nmid=50, seed=0, max_pairs=4 * 10 ** 7):
    """stream dual solutions: middle value = (energy, lag-1 PAF) of the left pair; yields (i0, i1, i2, i3)"""
    rng = np.random.default_rng(seed)
    P = [d['sig'] for d in D]
    key = lambda X: X[:, 0] * 10 ** 6 + X[:, 1]             # energy and lag class 1 (small ranges)
    k = [key(p) for p in P]
    tkey = None
    # candidate middle values: sample pairs
    n = [len(p) for p in P]
    i0 = rng.integers(n[0], size=200000); i1 = rng.integers(n[1], size=200000)
    L = k[0][i0] + k[1][i1]
    vals, cnt = np.unique(L, return_counts=True)
    order = vals[np.argsort(-cnt)][:nmid]
    by1 = defaultdict(list)
    for j, x in enumerate(k[1].tolist()): by1[x].append(j)
    by3 = defaultdict(list)
    for j, x in enumerate(k[3].tolist()): by3[x].append(j)
    return order, by1, by3, k


def main():
    a = sys.argv[1:]
    yi = int(a[a.index('--y') + 1]); tag = a[a.index('--tag') + 1] if '--tag' in a else f'y{yi}'
    nmid = int(a[a.index('--middles') + 1]) if '--middles' in a else 40
    seed = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    outdir = os.path.join(ROOT, 'logs/ab/p933'); os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, f'run_{tag}.jsonl')
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    S = Slice933(); sl = S.sl
    sols, _ = stage1(sl)
    Y = sols[yi]
    D = S.dual_lists(Y)
    log(dict(stage="setup", y=Y.tolist(), dual_list=[len(d['vecs']) for d in D]))
    order, by1, by3, k = wagner_dual(D, nmid=nmid, seed=seed)
    P = [d['sig'] for d in D]
    tk = int(S.Tq[S.qreps][0] * 10 ** 6 + S.Tq[S.qreps][1])
    Tsig = S.Tq[S.qreps]
    rng = np.random.default_rng(seed + 7); wv = rng.integers(1, 2 ** 62, size=P[0].shape[1], dtype=np.int64)
    gw = rng.integers(1, 2 ** 62, size=len(S.greps), dtype=np.int64)
    ncores = 0; tot_dual = 0; tot_cross_pairs = 0
    for mi, c in enumerate(order):
        t0 = time.time()
        # left pairs with key sum c, right pairs with key sum tk - c
        left = [(i, j) for i in range(len(k[0])) for j in by1.get(int(c - k[0][i]), [])]
        right = [(i, j) for i in range(len(k[2])) for j in by3.get(int(tk - c - k[2][i]), [])]
        if not left or not right:
            log(dict(mid=int(c), left=len(left), right=len(right), dual=0)); continue
        L = np.array(left); R = np.array(right)
        HL = P[0][L[:, 0]] @ wv + P[1][L[:, 1]] @ wv
        HR = np.int64(Tsig @ wv) - (P[2][R[:, 0]] @ wv + P[3][R[:, 1]] @ wv)
        o = np.argsort(HL); HLs = HL[o]
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        duals = []
        for r in np.nonzero(hi > lo)[0]:
            for p in o[lo[r]:hi[r]]:
                q4 = (L[p, 0], L[p, 1], R[r, 0], R[r, 1])
                if ((P[0][q4[0]] + P[1][q4[1]] + P[2][q4[2]] + P[3][q4[3]]) == Tsig).all(): duals.append(q4)
        tot_dual += len(duals)
        # cross level for each dual solution
        hits = 0; logm = []
        t1 = time.time()
        for q4 in duals:
            Xs = [S.fibre(D[kk]['alpha'], D[kk]['sigma'], D[kk]['vecs'][q4[kk]]) for kk in range(4)]
            logm.append(sum(np.log2(max(len(x), 1)) for x in Xs))
            if any(len(x) == 0 for x in Xs): continue
            Ps = [S.G.paf(x)[:, S.greps] for x in Xs]
            h = [p_ @ gw for p_ in Ps]
            hl = (h[0][:, None] + h[1][None, :]).ravel()
            hr = (np.int64(S.T @ gw) - (h[2][:, None] + h[3][None, :])).ravel()
            oo = np.argsort(hl); hs = hl[oo]
            lo2 = np.searchsorted(hs, hr, 'left'); hi2 = np.searchsorted(hs, hr, 'right')
            tot_cross_pairs += len(hl) + len(hr)
            n1, n3 = len(h[1]), len(h[3])
            for rr in np.nonzero(hi2 > lo2)[0]:
                c2, d2 = divmod(int(rr), n3)
                for pp in oo[lo2[rr]:hi2[rr]]:
                    a2, b2 = divmod(int(pp), n1)
                    Xc = np.array([Xs[0][a2], Xs[1][b2], Xs[2][c2], Xs[3][d2]])
                    if (S.G.paf(Xc).sum(0) == sl.T).all():
                        hits += 1; ncores += 1
                        np.save(os.path.join(outdir, f'core_{tag}_{ncores}.npy'), Xc)
                        log(dict(CORE=f'{tag}_{ncores}', verified_paf=True))
        log(dict(mid=int(c), left=len(left), right=len(right), dual=len(duals),
                 mean_log2_fibre_product=round(float(np.mean(logm)), 2) if logm else None, cores=hits,
                 dual_secs=round(t1 - t0, 1), cross_secs=round(time.time() - t1, 1)))
        if ncores and '--stop-on-core' in a: break
    log(dict(stage="end", total_dual=tot_dual, cores=ncores, cross_pairs=tot_cross_pairs))


if __name__ == "__main__":
    main()
