r"""b = 0 at v = 3q (q prime, 10 | q - 1) in the GP slice of index 10 on Z_q (G' = Z_3): exact staged solver [P, M].

Group G = Z_3 x Z_q (= Z_v by CRT).  Slice: x_k(g, 0) = alpha_k(g), x_k(g, z) = c_{k,j}(g) on the 10 cosets C_j of the
index-10 subgroup M of F_q^* (|M| = f = (q - 1)/10).  Levels (exact identities, theory/ab_b0_stage1.py, ab_gauss.py):
  stage 1  (compression along Z_q onto Z_3): y_k(g) = alpha_k(g) + f sigma_k(g), sum_k PAF_{Z_3}(y_k) = 4v delta;
  dual     (compression along Z_3 onto Z_q): xbar_k = (A_k at 0, T_{k,j} on C_j), A_k = sum_g alpha_k(g),
           T_{k,j} = sum_g c_{k,j}(g) in {+-1, +-3}, sum_k PAF_{Z_q}(xbar_k) = 4v delta (one equation per lag class of
           <M, -1> on Z_q, lag 0 included);
  cross    the remaining characters (phi != 1 on Z_3 and psi != 1 on Z_q): the full PAF over G at every lag class.
Given y, block k's class signs form a 3 x 10 sign matrix with row sums sigma_k(g) and column sums T_k; m(T) = number
of such matrices (the fibre of T).
Pipeline:
  1. per block, all dual vectors T with m(T) (DP over columns) under the spectral bound |xbar^(psi)|^2 <= 4v;
  2. per block and energy class keep the K vectors of largest m(T) (the yield of a dual solution is proportional to
     its fibre product m0 m1 m2 m3, the cost of its cross join to m0 m1 + m2 m3);
  3. dual join: all pairs (0,1) against all pairs (2,3), 64-bit linear hash of the dual signature, every match
     re-checked exactly; dual solutions sorted by fibre product;
  4. cross join per dual solution: list the four fibres, PAF over G at every lag class, 2+2 hashed join, every match
     re-checked exactly (sum_k PAF_G(x_k) = 4v delta); a core is saved, mapped to Z_v by CRT and completed by the
     Goethals-Seidel array; H H^T = 4v I is checked here and again by the independent scripts/verify_b0_gs.py.
Restart-safe: one JSONL line per dual solution processed (logs/ab/b0z3/run_<tag>.jsonl); done ones are skipped.
Dual solutions are processed in decreasing order of efficiency = log2(m0 m1 m2 m3) - log2(join cost + overhead).
usage: python3 theory/ab_b0z3.py --q 311 --y I [--K 300] [--max-dual N] [--min-score s] [--tag T] [--stop-on-core] [--stats]"""
import sys, os, json, time, itertools, hashlib
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_b0_stage1 import stage1
from ab_b0_dual import col_vectors
from ab_lib import AG, primitive_root
from ab_sweep import mgen

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_, S_ = 3, 10
OVH = 2 ** 15            # per-dual-solution overhead (fibre listing, PAFs) in join-entry units


def cfg_for(q, s=S_):
    return dict(v=3 * q, c=1, b=0, h=None, m=3 * q, Q=(), H1=(3,), q=q, s=s, f=(q - 1) // s, mods=(3, q), nQ=0, U=[])


class SliceZ3:
    def __init__(self, q):
        assert (q - 1) % S_ == 0
        self.q = q; self.f = (q - 1) // S_; self.v = 3 * q
        self.cfg = cfg_for(q); self.sl = W.build(self.cfg)
        g = mgen(q, S_); M = [1]
        while True:
            x = M[-1] * g % q
            if x == 1: break
            M.append(x)
        assert len(M) == self.f
        self.M = M
        gam = primitive_root(q)
        cls = -np.ones(q, np.int64); cls[0] = S_
        for j in range(S_):
            rep = pow(gam, j, q)
            for m_ in M: cls[rep * m_ % q] = j
        assert (cls >= 0).all()
        self.cls = cls
        self.Zq = AG((q,))
        self.qreps, _ = AGs.lag_reps(self.Zq, [self.Zq.act((m_,)) for m_ in M])
        assert self.qreps[0] == 0
        self.Tq = np.zeros(q, np.int64); self.Tq[0] = 4 * self.v
        zeta = np.exp(2j * np.pi * np.arange(q) / q)
        self.eta = np.array([[zeta[(pow(gam, i, q) * np.nonzero(cls == j)[0]) % q].sum() for j in range(S_)]
                             for i in range(S_)])
        self.G = self.sl.G
        self.greps, _ = AGs.lag_reps(self.G, self.sl.perms)
        assert self.greps[0] == 0
        self.T = self.sl.T[self.greps]
        self.gco = self.G.co

    def alpha_sigma(self, yk):
        f = self.f
        al = np.array([1 if (int(yk[b]) - 1) % (2 * f) == 0 else -1 for b in range(F_)])
        sig = (np.array(yk) - al) // f
        assert (al + f * sig == np.array(yk)).all() and (np.abs(sig) <= S_).all() and ((sig - S_) % 2 == 0).all()
        return al, sig

    def dual_lists(self, Y):
        out = []
        for k in range(4):
            al, sig = self.alpha_sigma(Y[k])
            cv = col_vectors([int(x) for x in sig], S_)
            vecs = np.array(list(cv.keys()), np.int64).reshape(-1, S_)
            mult = np.array(list(cv.values()), np.float64)
            A = int(al.sum())
            X = A + vecs.astype(float) @ self.eta.T
            ok = ((np.abs(X) ** 2) <= 4 * self.v + 1e-6).all(1)
            vecs, mult = vecs[ok], mult[ok]
            sig_ = np.zeros((len(vecs), len(self.qreps)), np.int64)
            for c0 in range(0, len(vecs), 20000):
                xb = np.zeros((min(20000, len(vecs) - c0), self.q), np.int64); xb[:, 0] = A
                for j in range(S_): xb[:, self.cls == j] = vecs[c0:c0 + len(xb), j:j + 1]
                sig_[c0:c0 + len(xb)] = self.Zq.paf(xb)[:, self.qreps]
            out.append(dict(alpha=al, sigma=sig, A=A, vecs=vecs, mult=mult, sig=sig_))
        return out

    def fibre_mats(self, sig, T):
        """all 3 x 10 sign matrices with row sums sig and column sums T (n x 3 x 10, int8)"""
        cols_opts = [[np.array(c) for c in itertools.product((1, -1), repeat=F_) if sum(c) == t] for t in T]
        mats = []; cur = []
        def rec(j, rows):
            if j == S_:
                if all(r == s_ for r, s_ in zip(rows, sig)): mats.append(np.array(cur).T.copy())
                return
            rem = S_ - j - 1
            for c in cols_opts[j]:
                r2 = [r + int(x) for r, x in zip(rows, c)]
                if all(abs(int(s_) - r) <= rem for s_, r in zip(sig, r2)):
                    cur.append(c); rec(j + 1, r2); cur.pop()
        rec(0, [0] * F_)
        return np.array(mats, np.int8).reshape(-1, F_, S_)

    def expand(self, Mx, al):
        """sign matrices (n x 3 x 10) and alpha -> x on G (n x |G|)"""
        g_ = self.gco[:, 0]; z_ = self.gco[:, 1]; zc = np.minimum(self.cls[z_], S_ - 1)
        X = Mx[:, g_, zc].astype(np.int8)
        X[:, z_ == 0] = np.asarray(al, np.int8)[g_[z_ == 0]][None, :]
        return X

    def fibre(self, al, sig, T):
        """all 3 x 10 sign matrices with row sums sig and column sums T, as x on G (rows: one per matrix)"""
        Mx = self.fibre_mats(sig, T)
        if len(Mx) == 0: return np.zeros((0, self.G.n), np.int8)
        return self.expand(Mx, al)

    def to_cyclic(self, Xc):
        """CRT: Z_v -> Z_3 x Z_q, n -> (n mod 3, n mod q)"""
        n = np.arange(self.v)
        return Xc[:, self.G.idx(np.stack([n % 3, n % self.q], 1))]


def gs_hadamard(Xcyc):
    """Goethals-Seidel array of four circulants (rows x[c - r]) with the back identity R"""
    v = Xcyc.shape[1]
    circ = lambda x: np.array([np.roll(x, j) for j in range(v)], np.int8)
    A, B, C, D = [circ(x) for x in Xcyc]
    R = lambda Mx: Mx[:, ::-1]
    H = np.block([[A, R(B), R(C), R(D)],
                  [-R(B), A, R(D.T), -R(C.T)],
                  [-R(C), -R(D.T), A, R(B.T)],
                  [-R(D), R(C.T), -R(B.T), A]]).astype(np.int8)
    Hf = H.astype(np.float64)
    return H, bool(np.array_equal(Hf @ Hf.T, 4 * v * np.eye(4 * v)))


def select(D, K):
    """per block and energy class, the K dual vectors of largest multiplicity"""
    out = []
    for d in D:
        E = d['sig'][:, 0]; keep = []
        for e in np.unique(E):
            ii = np.nonzero(E == e)[0]
            keep.append(ii[np.argsort(-d['mult'][ii], kind='stable')[:K]])
        out.append(np.sort(np.concatenate(keep)))
    return out


def dual_join(D, idx, Tsig, seed=5, max_out=2 * 10 ** 6, chunk=1 << 22):
    rng = np.random.default_rng(seed)
    wv = rng.integers(1, 2 ** 62, size=len(Tsig), dtype=np.int64)
    P = [D[k]['sig'][idx[k]] for k in range(4)]
    h = [p @ wv for p in P]
    n = [len(p) for p in P]
    HL = (h[0][:, None] + h[1][None, :]).ravel()
    o = np.argsort(HL, kind='stable'); HLs = HL[o]; del HL
    tw = np.int64(Tsig @ wv)
    out = []
    R = (h[2][:, None] + h[3][None, :]).ravel()
    for c0 in range(0, len(R), chunk):
        HR = tw - R[c0:c0 + chunk]
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        for r in np.nonzero(hi > lo)[0]:
            i2, i3 = divmod(int(c0 + r), n[3])
            for p in o[lo[r]:hi[r]]:
                i0, i1 = divmod(int(p), n[1])
                if (P[0][i0] + P[1][i1] + P[2][i2] + P[3][i3] == Tsig).all():
                    out.append((idx[0][i0], idx[1][i1], idx[2][i2], idx[3][i3]))
            if len(out) >= max_out: return out, True
    return out, False


def best_pairing(m):
    """the 2+2 split minimising the join cost |left| + |right| (left = smaller side, sorted; right streamed)"""
    best = None
    for (a, b), (c, d) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        L, R = m[a] * m[b], m[c] * m[d]
        if L > R: (a, b), (c, d), L, R = (c, d), (a, b), R, L
        if best is None or L + R < best[0]: best = (L + R, (a, b), (c, d))
    return best


def cross_join(S, Xs, seed=11, chunk=1 << 22):
    """all (a, b, c, d) with sum_k PAF_G(Xs[k][.]) = 4v delta: the smaller pair side is sorted, the larger streamed"""
    m = [len(x) for x in Xs]
    if min(m) == 0: return [], 0
    work, (a, b), (c, d) = best_pairing(m)
    rng = np.random.default_rng(seed)
    gw = rng.integers(1, 2 ** 62, size=len(S.greps), dtype=np.int64)
    h = {}
    for k in range(4):
        Pk = np.zeros((m[k], len(S.greps)), np.int64)
        for c0 in range(0, m[k], 4096): Pk[c0:c0 + 4096] = S.G.paf(Xs[k][c0:c0 + 4096].astype(np.int64))[:, S.greps]
        h[k] = Pk @ gw
    HL = (h[a][:, None] + h[b][None, :]).ravel()
    o = np.argsort(HL, kind='stable'); HLs = HL[o]; del HL
    tw = np.int64(S.T @ gw)
    hits = []
    for i0 in range(0, m[c], max(1, chunk // m[d])):
        HR = (tw - (h[c][i0:i0 + max(1, chunk // m[d])][:, None] + h[d][None, :])).ravel()
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        for r in np.nonzero(hi > lo)[0]:
            ic, id_ = divmod(int(r), m[d]); ic += i0
            for p in o[lo[r]:hi[r]]:
                ia, ib = divmod(int(p), m[b])
                sel = {a: ia, b: ib, c: ic, d: id_}
                Xc = np.array([Xs[k][sel[k]] for k in range(4)], np.int64)
                if (S.G.paf(Xc).sum(0) == S.sl.T).all(): hits.append(Xc)
    return hits, work


def cross_join_h(S, hs, row, gw, chunk=1 << 22):
    """as cross_join, from per-block hash vectors hs[k] (PAF_G at the lag classes, dotted with gw); row(k, i) gives
    the element i of block k on G for the exact re-check of every match"""
    m = [len(h) for h in hs]
    if min(m) == 0: return [], 0
    work, (a, b), (c, d) = best_pairing(m)
    HL = (hs[a][:, None] + hs[b][None, :]).ravel()
    o = np.argsort(HL, kind='stable'); HLs = HL[o]; del HL
    tw = np.int64(S.T @ gw); step = max(1, chunk // m[d]); hits = []
    for i0 in range(0, m[c], step):
        HR = (tw - (hs[c][i0:i0 + step][:, None] + hs[d][None, :])).ravel()
        lo = np.searchsorted(HLs, HR, 'left'); hi = np.searchsorted(HLs, HR, 'right')
        for r in np.nonzero(hi > lo)[0]:
            ic, id_ = divmod(int(r), m[d]); ic += i0
            for p in o[lo[r]:hi[r]]:
                ia, ib = divmod(int(p), m[b])
                sel = {a: ia, b: ib, c: ic, d: id_}
                Xc = np.array([row(k, sel[k]) for k in range(4)], np.int64)
                if (S.G.paf(Xc).sum(0) == S.sl.T).all(): hits.append(Xc)
    return hits, work


def main():
    a = sys.argv[1:]
    def opt(name, default=None):
        return a[a.index(name) + 1] if name in a else default
    q = int(opt('--q', '311')); yi = int(opt('--y', '0')); K = int(opt('--K', '300'))
    maxd = int(opt('--max-dual', '100000')); tag = opt('--tag', f'q{q}_y{yi}')
    minscore = float(opt('--min-score', '0'))
    outdir = os.path.join(ROOT, 'logs/ab/b0z3'); os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, f'run_{tag}.jsonl')
    done = set()
    if os.path.exists(logf):
        for l in open(logf):
            try: d = json.loads(l)
            except Exception: continue
            if 'dual_id' in d: done.add(d['dual_id'])
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    t0 = time.time()
    S = SliceZ3(q); sl = S.sl
    sols, _ = stage1(sl)
    Y = sols[yi]
    D = S.dual_lists(Y)
    stats = []
    for d in D:
        E = d['sig'][:, 0]
        stats.append({int(e): [int((E == e).sum()), round(float(np.log2(d['mult'][E == e].max())), 2)] for e in np.unique(E)})
    log(dict(stage="setup", q=q, v=S.v, f=S.f, y=Y.tolist(), stage1_solutions=len(sols), dual_list=[len(d['vecs']) for d in D],
             energy_classes=stats, secs=round(time.time() - t0, 1)))
    if '--stats' in a: return
    Tsig = S.Tq[S.qreps]
    idx = select(D, K)
    t1 = time.time()
    duals, capped = dual_join(D, idx, Tsig)
    score = np.array([sum(np.log2(D[k]['mult'][q4[k]]) for k in range(4)) for q4 in duals]) if duals else np.zeros(0)
    cost = np.array([np.log2(best_pairing([int(D[k]['mult'][q4[k]]) for k in range(4)])[0] + OVH) for q4 in duals]) if duals else np.zeros(0)
    eff = np.where(score >= minscore, score - cost, -np.inf)
    order = np.argsort(-eff, kind='stable')
    log(dict(stage="dual", K=K, kept=[len(i) for i in idx], dual_solutions=len(duals), capped=capped,
             log2_fibre_product_max=round(float(score.max()), 2) if len(score) else None,
             log2_fibre_product_median=round(float(np.median(score)), 2) if len(score) else None,
             top_efficiency=[round(float(eff[i]), 2) for i in order[:5]] if len(score) else None, secs=round(time.time() - t1, 1)))
    ncores = 0; nproc = 0
    gw = np.random.default_rng(11).integers(1, 2 ** 62, size=len(S.greps), dtype=np.int64)
    cache = {}
    def blockdata(k, i):
        key = (k, int(i))
        if key not in cache:
            if len(cache) > 4000: cache.clear()
            Mx = S.fibre_mats(D[k]['sigma'], D[k]['vecs'][i])
            assert len(Mx) == int(D[k]['mult'][i])
            h = np.zeros(len(Mx), np.int64)
            for c0 in range(0, len(Mx), 2048):
                h[c0:c0 + 2048] = S.G.paf(S.expand(Mx[c0:c0 + 2048], D[k]['alpha']).astype(np.int64))[:, S.greps] @ gw
            cache[key] = (Mx, h)
        return cache[key]
    for oi in order[:maxd]:
        if not np.isfinite(eff[oi]): break
        q4 = duals[oi]
        did = hashlib.sha1(np.array([D[k]['vecs'][q4[k]] for k in range(4)]).tobytes()).hexdigest()[:16]
        if did in done: continue
        t2 = time.time()
        bd = [blockdata(k, q4[k]) for k in range(4)]
        msz = [len(x[0]) for x in bd]
        row = lambda k, i: S.expand(bd[k][0][i:i + 1], D[k]['alpha'])[0]
        hits, work = cross_join_h(S, [x[1] for x in bd], row, gw)
        rec = dict(dual_id=did, rank=int(nproc), fibres=msz, log2_fibre_product=round(float(score[oi]), 2),
                   cross_pairs=int(work), cores=len(hits), secs=round(time.time() - t2, 1))
        for Xc in hits:
            ncores += 1
            name = f'{tag}_{ncores}'
            np.save(os.path.join(outdir, f'core_{name}.npy'), Xc.astype(np.int8))
            Xcyc = S.to_cyclic(Xc)
            P = np.rint(np.fft.ifft(np.abs(np.fft.fft(Xcyc.astype(float), axis=1)) ** 2, axis=1).real).astype(np.int64).sum(0)
            cyc_ok = bool(P[0] == 4 * S.v and (P[1:] == 0).all())
            with open(os.path.join(outdir, f'seq_{name}_Z{S.v}.txt'), 'w') as fh:
                for lab, x in zip('ABCD', Xcyc): fh.write(f"{lab}=" + ''.join('+' if t > 0 else '-' for t in x) + "\n")
            had = None
            if ncores == 1 or '--check-all' in a:
                _, had = gs_hadamard(Xcyc)
            log(dict(CORE=name, verified_paf_G=True, verified_paf_Zv=cyc_ok, hadamard=had, order=4 * S.v))
        log(rec)
        nproc += 1
        if ncores and '--stop-on-core' in a: break
    log(dict(stage="end", processed=nproc, cores=ncores, secs=round(time.time() - t0, 1)))


if __name__ == "__main__":
    main()
