"""Where do real solutions rank in the best-first bucket order of the slot k-tree (scripts/kt_join_slots*.py)?

The join visits tau-buckets (tau = PAF_A + PAF_B on k chunk classes, PAF_C + PAF_D = -tau) best first by the separable
score  s1(tau) = sum_q log2 P_AB,q(tau_q) + log2 P_CD,q(-tau_q)  (empirical one-class pair-sum frequencies, top 12
values per class).  s1 is the log of the pair count N_AB(tau) N_CD(-tau) only if the chunk classes are independent
of each other and of the rest; the expected number of solution tuples in the bucket is

    E(tau) = N_AB(tau) N_CD(-tau) P(rest of AB + rest of CD = 0 | chunk sums tau, -tau),

and under a Gaussian model of the class-PAF vectors the conditional rest probability is a Gaussian density in
delta(tau) = E[u_rest | u_c = tau] + E[w_rest | w_c = -tau] (u = PAF_A + PAF_B, w = PAF_C + PAF_D):
    log2 P_rest(tau) = c0 - (1/2) delta^T S^-1 delta / ln 2,   S = Cov(u_r | u_c) + Cov(w_r | w_c).
The corrected score s2 = s1 + log2 P_rest <= s1 + c0 is enumerated exactly best first with s1 as the bound.

For a known core (natural coordinates, 4 x v) every tuple in the pools' symmetry class is formed: all e rotations
x(i) -> x(g^s i) and every assignment of its blocks to slots with the profile's row sums (blocks normalised to positive
sum).  Prints, per tuple, its tau, whether tau is reachable (each value among the class's top 12), its rank under s1
and under s2, and log2 E(tau); then the predicted share of all solution tuples captured by the first N buckets under
both orders (from the Gaussian E(tau) summed over the enumerated buckets against the Gaussian total).

usage: python3 scripts/bucket_rank.py V D a,b,c,d POOLS_DIR CORE.npy [--k 8] [--seed 0] [--nb 20000] [--out J]
  POOLS_DIR holds slot_A..slot_D.npy (1-D uint64 or (N, W) masks, as for the joins); CORE.npy a (4, v) +-1 array."""
import os, sys, json, heapq, itertools, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kaggle'))
import kt_join_kaggle as KJ
sys.path.insert(0, HERE)
import flatw as FW


def class_pafs(R, X):
    """lag-class PAFs (one representative lag per class) of +-1 rows X (n, v), by FFT"""
    return R.paf(X)[:, R.lags]


def main():
    a = sys.argv[1:]
    v, d = int(a[0]), int(a[1]); prof = tuple(int(x) for x in a[2].split(',')); PD, CORE = a[3], a[4]
    k = int(a[a.index('--k') + 1]) if '--k' in a else 8
    seed = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    NB = int(a[a.index('--nb') + 1]) if '--nb' in a else 20000
    outp = a[a.index('--out') + 1] if '--out' in a else None
    outf = open(outp, 'w') if outp else None
    T0 = time.time()
    def log(**kw):
        kw['t'] = round(time.time() - T0, 1); s = json.dumps(kw); print(s, flush=True)
        if outf: outf.write(s + '\n'); outf.flush()
    R = KJ.Red(v, d); Nsp = KJ.sparse_N(R); e = (v - 1) // d; g = KJ.primroot(v)
    pools = {}
    for s in 'ABCD':
        x = np.load(os.path.join(PD, f'slot_{s}.npy'), mmap_mode='r').view(np.ndarray)
        pools[s] = x if x.ndim == 2 else x.reshape(-1, 1)
    paf = lambda s, idx: FW.paf_w(np.ascontiguousarray(pools[s][np.sort(idx)]), *Nsp, R.L).astype(np.int64)
    # the join's classes and its pair-sum statistics (same seeds and sample sizes as scripts/kt_join_slots*.py)
    rng_cls = np.random.default_rng(12345 + seed); rng_pd = np.random.default_rng(777 + seed)
    classes = rng_cls.permutation(R.L); chunk = classes[:k]; cc = classes[k]; rest = classes[k + 1:]
    ns = 2_000_000
    def pairdist(s1, s2):
        i = rng_pd.integers(0, len(pools[s1]), ns); j = rng_pd.integers(0, len(pools[s2]), ns)
        Pi = FW.paf_w(np.ascontiguousarray(pools[s1][i]), *Nsp, R.L)[:, chunk].astype(np.int16)
        Pj = FW.paf_w(np.ascontiguousarray(pools[s2][j]), *Nsp, R.L)[:, chunk].astype(np.int16)
        S = Pi + Pj
        return [dict(zip(*np.unique(S[:, q], return_counts=True))) for q in range(k)]
    dAB = pairdist('A', 'B'); dCD = pairdist('C', 'D'); opts = []
    for q in range(k):
        vals = sorted(set(dAB[q]) & {-x for x in dCD[q]}, key=lambda t: -(dAB[q][t] * dCD[q].get(-t, 0)))
        opts.append([(int(t), float(np.log2(dAB[q][t] / ns) + np.log2(dCD[q][-t] / ns))) for t in vals[:12]])
    optval = [{t: i for i, (t, _) in enumerate(o)} for o in opts]
    # Gaussian model on the independent classes: chunk + rest minus one (the row sum fixes sum_c PAF(c))
    rng = np.random.default_rng(99)
    mu, Cv = {}, {}
    ind = list(chunk) + [cc] + list(rest[:-1])
    for s in 'ABCD':
        P = paf(s, rng.choice(len(pools[s]), min(200000, len(pools[s])), replace=False))[:, ind].astype(float)
        mu[s] = P.mean(0); Cv[s] = np.cov(P.T)
    mAB, SAB = mu['A'] + mu['B'], Cv['A'] + Cv['B']; mCD, SCD = mu['C'] + mu['D'], Cv['C'] + Cv['D']
    c_ = np.arange(k); r_ = np.arange(k, len(ind))
    def cond(m, S):
        Kc = S[np.ix_(r_, c_)] @ np.linalg.inv(S[np.ix_(c_, c_)])
        return Kc, m[r_] - Kc @ m[c_], S[np.ix_(r_, r_)] - Kc @ S[np.ix_(c_, r_)]
    KAB, bAB, SrAB = cond(mAB, SAB); KCD, bCD, SrCD = cond(mCD, SCD)
    Sr = SrAB + SrCD; Sri = np.linalg.inv(Sr); nr = len(r_)
    sg, ld = np.linalg.slogdet(2 * np.pi * Sr)
    c0 = (nr * np.log(4) - 0.5 * ld) / np.log(2)                # log2 P_rest at delta = 0 (lattice spacing 4)
    def delta(tau):
        tau = np.asarray(tau, float); return (KAB @ tau + bAB) + (KCD @ (-tau) + bCD)
    def corr(tau):
        dl = delta(tau); return float(-0.5 * dl @ Sri @ dl / np.log(2))
    # Gaussian total of solution tuples (all independent classes), for the captured share
    Sz = SAB + SCD; mz = mAB + mCD; sg, ldz = np.linalg.slogdet(2 * np.pi * Sz)
    log2Pall = float((len(ind) * np.log(4) - 0.5 * ldz - 0.5 * mz @ np.linalg.solve(Sz, mz)) / np.log(2))
    log2L4 = float(sum(np.log2(len(pools[s])) for s in 'ABCD'))
    log(stage='model', v=v, d=d, profile=list(prof), k=k, chunk=[int(x) for x in chunk], hist=int(cc),
        n_prime_gauss=round(-log2Pall, 2), log2_tuples=round(log2L4, 2), log2_solution_tuples=round(log2L4 + log2Pall, 2),
        c0_rest=round(float(c0), 2), opts=[[t for t, _ in o] for o in opts])
    log2E = lambda tau, s1: float(log2L4 + s1 + c0 + corr(tau))    # N_AB N_CD = L^4 2^s1 (pairs with hist class free)
    # exact best-first enumeration under s1 and under s2 = s1 + corr (s1 is an upper bound of s2 - c0)
    def enum_s1(n):
        start = (0,) * k; sc = lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
        heap = [(-sc(start), start)]; seen = {start}; out = []
        while heap and len(out) < n:
            s_, ix = heapq.heappop(heap); out.append((ix, -s_))
            for q in range(k):
                if ix[q] + 1 < len(opts[q]):
                    nx = ix[:q] + (ix[q] + 1,) + ix[q + 1:]
                    if nx not in seen: seen.add(nx); heapq.heappush(heap, (-sc(nx), nx))
        return out
    tauof = lambda ix: [opts[q][ix[q]][0] for q in range(k)]
    E1 = enum_s1(NB * 10)                                        # s1 order; s2 order is re-sorted with the s1 bound
    rank1 = {ix: i for i, (ix, _) in enumerate(E1[:NB])}
    sec = []; order2 = []; ptr = 0
    while len(order2) < NB and ptr < len(E1):
        ix, s1 = E1[ptr]; ptr += 1; heapq.heappush(sec, (-(s1 + corr(tauof(ix))), ix, s1))
        bound = E1[ptr][1] if ptr < len(E1) else -1e300
        while sec and -sec[0][0] >= bound and len(order2) < NB:
            s2, ix2, s1b = heapq.heappop(sec); order2.append((ix2, -s2, s1b))
    exact2 = ptr < len(E1)
    rank2 = {ix: i for i, (ix, _, _) in enumerate(order2)}
    # s1 rank of any tau by meet in the middle over the 12^k product
    h1, h2 = list(range(k // 2)), list(range(k // 2, k))
    S1 = np.sort([sum(opts[q][i][1] for q, i in zip(h1, ix)) for ix in itertools.product(*[range(len(opts[q])) for q in h1])])
    S2 = np.array([sum(opts[q][i][1] for q, i in zip(h2, ix)) for ix in itertools.product(*[range(len(opts[q])) for q in h2])])
    def rank_s1(s):
        return int((len(S1) - np.searchsorted(S1, s - S2, 'right')).sum())
    # the core's tuples
    X0 = np.load(CORE).astype(np.int64)
    X0 = np.array([x if x.sum() > 0 else -x for x in X0])
    sums = [int(x.sum()) for x in X0]
    assigns = sorted({p for p in itertools.permutations(range(4)) if tuple(sums[i] for i in p) == prof})
    ii = np.arange(v); recs = []
    for s in range(e):
        gs = pow(g, s, v); Xs = X0[:, (gs * ii) % v]; P = class_pafs(R, Xs)
        assert np.all(P.sum(0) == 0), 'not a solution'
        for p in assigns:
            u = P[p[0]] + P[p[1]]; tau = [int(t) for t in u[chunk]]
            ix = tuple(optval[q].get(tau[q], -1) for q in range(k)); reach = all(i >= 0 for i in ix)
            s1 = sum(opts[q][ix[q]][1] for q in range(k)) if reach else None
            rec = dict(stage='tuple', rot=s, assign=list(p), tau=tau, hist=int(u[cc]), reachable=reach,
                       rank_s1=(rank1.get(ix, rank_s1(s1)) if reach else None), rank_s2=(rank2.get(ix, None) if reach else None),
                       log2E=(round(log2E(tau, s1), 2) if reach else None), corr=round(corr(tau), 2))
            recs.append(rec); log(**rec)
    r1 = np.array([r['rank_s1'] for r in recs if r['reachable']]); r2 = [r['rank_s2'] for r in recs if r['reachable']]
    r2a = np.array([x if x is not None else NB for x in r2])
    tot = 2.0 ** (log2L4 + log2Pall)
    cum1 = np.cumsum([2.0 ** log2E(tauof(ix), s1) for ix, s1 in E1[:NB]]) / tot
    cum2 = np.cumsum([2.0 ** (log2L4 + s2 + c0) for ix, s2, s1 in order2]) / tot
    marks = [n for n in (10, 100, 1000, 10000, 100000, 1000000) if n <= NB]
    log(stage='summary', tuples=len(recs), reachable=int(len(r1)),
        rank_s1_quantiles=[int(x) for x in np.quantile(r1, [0, .1, .25, .5])] if len(r1) else None,
        rank_s2_quantiles=[int(x) for x in np.quantile(r2a, [0, .1, .25, .5])] if len(r2a) else None,
        core_tuples_in_first={str(n): [int((r1 < n).sum()), int((r2a < n).sum())] for n in marks},
        predicted_share_first={str(n): [float('%.3g' % cum1[n - 1]), float('%.3g' % cum2[n - 1])] for n in marks},
        log2E_first_buckets=[round(log2E(tauof(ix), s1), 2) for ix, s1 in E1[:5]],
        log2E_first_buckets_s2=[round(log2L4 + s2 + c0, 2) for ix, s2, s1 in order2[:5]], s2_exact=exact2)


if __name__ == '__main__':
    main()
