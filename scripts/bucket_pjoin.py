"""In-bucket partial joins: measure the per-bucket solution rate of the slot k-tree directly.

For one tau-bucket of scripts/kt_join_slots*.py (same chunk classes, seeds and pair-sum statistics), enumerate the
AB pairs with PAF_A + PAF_B = tau and the CD pairs with PAF_C + PAF_D = -tau on the k chunk classes (A and C
subsampled by --frac), then count AB/CD pair-pairs that also match on the first j classes of the rest (hist class
first, as in the join) for each j in --js.  The measured bits per class against the Gaussian conditional prediction
(the model of scripts/bucket_rank.py) give the rate of the full match: E(bucket) = pairs_AB pairs_CD 2^-(bits for all
L - 1 - k independent rest classes).
usage: python3 scripts/bucket_pjoin.py V D a,b,c,d POOLS_DIR [--k 10] [--b 0 | --tau t1,..,tk] [--frac 0.125]
                                       [--js 2,4,6,8,10,12,14,16,18] [--seed 0] [--cache JOIN_OUT/cache]
  --cache: reuse the join's cached chunk-class PAFs (ch_B.npy, ch_D.npy) instead of recomputing them for B and D."""
import os, sys, json, time, heapq
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kaggle'))
import kt_join_kaggle as KJ
sys.path.insert(0, HERE)
import flatw as FW


def pairs_on_key(kA, kB):
    """all (i, j) with kA[i] == kB[j]"""
    oB = np.argsort(kB, kind='stable'); sB = kB[oB]
    lo = np.searchsorted(sB, kA, 'left'); hi = np.searchsorted(sB, kA, 'right'); n = hi - lo
    I = np.repeat(np.arange(len(kA)), n); off = np.repeat(lo - np.cumsum(n) + n, n) + np.arange(n.sum())
    return I, oB[off]


def main():
    a = sys.argv[1:]
    opt = lambda nm, df: type(df)(a[a.index(nm) + 1]) if nm in a else df
    v, d = int(a[0]), int(a[1]); prof = tuple(int(x) for x in a[2].split(',')); PD = a[3]
    cache = opt('--cache', ''); k = opt('--k', 10); frac = opt('--frac', 0.125); seed = opt('--seed', 0); bsel = opt('--b', 0)
    js = [int(x) for x in opt('--js', '2,4,6,8,10,12,14,16,18').split(',')]
    T0 = time.time(); log = lambda **kw: print(json.dumps(dict(kw, t=round(time.time() - T0, 1))), flush=True)
    R = KJ.Red(v, d); Nsp = KJ.sparse_N(R)
    pools = {}
    for s in 'ABCD':
        x = np.load(os.path.join(PD, f'slot_{s}.npy'), mmap_mode='r').view(np.ndarray)
        pools[s] = x if x.ndim == 2 else x.reshape(-1, 1)
    paf = lambda s, idx: FW.paf_w(np.ascontiguousarray(pools[s][idx]), *Nsp, R.L).astype(np.int16)
    rng_cls = np.random.default_rng(12345 + seed); rng_pd = np.random.default_rng(777 + seed)
    classes = rng_cls.permutation(R.L); chunk = classes[:k]; cc = classes[k]; rest = classes[k + 1:]
    order_rest = [cc] + list(rest)
    if '--tau' in a:
        tau = np.array([int(x) for x in opt('--tau', '').split(',')])
    else:                                        # the join's bucket b (same pair-sum statistics and best-first order)
        ns = 2_000_000
        def pairdist(s1, s2):
            i = rng_pd.integers(0, len(pools[s1]), ns); j = rng_pd.integers(0, len(pools[s2]), ns)
            S = paf(s1, i)[:, chunk] + paf(s2, j)[:, chunk]
            return [dict(zip(*np.unique(S[:, q], return_counts=True))) for q in range(k)]
        dAB = pairdist('A', 'B'); dCD = pairdist('C', 'D'); opts = []
        for q in range(k):
            vals = sorted(set(dAB[q]) & {-x for x in dCD[q]}, key=lambda t: -(dAB[q][t] * dCD[q].get(-t, 0)))
            opts.append([(t, np.log2(dAB[q][t] / ns) + np.log2(dCD[q][-t] / ns)) for t in vals[:12]])
        sc = lambda ix: sum(opts[q][ix[q]][1] for q in range(k)); start = (0,) * k; heap = [(-sc(start), start)]; seen = {start}
        for _ in range(bsel + 1):
            _, ix = heapq.heappop(heap)
            for q in range(k):
                if ix[q] + 1 < len(opts[q]):
                    nx = ix[:q] + (ix[q] + 1,) + ix[q + 1:]
                    if nx not in seen: seen.add(nx); heapq.heappush(heap, (-sc(nx), nx))
        tau = np.array([opts[q][ix[q]][0] for q in range(k)])
    rng = np.random.default_rng(5)
    # pair lists on the chunk (A and C subsampled), full PAF vectors of the participating blocks
    def side(s1, s2, tsign):
        i1 = np.sort(rng.choice(len(pools[s1]), max(1, int(frac * len(pools[s1]))), replace=False))
        P1 = paf(s1, i1)
        if cache: P2c = np.load(os.path.join(cache, f'ch_{s2}.npy')).astype(np.int16)
        else: P2c = np.concatenate([paf(s2, np.arange(s0, min(s0 + (1 << 21), len(pools[s2]))))[:, chunk]
                                    for s0 in range(0, len(pools[s2]), 1 << 21)])
        W = rng.integers(1, 2 ** 62, size=k).astype(np.uint64)
        key = lambda M: (M.astype(np.int64).astype(np.uint64) * W).sum(1)
        I, J = pairs_on_key(key(tsign * tau[None, :] - P1[:, chunk]), key(P2c)); del P2c
        uj, inv = np.unique(J, return_inverse=True); P2 = paf(s2, uj)
        U = (P1[I][:, order_rest].astype(np.int16) + P2[inv][:, order_rest].astype(np.int16))
        return U, len(i1)
    UAB, nA = side('A', 'B', 1); UCD, nC = side('C', 'D', -1)
    log(stage='bucket', tau=tau.tolist(), frac=frac, AB_pairs=len(UAB), CD_pairs=len(UCD),
        AB_pairs_full=round(len(UAB) / frac), CD_pairs_full=round(len(UCD) / frac))
    # Gaussian conditional prediction per j (pairs' own empirical covariance: the model inside this bucket)
    nr = len(order_rest) - 1                       # independent rest classes (the row sums fix one)
    mAB, mCD = UAB.mean(0), UCD.mean(0); SAB, SCD = np.cov(UAB.T.astype(float)), np.cov(UCD.T.astype(float))
    def gbits(j):
        idx = list(range(j)); S = SAB[np.ix_(idx, idx)] + SCD[np.ix_(idx, idx)]; m = mAB[idx] + mCD[idx]
        sg, ld = np.linalg.slogdet(2 * np.pi * S)
        return float(-(j * np.log(4) - 0.5 * ld - 0.5 * m @ np.linalg.solve(S, m)) / np.log(2))
    Wr = np.random.default_rng(9).integers(1, 2 ** 62, size=len(order_rest)).astype(np.uint64)
    for j in js:
        if j > nr: continue
        kab = (UAB[:, :j].astype(np.int64).astype(np.uint64) * Wr[:j]).sum(1)
        kcd = (np.uint64(0) - (UCD[:, :j].astype(np.int64).astype(np.uint64) * Wr[:j]).sum(1))
        kcd.sort(); cnt = int((np.searchsorted(kcd, kab, 'right') - np.searchsorted(kcd, kab, 'left')).sum())
        tot = np.log2(len(UAB)) + np.log2(len(UCD))
        log(stage='pjoin', j=j, matches=cnt, bits_measured=round(float(tot - np.log2(cnt)), 2) if cnt else None,
            bits_gaussian=round(gbits(j), 2))
    full = gbits(nr)
    log(stage='full', classes=nr, bits_gaussian=round(full, 2), per_class=round(full / nr, 3),
        log2E_bucket_gauss=round(float(np.log2(len(UAB) / frac) + np.log2(len(UCD) / frac) - full), 2))


if __name__ == '__main__':
    main()
