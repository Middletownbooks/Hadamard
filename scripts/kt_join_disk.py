"""Disk-capable slot k-tree join (pools of 2^26-2^28 blocks per slot), same join as scripts/kt_join_slots*.py.

Differences, all for memory:
- Prep (once, restart-safe, files under OUT/prep, memory-mapped in the join): per slot the k chunk-class PAFs (int8),
  the histogram class (int16), the 64-bit hash of the rest classes and the pool index (int64), all permuted into key
  order: A and C by partner key (tau - PAF, the order is the same for every tau [P]: each key field is
  (tau_q - PAF_q - v%4)/4 + offset with no wrap, monotone in -PAF_q), B and D by own key.  Slots that hold the same
  pool and the same order share the files (223: B and D).
- Keys are computed inside the merge kernels from the int8 chunk PAFs (no per-bucket key or sort arrays); fields of
  6 bits (k <= 10) or 5 bits (k <= 12; checked against the pools' PAF range).
- Stage 1 (pairs matching on the chunk, histogram class kept) is split into passes by histogram value under --mem,
  then the bloom match and exact verification of kaggle/kt_join_kaggle.py.
Buckets: best first by the separable score (the expected-count order, scripts/bucket_rank.py, scripts/bucket_pjoin.py).
Restart-safe: finished buckets (log) are skipped.
usage: python3 scripts/kt_join_disk.py V D a,b,c,d POOLS_DIR OUT [--k 12] [--mem 3e9] [--hours 100] [--only-tau t1,..,tk]"""
import os, sys, json, time, heapq
import numpy as np, numba as nb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kaggle'))
import kt_join_kaggle as KJ
sys.path.insert(0, HERE)
from pmatch import pmatch                          # exact radix-partitioned key match (faster than the byte bloom)
sys.path.insert(0, HERE)
import flatw as FW


@nb.njit(inline='always')
def key_at(ch, n, tau, sgn, r4, w, off):
    K = ch.shape[1]; key = np.uint64(0); m = (1 << w) - 1
    for i in range(K):
        val = (np.int64(tau[i]) + sgn * np.int64(ch[n, i]) - r4) // 4 + off
        key |= np.uint64(val & m) << np.uint64(w * (K - 1 - i))     # class 0 most significant
    return key


@nb.njit(parallel=True, cache=True)
def keys_all(ch, tau, sgn, r4, w, off):
    out = np.empty(ch.shape[0], np.uint64)
    for n in nb.prange(ch.shape[0]): out[n] = key_at(ch, n, tau, sgn, r4, w, off)
    return out


@nb.njit(cache=True)
def lower(chY, tauY, key, r4, w, off):
    lo = 0; hi = chY.shape[0]
    while lo < hi:
        mid = (lo + hi) // 2
        if key_at(chY, mid, tauY, 1, r4, w, off) < key: lo = mid + 1
        else: hi = mid
    return lo


@nb.njit(parallel=True, cache=True)
def stage(chX, cX, hX, chY, cY, hY, tau, r4, w, off, T, mask, mode, offs, K, I, J, G, outsign, nchunk):
    NX = chX.shape[0]; NY = chY.shape[0]; zero = np.zeros(tau.shape[0], np.int64)
    H = np.zeros((nchunk, 2 * T + 1), np.int64); cnt = np.zeros(nchunk, np.int64)
    b0 = np.linspace(0, NX, nchunk + 1).astype(np.int64)
    for b in nb.prange(nchunk):
        lo = b0[b]; hi = b0[b + 1]
        if lo >= hi: continue
        # do not split a run of equal X keys across chunks: start after the run that crosses lo (except chunk 0)
        if b > 0:
            kprev = key_at(chX, lo - 1, tau, -1, r4, w, off)
            while lo < hi and key_at(chX, lo, tau, -1, r4, w, off) == kprev: lo += 1
        if lo >= hi: continue
        while hi < NX and key_at(chX, hi, tau, -1, r4, w, off) == key_at(chX, hi - 1, tau, -1, r4, w, off): hi += 1
        i = lo; j = lower(chY, zero, key_at(chX, lo, tau, -1, r4, w, off), r4, w, off); n = offs[b] if mode == 2 else 0
        c = 0
        while i < hi and j < NY:
            a = key_at(chX, i, tau, -1, r4, w, off); bb = key_at(chY, j, zero, 1, r4, w, off)
            if a < bb: i += 1
            elif a > bb: j += 1
            else:
                j2 = j
                while j2 < NY and key_at(chY, j2, zero, 1, r4, w, off) == a: j2 += 1
                i2 = i
                while i2 < hi and key_at(chX, i2, tau, -1, r4, w, off) == a: i2 += 1
                for ii in range(i, i2):
                    cx = cX[ii]
                    for jj in range(j, j2):
                        t = cx + cY[jj]
                        if mode == 0: H[b, t + T] += 1
                        elif mask[t + T]:
                            if mode == 1: c += 1
                            else:
                                kk = hX[ii] + hY[jj] + np.uint64(np.int64(t)) * G
                                if outsign < 0: kk = np.uint64(0) - kk
                                K[n] = kk; I[n] = ii; J[n] = jj; n += 1
                i = i2; j = j2
        cnt[b] = c
    return H.sum(0), cnt


def main():
    a = sys.argv[1:]
    opt = lambda nm, df: type(df)(a[a.index(nm) + 1]) if nm in a else df
    v, d = int(a[0]), int(a[1]); prof = tuple(int(x) for x in a[2].split(',')); PD, OUT = a[3], a[4]
    k = opt('--k', 12); mem = opt('--mem', 3e9); hours = opt('--hours', 100.0); seed = opt('--seed', 0)
    os.makedirs(os.path.join(OUT, 'prep'), exist_ok=True); T0 = time.time()
    logf = open(os.path.join(OUT, 'log.jsonl'), 'a')
    def log(**kw):
        kw['t'] = round(time.time() - T0, 1); s = json.dumps(kw); print(s, flush=True); logf.write(s + '\n'); logf.flush()
    R = KJ.Red(v, d); Nsp = KJ.sparse_N(R); r4 = v % 4
    w = 6 if k <= 10 else 5; off = 1 << (w - 1); assert k <= 12
    raw = {}
    for s in 'ABCD':
        x = np.load(os.path.join(PD, f'slot_{s}.npy'), mmap_mode='r').view(np.ndarray)
        raw[s] = x if x.ndim == 2 else x.reshape(-1, 1)
    rng_cls = np.random.default_rng(12345 + seed); rng_pd = np.random.default_rng(777 + seed)
    classes = rng_cls.permutation(R.L); chunk = classes[:k]; cc = classes[k]; rest = classes[k + 1:]
    W = rng_cls.integers(1, 2 ** 63, size=len(rest), dtype=np.uint64) | np.uint64(1); G = np.uint64(0x9E3779B97F4A7C15)
    log(stage='config', v=v, d=d, profile=list(prof), k=k, field_bits=w, mem=mem,
        pools={s: len(x) for s, x in raw.items()}, log2_tuples=round(float(sum(np.log2(len(x)) for x in raw.values())), 2))
    # ---- prep: per (pool file, order) one set of permuted arrays
    pp = lambda nm: os.path.join(OUT, 'prep', nm + '.npy')
    src = {s: os.path.realpath(os.path.join(PD, f'slot_{s}.npy')) for s in 'ABCD'}
    order_of = {'A': -1, 'B': 1, 'C': -1, 'D': 1}
    tag = {}
    for s in 'ABCD':                                     # share files between slots with the same pool and order
        same = [t for t in tag if src[t] == src[s] and order_of[t] == order_of[s]]
        tag[s] = tag[same[0]] if same else s
    for s in 'ABCD':
        if tag[s] != s or os.path.exists(pp(f'idx_{s}')): continue
        x = raw[s]; N = len(x)
        ch = np.empty((N, k), np.int8); c_ = np.empty(N, np.int16); h = np.empty(N, np.uint64); amax = 0
        for s0 in range(0, N, 1 << 21):
            P = FW.paf_w(np.ascontiguousarray(x[s0:s0 + (1 << 21)]), *Nsp, R.L)
            ch[s0:s0 + len(P)] = P[:, chunk]; c_[s0:s0 + len(P)] = P[:, cc]; amax = max(amax, int(np.abs(P[:, chunk]).max()))
            h[s0:s0 + len(P)] = (P[:, rest].astype(np.int64).astype(np.uint64) * W[None, :]).sum(1)
        kk = keys_all(ch, np.zeros(k, np.int64), order_of[s], r4, w, off); o = np.argsort(kk, kind='stable'); del kk
        np.save(pp(f'ch_{s}'), ch[o]); np.save(pp(f'c_{s}'), c_[o]); np.save(pp(f'h_{s}'), h[o]); np.save(pp(f'idx_{s}'), o)
        log(stage='prep', slot=s, size=N, max_abs_paf=amax); del ch, c_, h, o
    D = {s: {nm: np.load(pp(f'{nm}_{tag[s]}'), mmap_mode='r').view(np.ndarray) for nm in ('ch', 'c', 'h', 'idx')} for s in 'ABCD'}
    # ---- bucket order (separable score, same statistics as the slot joins)
    ns = 2_000_000
    def pairdist(s1, s2):
        i = rng_pd.integers(0, len(raw[s1]), ns); j = rng_pd.integers(0, len(raw[s2]), ns)
        S = (FW.paf_w(np.ascontiguousarray(raw[s1][i]), *Nsp, R.L)[:, chunk].astype(np.int16)
             + FW.paf_w(np.ascontiguousarray(raw[s2][j]), *Nsp, R.L)[:, chunk].astype(np.int16))
        return [dict(zip(*np.unique(S[:, q], return_counts=True))) for q in range(k)]
    dAB = pairdist('A', 'B'); dCD = pairdist('C', 'D'); opts = []
    for q in range(k):
        vals = sorted(set(dAB[q]) & {-x for x in dCD[q]}, key=lambda t: -(dAB[q][t] * dCD[q].get(-t, 0)))
        opts.append([(int(t), float(np.log2(dAB[q][t] / ns) + np.log2(dCD[q][-t] / ns))) for t in vals[:12]])
    lim = (1 << (w - 1)) - 1
    amax = max(int(np.abs(D[s]['ch'][::97]).max()) for s in 'ABCD'); tmax = max(abs(t) for o in opts for t, _ in o)
    assert (amax + tmax + 3) // 4 < lim, f'key field overflow: |PAF| {amax} + |tau| {tmax} needs more than {w} bits'
    score = lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
    start = (0,) * k; heap = [(-score(start), start)]; seen = {start}; buckets = []
    while heap and len(buckets) < 200000:
        s_, ix = heapq.heappop(heap); buckets.append(ix)
        for q in range(k):
            if ix[q] + 1 < len(opts[q]):
                nx = ix[:q] + (ix[q] + 1,) + ix[q + 1:]
                if nx not in seen: seen.add(nx); heapq.heappush(heap, (-score(nx), nx))
    done = set()
    for l in open(os.path.join(OUT, 'log.jsonl')):
        try:
            z = json.loads(l)
            if z.get('stage') == 'bucket': done.add(z['b'])
        except Exception: pass
    if '--skip-log' in a:                         # skip tau-buckets already searched by another run (same chunk classes)
        skip = set()
        for l in open(a[a.index('--skip-log') + 1]):
            try:
                z = json.loads(l)
                if z.get('stage') == 'bucket': skip.add(tuple(z['tau']))
            except Exception: pass
        done |= {b for b, ix in enumerate(buckets) if tuple(int(opts[q][ix[q]][0]) for q in range(k)) in skip}
    if '--only-tau' in a:
        t_only = [int(x) for x in a[a.index('--only-tau') + 1].split(',')]
        opts = [[(t_only[q], 0.0)] for q in range(k)]; buckets = [(0,) * k]; done = set()
    log(stage='ready', buckets=len(buckets), done_before=len(done), chunk_classes=[int(c) for c in chunk])
    solf = open(os.path.join(OUT, f'solutions_v{v}.txt'), 'a'); nsol = 0; T = 300; NCH = 256
    A, B, C, Dd = D['A'], D['B'], D['C'], D['D']
    dummyK = np.zeros(1, np.uint64); dummyI = np.zeros(1, np.int64); m0 = np.zeros(2 * T + 1, np.bool_)
    for b, ix in enumerate(buckets):
        if b in done: continue
        if time.time() - T0 > hours * 3600: log(stage='time_limit'); break
        tau = np.array([opts[q][ix[q]][0] for q in range(k)], np.int64)
        hAB, _ = stage(A['ch'], A['c'], A['h'], B['ch'], B['c'], B['h'], tau, r4, w, off, T, m0, 0, np.zeros(NCH + 1, np.int64), dummyK, dummyI, dummyI, G, 1, NCH)
        hCD, _ = stage(C['ch'], C['c'], C['h'], Dd['ch'], Dd['c'], Dd['h'], -tau, r4, w, off, T, m0, 0, np.zeros(NCH + 1, np.int64), dummyK, dummyI, dummyI, G, -1, NCH)
        both = hAB * hCD[::-1]
        budget = mem / 24; passes = []; cur = []; load = 0
        for t in np.arange(-T, T + 1)[np.argsort(-both)]:
            if both[t + T] == 0: continue
            wt = hAB[t + T] + hCD[T - t]
            if load + wt > budget and cur: passes.append(cur); cur = []; load = 0
            cur.append(int(t)); load += wt
        if cur: passes.append(cur)
        found = 0; ncand = 0
        for pl in passes:
            mA = np.zeros(2 * T + 1, np.bool_); mC = np.zeros(2 * T + 1, np.bool_)
            for t in pl: mA[t + T] = True; mC[T - t] = True
            out = []
            for X, Y, tt, mk, sg in ((A, B, tau, mA, 1), (C, Dd, -tau, mC, -1)):
                _, cnt = stage(X['ch'], X['c'], X['h'], Y['ch'], Y['c'], Y['h'], tt, r4, w, off, T, mk, 1, np.zeros(NCH + 1, np.int64), dummyK, dummyI, dummyI, G, sg, NCH)
                offs = np.zeros(NCH + 1, np.int64); offs[1:] = np.cumsum(cnt)
                K = np.empty(offs[-1], np.uint64); I = np.empty(offs[-1], np.int64); J = np.empty(offs[-1], np.int64)
                stage(X['ch'], X['c'], X['h'], Y['ch'], Y['c'], Y['h'], tt, r4, w, off, T, mk, 2, offs, K, I, J, G, sg, NCH)
                out.append((K, I, J))
            (KA, IA, JA), (KC, IC, JC) = out
            ia, ic = pmatch(KA, KC, 10); ncand += len(ia)
            for q, qc in zip(ia, ic):
                bits = np.array([raw['A'][A['idx'][IA[q]]], raw['B'][B['idx'][JA[q]]], raw['C'][C['idx'][IC[qc]]], raw['D'][Dd['idx'][JC[qc]]]], np.uint64)
                X4 = R.expand(FW.unpack_w(bits.reshape(4, -1), R.r))
                if np.all(R.paf(X4).sum(0)[1:] == 0):
                    found += 1; nsol += 1
                    rec = dict(v=v, d=d, profile=list(prof), bucket=b, blocks={kk: ''.join('+' if t > 0 else '-' for t in x) for kk, x in zip('ABCD', X4)})
                    solf.write(json.dumps(rec) + '\n'); solf.flush(); log(stage='SOLUTION', b=b)
            del out, KA, IA, JA, KC, IC, JC
        log(stage='bucket', b=b, tau=tau.tolist(), AB=int(hAB.sum()), CD=int(hCD.sum()), passes=len(passes), cand=int(ncand),
            solutions=found, total=nsol)
    log(stage='end', solutions=nsol)


if __name__ == '__main__':
    main()
