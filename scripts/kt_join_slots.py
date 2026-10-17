"""Flat-pool k-tree join with one pool per slot (A, B, C, D), for rotation-closed pools with symmetry breaking.

Same join as kaggle/kt_join_kaggle.py (best-first tau-buckets on k chunk classes, a histogram class, a 64-bit hash
of the rest, bloom match, exact verification), but each slot reads its own pool, prepared by
scripts/build_slot_pools.py: slot A = one block per multiplier orbit (canonical rotation), slots B, C, D = uniform
samples of the rotation closures.  Every solution has a simultaneous rotation with A canonical [P], so nothing is
lost, and each hit is a distinct core up to that rotation.
Restart-safe: prep arrays cached under OUT/cache, finished buckets (log) skipped.
Keys: up to 8 chunk classes at 8 bits each (the kaggle key), or up to 10 at 6 bits (--k 9 or 10): every chunk-class
value tau_q + sgn PAF_q is = v (mod 4) [both PAF and tau = PAF_A + PAF_B have fixed residues], so (val - v % 4) / 4
fits in 6 bits for |val| <= 125, i.e. for flat blocks (|PAF| <= sqrt(F0)) at any rho used here.
usage: python3 scripts/kt_join_slots.py V D a,b,c,d POOLS_DIR OUT [--k 8] [--mem 3e9] [--bits 3.1] [--hours 100]
                                    [--only-tau t1,..,tk]   (planted tests: one bucket) [--skip-log LOG] [--no-fast]
  fast mode (default): when the previous bucket's pairs fit one pass with 30% room, the histogram merge is skipped and
  the bucket runs one pass over all histogram values (exp_hits is then logged as 0)."""
import os, sys, json, time, glob, heapq
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kaggle'))
import kt_join_kaggle as KJ
import numba as nb


@nb.njit(parallel=True, cache=True)
def keyof6(ch, tau, sgn, r4):
    """packed key of (tau + sgn ch) on the K <= 10 chunk classes, 6 bits per class after removing the residue mod 4"""
    N = ch.shape[0]; K = ch.shape[1]; out = np.empty(N, np.uint64)
    for n in nb.prange(N):
        w = np.uint64(0)
        for i in range(K):
            val = (np.int64(tau[i]) + sgn * np.int64(ch[n, i]) - r4) // 4 + 32
            w |= np.uint64(val & 63) << np.uint64(6 * i)
        out[n] = w
    return out


def main():
    a = sys.argv[1:]
    v, d = int(a[0]), int(a[1]); prof = tuple(int(x) for x in a[2].split(',')); PD, OUT = a[3], a[4]
    k = int(a[a.index('--k') + 1]) if '--k' in a else 8
    mem = float(a[a.index('--mem') + 1]) if '--mem' in a else 3e9
    bpc = float(a[a.index('--bits') + 1]) if '--bits' in a else 3.1
    hours = float(a[a.index('--hours') + 1]) if '--hours' in a else 100.0
    seed = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    os.makedirs(os.path.join(OUT, 'cache'), exist_ok=True); T0 = time.time()
    logf = open(os.path.join(OUT, 'log.jsonl'), 'a')
    def log(**kw):
        kw['t'] = round(time.time() - T0, 1); s = json.dumps(kw); print(s, flush=True); logf.write(s + '\n'); logf.flush()
    R = KJ.Red(v, d); Nsp = KJ.sparse_N(R)
    pools = {s: np.load(os.path.join(PD, f'slot_{s}.npy'), mmap_mode='r').view(np.ndarray) for s in 'ABCD'}
    log(stage='config', v=v, d=d, profile=list(prof), k=k, mem=mem, bits_per_class=bpc,
        pools={s: len(x) for s, x in pools.items()}, log2_tuples=round(float(sum(np.log2(len(x)) for x in pools.values())), 2))
    rng_cls = np.random.default_rng(12345 + seed); rng_pd = np.random.default_rng(777 + seed)
    classes = rng_cls.permutation(R.L); chunk = classes[:k]; cc = classes[k]; rest = classes[k + 1:]
    W = rng_cls.integers(1, 2 ** 63, size=len(rest), dtype=np.uint64) | np.uint64(1); G = np.uint64(0x9E3779B97F4A7C15)
    cp = lambda nm: os.path.join(OUT, 'cache', nm + '.npy')
    for s, x in pools.items():
        if os.path.exists(cp(f'h_{s}')): continue
        ch = np.empty((len(x), k), np.int8); c_ = np.empty(len(x), np.int16); h = np.empty(len(x), np.uint64)
        for s0 in range(0, len(x), 1 << 21):
            P = KJ.paf_from_bits(np.ascontiguousarray(x[s0:s0 + (1 << 21)]), *Nsp, R.L)
            ch[s0:s0 + len(P)] = P[:, chunk]; c_[s0:s0 + len(P)] = P[:, cc]
            h[s0:s0 + len(P)] = (P[:, rest].astype(np.int64).astype(np.uint64) * W[None, :]).sum(1)
        np.save(cp(f'ch_{s}'), ch); np.save(cp(f'c_{s}'), c_); np.save(cp(f'h_{s}'), h); del ch, c_, h
        log(stage='prep', slot=s)
    D = {s: dict(ch=np.load(cp(f'ch_{s}'), mmap_mode='r').view(np.ndarray), c=np.load(cp(f'c_{s}')),
                 h=np.load(cp(f'h_{s}'), mmap_mode='r').view(np.ndarray)) for s in 'ABCD'}
    if k <= 8:
        KE = k
        own = lambda ch: KJ.keyof_nb(np.ascontiguousarray(ch[:, :KE]), np.zeros(KE, np.int64), 1)
        partner = lambda ch, tau: KJ.keyof_nb(np.ascontiguousarray(ch[:, :KE]), np.asarray(tau, np.int64)[:KE], -1)
    else:
        KE = min(k, 10); r4 = v % 4
        own = lambda ch: keyof6(np.ascontiguousarray(ch[:, :KE]), np.zeros(KE, np.int64), 1, r4)
        partner = lambda ch, tau: keyof6(np.ascontiguousarray(ch[:, :KE]), np.asarray(tau, np.int64)[:KE], -1, r4)
    sB = own(D['B']['ch']); oB = np.argsort(sB, kind='stable'); sB = sB[oB]
    sD = own(D['D']['ch']); oD = np.argsort(sD, kind='stable'); sD = sD[oD]
    presort = k > 8                               # 6-bit key: every field of tau - ch is in range, so the partner-key order
    if presort:                                   # does not depend on tau: sort A and C once, not once per bucket
        oA0 = np.argsort(partner(D['A']['ch'], np.zeros(k, np.int64)), kind='stable'); chA0 = np.ascontiguousarray(D['A']['ch'][oA0])
        oC0 = np.argsort(partner(D['C']['ch'], np.zeros(k, np.int64)), kind='stable'); chC0 = np.ascontiguousarray(D['C']['ch'][oC0])
    ns = 2_000_000
    def pairdist(s1, s2):
        i = rng_pd.integers(0, len(pools[s1]), ns); j = rng_pd.integers(0, len(pools[s2]), ns)
        S = D[s1]['ch'][i].astype(np.int16) + D[s2]['ch'][j].astype(np.int16)
        return [dict(zip(*np.unique(S[:, q], return_counts=True))) for q in range(k)]
    dAB = pairdist('A', 'B'); dCD = pairdist('C', 'D'); opts = []
    for q in range(k):
        vals = sorted(set(dAB[q]) & {-x for x in dCD[q]}, key=lambda t: -(dAB[q][t] * dCD[q].get(-t, 0)))
        opts.append([(t, np.log2(dAB[q][t] / ns) + np.log2(dCD[q][-t] / ns)) for t in vals[:12]])
    score = lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
    start = tuple(0 for _ in range(k)); heap = [(-score(start), start)]; seen = {start}; buckets = []
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
    if '--only-tau' in a:                        # planted tests: process exactly one bucket (tau = PAF_A + PAF_B on the chunk classes)
        t_only = [int(x) for x in a[a.index('--only-tau') + 1].split(',')]
        opts = [[(t_only[q], 0.0)] for q in range(k)]; buckets = [tuple(0 for _ in range(k))]; done = set()
    log(stage='ready', buckets=len(buckets), done_before=len(done), chunk_classes=[int(c) for c in chunk], key_classes=KE)
    solf = open(os.path.join(OUT, f'solutions_v{v}.txt'), 'a'); nsol = 0; exp_cum = 0.0; T = 300
    fast = '--no-fast' not in a; prev_pairs = 0; npairs = [0, 0]
    A, B, Cc, Dd = D['A'], D['B'], D['C'], D['D']
    for b, ix in enumerate(buckets):
        if b in done: continue
        if time.time() - T0 > hours * 3600: log(stage='time_limit'); break
        tau = np.array([opts[q][ix[q]][0] for q in range(k)], np.int64)
        if presort:
            nA = partner(chA0, tau); oA = oA0; nC = partner(chC0, -tau); oC = oC0
            assert np.all(nA[1:] >= nA[:-1]) and np.all(nC[1:] >= nC[:-1]), 'partner keys not monotone in tau'
        else:
            nA = partner(A['ch'], tau); oA = np.argsort(nA); nA = nA[oA]
            nC = partner(Cc['ch'], -tau); oC = np.argsort(nC); nC = nC[oC]
        budget = mem / 24
        if fast and prev_pairs and 1.3 * prev_pairs < budget:      # the previous bucket fitted one pass with room: skip the
            hAB = hCD = None; exp_b = 0.0; passes = [list(range(-T, T + 1))]   # histogram merge, one pass over all t
        else:
            hAB = KJ.merge_thist(nA, oA, A['c'], sB, oB, B['c'], 64, T); hCD = KJ.merge_thist(nC, oC, Cc['c'], sD, oD, Dd['c'], 64, T)
            both = hAB * hCD[::-1]; exp_b = float(both.sum() * 2.0 ** (-bpc * (R.L - 1 - KE - 1)))
            passes = []; cur = []; load = 0
            for t in np.arange(-T, T + 1)[np.argsort(-both)]:
                if both[t + T] == 0: continue
                w = hAB[t + T] + hCD[T - t]
                if load + w > budget and cur: passes.append(cur); cur = []; load = 0
                cur.append(int(t)); load += w
            if cur: passes.append(cur)
        found = 0; ncand = 0
        for pl in passes:
            mA = np.zeros(2 * T + 1, np.bool_); mC = np.zeros(2 * T + 1, np.bool_)
            for t in pl: mA[t + T] = True; mC[T - t] = True
            KA, IA, JA = KJ.stage1(nA, oA, A['c'], A['h'], sB, oB, B['c'], B['h'], 1, G, 64, mA, T)
            KC, IC, JC = KJ.stage1(nC, oC, Cc['c'], Cc['h'], sD, oD, Dd['c'], Dd['h'], -1, G, 64, mC, T)
            ia, ic = KJ.bloom_match(KA, KC, 30); ncand += len(ia); npairs = [len(KA), len(KC)]
            for q, qc in zip(ia, ic):
                bits = [pools['A'][IA[q]], pools['B'][JA[q]], pools['C'][IC[qc]], pools['D'][JC[qc]]]
                Y = np.array([[1 if (int(bb) >> qq) & 1 else -1 for qq in range(R.r)] for bb in bits]); X = R.expand(Y)
                if np.all(R.paf(X).sum(0)[1:] == 0):
                    found += 1; nsol += 1
                    rec = dict(v=v, d=d, profile=list(prof), bucket=b, blocks={kk: ''.join('+' if t > 0 else '-' for t in x) for kk, x in zip('ABCD', X)})
                    solf.write(json.dumps(rec) + '\n'); solf.flush(); log(stage='SOLUTION', b=b)
            del KA, IA, JA, KC, IC, JC
        exp_cum += exp_b
        ABn, CDn = (int(hAB.sum()), int(hCD.sum())) if hAB is not None else npairs   # one-pass: pairs within |t| <= T
        prev_pairs = ABn + CDn if len(passes) == 1 else 0
        log(stage='bucket', b=b, tau=tau.tolist(), AB=ABn, CD=CDn, passes=len(passes), cand=int(ncand),
            solutions=found, total=nsol, exp_hits=round(exp_b, 5), exp_cum=round(exp_cum, 4))
    log(stage='end', solutions=nsol)


if __name__ == '__main__':
    main()
