r"""Exact Z_q-compression check of stage 2 with numpy reachable sets, for any b1 slice key and stage-1 list [M].

Same mathematics as theory/h4892_zqcheck.py (necessary condition sum_k PAF_{Z_q}(w_k)(t) = -4|G'|, t != 0, and
sum_k energy(w_k) = sum_k r_k^2 + (q-1) 4|G'|), shares no code with it.  Per block the reachable vectors
(w(0); w(C_0), .., w(C_{s-1})) are built orbit by orbit as int32 arrays, deduplicated (np.unique) after each orbit.
The PAF of the length-q expansion is a quadratic form in the s+1 values; it is evaluated on one lag per coset of M
(plus lag 0) for every reachable vector (exact int64).
Tests, in order, per stage-1 class:
  (1) range test: for each lag class, sum_k min_k <= target <= sum_k max_k (else no lift); then the window filter:
      a block value outside [T - sum_{others} max, T - sum_{others} min] in some coordinate cannot occur in a
      solution, so it is dropped, to a fixed point (outcome 'window' if a block becomes empty or a range fails);
  (2) exact 2+2 join over the distinct PAF vectors (left pairs bucketed by their lag-0 sum, right by the
      complementary lag-0 sum), when the pair counts are below --maxpairs; otherwise UNDECIDED for that class.
Planted: random reachable vectors of a random class, target = their PAF sum; must pass (1) and be found by (2).
usage: python3 theory/h4892_zqnp.py 'KEY' LIST [--planted N] [--maxpairs P] [--out FILE] [--only a:b]"""
import sys, json, itertools, time
import numpy as np


def main():
    np.seterr(over="ignore")
    a = sys.argv[1:]; key = json.loads(a[0]); L = json.load(open(a[1]))
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 3
    mp = int(float(a[a.index('--maxpairs') + 1])) if '--maxpairs' in a else 2 * 10 ** 8
    outp = a[a.index('--out') + 1] if '--out' in a else None
    fam, v, Q, H1, q, s, U = key; assert fam == 'b1'
    f = (q - 1) // s
    assert f >= 2, 'alpha is not determined by y when f = 1; not supported'
    primes = [p_ for p_ in range(2, q) if (q - 1) % p_ == 0 and all(p_ % d for d in range(2, p_))]
    g = next(x for x in range(2, q) if all(pow(x, (q - 1) // p_, q) != 1 for p_ in primes))
    Mg = sorted({pow(g, s * e, q) for e in range(f)})
    cos = [sorted({pow(g, j, q) * u % q for u in Mg}) for j in range(s)]
    cid = np.zeros(q, np.int64)                                   # position of z in the value vector (0 = w(0))
    for j, c in enumerate(cos): cid[c] = 1 + j
    lags = [0] + [c[0] for c in cos]
    # quadratic form: PAF(t) = sum_i u[cid[i]] u[cid[i+t]]
    Nq = np.zeros((len(lags), s + 1, s + 1), np.int64)
    for l, t in enumerate(lags):
        for i in range(q): Nq[l, cid[i], cid[(i + t) % q]] += 1
    sizes = [len(o) for o in L['orbits']]; nG = sum(sizes)
    pats = {sg: np.array([t for t in itertools.product((1, -1), repeat=s) if sum(t) == sg], np.int32)
            for sg in range(-s, s + 1, 2)}
    def decode(y):
        for al in (1, -1):
            if (y - al) % f == 0 and ((y - al) // f + s) % 2 == 0 and abs((y - al) // f) <= s: return al, (y - al) // f
        raise ValueError(y)
    cache = {}
    def reach(y):
        k = tuple(y)
        if k in cache: return cache[k]
        acc = np.zeros((1, s), np.int32); w0 = 0
        for val, m in sorted(zip(y, sizes), key=lambda t: -t[1]):
            al, sg = decode(val); w0 += m * al
            P = pats[sg] * m
            acc = np.unique((acc[:, None, :] + P[None, :, :]).reshape(-1, s), axis=0)
        U_ = np.concatenate([np.full((len(acc), 1), w0, np.int64), acc.astype(np.int64)], 1)
        PAF = np.stack([((U_ @ Nq[l]) * U_).sum(1) for l in range(len(lags))], 1)
        D = np.unique(PAF, axis=0)
        cache[k] = dict(n=len(acc), D=D, lo=D.min(0), hi=D.max(0)); return cache[k]
    T = np.array([4 + (q - 1) * 4 * nG] + [-4 * nG] * s, np.int64)
    rng = np.random.default_rng(q + s)
    hv = rng.integers(-2 ** 62, 2 ** 62, size=len(lags), dtype=np.int64)
    def decide(R, Tt):
        lo = sum(r['lo'] for r in R); hi = sum(r['hi'] for r in R)
        if (Tt < lo).any() or (Tt > hi).any(): return 'range', 0
        # window filter [P]: in a solution, block k has coordinate c within [Tt_c - sum_{others} hi_c,
        # Tt_c - sum_{others} lo_c]; iterate to a fixed point (applies to every lag class, energy first of all)
        R = [dict(D=r['D']) for r in R]
        while True:
            for r in R: r['lo'] = r['D'].min(0) if len(r['D']) else None; r['hi'] = r['D'].max(0) if len(r['D']) else None
            if any(r['lo'] is None for r in R): return 'window', 0
            lo = sum(r['lo'] for r in R); hi = sum(r['hi'] for r in R)
            if (Tt < lo).any() or (Tt > hi).any(): return 'window', 0
            ch = False
            for r in R:
                a_ = Tt - (hi - r['hi']); b_ = Tt - (lo - r['lo'])
                keep = ((r['D'] >= a_) & (r['D'] <= b_)).all(1)
                if not keep.all(): r['D'] = r['D'][keep]; ch = True
            if not ch: break
        # 2+2 join; pick the split with the fewest pairs
        nb = [len(r['D']) for r in R]
        (i, j), (k, l) = min([((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))],
                             key=lambda p: nb[p[0][0]] * nb[p[0][1]] + nb[p[1][0]] * nb[p[1][1]])
        if nb[i] * nb[j] > mp or nb[k] * nb[l] > mp: return 'undecided', None
        def ps(x, y):
            Hx = (R[x]['D'] * hv).sum(1); Hy = (R[y]['D'] * hv).sum(1)
            return (Hx[:, None] + Hy[None, :]).ravel()
        right = np.sort(ps(k, l)); need = (Tt * hv).sum() - ps(i, j)
        pos = np.minimum(np.searchsorted(right, need), len(right) - 1)
        hits = int((right[pos] == need).sum())
        if hits == 0: return 'join', 0
        # exact recheck of hash hits
        Hk = (R[k]['D'] * hv).sum(1); Hl = (R[l]['D'] * hv).sum(1); Hi = (R[i]['D'] * hv).sum(1); Hj = (R[j]['D'] * hv).sum(1)
        idx = np.nonzero(right[pos] == need)[0]; cnt = 0
        rd = {}
        for x in range(len(Hk)):
            for y in range(len(Hl)): rd.setdefault(int(Hk[x] + Hl[y]), []).append((x, y))
        th = (Tt * hv).sum()
        for h in idx:
            x, y = divmod(int(h), len(Hj))
            for (u, w) in rd.get(int(th - Hi[x] - Hj[y]), []):
                if (R[i]['D'][x] + R[j]['D'][y] + R[k]['D'][u] + R[l]['D'][w] == Tt).all(): cnt += 1
        return 'join', cnt
    sols = L['solutions']; t0 = time.time()
    pl = []
    for t in range(npl):
        Y = sols[int(rng.integers(len(sols)))]; R = [reach(y) for y in Y]
        Tp = sum(r['D'][int(rng.integers(len(r['D'])))] for r in R)
        pl.append(decide(R, Tp))
    res = []
    for Y in sols:
        res.append(decide([reach(y) for y in Y], T))
    from collections import Counter
    out = dict(key=key, script='theory/h4892_zqnp.py', classes=len(sols), planted=pl,
               outcome=Counter(f'{a_}:{b_}' for a_, b_ in res),
               reach_max=max(c['n'] for c in cache.values()), distinct_paf_max=max(len(c['D']) for c in cache.values()),
               secs=round(time.time() - t0, 1))
    print(json.dumps(out))
    if outp: json.dump(dict(out, per_class=res), open(outp, 'w'))


if __name__ == "__main__":
    main()
