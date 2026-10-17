r"""Exact Z_q-level screens of stage 2 over ALL stage-1 solutions of a staged c = 4 slice [M].

For a character psi of Q (the first two axes of G' = Q x H'; psi_0 trivial) put xbar_k^psi(t) = sum_{g'} psi(q(g'))
x_k(g', t), a K-invariant vector on Z_q.  Fourier-transforming sum_k PAF_G(x_k) = T over Q x H' at the characters
(psi, trivial on H') gives, for every core,                                                                   [P]
        sum_k PAF_{Z_q}(xbar_k^psi) = Tq^psi,   Tq^psi(t) = sum_{g in G, z(g) = t} psi(q(g)) T(g),
so each psi gives a necessary condition.  For a stage-1 solution y, block k's possible xbar_k^psi are sums over its
G'-orbits j of psi(q_j) times the contributions of the fibre patterns of the letter y_kj; the distribution depends
only on the multiset of those signed contribution sets (the block's psi-type).  Per type: an exact DP over its orbits
(states = xbar at the K-orbit representatives, exact integer counts) and the distribution of PAF signatures (P(0) and
every lag class but the last, which the fixed character sum zhat_k(psi) of y determines).  Per distinct quadruple of
types: the exact number of lifts meeting the condition (zero when the four minimal lag-0 energies exceed Tq^psi(0),
or the target lies outside every 4-sum; otherwise a 2+2 convolution pruned by the lag-0 budget, products of exact
counts summed in float64, so zero/non-zero is exact).  A y with count 0 for some psi has no lift: dead exactly.
usage: python3 theory/ab_s2_zq_screen.py CFG.json YBIN OUT.jsonl [--psi 0,1,2,3] [--keep-alive FILE]
  --keep-alive FILE: write the y's that pass every screen (raw int8 rows) to FILE"""
import sys, os, json, time, math
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s2_s4 as S4

CHARS = [(1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)]      # psi(q), q index = 2 a + b


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read())
    psis = [int(x) for x in a[a.index('--psi') + 1].split(',')] if '--psi' in a else [0, 1, 2, 3]
    run(W.build(cfg), a[1], a[2], psis, a[a.index('--keep-alive') + 1] if '--keep-alive' in a else None)


def run(sl, YB, OUT, psis, keep_alive=None):
    """the screens for the target sl.T (a test may replace it by a planted one)"""
    L1 = len(sl.breps); RB = 4 * L1; t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    assert tuple(sl.Gp.mods[:2]) == (2, 2)
    Z = S4.zq_setup(sl); kr = Z['kreps']; nl = len(Z['lreps'])
    qj = [int(sl.Gp.co[b, 0] * 2 + sl.Gp.co[b, 1]) for b in sl.breps]        # Q-coset of every G'-orbit
    qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]
    Tq = {}
    for p in psis:
        tt = np.zeros(Z['q'], np.int64); np.add.at(tt, sl.zco, np.array(CHARS[p])[qg] * np.asarray(sl.T)); Tq[p] = tt[Z['lreps']][:nl - 1]
    # per (psi, orbit, letter): the signed contribution set of the fibre patterns to xbar at the K-orbit reps
    csid = {}; lut = {p: -np.ones((L1, 256), np.int16) for p in psis}
    for j in range(L1):
        os_ = [o for o, _ in sl.parts[j]]
        for y in sl.letters[j]:
            C = AGs.fibre_patterns(sl, j, int(y))[0].astype(np.int64) @ Z['Mo'][os_][:, kr]
            for p in psis:
                key = tuple(sorted(map(tuple, (CHARS[p][qj[j]] * C).tolist())))
                lut[p][j, y % 256] = csid.setdefault(key, len(csid))
    cvec = {i: np.array(k, np.int64) for k, i in csid.items()}
    n_all = os.path.getsize(YB) // RB; tid = {}

    mult = np.random.default_rng(17).integers(1, 2 ** 62, L1).astype(np.int64) | 1

    def quad_keys(Y, p):
        C = np.sort(lut[p][np.arange(L1)[None, None, :], Y.view(np.uint8)], 2)  # (n, 4, L1) sorted set ids
        assert (C >= 0).all()
        h = np.zeros(C.shape[:2], np.int64)
        for j in range(L1): h += C[:, :, j].astype(np.int64) * mult[j]          # 64-bit key of the sorted ids
        u, first, inv = np.unique(h.ravel(), return_index=True, return_inverse=True); inv = inv.ravel()
        R = C.reshape(-1, L1)
        assert (R == R[first[inv]]).all(), "64-bit key collision"                 # keys are exact on this data
        tmap = np.array([tid.setdefault(tuple(r), len(tid)) for r in R[first].tolist()], np.int64)
        assert len(tid) < 1 << 15
        return (np.sort(tmap[inv].reshape(-1, 4), 1) << np.array([45, 30, 15, 0], np.int64)).sum(1)

    quads = {p: Counter() for p in psis}
    for s0 in range(0, n_all, 1 << 20):
        Y = np.fromfile(YB, np.int8, count=min(1 << 20, n_all - s0) * RB, offset=s0 * RB).reshape(-1, 4, L1)
        for p in psis:
            u, c = np.unique(quad_keys(Y, p), return_counts=True)
            for k_, c_ in zip(u.tolist(), c.tolist()): quads[p][k_] += c_
    log(stage='types', y_total=n_all, block_types=len(tid), type_quadruples={p: len(quads[p]) for p in psis})
    types = {i: t for t, i in tid.items()}; sig = {}

    def dist(i):
        if i in sig: return sig[i]
        st = np.zeros((1, len(kr)), np.int64); cnt = np.ones(1, np.int64); R = 2 * sl.Gp.n + 1
        for c_ in types[i]:
            V = cvec[c_]; nst = (st[:, None, :] + V[None, :, :]).reshape(-1, len(kr)); ncnt = np.repeat(cnt, len(V))
            key = ((nst + sl.Gp.n) * (R ** np.arange(len(kr)))).sum(1)
            u, inv = np.unique(key, return_inverse=True); inv = inv.ravel()
            c2 = np.zeros(len(u), np.int64); np.add.at(c2, inv, ncnt)
            first = np.zeros(len(u), np.int64); first[inv[::-1]] = np.arange(len(inv))[::-1]
            st, cnt = nst[first], c2
        Pq = S4.zq_paf(st[:, Z['kor']], Z['lreps'])[:, :nl - 1]
        u, inv = np.unique(Pq, axis=0, return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), cnt.astype(np.float64))
        sig[i] = (u, c, int(u[:, 0].min()), float(cnt.sum())); return sig[i]

    def count(tq, T):
        D = [dist(i) for i in tq]; mins = [d[2] for d in D]; T0 = int(T[0])
        if sum(mins) > T0: return 0.0, 'lag-0 energy'
        budget = [T0 - (sum(mins) - mins[k]) for k in range(4)]
        Dp = [(d[0][d[0][:, 0] <= budget[k]], d[1][d[0][:, 0] <= budget[k]]) for k, d in enumerate(D)]
        lo = np.min([d[0].min(0) for d in Dp], 0); hi = np.max([d[0].max(0) for d in Dp], 0); Rr = 4 * (hi - lo) + 1
        if ((T - 4 * lo) < 0).any() or ((T - 4 * lo) >= Rr).any(): return 0.0, 'signature range'
        assert float(np.prod(Rr.astype(float))) < 2 ** 62
        mul = np.cumprod(np.concatenate([[1], Rr[:-1]])).astype(np.int64)
        enc = lambda S: ((S - lo) * mul).sum(-1)

        def conv(x, y):
            k = (enc(x[0])[:, None] + enc(y[0])[None, :]).ravel(); v = (x[1][:, None] * y[1][None, :]).ravel()
            u, inv = np.unique(k, return_inverse=True); s = np.zeros(len(u)); np.add.at(s, inv.ravel(), v); return u, s
        A = conv(Dp[0], Dp[1]); B = conv(Dp[2], Dp[3])
        need = int(((T - 4 * lo) * mul).sum()) - B[0]
        pos = np.minimum(np.searchsorted(A[0], need), len(A[0]) - 1); hit = A[0][pos] == need
        return float((A[1][pos[hit]] * B[1][hit]).sum()), 'exact convolution'

    done = {}
    if os.path.exists(OUT):
        for l in open(OUT):
            z = json.loads(l)
            if 'key' in z: done[(z['psi'], z['key'])] = z
    fo = open(OUT, 'a'); t_last = 0.0; nq = 0
    for p in psis:
        for key, ny in sorted(quads[p].items(), key=lambda kv: -kv[1]):
            nq += 1
            if (p, key) in done: continue
            tq = [(key >> s) & 0x7FFF for s in (45, 30, 15, 0)]
            ex, by = count(tq, Tq[p])
            rec = dict(psi=p, key=key, types=tq, y=ny, min_P0=[dist(i)[2] for i in tq], count=ex,
                       log2_count=round(math.log2(ex), 3) if ex > 0 else None, by=by)
            fo.write(json.dumps(rec) + "\n"); fo.flush(); done[(p, key)] = rec
            if time.time() - t_last > 60: log(stage='progress', quads_done=nq, of=sum(len(quads[q]) for q in psis)); t_last = time.time()
    per_psi = {p: dict(y_dead=sum(r['y'] for (q, k), r in done.items() if q == p and r['count'] == 0)) for p in psis}
    # pass 2: the y's alive under every character
    alive = 0; fk = open(keep_alive, 'wb') if keep_alive else None
    for s0 in range(0, n_all, 1 << 20):
        Y = np.fromfile(YB, np.int8, count=min(1 << 20, n_all - s0) * RB, offset=s0 * RB).reshape(-1, 4, L1)
        ok = np.ones(len(Y), bool)
        for p in psis:
            ks = quad_keys(Y, p); u, inv = np.unique(ks, return_inverse=True)
            ok &= np.array([done[(p, k_)]['count'] > 0 for k_ in u.tolist()], bool)[inv.ravel()]
        alive += int(ok.sum())
        if fk is not None and ok.any(): fk.write(np.ascontiguousarray(Y[ok]).tobytes())
    if fk is not None: fk.close()
    log(stage='done', y_total=n_all, per_psi=per_psi, y_alive_all=alive, y_dead=n_all - alive)
    return dict(per_psi=per_psi, alive=alive, done=done)


if __name__ == "__main__":
    main()
