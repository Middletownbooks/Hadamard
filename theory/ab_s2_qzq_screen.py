r"""Exact Q x Z_q-level screen of stage 2 (a staged c = 4 slice, Q = Z_2^2) [P, M].

Compress a core along H' only: xbar_k(q, t) = sum over h in H' of x_k(q, h, t), a K-invariant function on Q x Z_q.
Summing sum_k PAF_G(x_k) = T over H' gives, for every core,                                                    [P]
        sum_k PAF_{Q x Z_q}(xbar_k) = T_Q,    T_Q(delta, t) = sum of T over the points of G projecting to (delta, t).
This is all four Q-characters of theory/ab_s2_zq_screen.py taken jointly (it contains every Z_2 x Z_q level of
theory/ab_s2_z2zq_screen.py).  For a stage-1 solution y, xbar_k(q, .) is a sum over the G'-orbits j of the coset q
of the contributions of the fibre patterns of the letter y_kj (weighted by |j|): each coset of each block has an
exact DP distribution of K-invariant vectors (depending only on the multiset of its orbits' contribution sets), and
the block's joint states are the 4-tuples of coset states with lag-0 energy within the budget.  Per distinct
quadruple of block types: the number of lifts by a 2+2 join of the PAF signatures on Q x Z_q (every lag class but
one, which the fixed row sum of y determines) by a 64-bit linear hash, in energy bands: a zero count is exact (every
true solution matches its hash), a positive one is the hashed count.  Cached by content (restart-safe), shardable.
usage: python3 theory/ab_s2_qzq_screen.py CFG.json YBIN OUT.jsonl [--shard i/n] [--keep-alive ALIVE.ybin]"""
import sys, os, json, time, math, glob
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s2_s4 as S4


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); YB, OUT = a[1], a[2]
    keep = a[a.index('--keep-alive') + 1] if '--keep-alive' in a else None
    shard = tuple(int(x) for x in a[a.index('--shard') + 1].split('/')) if '--shard' in a else None
    run(W.build(cfg), YB, OUT, keep, shard)


def run(sl, YB, OUT, keep=None, shard=None, maxjoint=1 << 24):
    L1 = len(sl.breps); RB = 4 * L1; t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    assert tuple(sl.Gp.mods[:2]) == (2, 2)
    Z = S4.zq_setup(sl); kr = Z['kreps']; nk = len(kr); kor = Z['kor']; lreps = Z['lreps']; nl = len(lreps)
    ksz = np.bincount(kor, minlength=nk)
    qj = np.array([int(sl.Gp.co[b, 0] * 2 + sl.Gp.co[b, 1]) for b in sl.breps])
    qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]
    tt = np.zeros((4, Z['q']), np.int64); np.add.at(tt, (qg, sl.zco), np.asarray(sl.T, np.int64))
    Tt = np.concatenate([tt[dl][lreps] for dl in range(4)])[:-1]; T0 = int(Tt[0]); d = len(Tt)
    cont = {}
    for j in range(L1):
        os_ = [o for o, _ in sl.parts[j]]
        for y in sl.letters[j]:
            C = AGs.fibre_patterns(sl, j, int(y))[0].astype(np.int64) @ Z['Mo'][os_][:, kr]
            cont[(j, int(y))] = tuple(sorted(map(tuple, C.tolist())))
    csid = {k: i for i, k in enumerate(sorted(set(cont.values())))}
    cvec = {i: np.array(k, np.int64) for k, i in csid.items()}
    lut = -np.ones((L1, 256), np.int16)
    for (j, y), k in cont.items(): lut[j, y % 256] = csid[k]
    jq = [np.nonzero(qj == q)[0] for q in range(4)]
    mult = np.random.default_rng(31).integers(1, 2 ** 62, L1).astype(np.int64) | 1
    tid = {}

    def quad_keys(Y):
        C = lut[np.arange(L1)[None, None, :], Y.view(np.uint8)]
        C = np.concatenate([np.sort(C[:, :, jq[q]], 2) for q in range(4)], 2)   # block type: the four coset types
        h = np.zeros(C.shape[:2], np.int64)
        for j in range(L1): h += C[:, :, j].astype(np.int64) * mult[j]
        u, first, inv = np.unique(h.ravel(), return_index=True, return_inverse=True); inv = inv.ravel()
        R = C.reshape(-1, L1); assert (R == R[first[inv]]).all(), "64-bit key collision"
        tmap = np.array([tid.setdefault(tuple(x), len(tid)) for x in R[first].tolist()], np.int64)
        assert len(tid) < 1 << 15
        return (np.sort(tmap[inv].reshape(-1, 4), 1) << np.array([45, 30, 15, 0], np.int64)).sum(1)

    n_all = os.path.getsize(YB) // RB; quads = Counter()
    for s0 in range(0, n_all, 1 << 20):
        Y = np.fromfile(YB, np.int8, count=min(1 << 20, n_all - s0) * RB, offset=s0 * RB).reshape(-1, 4, L1)
        u, c = np.unique(quad_keys(Y), return_counts=True)
        for k_, c_ in zip(u.tolist(), c.tolist()): quads[k_] += c_
    log(stage='types', y_total=n_all, block_types=len(tid), type_quadruples=len(quads))
    types = {i: t for t, i in tid.items()}; ccache = {}
    cut = np.cumsum([0] + [len(x) for x in jq])

    def coset_dist(ids):
        if ids in ccache: return ccache[ids]
        st = np.zeros((1, nk), np.int64); cnt = np.ones(1, np.int64); R = 2 * sl.Gp.n + 1
        for c_ in ids:
            V = cvec[c_]; nst = (st[:, None, :] + V[None, :, :]).reshape(-1, nk); ncnt = np.repeat(cnt, len(V))
            key = ((nst + sl.Gp.n) * (R ** np.arange(nk))).sum(1)
            u, inv = np.unique(key, return_inverse=True); inv = inv.ravel()
            c2 = np.zeros(len(u), np.int64); np.add.at(c2, inv, ncnt)
            first = np.zeros(len(u), np.int64); first[inv[::-1]] = np.arange(len(inv))[::-1]
            st, cnt = nst[first], c2
        e = (st ** 2 * ksz).sum(1); o = np.argsort(e, kind='stable')
        ccache[ids] = (st[o], cnt[o].astype(np.float64), e[o]); return ccache[ids]

    def cosets_of(i):
        t = types[i]; return [coset_dist(t[cut[q]:cut[q + 1]]) for q in range(4)]

    def block_min(i): return int(sum(c[2][0] for c in cosets_of(i)))

    Wh = np.random.default_rng(37).integers(1, 2 ** 62, d).astype(np.int64); bsig = {}; gmin = None

    def block_sigs(i):
        """the joint states of type i with lag-0 energy <= T0 - 3 gmin, aggregated by (energy, signature hash):
        energies (sorted), multiplicities, hashes"""
        if i in bsig: return bsig[i]
        cs = cosets_of(i); budget = T0 - 3 * gmin; mins = [int(c[2][0]) for c in cs]
        idx = np.zeros((1, 0), np.int64); en = np.zeros(1, np.int64); ct = np.ones(1)
        for q in range(4):
            S, c, e = cs[q]; room = budget - sum(mins[q + 1:])
            n_ = np.searchsorted(e, room - en, 'right')                         # per partial state: coset states that fit
            rep = np.repeat(np.arange(len(en)), n_); sub = np.arange(int(n_.sum())) - np.repeat(np.cumsum(n_) - n_, n_)
            if len(rep) > maxjoint: raise MemoryError(f'block type {i}: {len(rep)} joint states')
            idx = np.concatenate([idx[rep], sub[:, None]], 1); en = en[rep] + e[sub]; ct = ct[rep] * c[sub]
            if len(en) == 0: break
        if len(en) == 0:
            bsig[i] = (np.zeros(0, np.int64), np.zeros(0), np.zeros(0, np.int64)); return bsig[i]
        X = [cs[q][0][idx[:, q]][:, kor] for q in range(4)]                     # (n, |Z_q|) per coset
        P = []
        for dl in range(4):
            for t_ in lreps:
                P.append(sum((X[q] * np.roll(X[q ^ dl], -t_, 1)).sum(1) for q in range(4)))
        P = np.stack(P, 1)[:, :-1]
        assert (P[:, 0] == en).all()
        h = P @ Wh                                                              # aggregate by (energy, hash): states merged
        u, inv = np.unique(np.stack([P[:, 0], h], 1), axis=0, return_inverse=True)  # always share their energy band, so a
        c = np.zeros(len(u)); np.add.at(c, inv.ravel(), ct)                     # hash collision can only add, never hide
        o = np.argsort(u[:, 0], kind='stable')
        bsig[i] = (u[o, 0].copy(), c[o], u[o, 1].copy()); return bsig[i]

    hT = np.int64(Tt @ Wh)

    def pairs(x, y, lo, hi):
        a0 = np.searchsorted(y[0], lo - x[0], 'left'); a1 = np.searchsorted(y[0], hi - x[0], 'left')
        cnt = np.maximum(a1 - a0, 0)
        if cnt.sum() == 0: return np.zeros(0, np.int64), np.zeros(0)
        rep = np.repeat(np.arange(len(cnt)), cnt); ix = np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt) + np.repeat(a0, cnt)
        return x[2][rep] + y[2][ix], x[1][rep] * y[1][ix]                        # pair hashes and multiplicities

    def npairs(x, y, e): return int(np.maximum(np.searchsorted(y[0], e - x[0], 'left'), 0).sum())

    def count(tq, maxpairs=1 << 23):
        mins = [block_min(i) for i in tq]
        if sum(mins) > T0: return 0.0, 'lag-0 energy'
        D = []
        for k, i in enumerate(tq):
            e, c, h = block_sigs(i); m = np.searchsorted(e, T0 - (sum(mins) - mins[k]), 'right')
            if m == 0: return 0.0, 'lag-0 energy'
            D.append((e[:m], c[:m], h[:m]))
        sz = [len(x[0]) for x in D]                                              # the 2+2 split with the smallest pair lists
        (a_, b_), (c_, d_) = min((((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))),
                                 key=lambda sp: max(sz[sp[0][0]] * sz[sp[0][1]], sz[sp[1][0]] * sz[sp[1][1]]))
        D = [D[a_], D[b_], D[c_], D[d_]]
        lo_all, hi_all = int(D[0][0][0] + D[1][0][0]), int(T0 - D[2][0][0] - D[3][0][0] + 1)
        tot = 0.0; e0 = lo_all
        while e0 < hi_all:
            e1 = hi_all; base = npairs(D[0], D[1], e0)
            while npairs(D[0], D[1], e1) - base > maxpairs and e1 > e0 + 1: e1 = (e0 + e1 + 1) // 2
            hA, cA = pairs(D[0], D[1], e0, e1)
            if len(hA):
                hB, cB = pairs(D[2], D[3], T0 - e1 + 1, T0 - e0 + 1)
                if len(hB):                                                     # sorted A, prefix sums of multiplicities
                    o = np.argsort(hA, kind='stable'); hs = hA[o]; cs = np.concatenate([[0.0], np.cumsum(cA[o])])
                    q_ = hT - hB; lo = np.searchsorted(hs, q_, 'left'); hi = np.searchsorted(hs, q_, 'right'); m = hi > lo
                    if m.any(): tot += float((cB[m] * (cs[hi[m]] - cs[lo[m]])).sum())
            e0 = e1
        return tot, 'hashed convolution'

    gmin = min(block_min(i) for i in types)
    done = {}
    for fn in [OUT] + sorted(glob.glob(OUT + '.s*')):
        if not os.path.exists(fn): continue
        for l in open(fn):
            try: z = json.loads(l)
            except Exception: continue
            if 'types_content' in z: done[tuple(sorted(tuple(t) for t in z['types_content']))] = z
    fo = open(OUT if shard is None else OUT + '.s%d' % shard[0], 'a'); t_last = 0.0; cnt_of = {}; dead_y = 0
    for qi, (key, ny) in enumerate(sorted(quads.items(), key=lambda kv: -kv[1])):
        tq = [(key >> s) & 0x7FFF for s in (45, 30, 15, 0)]; ck = tuple(sorted(types[i] for i in tq))
        if shard is not None and qi % shard[1] != shard[0]: continue
        if ck not in done:
            ex, by = count(tq)
            rec = dict(types_content=[list(t) for t in ck], count=ex, log2_count=round(math.log2(ex), 3) if ex > 0 else None,
                       by=by, y_first_run=ny)
            fo.write(json.dumps(rec) + "\n"); fo.flush(); done[ck] = rec
        cnt_of[key] = done[ck]['count']; dead_y += ny if done[ck]['count'] == 0 else 0
        if time.time() - t_last > 60: log(stage='progress', quads_done=qi + 1, of=len(quads), y_dead_so_far=dead_y); t_last = time.time()
    if shard is not None:
        log(stage='shard-done', shard=list(shard)); return dict(done=done)
    alive = 0; fk = open(keep, 'wb') if keep else None
    for s0 in range(0, n_all, 1 << 20):
        Y = np.fromfile(YB, np.int8, count=min(1 << 20, n_all - s0) * RB, offset=s0 * RB).reshape(-1, 4, L1)
        ks = quad_keys(Y); u, inv = np.unique(ks, return_inverse=True)
        ok = np.array([cnt_of[k_] > 0 for k_ in u.tolist()], bool)[inv.ravel()]
        alive += int(ok.sum())
        if fk is not None and ok.any(): fk.write(np.ascontiguousarray(Y[ok]).tobytes())
    if fk is not None: fk.close()
    log(stage='done', y_total=n_all, type_quadruples=len(quads), y_dead=n_all - alive, y_alive=alive)
    return dict(alive=alive, done=done)


if __name__ == "__main__":
    main()
