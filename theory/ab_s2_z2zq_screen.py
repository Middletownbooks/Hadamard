r"""Exact Z_2 x Z_q-level screens of stage 2 (a staged c = 4 slice, Q = Z_2^2) [P, M].

For an element r of order 2 of Q, compress a core along <r> x H': xbar_k(a, t) = sum over q in the coset a of <r>
and h in H' of x_k(q, h, t), a K-invariant function on Z_2 x Z_q (a = 0 for the coset of 0).  Summing
sum_k PAF_G(x_k) = T over the same subgroup gives, for every core,                                           [P]
        sum_k PAF_{Z_2 x Z_q}(xbar_k) = T_r,    T_r(a, t) = sum of T over the points of G projecting to (a, t).
This is the pair of characters {psi_0, psi_r} of theory/ab_s2_zq_screen.py taken jointly, so it is stronger than both.
For a stage-1 solution y, xbar_k(a, .) is a sum over the G'-orbits j in half a of the contributions (fibre
patterns of the letter y_kj, weighted by |j|); so each half of each block has an exact DP distribution over
K-invariant vectors (states at the K-orbit representatives), depending only on the multiset of its orbits'
contribution sets (the half's type), and the block's joint states are pairs of half states with lag-0 energy within
the budget.  Per distinct quadruple of block types: the exact number of lifts by a 2+2 convolution of the PAF
signatures on Z_2 x Z_q (all lag classes but one, which the fixed row sum of y determines), products of exact counts
summed in float64 (zero/non-zero exact).  Cached by content (restart-safe).
usage: python3 theory/ab_s2_z2zq_screen.py CFG.json YBIN OUT.jsonl R --keep-alive ALIVE.ybin
       (R = 1, 2 or 3: the element r = (b, a) of Q with index 2 a + b, as in the Q-character order)"""
import sys, os, json, time, math
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s2_s4 as S4


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); YB, OUT, r = a[1], a[2], int(a[3])
    keep = a[a.index('--keep-alive') + 1] if '--keep-alive' in a else None
    run(W.build(cfg), YB, OUT, r, keep)


def run(sl, YB, OUT, r, keep=None):
    L1 = len(sl.breps); RB = 4 * L1; t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    assert tuple(sl.Gp.mods[:2]) == (2, 2) and r in (1, 2, 3)
    Z = S4.zq_setup(sl); kr = Z['kreps']; nk = len(kr); kor = Z['kor']; lreps = Z['lreps']; nl = len(lreps)
    ksz = np.bincount(kor, minlength=nk)
    qj = np.array([int(sl.Gp.co[b, 0] * 2 + sl.Gp.co[b, 1]) for b in sl.breps])
    half = np.where((qj == 0) | (qj == r), 0, 1)                                 # half of every G'-orbit
    # target on Z_2 x Z_q at the lag classes (a, tau-class)
    qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]; ag = np.where((qg == 0) | (qg == r), 0, 1)
    tt = np.zeros((2, Z['q']), np.int64); np.add.at(tt, (ag, sl.zco), np.asarray(sl.T, np.int64))
    Tt = np.concatenate([tt[0][lreps], tt[1][lreps]])[:-1]; T0 = int(Tt[0]); d = len(Tt)
    # contribution sets (unsigned), ids = ranks of their contents
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
    j0 = np.nonzero(half == 0)[0]; j1 = np.nonzero(half == 1)[0]
    mult = np.random.default_rng(23).integers(1, 2 ** 62, L1).astype(np.int64) | 1
    tid = {}

    def quad_keys(Y):
        C = lut[np.arange(L1)[None, None, :], Y.view(np.uint8)]
        C = np.concatenate([np.sort(C[:, :, j0], 2), np.sort(C[:, :, j1], 2)], 2)  # block type: (half 0, half 1)
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
    log(stage='types', r=r, y_total=n_all, block_types=len(tid), type_quadruples=len(quads))
    types = {i: t for t, i in tid.items()}; hcache = {}; bcache = {}

    def half_dist(ids):
        if ids in hcache: return hcache[ids]
        st = np.zeros((1, nk), np.int64); cnt = np.ones(1, np.int64); R = 2 * sl.Gp.n + 1
        for c_ in ids:
            V = cvec[c_]; nst = (st[:, None, :] + V[None, :, :]).reshape(-1, nk); ncnt = np.repeat(cnt, len(V))
            key = ((nst + sl.Gp.n) * (R ** np.arange(nk))).sum(1)
            u, inv = np.unique(key, return_inverse=True); inv = inv.ravel()
            c2 = np.zeros(len(u), np.int64); np.add.at(c2, inv, ncnt)
            first = np.zeros(len(u), np.int64); first[inv[::-1]] = np.arange(len(inv))[::-1]
            st, cnt = nst[first], c2
        e = (st ** 2 * ksz).sum(1); o = np.argsort(e)
        hcache[ids] = (st[o], cnt[o].astype(np.float64), e[o]); return hcache[ids]

    def block_min(i):
        t = types[i]; h0 = half_dist(t[:len(j0)]); h1 = half_dist(t[len(j0):]); return int(h0[2][0] + h1[2][0])

    def block_sigs(i, budget):
        t = types[i]; A0, c0, e0 = half_dist(t[:len(j0)]); A1, c1, e1 = half_dist(t[len(j0):])
        n0 = np.searchsorted(e0, budget - e1[0], 'right')
        ia, ib = [], []
        for x in range(n0):
            m = np.searchsorted(e1, budget - e0[x], 'right'); ia.append(np.full(m, x)); ib.append(np.arange(m))
        ia = np.concatenate(ia) if ia else np.zeros(0, np.int64); ib = np.concatenate(ib) if ib else np.zeros(0, np.int64)
        X0 = A0[ia][:, kor]; X1 = A1[ib][:, kor]
        P0 = np.stack([(X0 * np.roll(X0, -t_, 1)).sum(1) + (X1 * np.roll(X1, -t_, 1)).sum(1) for t_ in lreps], 1)
        P1 = np.stack([(X0 * np.roll(X1, -t_, 1)).sum(1) + (X1 * np.roll(X0, -t_, 1)).sum(1) for t_ in lreps], 1)
        P = np.concatenate([P0, P1], 1)[:, :-1]
        u, inv = np.unique(P, axis=0, return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), c0[ia] * c1[ib])
        return u, c

    dt = np.dtype([('f%d' % i_, np.int64) for i_ in range(d)])
    sv = lambda R: np.ascontiguousarray(R, np.int64).view(dt).ravel()

    def conv(x, y, chunk=1 << 22):
        us, cs = [], []; step = max(1, chunk // max(1, len(y[0])))
        for s0 in range(0, len(x[0]), step):
            R = (x[0][s0:s0 + step][:, None, :] + y[0][None, :, :]).reshape(-1, d); v = (x[1][s0:s0 + step][:, None] * y[1][None, :]).ravel()
            u, inv = np.unique(sv(R), return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), v); us.append(u); cs.append(c)
        u, inv = np.unique(np.concatenate(us), return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), np.concatenate(cs))
        return u, c

    def count(tq):
        mins = [block_min(i) for i in tq]
        if sum(mins) > T0: return 0.0, 'lag-0 energy'
        D = [block_sigs(i, T0 - (sum(mins) - mins[k])) for k, i in enumerate(tq)]
        if any(len(x[0]) == 0 for x in D): return 0.0, 'lag-0 energy'
        A_ = conv(D[0], D[1]); B_ = conv(D[2], D[3])
        q_ = sv(Tt[None, :] - B_[0].view(np.int64).reshape(-1, d))
        pos = np.minimum(np.searchsorted(A_[0], q_), len(A_[0]) - 1); hit = A_[0][pos] == q_
        return float((A_[1][pos[hit]] * B_[1][hit]).sum()), 'exact convolution'

    done = {}
    if os.path.exists(OUT):
        for l in open(OUT):
            z = json.loads(l)
            if z.get('r') == r and 'types_content' in z: done[tuple(sorted(tuple(t) for t in z['types_content']))] = z
    fo = open(OUT, 'a'); t_last = 0.0; cnt_of = {}; dead_y = 0
    for qi, (key, ny) in enumerate(sorted(quads.items(), key=lambda kv: -kv[1])):
        tq = [(key >> s) & 0x7FFF for s in (45, 30, 15, 0)]; ck = tuple(sorted(types[i] for i in tq))
        if ck not in done:
            ex, by = count(tq)
            rec = dict(r=r, types_content=[list(t) for t in ck], count=ex, log2_count=round(math.log2(ex), 3) if ex > 0 else None,
                       by=by, y_first_run=ny)
            fo.write(json.dumps(rec) + "\n"); fo.flush(); done[ck] = rec
        cnt_of[key] = done[ck]['count']; dead_y += ny if done[ck]['count'] == 0 else 0
        if time.time() - t_last > 60: log(stage='progress', quads_done=qi + 1, of=len(quads), y_dead_so_far=dead_y); t_last = time.time()
    alive = 0; fk = open(keep, 'wb') if keep else None
    for s0 in range(0, n_all, 1 << 20):
        Y = np.fromfile(YB, np.int8, count=min(1 << 20, n_all - s0) * RB, offset=s0 * RB).reshape(-1, 4, L1)
        ks = quad_keys(Y); u, inv = np.unique(ks, return_inverse=True)
        ok = np.array([cnt_of[k_] > 0 for k_ in u.tolist()], bool)[inv.ravel()]
        alive += int(ok.sum())
        if fk is not None and ok.any(): fk.write(np.ascontiguousarray(Y[ok]).tobytes())
    if fk is not None: fk.close()
    log(stage='done', r=r, y_total=n_all, type_quadruples=len(quads), y_dead=n_all - alive, y_alive=alive)
    return dict(alive=alive, done=done)


if __name__ == "__main__":
    main()
