r"""Exact cross level of 647 x353 given a joint stage-1 solution (y_A, y_B) [M]; plants and calibration data.

Setting (theory/nd_cross.py).  Per block k and Z_2-layer e the labels of the 35 M-orbits of that layer form a 5 x 7
+-1 matrix with prescribed weighted row sums (y_B) and column sums (y_A); nd_cross.Cross.layer_labellings lists them
exactly (L0: n0 rows, L1: n1 rows).  With a0, a1 the two layers at the 35 orbit representatives of Z_323 (Z_646 =
Z_2 x Z_323), sigma = (a0 + a1)/2 and delta = (a0 - a1)/2, and the core equations are [P]
    sum_k PAF_323(sigma_k) = TS = (645, -2, ..., -2),   sum_k PAF_323(delta_k) = TD = (647, 0, ..., 0)
at the 20 <M, -1>-lag classes.  Two identities make the pipeline cheap [P]:
  * code(sigma) = sum_b (sigma_b + 1) 3^b = c0(a0) + c1(a1) with c_e(a) = sum_b ((a_b + 1)/2) 3^b, so the distinct
    sigma of a block are the distinct values of c0[i] + c1[j] (one int64 unique, no row sorting);
  * PAF(delta) = (PAF(a0) + PAF(a1))/2 - PAF(sigma), so given an S-solution the D-equations read
    sum_k (P0_k[i_k] + P1_k[j_k]) = TE = 2 (TS + TD) = (2584, -4, ..., -4)   (P_e = PAF of the layer),
    i.e. the even-lag PAF of x on Z_646; no per-pair D-signature is stored.
Pipeline: distinct sigma per block (S-signature = PAF at the 20 lags) -> exact S-join (theory/p647_join.py) ->
D-check of every S-solution over the product of its sigma-preimages (compiled 2+2 on 64-bit hashes, every match
re-checked exactly) -> x_k on Z_646 -> full cyclic PAF check (sum_k PAF_k = (2584, -4, ..., -4)).
usage: python3 theory/p647_cross.py --plant N [--full]      (plants on random joint solutions of y_B #0)
       python3 theory/p647_cross.py --brute N                (reduced instances against brute force)"""
import sys, os, json, math, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
from p647_join import join4, join4p
from nd_cross import orbit_maps

PW3 = 3 ** np.arange(35, dtype=np.int64)          # 647: 35 orbit reps on Z_323 (3^35 < 2^63)

def pw3(n):
    assert n <= 39, "sigma code needs 3^n < 2^63"
    return 3 ** np.arange(n, dtype=np.int64)

def lag_tensor(X, oi):
    """NH[h, o, o'] = #{b in Z_323: b in orbit o, b + lag_h in orbit o'}, so PAF_h(v) = v^T NH[h] v"""
    mp = X.mp; no = int(oi.max()) + 1
    NH = np.zeros((len(X.lags), no, no), np.float64)
    b = np.arange(mp)
    for t, h in enumerate(X.lags):
        np.add.at(NH[t], (oi[b], oi[(b + h) % mp]), 1.0)
    return NH

def paf_orb(V, NH, chunk=1 << 16):
    """PAF on Z_323 at the lag classes of vectors given at the orbit reps (exact in float64, returned int64)"""
    V = np.asarray(V); out = np.empty((len(V), NH.shape[0]), np.int64)
    for s in range(0, len(V), chunk):
        F = V[s:s + chunk].astype(np.float64)
        for h in range(NH.shape[0]):
            out[s:s + chunk, h] = np.rint(((F @ NH[h]) * F).sum(1)).astype(np.int64)
    return out

class Block:
    """one block of one joint solution: layer lists, distinct sigma, preimage groups, S-signatures, layer PAFs"""
    def __init__(self, L0, L1, maps, NH):
        A0 = L0[:, maps[0]].astype(np.int8); A1 = L1[:, maps[1]].astype(np.int8)
        self.A0, self.A1 = A0, A1; self.n0, self.n1 = len(A0), len(A1)
        W = pw3(A0.shape[1]); self.W = W
        c0 = ((A0.astype(np.int64) + 1) // 2) @ W; c1 = ((A1.astype(np.int64) + 1) // 2) @ W
        codes = (c0[:, None] + c1[None, :]).ravel()
        U, inv = np.unique(codes, return_inverse=True)
        o = np.argsort(inv, kind='stable')
        self.pairs = o.astype(np.int64)                                     # pair index = i * n1 + j, grouped by sigma
        self.bd = np.searchsorted(inv[o], np.arange(len(U) + 1)).astype(np.int64)
        self.mult = np.diff(self.bd)
        self.U = U
        sig = ((U[:, None] // W[None, :]) % 3 - 1).astype(np.int8)          # distinct sigma at the orbit reps
        self.sig = sig
        self.Ssig = paf_orb(sig, NH)
        self.P0 = paf_orb(A0, NH); self.P1 = paf_orb(A1, NH)
        npts = int(NH[0].sum()); assert (self.P0[:, 0] == npts).all() and (self.P1[:, 0] == npts).all()
        self.N = len(U)

    def sigma_of_pair(self, p):
        i, j = divmod(int(p), self.n1)
        return (self.A0[i].astype(np.int64) + self.A1[j]) // 2

def build(X, yb, ya, maps=None, NH=None, oi=None):
    if maps is None:
        reps, oi, maps = orbit_maps(X)
    if NH is None: NH = lag_tensor(X, oi)
    Ls = X.lists(yb, ya)
    if min(min(len(a), len(b)) for a, b in Ls) == 0: return None
    return [Block(L0, L1, maps, NH) for L0, L1 in Ls]

# ---------------------------------------------------------------- D-check
def dprep(blocks, TE, seed=5):
    """concatenated per-block arrays for the compiled D-check"""
    L = len(TE); w = np.random.default_rng(seed).integers(1, 2 ** 62, L, dtype=np.uint64) * np.uint64(2) + np.uint64(1)
    d = dict(w=w, hTE=np.uint64((np.asarray(TE).astype(np.uint64) * w).sum(dtype=np.uint64)), TE=np.asarray(TE, np.int64))
    for k, B in enumerate(blocks):
        d[f'h0_{k}'] = (B.P0.astype(np.uint64) * w).sum(1, dtype=np.uint64)
        d[f'h1_{k}'] = (B.P1.astype(np.uint64) * w).sum(1, dtype=np.uint64)
    return d

@njit(cache=True)
def _dcheck(sols, pr0, bd0, n10, pr1, bd1, n11, pr2, bd2, n12, pr3, bd3, n13,
            h00, h10, h01, h11, h02, h12, h03, h13, P00, P10, P01, P11, P02, P12, P03, P13, hTE, TE, out, maxout, cand,
            tk, tv, lp):
    """for every S-solution: all preimage choices (p0..p3) with sum_k P0_k[i_k] + P1_k[j_k] = TE (2+2 on hashes,
    exact re-check); out rows = (solution index, p0, p1, p2, p3); cand[s] = number of D-candidates"""
    nout = 0; L = TE.shape[0]
    for s in range(sols.shape[0]):
        a0 = bd0[sols[s, 0]]; a1 = bd0[sols[s, 0] + 1]; b0 = bd1[sols[s, 1]]; b1 = bd1[sols[s, 1] + 1]
        c0 = bd2[sols[s, 2]]; c1 = bd2[sols[s, 2] + 1]; e0 = bd3[sols[s, 3]]; e1 = bd3[sols[s, 3] + 1]
        nl = (a1 - a0) * (b1 - b0); nr = (c1 - c0) * (e1 - e0)
        cand[s] = nl * nr
        tb = 4
        while (1 << tb) < 2 * nl: tb += 1
        msk = (1 << tb) - 1
        q = 0
        for x in range(a0, a1):
            p = pr0[x]; i = p // n10; j = p % n10
            hx = h00[i] + h10[j]
            for y in range(b0, b1):
                pp = pr1[y]; ii = pp // n11; jj = pp % n11
                key = hx + h01[ii] + h11[jj]
                slot = np.int64((key * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(64 - tb)) & msk
                while tv[slot] >= 0: slot = (slot + 1) & msk
                tk[slot] = key; tv[slot] = q; lp[q, 0] = x; lp[q, 1] = y; q += 1
        for x in range(c0, c1):
            p = pr2[x]; i = p // n12; j = p % n12
            hx = hTE - h02[i] - h12[j]
            for y in range(e0, e1):
                pp = pr3[y]; ii = pp // n13; jj = pp % n13
                key = hx - h03[ii] - h13[jj]
                slot = np.int64((key * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(64 - tb)) & msk
                while tv[slot] >= 0:
                    if tk[slot] == key:
                        qq = tv[slot]; xa = lp[qq, 0]; ya = lp[qq, 1]
                        pa = pr0[xa]; pb = pr1[ya]; pc = pr2[x]; pd = pr3[y]
                        good = True
                        for t in range(L):
                            v = (P00[pa // n10, t] + P10[pa % n10, t] + P01[pb // n11, t] + P11[pb % n11, t]
                                 + P02[pc // n12, t] + P12[pc % n12, t] + P03[pd // n13, t] + P13[pd % n13, t])
                            if v != TE[t]:
                                good = False; break
                        if good:
                            if nout < maxout:
                                out[nout, 0] = s; out[nout, 1] = pa; out[nout, 2] = pb; out[nout, 3] = pc; out[nout, 4] = pd
                            nout += 1
                    slot = (slot + 1) & msk
        # clear
        for x in range(a0, a1):
            p = pr0[x]; i = p // n10; j = p % n10
            hx = h00[i] + h10[j]
            for y in range(b0, b1):
                pp = pr1[y]; ii = pp // n11; jj = pp % n11
                key = hx + h01[ii] + h11[jj]
                slot = np.int64((key * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(64 - tb)) & msk
                while tv[slot] >= 0:
                    tv[slot] = -1; slot = (slot + 1) & msk
    return nout

def dcheck(blocks, sols, TE, maxout=1 << 16):
    """D-check of S-solutions (rows of 4 distinct-sigma indices); returns (hits as (s, p0..p3), candidates per sol)"""
    sols = np.ascontiguousarray(np.asarray(sols, np.int64).reshape(-1, 4))
    d = dprep(blocks, TE)
    out = np.zeros((maxout, 5), np.int64); cand = np.zeros(len(sols), np.int64)
    B = blocks
    nlmax = int((B[0].mult[sols[:, 0]] * B[1].mult[sols[:, 1]]).max()) if len(sols) else 1
    tb = 4
    while (1 << tb) < 2 * nlmax: tb += 1
    tk = np.zeros(1 << tb, np.uint64); tv = -np.ones(1 << tb, np.int64); lp = np.empty((nlmax, 2), np.int64)
    n = _dcheck(sols, B[0].pairs, B[0].bd, B[0].n1, B[1].pairs, B[1].bd, B[1].n1, B[2].pairs, B[2].bd, B[2].n1,
                B[3].pairs, B[3].bd, B[3].n1, d['h0_0'], d['h1_0'], d['h0_1'], d['h1_1'], d['h0_2'], d['h1_2'],
                d['h0_3'], d['h1_3'], B[0].P0, B[0].P1, B[1].P0, B[1].P1, B[2].P0, B[2].P1, B[3].P0, B[3].P1,
                d['hTE'], d['TE'], out, maxout, cand, tk, tv, lp)
    if n > maxout: return dcheck(blocks, sols, TE, maxout=2 * n)
    return out[:n], cand

def assemble(X, blocks, pr):
    """x_k on Z_646 from one pair index per block"""
    reps, oi, maps = orbit_maps(X)
    xs = []
    for B, p in zip(blocks, pr):
        i, j = divmod(int(p), B.n1)
        x = np.zeros(X.m, np.int64); x[X.E[0]] = B.A0[i].astype(np.int64)[oi]; x[X.E[1]] = B.A1[j].astype(np.int64)[oi]
        xs.append(x)
    return np.array(xs)

def full_paf(xs):
    return sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x.astype(float))) ** 2).real).astype(np.int64) for x in xs)

def is_core(xs, m=None):
    m = xs.shape[1] if m is None else m
    P = full_paf(xs); return bool(P[0] == 4 * m and (P[1:] == -4).all())

def solve(X, blocks, TS=None, TE=None, log=None, **kw):
    """exact cross level: S-join, D-check, cores (as pair indices); stats"""
    TS = X.TS if TS is None else TS
    TE = 2 * (X.TS + X.TD) if TE is None else TE
    t0 = time.time()
    Ssols, st = join4([B.Ssig for B in blocks], TS, log=log, **kw)
    t1 = time.time()
    hits, cand = dcheck(blocks, Ssols, TE) if Ssols else (np.zeros((0, 5), np.int64), np.zeros(0, np.int64))
    rec = dict(N=[B.N for B in blocks], Hs=round(sum(math.log2(B.N) for B in blocks), 3),
               H=round(sum(math.log2(B.n0 * B.n1) for B in blocks), 3), S_solutions=len(Ssols),
               S_pairs=int(cand.sum()), join=st, dcheck_secs=round(time.time() - t1, 2), D_hits=len(hits),
               secs=round(time.time() - t0, 1))
    return rec, Ssols, hits, cand

# ---------------------------------------------------------------- plants
def plant(X, blocks, rng):
    """one random pair per block; targets = its own S and even-lag sums"""
    pr = []; TS = 0; TE = 0
    for B in blocks:
        i = int(rng.integers(B.n0)); j = int(rng.integers(B.n1)); p = i * B.n1 + j
        s = int(np.searchsorted(B.bd, np.nonzero(B.pairs == p)[0][0], 'right') - 1)
        pr.append((p, s)); TS = TS + B.Ssig[s]; TE = TE + B.P0[i] + B.P1[j]
    return pr, np.asarray(TS, np.int64), np.asarray(TE, np.int64)

def run_plants(n, full=False, seed=0, log=print):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    yb = ys[0]; rng = np.random.default_rng(seed)
    reps, oi, maps = orbit_maps(X); NH = lag_tensor(X, oi)
    res = []
    for t in range(n):
        j = int(rng.integers(len(sols)))
        ya = X.J.ya_full(sols[j]['layers'])
        blocks = build(X, yb, ya, maps, NH, oi)
        pr, TS, TE = plant(X, blocks, rng)
        sidx = [s for p, s in pr]
        # energy class of the plant under the pairing join4 picks (or all energies with --full)
        Ssig = [B.Ssig for B in blocks]
        _, st0 = join4(Ssig, TS, stats_only=True)
        a, b = st0['order'][:2]
        els = None if full else [int(blocks[a].Ssig[sidx[a], 0] + blocks[b].Ssig[sidx[b], 0])]
        r, Ssols, hits, cand = solve(X, blocks, TS=TS, TE=TE, els=els)
        s_found = tuple(sidx) in set(map(tuple, Ssols))
        cores = [tuple(int(v) for v in h[1:]) for h in hits]
        p_found = tuple(p for p, s in pr) in set(cores)
        xs = assemble(X, blocks, [p for p, s in pr])
        P = full_paf(xs)
        out = dict(plant=t, joint_index=j, els=els, S_solutions=r['S_solutions'], S_pairs=r['S_pairs'], D_hits=r['D_hits'],
                   plant_S_found=s_found, plant_found=p_found, work_log2=r['join']['work_log2'],
                   left=r['join']['left'], right=r['join']['right'], ns_per_pair=r['join']['ns_per_pair'], secs=r['secs'])
        log(json.dumps(out)); res.append(out)
    log(json.dumps(dict(plants=n, full=full, found=sum(o['plant_found'] for o in res))))
    return res

if __name__ == "__main__":
    a = sys.argv[1:]
    if "--plant" in a:
        n = int(a[a.index("--plant") + 1]); seed = int(a[a.index("--seed") + 1]) if "--seed" in a else 0
        run_plants(n, full="--full" in a, seed=seed, log=lambda s: print(s, flush=True))

# ---------------------------------------------------------------- brute force on reduced instances
def brute_reduced(X, blocks, rng, keep=40):
    """plant one pair per block; keep the plant's sigma and keep-1 random other sigma per block (all their
    preimage pairs).  Pipeline (S-join on the kept sigma, D-check) against a brute force that never uses the
    pipeline's signatures: every kept pair -> its full sequence on Z_646 -> full cyclic PAF (FFT, 646 lags); all
    quadruples whose PAF sum equals the plant's, by a 2+2 dictionary on an exact 646-dim key."""
    pr, TS, TE = plant(X, blocks, rng)
    sub = []
    for B, (p, s) in zip(blocks, pr):
        others = rng.choice(B.N, size=min(keep - 1, B.N - 1), replace=False)
        idx = np.unique(np.concatenate([[s], others[others != s]]))
        sub.append(idx)
    # pipeline on the reduced lists
    Ssig = [B.Ssig[idx] for B, idx in zip(blocks, sub)]
    Ss, st = join4p(Ssig, TS, rbits=2, budget=1 << 12)
    Ss = np.array([[sub[k][t[k]] for k in range(4)] for t in Ss], np.int64).reshape(-1, 4)
    hits, cand = dcheck(blocks, Ss, TE)
    pipe = set(tuple(int(v) for v in h[1:]) for h in hits)
    # brute force from the raw sequences
    reps, oi, maps = orbit_maps(X)
    lists = []
    for k, (B, idx) in enumerate(zip(blocks, sub)):
        L = []
        for s in idx:
            for p in B.pairs[B.bd[s]:B.bd[s + 1]]:
                i, j = divmod(int(p), B.n1)
                x = np.zeros(X.m); x[X.E[0]] = B.A0[i][oi]; x[X.E[1]] = B.A1[j][oi]
                L.append((int(p), np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64)))
        lists.append(L)
    xs = assemble(X, blocks, [p for p, s in pr]); Tfull = full_paf(xs)
    left = {}
    for p0, v0 in lists[0]:
        for p1, v1 in lists[1]:
            left.setdefault((v0 + v1).tobytes(), []).append((p0, p1))
    brute = set()
    for p2, v2 in lists[2]:
        for p3, v3 in lists[3]:
            for (p0, p1) in left.get((Tfull - v2 - v3).tobytes(), []):
                brute.add((p0, p1, p2, p3))
    planted = tuple(p for p, s in pr)
    return dict(kept=[len(i) for i in sub], pairs=[len(L) for L in lists], S_solutions=len(Ss), D_candidates=int(cand.sum()),
                pipeline=len(pipe), brute=len(brute), agree=pipe == brute, plant_in=planted in pipe)

def run_brute(n, seed=0, log=print):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = lag_tensor(X, oi)
    rng = np.random.default_rng(seed); ok = 0
    for t in range(n):
        j = int(rng.integers(len(sols)))
        blocks = build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
        r = brute_reduced(X, blocks, rng)
        ok += r['agree'] and r['plant_in']
        log(json.dumps(dict(trial=t, joint_index=j, **r)))
    log(json.dumps(dict(brute_trials=n, agree_and_plant=ok)))

if __name__ == "__main__" and "--brute" in sys.argv:
    run_brute(int(sys.argv[sys.argv.index("--brute") + 1]), log=lambda s: print(s, flush=True))

def brute_real(X, blocks, S_full, rng, nsol=25, extra=40, nd=3):
    """validation with the real target on reduced lists built around real S-solutions S_full (complete set of the
    full lists): (1) the pipeline join on the reduced lists returns exactly the members of S_full inside them;
    (2) an S-level brute force from scratch (sigma on Z_323, PAF by direct shifts, 2+2 dictionary) agrees;
    (3) for nd of them the D-check agrees with the full 646-lag PAF of every candidate sequence quadruple."""
    reps, oi, maps = orbit_maps(X)
    pick = S_full[rng.choice(len(S_full), size=min(nsol, len(S_full)), replace=False)]
    sub = []
    for k, B in enumerate(blocks):
        idx = np.unique(np.concatenate([pick[:, k], rng.choice(B.N, size=extra, replace=False)]))
        sub.append(idx)
    keep = [np.zeros(B.N, bool) for B in blocks]
    for k in range(4): keep[k][sub[k]] = True
    inside = set(tuple(int(v) for v in t) for t in S_full if all(keep[k][t[k]] for k in range(4)))
    Ss, st = join4p([B.Ssig[i] for B, i in zip(blocks, sub)], X.TS, rbits=2, budget=1 << 14)
    pipe = set(tuple(int(sub[k][t[k]]) for k in range(4)) for t in Ss)
    # independent S signatures: sigma on Z_323, PAF at the lag reps by direct shifts
    def sig_scratch(B, s):
        v = B.sig[s].astype(np.int64)[oi]
        return np.array([int(np.dot(v, np.roll(v, -h))) for h in X.lags], np.int64)
    sigs = [{int(s): sig_scratch(B, s) for s in idx} for B, idx in zip(blocks, sub)]
    left = {}
    for s0, v0 in sigs[0].items():
        for s1, v1 in sigs[1].items(): left.setdefault((v0 + v1).tobytes(), []).append((s0, s1))
    brute = set()
    for s2, v2 in sigs[2].items():
        for s3, v3 in sigs[3].items():
            for (s0, s1) in left.get((X.TS - v2 - v3).tobytes(), []): brute.add((s0, s1, s2, s3))
    # D-level on a few S-solutions: every candidate, full PAF on Z_646 from the raw sequences
    dagree = 0; dcands = 0
    TE = 2 * (X.TS + X.TD)
    for t in list(inside)[:nd]:
        hits, cand = dcheck(blocks, np.array([t], np.int64), TE)
        pipe_hits = set(tuple(int(v) for v in h[1:]) for h in hits)
        seqs = []
        for k, B in enumerate(blocks):
            L = []
            for p in B.pairs[B.bd[t[k]]:B.bd[t[k] + 1]]:
                i, j = divmod(int(p), B.n1)
                x = np.zeros(X.m); x[X.E[0]] = B.A0[i][oi]; x[X.E[1]] = B.A1[j][oi]
                L.append((int(p), np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64)))
            seqs.append(L)
        tgt = np.full(X.m, -4, np.int64); tgt[0] = 4 * X.m
        lf = {}
        for p0, v0 in seqs[0]:
            for p1, v1 in seqs[1]: lf.setdefault((v0 + v1).tobytes(), []).append((p0, p1))
        bh = set()
        for p2, v2 in seqs[2]:
            for p3, v3 in seqs[3]:
                for (p0, p1) in lf.get((tgt - v2 - v3).tobytes(), []): bh.add((p0, p1, p2, p3))
        dcands += int(cand.sum()); dagree += (bh == pipe_hits)
    return dict(kept=[len(i) for i in sub], S_inside=len(inside), S_pipeline=len(pipe), S_brute=len(brute),
                S_agree=(pipe == inside == brute), D_checked=min(nd, len(inside)), D_candidates=dcands, D_agree=dagree)

if __name__ == "__main__" and "--brute-real" in sys.argv:
    j = int(sys.argv[sys.argv.index("--brute-real") + 1]); n = int(sys.argv[sys.argv.index("--brute-real") + 2])
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = lag_tensor(X, oi)
    blocks = build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
    S_full = np.load(os.path.join(PC.LOG, "calib", f"j{j}_ssols.npy"))
    rng = np.random.default_rng(j); ok = 0
    for t in range(n):
        r = brute_real(X, blocks, S_full, rng)
        ok += r['S_agree'] and r['D_agree'] == r['D_checked']
        print(json.dumps(dict(trial=t, j=j, **r)), flush=True)
    print(json.dumps(dict(brute_real_trials=n, all_agree=ok)), flush=True)
