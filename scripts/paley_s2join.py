"""Fast exact stage-2 join for the Paley(p) + multiplier slice (b = 1), report/paley_design_answer.md section 7.5.

Given a stage-1 solution (alpha A, constant-coset signs W; beta = A - W), find u (+-1 on the twisted multiplier orbits,
0 on constant cosets) with
    sum_k [PAF(beta_k) + p PAF(u_k)] = 4v delta        (S2 real)
    sum_k [CC(beta_k, u_k) - CC(u_k, beta_k)] = 0      (S2 imaginary, "X = 0")
Per block k the candidates are all 2^L_k sign vectors on its L_k twisted orbits.  Signatures are reduced to independent
coordinates (one lag per orbit of <mu, -1> for the PAF part, one per <mu>-orbit pair {h, -h} for X); constant
coordinates are checked against the target up front (early reject).
Join: 2+2 meet in the middle on a 64-bit linear hash of the signature.  Lists are paired (largest with smallest) and the
(A, B) pair space is cut into M residue classes of the hash mod M, so the hash table holds about |A||B|/M entries
(stamped open addressing, never cleared).  Work = |A||B| + |C||D| probes per y; every hash hit is re-checked exactly.
--dcomp d (d | n'): first restrict to class pairs whose Z_d-compressed S2 equations (additive over blocks) can be
completed by the other side; at 573 (d = 4) this keeps about 2^23-2^25 of the 2^34 pairs per side.
A solution is expanded, the core verified by 2-D FFT and completed to a bordered GS Hadamard matrix
(theory/paley_b1pipe.finish, exact H H^T = N I), which writes data/core_v{v}_b1_paley.txt.

usage: paley_s2join.py npz v p n mu [--first i] [--count c] [--table-log2 t] [--seed s] [--dcomp d]   (recorded y's)
       paley_s2join.py check419 [n_y] [d]   (calibration: y's from the published 419 gadget; compare with the numpy join)"""
import sys, os, json, time, math
import numpy as np
import numba
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))

def orbits(n, mu):
    seen, out = set(), []
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * mu) % n
        seen.update(o); out.append(sorted(o))
    return out

def pafc(x):
    x = np.atleast_2d(x); n = x.shape[1]
    return np.stack([(x * np.roll(x, -h, axis=1)).sum(1) for h in range(n)], 1)

def ccx(f, g):
    """CC(f, g)(h) = sum_x f(x) g(x + h) for each row g"""
    g = np.atleast_2d(g); n = g.shape[1]
    return np.stack([(f[None, :] * np.roll(g, -h, axis=1)).sum(1) for h in range(n)], 1)

def ccx_rev(g, f):
    """CC(g, f)(h) = sum_x g(x) f(x + h) for each row g"""
    g = np.atleast_2d(g); n = g.shape[1]
    return np.stack([(g * np.roll(f, -h)[None, :]).sum(1) for h in range(n)], 1)

def coord_reps(n, mu):
    """PAF coordinates: one lag per orbit of <mu, -1> (h != 0; lag 0 is constant per block);
    X coordinates: one lag per <mu>-orbit, keeping only one of each pair {O, -O} (X(-h) = -X(h); X = 0 on O = -O)"""
    seen, paf_reps = {0}, []
    for h in range(1, n):
        if h in seen: continue
        o = set(); stack = [h]
        while stack:
            x = stack.pop()
            if x in o: continue
            o.add(x); stack += [(x * mu) % n, (-x) % n]
        seen |= o; paf_reps.append(h)
    seen, x_reps = {0}, []
    for h in range(1, n):
        if h in seen: continue
        o = set(); x = h
        while x not in o: o.add(x); x = (x * mu) % n
        neg = {(-y) % n for y in o}
        seen |= o | neg
        if neg != o: x_reps.append(h)
    return paf_reps, x_reps

def block_lists(A, W, p, n, orbs, paf_reps, x_reps):
    """per block: U (2^L, n) int8 and reduced signature (2^L, D) int64 (FFT-based); plus the lag-0 value"""
    out = []
    for k in range(4):
        beta = (A[k] - W[k]).astype(np.float64)
        tw = [o for o in orbs if W[k, o[0]] == 0]; L = len(tw)
        idx = np.arange(2 ** L, dtype=np.int64)
        S = (1 - 2 * ((idx[:, None] >> np.arange(L)) & 1)).astype(np.int8)
        U = np.zeros((2 ** L, n), np.int8)
        for j, o in enumerate(tw): U[:, o] = S[:, [j]]
        Uh = np.fft.rfft(U.astype(np.float32), axis=1); bh = np.fft.rfft(beta)
        Pu = np.rint(np.fft.irfft(np.abs(Uh) ** 2, n=n, axis=1)).astype(np.int64)
        Pb = np.rint(np.fft.irfft(np.abs(bh) ** 2, n=n)).astype(np.int64)
        # CC(beta, u)(h) = sum_x beta(x) u(x+h) = ifft(conj(bh) Uh);  CC(u, beta)(h) = ifft(conj(Uh) bh)
        X = np.rint(np.fft.irfft(np.conj(bh)[None, :] * Uh - np.conj(Uh) * bh[None, :], n=n, axis=1)).astype(np.int64)
        P = Pb[None, :] + p * Pu
        sig = np.concatenate([P[:, paf_reps], X[:, x_reps]], 1)
        del Uh, Pu, X
        out.append((U, sig, int(P[0, 0])))
    return out

@numba.njit(cache=True)
def _join(kA, kB, kC, kD, kT, M, bstart, border, dstart, dorder, tbits, maxhits):
    """all (a, b, c, d) with kA[a] + kB[b] + kC[c] + kD[d] == kT (uint64 wraparound), residue-bucketed on mod M.
    bstart/border: B indices grouped by kB mod M (CSR); same for D.  Returns up to maxhits candidate quadruples."""
    size = 1 << tbits; mask = size - 1
    tkey = np.zeros(size, np.uint64); ta = np.zeros(size, np.int32); tb = np.zeros(size, np.int32)
    stamp = np.zeros(size, np.int32)
    hits = np.zeros((maxhits, 4), np.int64); nh = 0
    overflow = 0; probes = 0
    Mu = np.uint64(M)
    for c in range(M):
        st = c + 1; cu = np.uint64(c); filled = 0
        # insert all (a, b) with (kA + kB) mod M == c
        for a in range(kA.shape[0]):
            rb = (cu - kA[a] % Mu) % Mu
            for t in range(bstart[rb], bstart[rb + 1]):
                b = border[t]; key = kA[a] + kB[b]
                slot = int(key >> np.uint64(64 - tbits)) & mask
                dup = False
                while stamp[slot] == st:
                    if tkey[slot] == key:
                        dup = True; break
                    slot = (slot + 1) & mask
                if dup: continue
                stamp[slot] = st; tkey[slot] = key; ta[slot] = a; tb[slot] = b
                filled += 1
                if filled * 10 > size * 7: overflow += 1
        probes += filled
        # probe all (c_, d) with (kC + kD) mod M == (kT - c) mod M
        cc = (kT % Mu - cu + Mu) % Mu
        for ci in range(kC.shape[0]):
            rd = (cc - kC[ci] % Mu + Mu) % Mu
            for t in range(dstart[rd], dstart[rd + 1]):
                d = dorder[t]; need = kT - kC[ci] - kD[d]
                slot = int(need >> np.uint64(64 - tbits)) & mask
                while stamp[slot] == st:
                    if tkey[slot] == need:
                        if nh < maxhits:
                            hits[nh, 0] = ta[slot]; hits[nh, 1] = tb[slot]; hits[nh, 2] = ci; hits[nh, 3] = d
                        nh += 1; break
                    slot = (slot + 1) & mask
                probes += 1
    return hits[:min(nh, maxhits)], nh, overflow, probes

def csr(keys, M):
    r = (keys % np.uint64(M)).astype(np.int64)
    order = np.argsort(r, kind='stable').astype(np.int64)
    start = np.zeros(M + 1, np.int64); np.add.at(start, r + 1, 1); start = np.cumsum(start)
    return start, order

def compress_features(U, beta, p, d):
    """Z_d-level statistic per candidate (additive over blocks): compressions ub[r] = sum_{g = r mod d} u(g);
    F = (PAF_d(betabar) + p PAF_d(ubar) at lags 0..d//2,  Xbar_d at lags 1..(d-1)//2).  Target: (4v, 0, ..., 0)."""
    U = np.asarray(U, np.int64)
    ub = np.stack([U[:, r::d].sum(1) for r in range(d)], 1); bb = np.array([beta[r::d].sum() for r in range(d)])
    cols = []
    for h in range(d // 2 + 1):
        cols.append(p * (ub * np.roll(ub, -h, axis=1)).sum(1) + int((bb * np.roll(bb, -h)).sum()))
    for h in range(1, (d - 1) // 2 + 1):
        cols.append((bb[None, :] * np.roll(ub, -h, axis=1)).sum(1) - (ub * np.roll(bb, -h)[None, :]).sum(1))
    return np.stack(cols, 1)

@numba.njit(cache=True)
def _join_classes(kA, kB, kC, kD, rA, rB, rC, rD, P, kT, A_st, A_ix, B_st, B_ix, C_st, C_ix, D_st, D_ix, ABp, CDp,
                  tbits, maxhits):
    """all (a, b, c, d) with kA[a] + kB[b] + kC[c] + kD[d] == kT (uint64 wraparound), restricted to the class pairs
    ABp (rows (iA, iB)) and CDp.  Residue passes: r = (rX + rY) mod P (a separate additive hash mod the prime P);
    pass c holds the AB pairs with r == c and probes the CD pairs with r == -c (mod P), the target residue being 0.
    Duplicate AB keys are stored once.  Returns (hits, n_hits, overflow, probes); overflow > 0 means a pass filled
    the table past 85% and the result is INCOMPLETE (caller retries with a larger P)."""
    size = 1 << tbits; mask = size - 1; cap = (size * 85) // 100
    tkey = np.zeros(size, np.uint64); ta = np.zeros(size, np.int32); tb = np.zeros(size, np.int32)
    stamp = np.zeros(size, np.int32)
    hits = np.zeros((maxhits, 4), np.int64); nh = 0; overflow = 0; probes = 0
    for c in range(P):
        st = c + 1; filled = 0
        for q in range(ABp.shape[0]):
            ia = ABp[q, 0]; ib = ABp[q, 1]
            for s in range(A_st[ia], A_st[ia + 1]):
                a = A_ix[s]; need_r = (c - rA[a]) % P
                for t in range(B_st[ib], B_st[ib + 1]):
                    b = B_ix[t]
                    if rB[b] != need_r: continue
                    if filled >= cap:
                        overflow += 1; continue
                    key = kA[a] + kB[b]
                    slot = int(key >> np.uint64(64 - tbits)) & mask
                    dup = False
                    while stamp[slot] == st:
                        if tkey[slot] == key:
                            dup = True; break
                        slot = (slot + 1) & mask
                    if dup: continue
                    stamp[slot] = st; tkey[slot] = key; ta[slot] = a; tb[slot] = b
                    filled += 1
        probes += filled
        cc = (P - c) % P
        for q in range(CDp.shape[0]):
            ic = CDp[q, 0]; idd = CDp[q, 1]
            for s in range(C_st[ic], C_st[ic + 1]):
                ci = C_ix[s]; need_r = (cc - rC[ci]) % P
                for t in range(D_st[idd], D_st[idd + 1]):
                    d = D_ix[t]
                    if rD[d] != need_r: continue
                    need = kT - kC[ci] - kD[d]
                    slot = int(need >> np.uint64(64 - tbits)) & mask
                    while stamp[slot] == st:
                        if tkey[slot] == need:
                            if nh < maxhits:
                                hits[nh, 0] = ta[slot]; hits[nh, 1] = tb[slot]; hits[nh, 2] = ci; hits[nh, 3] = d
                            nh += 1; break
                        slot = (slot + 1) & mask
                    probes += 1
    return hits[:min(nh, maxhits)], nh, overflow, probes

def next_prime(x):
    x = max(2, int(x))
    while any(x % q == 0 for q in range(2, int(x ** 0.5) + 1)): x += 1
    return x

def classes(F):
    """unique rows of F -> (class feature array, CSR start, element index order, class id per element)"""
    uniq, inv = np.unique(F, axis=0, return_inverse=True); inv = inv.ravel()
    order = np.argsort(inv, kind='stable').astype(np.int64)
    st = np.zeros(len(uniq) + 1, np.int64); np.add.at(st, inv + 1, 1)
    return uniq, np.cumsum(st), order

def compatible_pairs(FX, FY, FZ, FW, T):
    """class pairs (i, j) of lists X, Y whose feature sum s has a class pair of Z, W with sum T - s, and vice versa
    (features encoded as mixed-radix int64 codes, so sums of codes are codes of sums)"""
    allF = np.concatenate([FX, FY, FZ, FW]); lo = allF.min(0); span = allF.max(0) - lo
    # pair-sum digits lie in [0, 2 span]; (T - 4 lo) - digits differs from any pair digit by at most
    # 4 span + |T - 4 lo|, so this radix makes integer equality of codes equivalent to digit-wise equality
    rad = 4 * span + 1 + np.abs(T - 4 * lo)
    mult = np.concatenate([[1], np.cumprod(rad[:-1])]).astype(np.int64)
    assert np.log2(float(np.prod(rad.astype(float)))) < 62
    code = lambda F: ((F - lo) * mult).sum(1)
    tcode = int(((T - 4 * lo) * mult).sum())
    xy = (code(FX)[:, None] + code(FY)[None, :]); zw = (code(FZ)[:, None] + code(FW)[None, :])
    def member(q, ref):
        u = np.unique(ref.ravel()); i = np.searchsorted(u, q.ravel()).clip(0, len(u) - 1)
        return (u[i] == q.ravel()).reshape(q.shape)
    okxy = member(tcode - xy, zw); okzw = member(tcode - zw, xy)
    return np.argwhere(okxy).astype(np.int64), np.argwhere(okzw).astype(np.int64)

_CACHE = {}
_CACHE_MAX = 12

def block_data(Ak, Wk, v, p, n, mu, seed, dcomp):
    """cached per-block data: U (int8), 64-bit keys, signature (int32), lag-0 value, and Z_d classes"""
    key = (Ak.tobytes(), Wk.tobytes(), p, n, mu, seed, dcomp)
    if key in _CACHE:
        _CACHE[key] = _CACHE.pop(key); return _CACHE[key]          # LRU refresh
    orbs = orbits(n, mu); paf_reps, x_reps = coord_reps(n, mu)
    A4 = np.zeros((4, n), np.int64); W4 = np.zeros((4, n), np.int64); A4[0] = Ak; W4[0] = Wk
    beta = (Ak - Wk).astype(np.float64)
    tw = [o for o in orbs if Wk[o[0]] == 0]; L = len(tw)
    idx = np.arange(2 ** L, dtype=np.int64)
    S = (1 - 2 * ((idx[:, None] >> np.arange(L)) & 1)).astype(np.int8)
    U = np.zeros((2 ** L, n), np.int8)
    for j, o in enumerate(tw): U[:, o] = S[:, [j]]
    Uh = np.fft.rfft(U.astype(np.float32), axis=1); bh = np.fft.rfft(beta)
    Pu = np.rint(np.fft.irfft(np.abs(Uh) ** 2, n=n, axis=1)).astype(np.int64)
    Pb = np.rint(np.fft.irfft(np.abs(bh) ** 2, n=n)).astype(np.int64)
    X = np.rint(np.fft.irfft(np.conj(bh)[None, :] * Uh - np.conj(Uh) * bh[None, :], n=n, axis=1)).astype(np.int64)
    del Uh
    lag0 = int(Pb[0] + p * Pu[0, 0])
    sig = np.concatenate([(Pb[None, paf_reps] + p * Pu[:, paf_reps]), X[:, x_reps]], 1).astype(np.int32)
    del Pu, X
    wts = np.random.default_rng(seed).integers(1, 2 ** 63, sig.shape[1], dtype=np.int64).astype(np.uint64) | np.uint64(1)
    keys = (sig.astype(np.int64).astype(np.uint64) * wts[None, :]).sum(1)
    # equal signatures complete identically: keep one representative per distinct signature (hash key; re-checked
    # exactly on every hit).  Without this, repeated pair sums form long probe chains in the hash table.
    nfull = len(keys)
    keys, first = np.unique(keys, return_index=True)
    U = U[first]; sig = sig[first]
    if dcomp: cl = classes(compress_features(U, Ak - Wk, p, dcomp))
    else: cl = (np.zeros((1, 1), np.int64), np.array([0, len(U)], np.int64), np.arange(len(U), dtype=np.int64))
    rv = np.random.default_rng(seed + 1).integers(0, 1 << 30, sig.shape[1], dtype=np.int64)
    rec = dict(U=U, sig=sig, keys=keys, lag0=lag0, cl=cl, nfull=nfull, rv=rv, res={})
    _CACHE[key] = rec
    while len(_CACHE) > _CACHE_MAX: _CACHE.pop(next(iter(_CACHE)))
    return rec

def solve_y(A, W, v, p, n, mu, seed=0, table_log2=24, maxhits=64, finish=True, dcomp=0):
    """stage 2 for one y.  Returns (record, solutions)"""
    t0 = time.time()
    B = [block_data(A[k], W[k], v, p, n, mu, seed, dcomp) for k in range(4)]
    sizes = [len(b['U']) for b in B]
    rec = dict(sizes=sizes, full_sizes=[b['nfull'] for b in B], setup_secs=round(time.time() - t0, 2))
    if sum(b['lag0'] for b in B) != 4 * v:
        rec.update(reject='lag0', solutions=0, secs=round(time.time() - t0, 2)); return rec, []
    order = sorted(range(4), key=lambda k: sizes[k])
    perm = [order[3], order[0], order[2], order[1]]
    kA, kB, kC, kD = [B[k]['keys'] for k in perm]; kT = np.uint64(0)
    cl = [B[k]['cl'] for k in perm]
    Td = np.zeros(cl[0][0].shape[1], np.int64); Td[0] = 4 * v if dcomp else 0
    if dcomp:
        ABp, CDp = compatible_pairs(cl[0][0], cl[1][0], cl[2][0], cl[3][0], Td)
    else:
        ABp = CDp = np.zeros((1, 2), np.int64)
    cA, cB, cC, cD = [np.diff(c[1]) for c in cl]
    nAB = int((cA[ABp[:, 0]] * cB[ABp[:, 1]]).sum()) if len(ABp) else 0
    nCD = int((cC[CDp[:, 0]] * cD[CDp[:, 1]]).sum()) if len(CDp) else 0
    rec.update(classes=[len(c[0]) for c in cl], kept_log2=[round(math.log2(max(1, nAB)), 1), round(math.log2(max(1, nCD)), 1)])
    if nAB == 0 or nCD == 0:
        rec.update(solutions=0, reject='Z%d-level' % dcomp, secs=round(time.time() - t0, 2)); return rec, []
    P = next_prime(max(1, nAB / (1 << (table_log2 - 2))))
    while True:
        def resid(k):
            b_ = B[k]
            if P not in b_['res']: b_['res'][P] = ((b_['sig'].astype(np.int64) % P) @ (b_['rv'] % P)) % P
            return b_['res'][P]
        rA, rB, rC, rD = [resid(k) for k in perm]
        hits, nh, overflow, probes = _join_classes(kA, kB, kC, kD, rA, rB, rC, rD, P, kT, cl[0][1], cl[0][2], cl[1][1],
                                                   cl[1][2], cl[2][1], cl[2][2], cl[3][1], cl[3][2], ABp, CDp, table_log2,
                                                   maxhits)
        if not overflow: break
        P = next_prime(2 * P)                                        # incomplete pass: retry with more buckets
    M = P
    rec.update(P=M, hash_hits=int(nh), overflow=int(overflow), probes_log2=round(math.log2(max(1, probes)), 1),
               perm=perm, join_secs=round(time.time() - t0, 2))
    sols = []
    for hit in hits:
        idx = dict(zip(perm, hit))
        full = sum(B[k]['sig'][idx[k]].astype(np.int64) for k in range(4))
        if (full == 0).all():
            sols.append(np.array([B[k]['U'][idx[k]] for k in range(4)], np.int64))
    rec['solutions'] = len(sols)
    if sols and finish:
        import paley_b1pipe as BP
        fin = BP.finish(A, sols[0], W, v, p, n)
        rec.update(core_ok=fin['core_ok'], hadamard=fin['hadamard'])
    rec['secs'] = round(time.time() - t0, 2)
    return rec, sols

def check419(n_y=40, dcomp=2):
    """calibration: stage-1 solutions from the published 419 gadget; compare existence of stage-2 solutions with the
    numpy reference join (theory/paley_slicepipe.stage2_slice)"""
    import paley_slicepipe as SP
    v, p, n, mu = 419, 11, 38, 7; orbs = orbits(n, mu)
    W, Jp = SP.published419_gadget()
    class Capped(list):
        def append(self, x):
            super().append(x)
            if len(self) >= n_y: raise StopIteration
    rec = Capped(); stats = {}
    try: SP.run_gadget(W, Jp, v, p, n, orbs, stats=stats, stage2=False, record=rec)
    except StopIteration: pass
    print(json.dumps(dict(stage1_ys=len(rec), **{k: stats[k] for k in stats if k != 'alpha_lists'})), flush=True)
    agree = 0
    for i, (A, W_, Jp_) in enumerate(rec[:n_y]):
        t = time.time(); ref = SP.stage2_slice(A, W_, p, v, n, orbs, max_sols=10 ** 6); tref = time.time() - t
        r, sols = solve_y(A, W_, v, p, n, mu, maxhits=10 ** 6, finish=(i == 0 or len(ref) > 0), dcomp=dcomp)
        same = (len(ref) > 0) == (r['solutions'] > 0)
        agree += same
        print(json.dumps(dict(y=i, ref_solutions=len(ref), new=r, ref_secs=round(tref, 2), agree=same)), flush=True)
    print(json.dumps(dict(summary=True, compared=min(n_y, len(rec)), agree=agree)), flush=True)

if __name__ == "__main__":
    if sys.argv[1] == "check419":
        check419(int(sys.argv[2]) if len(sys.argv) > 2 else 40, int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    else:
        a = sys.argv[1:]
        opt = {"--first": "0", "--count": "1000000", "--table-log2": "24", "--seed": "0", "--dcomp": "0"}
        pos = []; i = 0
        while i < len(a):
            if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
            else: pos.append(a[i]); i += 1
        path = pos[0]; v, p, n, mu = map(int, pos[1:5])
        d = np.load(path); As, Ws = d['A'], d['W']
        f0 = int(opt["--first"]); cnt = int(opt["--count"])
        for j in range(f0, min(len(As), f0 + cnt)):
            r, sols = solve_y(As[j], Ws[j], v, p, n, mu, seed=int(opt["--seed"]), table_log2=int(opt["--table-log2"]),
                              dcomp=int(opt["--dcomp"]))
            print(json.dumps(dict(file=os.path.basename(path), y=j, **r)), flush=True)
            if sols and r.get('hadamard'):
                print(json.dumps(dict(FOUND=True, v=v, file=path, y=j)), flush=True)
                break
