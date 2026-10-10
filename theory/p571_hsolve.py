r"""Hierarchical exact solver for b = 1 GS cores in a multiplier slice fixing Z6 pointwise (571 x27 slice and its
same-skeleton control v = 211), built on the exact Z6-level decomposition (theory/p571_levels.py):

  stage 1 (y = Z3-compression):  L1 class X1 = y0 + y1  ->  L2 refinement Xm = y0 - y1   (S1 = L1 + L2)
  stage 2 (Z' = y w^{z*}):       L3 class S' = Z'0 + Z'1 ->  L4 refinement D' = Z'0 - Z'1 (S2 = L3 + L4)

Each level is a 4-block condition on Z_q (q = m/6) that is additive over blocks; signatures are exact integer PAF
values at one lag per orbit of <M, -1> (plus lag 0), computed from orbit structure constants
N[h, i, j] = #{t in O_i : t + h in O_j}  (P(h) = sum_ij N[h,i,j] c_i conj c_j).  Eisenstein values a + b w are integer
pairs.  Joins: numba 2+2 meet in the middle on a random linear 64-bit hash with residue passes (as in
scripts/paley_s2join.py), every hit re-checked exactly.
usage: python3 theory/p571_hsolve.py stage1 m gen [--L 8192] [--rounds 1] [--seed 0] [--out file.npz]"""
import sys, os, json, math, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_levels import Skel

# ---------------------------------------------------------------- instance
class Inst:
    def __init__(self, m, gens):
        S = Skel(m, gens); self.S = S
        self.m, self.q, self.v, self.R = S.m, S.q, S.v, S.R
        self.sizes = S.sizes.astype(np.int64)
        self.lags = [0] + list(S.reps); H = len(self.lags)
        N = np.zeros((H, self.R, self.R), np.int64)
        for a, h in enumerate(self.lags):
            for t in range(self.q): N[a, S.oid[t], S.oid[(t + h) % self.q]] += 1
        self.N = N
        self.t1 = np.array([4 * self.v - 24] + [-24] * (H - 1), np.int64)     # L1 target (sum over blocks)
        self.t0 = np.array([4 * self.v] + [0] * (H - 1), np.int64)            # L2 target
        assert (4 * self.v - 4) // self.q == 24

# ---------------------------------------------------------------- signatures
@numba.njit(cache=True)
def sig_real(V, N):
    """V (L, R) int64 orbit values -> (L, H) PAF signatures"""
    L, R = V.shape; H = N.shape[0]
    out = np.zeros((L, H), np.int64)
    for l in range(L):
        for h in range(H):
            s = 0
            for i in range(R):
                vi = V[l, i]
                if vi == 0: continue
                acc = 0
                for j in range(R): acc += N[h, i, j] * V[l, j]
                s += vi * acc
            out[l, h] = s
    return out

@numba.njit(cache=True)
def sig_eis(A, B, N):
    """Eisenstein orbit values c = A + B w (L, R) -> signature (L, 1 + 2(H-1)): P(0) (real), then (a, b) of
    P(h) = sum_ij N[h,i,j] c_i conj(c_j) = a + b w for h = lags[1:]"""
    L, R = A.shape; H = N.shape[0]
    out = np.zeros((L, 2 * H - 1), np.int64)
    for l in range(L):
        for h in range(H):
            sa = 0; sb = 0
            for i in range(R):
                ai = A[l, i]; bi = B[l, i]
                if ai == 0 and bi == 0: continue
                for j in range(R):
                    n = N[h, i, j]
                    if n == 0: continue
                    e = A[l, j] - B[l, j]; f = -B[l, j]          # conj(c_j) = e + f w
                    sa += n * (ai * e - bi * f); sb += n * (ai * f + bi * e - bi * f)
            if h == 0: out[l, 0] = sa
            else:
                out[l, 2 * h - 1] = sa; out[l, 2 * h] = sb
    return out

# ---------------------------------------------------------------- generic 2+2 join (residue passes)
@numba.njit(cache=True)
def join4(kA, kB, kC, kD, rA, rB, rC, rD, rT, P, bstart, border, dstart, dorder, kT, tbits, maxhits):
    """all (a, b, c, d) with kA[a] + kB[b] + kC[c] + kD[d] == kT (uint64 wrap).  Pass c holds the AB pairs with
    (rA + rB) mod P == c, where r = an additive hash mod the prime P computed exactly from the integer signature
    (so residues add without wraparound); B is grouped by rB (CSR bstart/border), D by rD.  Duplicate AB keys are
    chained, so every combination is reported (up to maxhits)."""
    size = 1 << tbits; mask = size - 1
    tkey = np.zeros(size, np.uint64); ta = np.zeros(size, np.int32); tb = np.zeros(size, np.int32)
    stamp = np.zeros(size, np.int32)
    hits = np.zeros((maxhits, 4), np.int64); nh = 0; overflow = 0; probes = 0
    cap = (size * 85) // 100
    for c in range(P):
        st = c + 1; filled = 0
        for a in range(kA.shape[0]):
            rb = (c - rA[a]) % P
            for t in range(bstart[rb], bstart[rb + 1]):
                b = border[t]; key = kA[a] + kB[b]
                if filled >= cap:
                    overflow += 1; continue
                slot = int(key >> np.uint64(64 - tbits)) & mask
                while stamp[slot] == st: slot = (slot + 1) & mask
                stamp[slot] = st; tkey[slot] = key; ta[slot] = a; tb[slot] = b; filled += 1
        probes += filled
        cc = (rT - c) % P
        for ci in range(kC.shape[0]):
            rd = (cc - rC[ci]) % P
            for t in range(dstart[rd], dstart[rd + 1]):
                d = dorder[t]; need = kT - kC[ci] - kD[d]
                slot = int(need >> np.uint64(64 - tbits)) & mask
                while stamp[slot] == st:
                    if tkey[slot] == need:
                        if nh < maxhits:
                            hits[nh, 0] = ta[slot]; hits[nh, 1] = tb[slot]; hits[nh, 2] = ci; hits[nh, 3] = d
                        nh += 1
                    slot = (slot + 1) & mask
                probes += 1
    return hits[:min(nh, maxhits)], nh, overflow, probes

def next_prime(x):
    x = max(2, int(x))
    while any(x % p == 0 for p in range(2, int(x ** 0.5) + 1)): x += 1
    return x

def csr_res(r, P):
    order = np.argsort(r, kind='stable').astype(np.int64)
    start = np.zeros(P + 1, np.int64); np.add.at(start, r + 1, 1)
    return np.cumsum(start), order

def hash_keys(sig, wts):
    return (sig.astype(np.int64).astype(np.uint64) * wts[None, :]).sum(1)

def residues(sig, rho, P):
    return ((sig.astype(np.int64) % P) * (rho % P)[None, :] % P).sum(1) % P

def join_all(sigs, target, seed=0, tbits=22, maxhits=1 << 20):
    """exact 2+2 join of four signature arrays; returns verified index quadruples (n, 4)"""
    D = sigs[0].shape[1]; rng = np.random.default_rng(seed)
    wts = rng.integers(1, 2 ** 63, D, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    rho = rng.integers(1, 1 << 30, D, dtype=np.int64)
    keys = [hash_keys(s, wts) for s in sigs]; kT = hash_keys(target[None, :], wts)[0]
    nAB = len(keys[0]) * len(keys[1])
    P = next_prime(max(2, nAB / ((1 << tbits) * 0.5)))
    while True:
        rs = [residues(s, rho, P) for s in sigs]; rT = int(residues(target[None, :], rho, P)[0])
        bs, bo = csr_res(rs[1], P); ds, do = csr_res(rs[3], P)
        hits, nh, ov, pr = join4(keys[0], keys[1], keys[2], keys[3], rs[0], rs[1], rs[2], rs[3], rT, P, bs, bo, ds, do,
                                 kT, tbits, maxhits)
        if ov == 0: break
        P = next_prime(2 * P)
    st = dict(probes=int(pr), passes=P, hash_hits=int(nh), truncated=bool(nh > maxhits))
    if len(hits) == 0: return np.zeros((0, 4), np.int64), st
    tot = sigs[0][hits[:, 0]] + sigs[1][hits[:, 1]] + sigs[2][hits[:, 2]] + sigs[3][hits[:, 3]]
    ok = (tot == target[None, :]).all(1)
    return hits[ok], st

# ---------------------------------------------------------------- stage 1
YPAIRS = np.array([(a, b) for a in (1, -1, 3, -3) for b in (1, -1, 3, -3)], np.int64)   # 16 (y0, y1)

def x1_sampler(inst, allowed):
    """allowed (R, 16) bool: permitted (y0, y1) per orbit.  Returns (values, weights) per orbit for X1 = y0 + y1,
    weighted by x-multiplicity (so X1 is distributed as under uniform x: E X1^2 = 6, matching the L1 lag-0 target)"""
    vals = []; wts = []
    for i in range(inst.R):
        d = {}
        for p in np.nonzero(allowed[i])[0]:
            # weight = number of x-configurations of this (y0, y1): 3 per side with |y| = 1 (uniform-x pushforward)
            x1 = int(YPAIRS[p].sum()); d[x1] = d.get(x1, 0) + 3 ** int((np.abs(YPAIRS[p]) == 1).sum())
        vals.append(np.array(sorted(d), np.int64)); wts.append(np.array([d[x] for x in sorted(d)], float))
    return vals, wts

def sample_x1(inst, vals, wts, r, L, rng):
    """L distinct X1 orbit vectors with sum |O_i| X1_i = r, drawn with probability prop. to prod weights"""
    R = inst.R; sz = inst.sizes
    # backward DP over partial sums (dict per orbit)
    cnt = [None] * (R + 1); cnt[R] = {0: 1.0}
    for i in range(R - 1, -1, -1):
        E = {}
        for t, c in cnt[i + 1].items():
            for x, w in zip(vals[i], wts[i]):
                E[t + sz[i] * x] = E.get(t + sz[i] * x, 0.0) + c * w
        cnt[i] = E
    if r not in cnt[0]: return np.zeros((0, R), np.int64)
    out = set(); tries = 0
    last = 0
    while len(out) < L and tries < 20 * L and tries - last < max(4000, 4 * len(out)):
        tries += 1; t = r; row = []
        for i in range(R):
            ps = np.array([w * cnt[i + 1].get(t - sz[i] * x, 0.0) for x, w in zip(vals[i], wts[i])])
            j = rng.choice(len(ps), p=ps / ps.sum()); row.append(int(vals[i][j])); t -= sz[i] * vals[i][j]
        n0 = len(out); out.add(tuple(row))
        if len(out) > n0: last = tries
    return np.array(sorted(out), np.int64)

def refinements(inst, x1row, allowed):
    """all (y0, y1) orbit configurations with y0 + y1 = X1 (allowed pairs): returns Y0, Y1 (n, R)"""
    opts = []
    for i in range(inst.R):
        ps = [p for p in np.nonzero(allowed[i])[0] if YPAIRS[p].sum() == x1row[i]]
        opts.append(np.array(ps, np.int64))
    n = int(np.prod([len(o) for o in opts]))
    idx = np.zeros((n, inst.R), np.int64); rep = n
    for i, o in enumerate(opts):
        rep //= len(o); idx[:, i] = np.tile(np.repeat(o, rep), n // (rep * len(o)))
    return YPAIRS[idx, 0], YPAIRS[idx, 1]

def stage1(inst, L=8192, rounds=1, seed=0, allowed=None, rows=(2, 0, 0, 0), max_y=10 ** 6, log=None):
    """returns list of y = (4, 2, R) arrays (S1 solutions, verified exactly)"""
    rng = np.random.default_rng(seed)
    if allowed is None: allowed = np.ones((inst.R, 16), bool)
    vals, wts = x1_sampler(inst, allowed)
    found = []; stats = []
    for rd in range(rounds):
        t0 = time.time()
        X = []
        for r in rows:
            ex = enumerate_x1(inst, vals, wts, r, limit=L)
            X.append(ex if ex is not None else sample_x1(inst, vals, wts, r, L, rng))
        sig = [sig_real(x, inst.N) for x in X]
        l1, st1 = join_all(sig, inst.t1, seed=seed * 1000 + rd)
        n_y = 0; l2work = 0
        for q in l1:
            Ys = [refinements(inst, X[k][q[k]], allowed) for k in range(4)]
            sig2 = [sig_real(y0 - y1, inst.N) for (y0, y1) in Ys]
            l2work += len(sig2[0]) * len(sig2[1]) + len(sig2[2]) * len(sig2[3])
            l2, _ = join_all(sig2, inst.t0, seed=seed * 1000 + rd + 7, tbits=18)
            for h in l2:
                y = np.stack([np.stack([Ys[k][0][h[k]], Ys[k][1][h[k]]]) for k in range(4)])
                found.append(y); n_y += 1
            if len(found) >= max_y: break
        rec = dict(round=rd, list_sizes=[len(x) for x in X], l1_solutions=int(len(l1)), l1_probes_log2=round(math.log2(max(1, st1['probes'])), 1),
                   l2_pairs_log2=round(math.log2(max(1, l2work)), 1), ys=n_y, secs=round(time.time() - t0, 1))
        stats.append(rec)
        if log: print(json.dumps(rec), flush=True)
        if len(found) >= max_y: break
    return found, stats

def y_to_x_check(inst, y):
    """exact S1 check of y (4, 2, R): sum_k PAF(y_k) on Z_2q = T_J (4m - 8 at 0, -12 elsewhere), rows (2,0,0,0)"""
    q = inst.q; S = inst.S
    Y = np.zeros((4, 2 * q), np.int64); i = np.arange(2 * q)
    Y[:, :] = y[:, i % 2, S.oid[i % q]]
    P = np.rint(np.fft.ifft(np.abs(np.fft.fft(Y, axis=1)) ** 2, axis=1).real).astype(np.int64).sum(0)
    T = np.full(2 * q, -12); T[0] = 4 * inst.m - 8
    rows = sorted(np.abs(Y.sum(1)).tolist())
    return bool((P == T).all()) and rows == [0, 0, 0, 2]

# ---------------------------------------------------------------- stage 2
EIS = np.array([(1, 0), (0, 1), (-1, -1)], np.int64)          # w^j as a + b w, j = 0, 1, 2

def zopts(ya):
    """options for Z'_a on one orbit: list of (eisenstein (a, b), j); |y| = 3 -> [(0, 0), -1]"""
    if abs(ya) == 3: return [((0, 0), -1)]
    return [((ya * EIS[j, 0], ya * EIS[j, 1]), j) for j in range(3)]

class Stage2Block:
    """per block: S-classes (product of per-orbit S' options) and per-orbit D'-option tables"""
    def __init__(self, inst, yk):
        R = inst.R
        self.nS = np.zeros(R, np.int64); Sval = []; Dtab = []
        for i in range(R):
            d = {}
            for (z0, j0) in zopts(int(yk[0, i])):
                for (z1, j1) in zopts(int(yk[1, i])):
                    s = (z0[0] + z1[0], z0[1] + z1[1]); dd = (z0[0] - z1[0], z0[1] - z1[1])
                    d.setdefault(s, []).append((dd, j0, j1))
            keys = sorted(d); self.nS[i] = len(keys); Sval.append(keys); Dtab.append([d[s] for s in keys])
        self.Sval, self.Dtab = Sval, Dtab
        self.n_classes = int(np.prod(self.nS)); self.log2n = float(np.log2(self.nS).sum())
        # D-option tables as arrays: dA, dB, dj0, dj1 (R, maxS, 3), dcnt (R, maxS)
        mS = int(self.nS.max())
        self.dA = np.zeros((R, mS, 3), np.int64); self.dB = np.zeros((R, mS, 3), np.int64)
        self.dj0 = np.full((R, mS, 3), -1, np.int64); self.dj1 = np.full((R, mS, 3), -1, np.int64)
        self.dcnt = np.zeros((R, mS), np.int64)
        for i in range(R):
            for s, lst in enumerate(Dtab[i]):
                self.dcnt[i, s] = len(lst)
                for t, (dd, j0, j1) in enumerate(lst):
                    self.dA[i, s, t], self.dB[i, s, t] = dd; self.dj0[i, s, t] = j0; self.dj1[i, s, t] = j1
        self.sA = np.zeros((R, mS), np.int64); self.sB = np.zeros((R, mS), np.int64)
        for i in range(R):
            for s, val in enumerate(Sval[i]): self.sA[i, s], self.sB[i, s] = val
    def classes(self, idx=None):
        """S' orbit vectors (A, B coords) for class indices idx (mixed radix over orbits; all if None)"""
        n = self.n_classes if idx is None else len(idx)
        if idx is None: idx = np.arange(n, dtype=np.int64)
        digits = np.zeros((n, len(self.nS)), np.int64); rem = idx.copy()
        for i in range(len(self.nS) - 1, -1, -1):
            digits[:, i] = rem % self.nS[i]; rem //= self.nS[i]
        R = len(self.nS)
        return self.sA[np.arange(R)[None, :], digits], self.sB[np.arange(R)[None, :], digits], digits

@numba.njit(cache=True)
def _dsigs(digits_row, dA, dB, dcnt, N):
    """all D' orbit vectors for one S-class (digits_row (R,)), their signatures and option indices"""
    R = digits_row.shape[0]; H = N.shape[0]
    n = 1
    for i in range(R): n *= dcnt[i, digits_row[i]]
    A = np.zeros((n, R), np.int64); B = np.zeros((n, R), np.int64); T = np.zeros((n, R), np.int64)
    for l in range(n):
        rem = l
        for i in range(R - 1, -1, -1):
            c = dcnt[i, digits_row[i]]; t = rem % c; rem //= c
            A[l, i] = dA[i, digits_row[i], t]; B[l, i] = dB[i, digits_row[i], t]; T[l, i] = t
    return A, B, T

@numba.njit(cache=True)
def _sig_eis1(A, B, N):
    L, R = A.shape; H = N.shape[0]
    out = np.zeros((L, 2 * H - 1), np.int64)
    for l in range(L):
        for h in range(H):
            sa = 0; sb = 0
            for i in range(R):
                ai = A[l, i]; bi = B[l, i]
                if ai == 0 and bi == 0: continue
                for j in range(R):
                    n = N[h, i, j]
                    if n == 0: continue
                    e = A[l, j] - B[l, j]; f = -B[l, j]
                    sa += n * (ai * e - bi * f); sb += n * (ai * f + bi * e - bi * f)
            if h == 0: out[l, 0] = sa
            else:
                out[l, 2 * h - 1] = sa; out[l, 2 * h] = sb
    return out

def dcheck(blocks, digits4, inst, target):
    """D-level for one S-solution: returns list of (Topt per block) for D-solutions"""
    D = []
    for k in range(4):
        A, B, T = _dsigs(digits4[k], blocks[k].dA, blocks[k].dB, blocks[k].dcnt, inst.N)
        D.append((_sig_eis1(A, B, inst.N), T))
    # brute 2+2: AB sums dict
    ab = {}
    for a in range(len(D[0][0])):
        for b in range(len(D[1][0])):
            ab.setdefault((D[0][0][a] + D[1][0][b]).tobytes(), []).append((a, b))
    out = []
    for c in range(len(D[2][0])):
        for d in range(len(D[3][0])):
            need = (target - D[2][0][c] - D[3][0][d]).tobytes()
            for (a, b) in ab.get(need, []):
                out.append([D[0][1][a], D[1][1][b], D[2][1][c], D[3][1][d]])
    return out

def assemble(inst, y, blocks, digits4, Topt):
    """x (4, m) from y and the chosen (S-option digit, D-option index) per orbit"""
    R = inst.R; sym = np.zeros((4, 2, 3, R), np.int64)
    for k in range(4):
        for i in range(R):
            s = digits4[k][i]; t = Topt[k][i]
            j0 = blocks[k].dj0[i, s, t]; j1 = blocks[k].dj1[i, s, t]
            for a, j in ((0, j0), (1, j1)):
                ya = int(y[k, a, i])
                if abs(ya) == 3: sym[k, a, :, i] = ya // 3
                else:
                    sym[k, a, :, i] = ya; sym[k, a, j, i] = -ya
    return inst.S.from_orbits(sym)

def core_ok(x, m):
    P = np.rint(np.fft.ifft(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=1).real).astype(np.int64).sum(0)
    return bool(P[0] == 4 * m and (P[1:] == -4).all()) and sorted(np.abs(x.sum(1)).tolist()) == [0, 0, 0, 2]

def stage2(inst, y, sub=None, seed=0, tbits=22, maxhits=1 << 22, log=False, force=None):
    """exact stage 2 for one y: S-level join on (optionally subsampled) S-class lists, then the D-level.
    sub: per-block cap on the number of S-classes (random subset).  Returns (record, list of cores x)."""
    t0 = time.time(); rng = np.random.default_rng(seed)
    blocks = [Stage2Block(inst, y[k]) for k in range(4)]
    rec = dict(sclass_log2=[round(b.log2n, 1) for b in blocks], configs_log2=None)
    idx = []
    for b in blocks:
        n = b.n_classes
        if sub is not None and n > sub: ix = rng.choice(n, size=sub, replace=False).astype(np.int64)
        else: ix = np.arange(n, dtype=np.int64)
        idx.append(ix)
    if force is not None:                                    # test hook: make sure given class indices are present
        idx = [np.unique(np.concatenate([ix, np.array([f], np.int64)])) for ix, f in zip(idx, force)]
    sigs = []; digs = []
    for b, ix in zip(blocks, idx):
        A, B, dg = b.classes(ix); sigs.append(_sig_eis1(A, B, inst.N)); digs.append(dg)
    tS = np.zeros(sigs[0].shape[1], np.int64); tS[0] = inst.v
    order = sorted(range(4), key=lambda k: len(sigs[k]))
    perm = [order[3], order[0], order[2], order[1]]                  # balance the two sides
    hits, st = join_all([sigs[k] for k in perm], tS, seed=seed + 5, tbits=tbits, maxhits=maxhits)
    rec.update(used_log2=[round(math.log2(len(s)), 1) for s in sigs], s_solutions=int(len(hits)), s_join=st,
               s_secs=round(time.time() - t0, 1))
    cores = []; nd = 0
    for h in hits:
        dig4 = [None] * 4
        for pos, k in enumerate(perm): dig4[k] = digs[k][h[pos]]
        for Topt in dcheck(blocks, dig4, inst, tS):
            nd += 1
            x = assemble(inst, y, blocks, dig4, Topt)
            if core_ok(x, inst.m): cores.append(x)
    rec.update(d_solutions=nd, cores=len(cores), secs=round(time.time() - t0, 1))
    if log: print(json.dumps(rec), flush=True)
    return rec, cores


# ---------------------------------------------------------------- stage 1 with the theta_5 compression (exact)
def c0_mask(inst):
    """per orbit of Z_q (q = 5 q'): number of its elements with Z_5-coordinate 0 (t = 0 mod 5)"""
    assert inst.q % 5 == 0
    return np.array([sum(1 for t in o if t % 5 == 0) for o in inst.S.orbs], np.int64)

def l1_profiles(inst, rows=(2, 0, 0, 0)):
    """[P] all (A_0..A_3) with sum_k ((5 A_k - r_k)/4)^2 = 4v and A_k = r_k mod 4 (exact): A_k = sum of X1 over the
    c = 0 positions; theta_5 (rational: M maps onto Z_5^*) of the L1 level equals (5 A - r)/4, B = (r - A)/4."""
    c0 = c0_mask(inst); amax = 6 * int(c0.sum())
    T = 4 * inst.v; out = []
    sq = [[(A, ((5 * A - r) // 4) ** 2) for A in range(-amax, amax + 1) if (A - r) % 4 == 0] for r in rows]
    for A0, s0 in sq[0]:
        if s0 > T: continue
        for A1, s1 in sq[1]:
            if s0 + s1 > T: continue
            for A2, s2 in sq[2]:
                if s0 + s1 + s2 > T: continue
                rem = T - s0 - s1 - s2
                for A3, s3 in sq[3]:
                    if s3 == rem: out.append((A0, A1, A2, A3))
    return out

class CondSampler:
    """X1 vectors with fixed row sum r and c = 0 sum A (x-multiplicity weights); backward DP over (row, A)"""
    def __init__(self, inst, vals, wts):
        R = inst.R; sz = inst.sizes; c0 = c0_mask(inst)
        cnt = [None] * (R + 1); cnt[R] = {(0, 0): 1.0}
        for i in range(R - 1, -1, -1):
            E = {}
            for (t, u), c in cnt[i + 1].items():
                for x, w in zip(vals[i], wts[i]):
                    key = (t + sz[i] * x, u + c0[i] * x); E[key] = E.get(key, 0.0) + c * w
            cnt[i] = E
        self.cnt, self.vals, self.wts, self.sz, self.c0, self.R = cnt, vals, wts, sz, c0, R
    def total(self, r, A): return self.cnt[0].get((r, A), 0.0)
    def sample(self, r, A, L, rng):
        if self.total(r, A) == 0: return np.zeros((0, self.R), np.int64)
        out = set(); tries = 0
        last = 0
        while len(out) < L and tries < 20 * L and tries - last < max(4000, 4 * len(out)):
            tries += 1; t, u = r, A; row = []
            for i in range(self.R):
                ps = np.array([w * self.cnt[i + 1].get((t - self.sz[i] * x, u - self.c0[i] * x), 0.0)
                               for x, w in zip(self.vals[i], self.wts[i])])
                j = rng.choice(len(ps), p=ps / ps.sum()); x = int(self.vals[i][j])
                row.append(x); t -= self.sz[i] * x; u -= self.c0[i] * x
            n0 = len(out); out.add(tuple(row))
            if len(out) > n0: last = tries
        return np.array(sorted(out), np.int64)

def stage1c(inst, L=4096, n_profiles=10, seed=0, allowed=None, rows=(2, 0, 0, 0), log=False, max_y=10 ** 6, skip=0):
    """stage 1 over theta_5 profiles (largest weight first): conditioned L1 lists, L1 join, L2 refinement join"""
    rng = np.random.default_rng(seed)
    if allowed is None: allowed = np.ones((inst.R, 16), bool)
    vals, wts = x1_sampler(inst, allowed); cs = CondSampler(inst, vals, wts)
    profs = l1_profiles(inst, rows)
    tots = np.array([np.prod([cs.total(rows[k], P[k]) for k in range(4)]) for P in profs])
    order = [i for i in np.argsort(-tots) if tots[i] > 0]
    if log: print(json.dumps(dict(profiles=len(profs), nonzero=len(order),
                                  top=[(list(profs[i]), round(math.log2(tots[i]), 1)) for i in order[:5]])), flush=True)
    found = []
    for pi in order[skip:skip + n_profiles]:
        P = profs[pi]; t0 = time.time()
        X = [cs.sample(rows[k], P[k], L, rng) for k in range(4)]
        sig = [sig_real(x, inst.N) for x in X]
        l1, st1 = join_all(sig, inst.t1, seed=seed * 1000 + int(pi))
        n_y = 0
        for q in l1:
            Ys = [refinements(inst, X[k][q[k]], allowed) for k in range(4)]
            sig2 = [sig_real(y0 - y1, inst.N) for (y0, y1) in Ys]
            l2, _ = join_all(sig2, inst.t0, seed=seed * 1000 + int(pi) + 7, tbits=18)
            for h in l2:
                found.append(np.stack([np.stack([Ys[k][0][h[k]], Ys[k][1][h[k]]]) for k in range(4)])); n_y += 1
        if log: print(json.dumps(dict(profile=list(P), log2_weight=round(math.log2(tots[pi]), 1), lists=[len(x) for x in X],
                                      l1=int(len(l1)), ys=n_y, secs=round(time.time() - t0, 1))), flush=True)
        if len(found) >= max_y: break
    return found

# ---------------------------------------------------------------- stage 1, L2 conditioned on (R, A') profiles [P]
# For X = Xm (L2): R = row sum of Xm (alternating row sum of x), A' = sum over c = 0 positions; exact equations
#   theta_0: sum_k R_k^2 = 4v,   theta_5: sum_k ((5 A'_k - R_k)/4)^2 = 4v,   A'_k = R_k (mod 8),  R_k = 2 - r_k (mod 4).
class XmCond:
    """for one block with fixed X1 (orbit vector): weighted Xm options per orbit and a DP over (R, A')"""
    def __init__(self, inst, x1row, allowed):
        self.R = inst.R; sz = inst.sizes; c0 = c0_mask(inst) if inst.q % 5 == 0 else np.zeros(inst.R, np.int64)
        self.opts = []
        for i in range(inst.R):
            d = {}
            for p in np.nonzero(allowed[i])[0]:
                if YPAIRS[p].sum() != x1row[i]: continue
                xm = int(YPAIRS[p, 0] - YPAIRS[p, 1]); d[xm] = d.get(xm, 0) + 3 ** int((np.abs(YPAIRS[p]) == 1).sum())
            self.opts.append((np.array(sorted(d), np.int64), np.array([d[k] for k in sorted(d)], float)))
        cnt = [None] * (self.R + 1); cnt[self.R] = {(0, 0): 1.0}
        for i in range(self.R - 1, -1, -1):
            E = {}
            vals, wts = self.opts[i]
            for (t, u), c in cnt[i + 1].items():
                for x, w in zip(vals, wts):
                    key = (t + sz[i] * x, u + c0[i] * x); E[key] = E.get(key, 0.0) + c * w
            cnt[i] = E
        self.cnt, self.sz, self.c0 = cnt, sz, c0
    def profiles(self): return self.cnt[0]
    def sample(self, Rk, Ak, L, rng):
        if self.cnt[0].get((Rk, Ak), 0) == 0: return np.zeros((0, self.R), np.int64)
        out = set(); tries = 0
        last = 0
        while len(out) < L and tries < 20 * L and tries - last < max(4000, 4 * len(out)):
            tries += 1; t, u = Rk, Ak; row = []
            for i in range(self.R):
                vals, wts = self.opts[i]
                ps = np.array([w * self.cnt[i + 1].get((t - self.sz[i] * x, u - self.c0[i] * x), 0.0) for x, w in zip(vals, wts)])
                j = rng.choice(len(ps), p=ps / ps.sum()); x = int(vals[j]); row.append(x)
                t -= self.sz[i] * x; u -= self.c0[i] * x
            n0 = len(out); out.add(tuple(row))
            if len(out) > n0: last = tries
        return np.array(sorted(out), np.int64)

def l2_profiles(inst, conds, top=50):
    """quadruples ((R_k, A'_k))_k satisfying both exact equations, ranked by product of DP weights"""
    T = 4 * inst.v; per = []
    for c in conds:
        P = c.profiles()
        if inst.q % 5 == 0: per.append([(R, A, w, R * R, ((5 * A - R) // 4) ** 2) for (R, A), w in P.items() if (A - R) % 8 == 0 and (5 * A - R) % 4 == 0])
        else: per.append([(R, A, w, R * R, None) for (R, A), w in P.items()])
    has5 = inst.q % 5 == 0
    best = []
    # index block 3 by (s0, s5)
    idx3 = {}
    for e in per[3]: idx3.setdefault((e[3], e[4]), []).append(e)
    for e0 in per[0]:
        for e1 in per[1]:
            if e0[3] + e1[3] > T or (has5 and e0[4] + e1[4] > T): continue
            for e2 in per[2]:
                s0 = T - e0[3] - e1[3] - e2[3]
                if s0 < 0: continue
                s5 = (T - e0[4] - e1[4] - e2[4]) if has5 else None
                if has5 and s5 < 0: continue
                for e3 in idx3.get((s0, s5), []):
                    best.append((e0[2] * e1[2] * e2[2] * e3[2], [(e[0], e[1]) for e in (e0, e1, e2, e3)]))
    best.sort(key=lambda t: -t[0])
    return best[:top], len(best)

def stage1d(inst, L1=4096, L2=8192, n_profiles=10, n_l2prof=4, seed=0, allowed=None, rows=(2, 0, 0, 0), log=False,
            max_y=200, skip=0, max_l1=50):
    """stage 1 with exact compression profiles at both levels: theta_5 profile for L1 (when 5 | q), (R, A') profiles
    for L2; lists sampled conditioned on the profile; exact joins on the full level signatures"""
    rng = np.random.default_rng(seed)
    if allowed is None: allowed = np.ones((inst.R, 16), bool)
    vals, wts = x1_sampler(inst, allowed)
    has5 = inst.q % 5 == 0 and len({t % 5 for t in inst.S.Gq}) == 4
    if has5:
        cs = CondSampler(inst, vals, wts); profs = l1_profiles(inst, rows)
        tots = np.array([np.prod([cs.total(rows[k], P[k]) for k in range(4)]) for P in profs])
        order = [i for i in np.argsort(-tots) if tots[i] > 0]
    else:
        profs = [None]; order = [0]; tots = np.array([1.0])
    found = []
    for pi in order[skip:skip + n_profiles]:
        P = profs[pi]; t0 = time.time()
        if has5: X = [cs.sample(rows[k], P[k], L1, rng) for k in range(4)]
        else: X = [sample_x1(inst, vals, wts, rows[k], L1, rng) for k in range(4)]
        sig = [sig_real(x, inst.N) for x in X]
        l1, st1 = join_all(sig, inst.t1, seed=seed * 1000 + int(pi))
        n_y = 0; l2stats = []
        for q in l1[:max_l1]:
            conds = [XmCond(inst, X[k][q[k]], allowed) for k in range(4)]
            best, nprof = l2_profiles(inst, conds)
            for w, prof in best[:n_l2prof]:
                M_ = [conds[k].sample(prof[k][0], prof[k][1], L2, rng) for k in range(4)]
                if any(len(mm) == 0 for mm in M_): continue
                sig2 = [sig_real(mm, inst.N) for mm in M_]
                l2, _ = join_all(sig2, inst.t0, seed=seed * 1000 + int(pi) + 7, tbits=20)
                l2stats.append((nprof, [len(mm) for mm in M_], int(len(l2))))
                for h in l2:
                    y = np.zeros((4, 2, inst.R), np.int64)
                    for k in range(4):
                        x1 = X[k][q[k]]; xm = M_[k][h[k]]; y[k, 0] = (x1 + xm) // 2; y[k, 1] = (x1 - xm) // 2
                    found.append(y); n_y += 1
                if len(found) >= max_y: break
            if len(found) >= max_y: break
        if log: print(json.dumps(dict(profile=(list(P) if P else None), lists=[len(x) for x in X], l1=int(len(l1)),
                                      l2=l2stats[:6], ys=n_y, secs=round(time.time() - t0, 1))), flush=True)
        if len(found) >= max_y: break
    return found

def enumerate_x1(inst, vals, wts, r, limit=1 << 20):
    """all X1 orbit vectors with row sum r (exact DFS over the reachability DP); None if more than limit"""
    R = inst.R; sz = inst.sizes
    reach = [None] * (R + 1); reach[R] = {0}
    for i in range(R - 1, -1, -1): reach[i] = {t + sz[i] * x for t in reach[i + 1] for x in vals[i]}
    if r not in reach[0]: return np.zeros((0, R), np.int64)
    out = []; row = [0] * R
    def rec(i, t):
        if len(out) > limit: return
        if i == R:
            if t == 0: out.append(list(row))
            return
        for x in vals[i]:
            if t - sz[i] * x in reach[i + 1]:
                row[i] = int(x); rec(i + 1, t - sz[i] * x)
    rec(0, r)
    return None if len(out) > limit else np.array(out, np.int64)


# ---------------------------------------------------------------- stage 1 on exact Z10-compression profiles
YV = np.array([1, -1, 3, -3], np.int64)
YW = np.array([3, 3, 1, 1], np.float64)                     # x-multiplicity of each y value

def group_lists(weights, target):
    """all assignments of y in {+-1, +-3} to the orbits with nonzero weight w (sum w*y = target): (configs, weights)"""
    idx = [i for i, w in enumerate(weights) if w]
    if not idx: return (np.zeros((1, 0), np.int64), np.ones(1)) if target == 0 else (np.zeros((0, 0), np.int64), np.zeros(0))
    combos = np.array(np.meshgrid(*[np.arange(4)] * len(idx), indexing='ij')).reshape(len(idx), -1).T
    vals = YV[combos]; s = (vals * np.array([weights[i] for i in idx])[None, :]).sum(1)
    keep = s == target
    return vals[keep], YW[combos[keep]].prod(1)

class ProfileBlock:
    """per block and Z10 profile (Y00, Y0s, Y10, Y1s): the four exact factor lists (a in {0,1}) x (c = 0, c != 0)"""
    def __init__(self, inst, w0, w1, prof):
        Y00, Y0s, Y10, Y1s = prof
        self.R = inst.R
        self.g0 = [i for i in range(self.R) if w0[i]]; self.g1 = [i for i in range(self.R) if w1[i]]
        self.L = {(0, 0): group_lists(w0, Y00), (0, 1): group_lists(w1, Y0s),
                  (1, 0): group_lists(w0, Y10), (1, 1): group_lists(w1, Y1s)}
        self.ok = all(len(self.L[k][0]) for k in self.L)
        self.prof = prof; self.w0, self.w1 = w0, w1
    def sample(self, N, rng):
        """N weighted random y-configs (2, R)"""
        out = np.zeros((N, 2, self.R), np.int64)
        for a in (0, 1):
            for g, orb in ((0, self.g0), (1, self.g1)):
                V, W = self.L[(a, g)]
                pick = rng.choice(len(V), size=N, p=W / W.sum())
                out[:, a, orb] = V[pick]
        return out
    def refinements(self, x1):
        """all y (2, R) in this profile with y0 + y1 = x1: per group, y0 options from the a = 0 list filtered so that
        y1 = x1 - y0 is valid and lies in the a = 1 list"""
        parts = []
        for g, orb in ((0, self.g0), (1, self.g1)):
            V0, _ = self.L[(0, g)]; V1, _ = self.L[(1, g)]
            Y1 = x1[orb][None, :] - V0
            valid = np.isin(np.abs(Y1), (1, 3)).all(1)
            if not valid.any(): return np.zeros((0, 2, self.R), np.int64)
            # y1 must be one of the a = 1 configs (its group sum is then automatically right)
            set1 = {tuple(r) for r in V1}
            sel = [j for j in np.nonzero(valid)[0] if tuple(Y1[j]) in set1]
            parts.append((orb, V0[sel], Y1[sel]))
        n0, n1 = len(parts[0][1]), len(parts[1][1])
        out = np.zeros((n0 * n1, 2, self.R), np.int64)
        for (orb, A0, A1), rep, til in ((parts[0], n1, 1), (parts[1], 1, n0)):
            out[:, 0, orb] = np.tile(np.repeat(A0, rep, axis=0), (til, 1))
            out[:, 1, orb] = np.tile(np.repeat(A1, rep, axis=0), (til, 1))
        return out

def stage1z(inst, profs, L1=4096, n_profiles=4, seed=0, log=False, max_y=1000, max_l1=10 ** 6):
    """exact stage 1 inside Z10-compression profiles (theory/p571_z10profiles.py): L1 join on X1 classes of
    profile-consistent y's, then the L2 join over all profile-consistent refinements of each L1 solution"""
    from p571_z10profiles import fibre_weights
    w0, w1 = fibre_weights(inst); rng = np.random.default_rng(seed); found = []
    for P in profs[:n_profiles]:
        t0 = time.time()
        blocks = [ProfileBlock(inst, w0, w1, P[k]) for k in range(4)]
        if not all(b.ok for b in blocks): continue
        X = []
        for b in blocks:
            Ys = b.sample(4 * L1, rng); x1 = np.unique(Ys[:, 0] + Ys[:, 1], axis=0)
            X.append(x1[rng.permutation(len(x1))[:L1]])
        sig = [sig_real(x, inst.N) for x in X]
        l1, st1 = join_all(sig, inst.t1, seed=seed * 7919 + len(found))
        n_y = 0; l2pairs = 0
        for q in l1[:max_l1]:
            refs = [blocks[k].refinements(X[k][q[k]]) for k in range(4)]
            if any(len(r) == 0 for r in refs): continue
            sig2 = [sig_real(r[:, 0] - r[:, 1], inst.N) for r in refs]
            l2pairs += len(refs[0]) * len(refs[1]) + len(refs[2]) * len(refs[3])
            l2, _ = join_all(sig2, inst.t0, seed=seed + 11, tbits=18)
            for h in l2:
                found.append(np.stack([refs[k][h[k]] for k in range(4)])); n_y += 1
            if len(found) >= max_y: break
        if log: print(json.dumps(dict(profile=[list(p) for p in P], lists=[len(x) for x in X], l1=int(len(l1)),
                                      l2_pairs_log2=round(math.log2(max(1, l2pairs)), 1), ys=n_y, secs=round(time.time() - t0, 1))), flush=True)
        if len(found) >= max_y: break
    return found

if __name__ == "__main__" and sys.argv[1] in ("stage1",):
    cmd = sys.argv[1]; a = sys.argv[2:]
    opt = {"--L": "8192", "--rounds": "1", "--seed": "0", "--out": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    inst = Inst(pos[0], pos[1:])
    if cmd == 'stage1':
        ys, st = stage1(inst, L=int(opt['--L']), rounds=int(opt['--rounds']), seed=int(opt['--seed']), log=True)
        ok = [y_to_x_check(inst, y) for y in ys]
        print(json.dumps(dict(total_ys=len(ys), all_S1_exact=bool(all(ok)) if ys else None)), flush=True)
        if opt['--out'] and ys: np.savez_compressed(opt['--out'], Y=np.array(ys))

if __name__ == "__main__" and sys.argv[1] in ('stage1c', 'stage1d', 'stage1z'):
    a = sys.argv[2:]; opt = {"--L": "4096", "--np": "10", "--seed": "0", "--out": "", "--skip": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    inst = Inst(pos[0], pos[1:])
    if sys.argv[1] == 'stage1z':
        from p571_z10profiles import profiles as zprof
        sols, _ = zprof(inst); sols.sort(key=lambda t: -t[1]); rng0 = np.random.default_rng(int(opt['--seed']))
        sk = int(opt['--skip']); profs = [p for p, w in sols[sk:]]
        ys = stage1z(inst, profs, L1=int(opt['--L']), n_profiles=int(opt['--np']), seed=int(opt['--seed']), log=True)
        ok = [y_to_x_check(inst, y) for y in ys]
        print(json.dumps(dict(total_ys=len(ys), all_S1_exact=bool(all(ok)) if ys else None)), flush=True)
        if opt['--out'] and ys: np.savez_compressed(opt['--out'], Y=np.array(ys))
    if sys.argv[1] == 'stage1d':
        ys = stage1d(inst, L1=int(opt['--L']), n_profiles=int(opt['--np']), seed=int(opt['--seed']), log=True, skip=int(opt['--skip']))
        ok = [y_to_x_check(inst, y) for y in ys]
        print(json.dumps(dict(total_ys=len(ys), all_S1_exact=bool(all(ok)) if ys else None)), flush=True)
        if opt['--out'] and ys: np.savez_compressed(opt['--out'], Y=np.array(ys))
    if sys.argv[1] == 'stage1c':
        ys = stage1c(inst, L=int(opt['--L']), n_profiles=int(opt['--np']), seed=int(opt['--seed']), log=True,
                     skip=int(opt['--skip']))
        ok = [y_to_x_check(inst, y) for y in ys]
        print(json.dumps(dict(total_ys=len(ys), all_S1_exact=bool(all(ok)) if ys else None)), flush=True)
        if opt['--out'] and ys: np.savez_compressed(opt['--out'], Y=np.array(ys))
