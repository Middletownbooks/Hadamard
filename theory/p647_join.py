r"""Exact 4-list join, energy split + residue passes + radix partitions (the 647 x353 S-join kernel) [M tool].

Problem (as theory/nd_ejoin.py): four integer signature matrices S_k (N_k x L, column 0 an 'energy' >= 0) and a
target T; find every (i0, i1, i2, i3) with S_0[i0] + S_1[i1] + S_2[i2] + S_3[i3] = T on all L coordinates.

Method.  Pairing (a, b | c, d).  Every row gets a 64-bit hash h (random odd weights on columns 1..L-1, wrapping) and
a residue r = (random weights . columns 1..L-1) mod P, P prime; both are additive, so a solution has
    e_a + e_b = EL, e_c + e_d = T0 - EL,  r_a + r_b = r_T - r_c - r_d = p (mod P),  h_a + h_b = h_T - h_c - h_d.
For every energy EL and every group of G consecutive residues p the left pairs (a, b) and right pairs (c, d) of
that class are emitted (the inner lists b, d are sorted by (energy, residue) with a CSR index, so the partners of a
are one or two contiguous ranges), scattered into 2^rbits partitions by the top key bits, and every partition is
joined in cache (open-addressing table on the smaller side).  The exact pair count of every (EL, p) is precomputed
by a 2-D FFT (linear in energy, cyclic in residue), so passes have a known size.  Every key match is re-checked on
the exact integer signature.  Work = number of emitted left + right pairs (each pair once).
Tests: python3 theory/p647_join.py --selftest  (random instances against brute force and nd_ejoin.ejoin4)."""
import sys, os, json, math, time
import numpy as np
from numba import njit

MIX = np.uint64(0x9E3779B97F4A7C15)

def hash_w(L, seed):
    return np.random.default_rng(seed).integers(1, 2 ** 62, L - 1, dtype=np.uint64) * np.uint64(2) + np.uint64(1)

def res_w(L, P, seed):
    return np.random.default_rng(seed + 1000).integers(0, P, L - 1).astype(np.int64)

def next_prime(n):
    def isp(k):
        if k < 2: return False
        i = 2
        while i * i <= k:
            if k % i == 0: return False
            i += 1
        return True
    while not isp(n): n += 1
    return n

class SL:
    """one list: rows with energy <= E, sorted by (energy, residue); hashes, residues, energy starts, CSR (inner)"""
    def __init__(self, S, E, w, u, P):
        S = np.asarray(S, np.int64)
        idx = np.nonzero(S[:, 0] <= E)[0]
        e = S[idx, 0]; r = (S[idx, 1:] * u).sum(1) % P
        o = np.lexsort((r, e)); idx = idx[o]; e = e[o]; r = r[o]
        self.idx = idx; self.S = S[idx]; self.e = e; self.r = r.astype(np.int64)
        self.h = (self.S[:, 1:].astype(np.uint64) * w).sum(1, dtype=np.uint64)
        self.E = E; self.P = P
        comb = e * P + r
        self.csr = np.searchsorted(comb, np.arange((E + 1) * P + 1)).astype(np.int64)   # flat: start of (e, r)
        self.st = self.csr[np.arange(E + 1) * P]; self.en = self.csr[np.arange(1, E + 2) * P]
        self.hist = np.zeros((E + 1, P), np.float64)
        np.add.at(self.hist, (e, r), 1.0)
        self.emax = int(e.max()) if len(e) else -1

def pair_counts(X, Y, E, P):
    """cnt[EL, p] = #{(x, y): e_x + e_y = EL, r_x + r_y = p mod P}, EL = 0..E (exact, by FFT)"""
    n = 1 << int(math.ceil(math.log2(2 * (E + 1))))
    FX = np.fft.fft2(X.hist, s=(n, P)); FY = np.fft.fft2(Y.hist, s=(n, P))
    c = np.fft.ifft2(FX * FY).real[:E + 1]
    return np.rint(c).astype(np.int64)

@njit(cache=True)
def _emit(hX, rX, stX, enX, hY, csrY, EL, emaxX, emaxY, P, p0, G, right, hT, rT, sh, pos, endp, outk, outi):
    ok = True
    lo = EL - emaxY
    if lo < 0: lo = 0
    hi = EL if EL < emaxX else emaxX
    for ex in range(lo, hi + 1):
        ey = EL - ex
        if stX[ex] == enX[ex]: continue
        base = ey * P
        if csrY[base] == csrY[base + P]: continue
        for x in range(stX[ex], enX[ex]):
            if right:
                s0 = (rT - rX[x] - p0 - G + 1) % P
                if s0 < 0: s0 += P
                kx = hT - hX[x]
            else:
                s0 = (p0 - rX[x]) % P
                if s0 < 0: s0 += P
                kx = hX[x]
            s1 = s0 + G
            if s1 <= P:
                a0 = csrY[base + s0]; a1 = csrY[base + s1]; b0 = 0; b1 = 0
            else:
                a0 = csrY[base + s0]; a1 = csrY[base + P]; b0 = csrY[base]; b1 = csrY[base + s1 - P]
            xi = np.int64(x) << 32
            for t in range(a0, a1):
                key = kx - hY[t] if right else kx + hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    outk[i] = key; outi[i] = xi | np.int64(t); pos[part] = i + 1
                else:
                    ok = False
            for t in range(b0, b1):
                key = kx - hY[t] if right else kx + hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    outk[i] = key; outi[i] = xi | np.int64(t); pos[part] = i + 1
                else:
                    ok = False
    return ok

@njit(cache=True)
def _emit_count(hX, rX, stX, enX, hY, csrY, EL, emaxX, emaxY, P, p0, G, right, hT, rT, sh, cnt):
    lo = EL - emaxY
    if lo < 0: lo = 0
    hi = EL if EL < emaxX else emaxX
    for ex in range(lo, hi + 1):
        ey = EL - ex
        if stX[ex] == enX[ex]: continue
        base = ey * P
        if csrY[base] == csrY[base + P]: continue
        for x in range(stX[ex], enX[ex]):
            if right:
                s0 = (rT - rX[x] - p0 - G + 1) % P
                if s0 < 0: s0 += P
                kx = hT - hX[x]
            else:
                s0 = (p0 - rX[x]) % P
                if s0 < 0: s0 += P
                kx = hX[x]
            s1 = s0 + G
            if s1 <= P:
                a0 = csrY[base + s0]; a1 = csrY[base + s1]; b0 = 0; b1 = 0
            else:
                a0 = csrY[base + s0]; a1 = csrY[base + P]; b0 = csrY[base]; b1 = csrY[base + s1 - P]
            for t in range(a0, a1):
                key = kx - hY[t] if right else kx + hY[t]
                cnt[np.int64(key >> sh)] += 1
            for t in range(b0, b1):
                key = kx - hY[t] if right else kx + hY[t]
                cnt[np.int64(key >> sh)] += 1

@njit(cache=True)
def _join_parts(Lk, Li, stL, cntL, Rk, Ri, stR, cntR, nparts, tk, tv, tbits, hits, nh, maxhits):
    tsize = np.int64(1) << tbits; tmask = tsize - 1; shs = np.uint64(64 - tbits)
    used = np.empty(tsize, np.int64)
    for i in range(nparts):
        nl = cntL[i]; nr = cntR[i]
        if nl == 0 or nr == 0: continue
        if nl <= nr:
            Sk = Lk; Si = Li; s0 = stL[i]; s1 = stL[i] + nl; Bk = Rk; Bi = Ri; b0 = stR[i]; b1 = stR[i] + nr; sl = True
        else:
            Sk = Rk; Si = Ri; s0 = stR[i]; s1 = stR[i] + nr; Bk = Lk; Bi = Li; b0 = stL[i]; b1 = stL[i] + nl; sl = False
        nu = 0
        for j in range(s0, s1):
            key = Sk[j]; slot = np.int64((key * MIX) >> shs)
            while tv[slot] >= 0: slot = (slot + 1) & tmask
            tk[slot] = key; tv[slot] = j; used[nu] = slot; nu += 1
        for j in range(b0, b1):
            key = Bk[j]; slot = np.int64((key * MIX) >> shs)
            while tv[slot] >= 0:
                if tk[slot] == key:
                    if nh < maxhits:
                        if sl:
                            hits[nh, 0] = Si[tv[slot]]; hits[nh, 1] = Bi[j]
                        else:
                            hits[nh, 0] = Bi[j]; hits[nh, 1] = Si[tv[slot]]
                    nh += 1
                slot = (slot + 1) & tmask
        for u in range(nu): tv[used[u]] = -1
    return nh

@njit(cache=True)
def _verify(hits, n, SA, SB, SC, SD, T, out):
    """unpack (left packed, right packed) hits and keep those with an exact signature match; returns count"""
    L = T.shape[0]; m = 0
    for q in range(n):
        x = hits[q, 0] >> 32; y = hits[q, 0] & 0xFFFFFFFF; c = hits[q, 1] >> 32; d = hits[q, 1] & 0xFFFFFFFF
        good = True
        for j in range(L):
            if SA[x, j] + SB[y, j] + SC[c, j] + SD[d, j] != T[j]:
                good = False; break
        if good:
            out[m, 0] = x; out[m, 1] = y; out[m, 2] = c; out[m, 3] = d; m += 1
    return m

PAIRINGS = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))

def plan_work(S, T, P=4093, seed=11):
    """exact left/right pair counts (energy compatible) for the three pairings, log2"""
    T = np.asarray(T, np.int64); E = int(T[0]); w = hash_w(len(T), seed); u = res_w(len(T), P, seed)
    Ls = [SL(s_, E, w, u, P) for s_ in S]
    out = {}
    for a, b, c, d in PAIRINGS:
        cL = pair_counts(Ls[a], Ls[b], E, P).sum(1); cR = pair_counts(Ls[c], Ls[d], E, P).sum(1)
        both = (cL > 0) & (cR[::-1] > 0)
        out[(a, b, c, d)] = (float(cL[both].sum()), float(cR[::-1][both].sum()))
    return out

def join4(S, T, order=None, P=4093, rbits=9, budget=1 << 23, seed=11, log=None, els=None, done_cb=None,
          maxhits=1 << 16, lists=None, stats_only=False):
    """all (i0, i1, i2, i3) with sum_k S_k[i_k] = T, in the caller's numbering (sorted, exact); and stats.
    order: (a, b, c, d) (default: smallest emitted work).  els: subset of left energies (restart support);
    done_cb(EL, sols_of_EL) after each finished EL."""
    t0 = time.time()
    T = np.asarray(T, np.int64); E = int(T[0]); L = len(T)
    w = hash_w(L, seed); u = res_w(L, P, seed)
    hT = np.uint64((T[1:].astype(np.uint64) * w).sum(dtype=np.uint64)); rT = int((T[1:] * u).sum() % P)
    Ls = lists if lists is not None else [SL(s_, E, w, u, P) for s_ in S]
    if order is None:
        best = None
        for o4 in PAIRINGS:
            a, b, c, d = o4
            cL = pair_counts(Ls[a], Ls[b], E, P).sum(1); cR = pair_counts(Ls[c], Ls[d], E, P).sum(1)
            both = (cL > 0) & (cR[::-1] > 0); wk = cL[both].sum() + cR[::-1][both].sum()
            if best is None or wk < best[0]: best = (wk, o4)
        order = best[1]
    a, b, c, d = order
    assert rbits >= 1                                  # key >> 64 is undefined
    A, B, C, D = Ls[a], Ls[b], Ls[c], Ls[d]
    cL = pair_counts(A, B, E, P)                       # [EL, p]
    cRq = pair_counts(C, D, E, P)                      # [ER, q], q = r_c + r_d
    # right residue p = rT - q  ->  cR[ER, p] = cRq[ER, (rT - p) mod P]
    cR = cRq[:, (rT - np.arange(P)) % P]
    nparts = 1 << rbits; sh = np.uint64(64 - rbits)
    tbits = max(8, int(math.ceil(math.log2(max(16, 2.6 * budget / nparts)))))
    tk = np.zeros(1 << tbits, np.uint64); tv = -np.ones(1 << tbits, np.int64)
    cap = int(budget * 1.25) + nparts * 64
    Lk = np.zeros(cap, np.uint64); Li = np.zeros(cap, np.int64); Rk = np.zeros(cap, np.uint64); Ri = np.zeros(cap, np.int64)
    hits = np.zeros((maxhits, 2), np.int64)
    stats = dict(order=list(order), P=P, rbits=rbits, budget=budget, passes=0, left=0, right=0, key_matches=0, overflow=0)
    wL = 0; wR = 0
    for EL in range(E + 1):
        if cL[EL].sum() and cR[E - EL].sum(): wL += int(cL[EL].sum()); wR += int(cR[E - EL].sum())
    stats['work_log2'] = round(math.log2(wL + wR), 3) if wL + wR else None
    if stats_only: return [], stats
    sols = []
    ELs = range(E + 1) if els is None else els
    for EL in ELs:
        ER = E - EL
        rowL = cL[EL]; rowR = cR[ER]
        if rowL.sum() == 0 or rowR.sum() == 0:
            if done_cb: done_cb(EL, [])
            continue
        found = []
        p = 0
        mx = np.maximum(rowL, rowR); csum = np.concatenate([[0], np.cumsum(mx)])
        while p < P:
            # largest G with max(left, right) of the group <= budget (at least 1)
            g = int(np.searchsorted(csum, csum[p] + budget, 'right')) - 1 - p
            g = max(1, min(g, P - p))
            nL = int(rowL[p:p + g].sum()); nR = int(rowR[p:p + g].sum())
            if nL == 0 or nR == 0: p += g; continue
            need = max(nL, nR)
            if need * 1.25 + nparts * 64 > len(Lk):
                cap = int(need * 1.25) + nparts * 64
                Lk = np.zeros(cap, np.uint64); Li = np.zeros(cap, np.int64); Rk = np.zeros(cap, np.uint64); Ri = np.zeros(cap, np.int64)
            slL = int(nL / nparts * 1.15) + 64; slR = int(nR / nparts * 1.15) + 64
            posL = np.arange(nparts, dtype=np.int64) * slL; endL = posL + slL
            posR = np.arange(nparts, dtype=np.int64) * slR; endR = posR + slR
            stL = posL.copy(); stR = posR.copy()
            okL = _emit(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, sh, posL, endL, Lk, Li)
            okR = _emit(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, sh, posR, endR, Rk, Ri)
            if not (okL and okR):
                stats['overflow'] += 1
                cntA = np.zeros(nparts, np.int64); cntB = np.zeros(nparts, np.int64)
                _emit_count(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, sh, cntA)
                _emit_count(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, sh, cntB)
                stL = np.concatenate([[0], np.cumsum(cntA)[:-1]]).astype(np.int64); stR = np.concatenate([[0], np.cumsum(cntB)[:-1]]).astype(np.int64)
                posL = stL.copy(); endL = stL + cntA; posR = stR.copy(); endR = stR + cntB
                _emit(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, sh, posL, endL, Lk, Li)
                _emit(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, sh, posR, endR, Rk, Ri)
            cntL_ = posL - stL; cntR_ = posR - stR
            assert int(cntL_.sum()) == nL and int(cntR_.sum()) == nR, (int(cntL_.sum()), nL, int(cntR_.sum()), nR)
            small = int(np.minimum(cntL_, cntR_).max())
            if 2 * small + 2 > len(tk):              # the table must hold the smaller side of every partition
                tbits = int(math.ceil(math.log2(2.6 * small + 4)))
                tk = np.zeros(1 << tbits, np.uint64); tv = -np.ones(1 << tbits, np.int64)
            while True:
                nh = _join_parts(Lk, Li, stL, cntL_, Rk, Ri, stR, cntR_, nparts, tk, tv, tbits, hits, 0, len(hits))
                if nh <= len(hits): break
                hits = np.zeros((2 * nh, 2), np.int64)
            stats['passes'] += 1; stats['left'] += nL; stats['right'] += nR; stats['key_matches'] += int(nh)
            if nh:
                out = np.zeros((nh, 4), np.int64)
                m = _verify(hits, nh, A.S, B.S, C.S, D.S, T, out)
                for q in range(m):
                    t = [0] * 4
                    t[a] = int(A.idx[out[q, 0]]); t[b] = int(B.idx[out[q, 1]]); t[c] = int(C.idx[out[q, 2]]); t[d] = int(D.idx[out[q, 3]])
                    found.append(tuple(t))
            p += g
        sols.extend(found)
        if done_cb: done_cb(EL, sorted(set(found)))
        if log: log(json.dumps(dict(EL=EL, left=int(rowL.sum()), right=int(rowR.sum()), sols=len(sols), secs=round(time.time() - t0, 1))))
    stats['secs'] = round(time.time() - t0, 2); stats['solutions'] = len(set(sols))
    tot = stats['left'] + stats['right']
    stats['ns_per_pair'] = round(1e9 * stats['secs'] / tot, 2) if tot else None
    return sorted(set(sols)), stats

def _brute(S, T):
    N = [len(s_) for s_ in S]
    AB = (S[0][:, None, :] + S[1][None, :, :]).reshape(-1, len(T))
    CD = (S[2][:, None, :] + S[3][None, :, :]).reshape(-1, len(T))
    d = {}
    for i, v in enumerate(map(tuple, AB)): d.setdefault(v, []).append(i)
    out = []
    for j, v in enumerate(map(tuple, T[None, :] - CD)):
        for i in d.get(v, []):
            out.append((i // N[1], i % N[1], j // N[3], j % N[3]))
    return sorted(out)

def selftest(trials=40):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from nd_ejoin import ejoin4
    rng = np.random.default_rng(3); ok = 0; okE = 0
    for t in range(trials):
        L = int(rng.integers(2, 6)); N = [int(rng.integers(5, 120)) for _ in range(4)]
        S = []
        for n in N:
            s_ = rng.integers(-2, 3, (n, L)); s_[:, 0] = rng.integers(0, 15, n); S.append(s_)
        pick = [int(rng.integers(n)) for n in N]
        T = sum(S[k][pick[k]] for k in range(4))
        P = int(rng.choice([2, 3, 7, 31, 101]))
        sols, st = join4(S, T, P=P, rbits=int(rng.integers(1, 6)), budget=int(rng.choice([16, 256, 4096])), seed=t)
        br = _brute(S, T)
        ok += (sols == br and tuple(pick) in set(sols))
        e4, _ = ejoin4(S, T, max_table_bits=8)
        okE += (sorted(e4) == br)
    return dict(selftest=f"{ok}/{trials}", ejoin4_agree=f"{okE}/{trials}")


# ---------------------------------------------------------------- packed 8-byte variant (production)
# entry = (fingerprint << 40) | (x << 20) | y; fingerprint = key bits just below the partition bits (24 bits); a
# fingerprint match inside a partition is confirmed on the full 64-bit key (h arrays), then on the exact signature.
IB = 20
IMASK = (1 << IB) - 1

@njit(cache=True)
def _emit8(hX, rX, stX, enX, hY, csrY, EL, emaxX, emaxY, P, p0, G, right, hT, rT, rbits, pos, endp, out):
    ok = True
    sh = np.uint64(64 - rbits); shf = np.uint64(40 - rbits)
    lo = EL - emaxY
    if lo < 0: lo = 0
    hi = EL if EL < emaxX else emaxX
    for ex in range(lo, hi + 1):
        ey = EL - ex
        if stX[ex] == enX[ex]: continue
        base = ey * P
        if csrY[base] == csrY[base + P]: continue
        for x in range(stX[ex], enX[ex]):
            if right:
                s0 = (rT - rX[x] - p0 - G + 1) % P
                if s0 < 0: s0 += P
                kx = hT - hX[x]
            else:
                s0 = (p0 - rX[x]) % P
                if s0 < 0: s0 += P
                kx = hX[x]
            s1 = s0 + G
            if s1 <= P:
                a0 = csrY[base + s0]; a1 = csrY[base + s1]; b0 = 0; b1 = 0
            else:
                a0 = csrY[base + s0]; a1 = csrY[base + P]; b0 = csrY[base]; b1 = csrY[base + s1 - P]
            xi = np.uint64(x) << np.uint64(IB)
            for t in range(a0, a1):
                key = kx - hY[t] if right else kx + hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    out[i] = ((key << np.uint64(rbits)) >> np.uint64(40)) << np.uint64(40) | xi | np.uint64(t); pos[part] = i + 1
                else:
                    ok = False
            for t in range(b0, b1):
                key = kx - hY[t] if right else kx + hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    out[i] = ((key << np.uint64(rbits)) >> np.uint64(40)) << np.uint64(40) | xi | np.uint64(t); pos[part] = i + 1
                else:
                    ok = False
    return ok

@njit(cache=True)
def _join8(Lo, stL, cntL, Ro, stR, cntR, nparts, tab, tbits, hA, hB, hC, hD, hT, hits, nh, maxhits):
    """exact 64-bit key matches between left entries (a, b) and right entries (c, d) of every partition"""
    tsize = np.int64(1) << tbits; tmask = np.uint64(tsize - 1)
    EMPTY = np.uint64(0xFFFFFFFFFFFFFFFF); s40 = np.uint64(40); im = np.uint64(IMASK); ib = np.uint64(IB)
    nfp = 0
    for i in range(nparts):
        nl = cntL[i]; nr = cntR[i]
        if nl == 0 or nr == 0: continue
        if nl <= nr:
            S = Lo; s0 = stL[i]; s1 = stL[i] + nl; B = Ro; b0 = stR[i]; b1 = stR[i] + nr; sl = True
        else:
            S = Ro; s0 = stR[i]; s1 = stR[i] + nr; B = Lo; b0 = stL[i]; b1 = stL[i] + nl; sl = False
        for j in range(s0, s1):
            e = S[j]; slot = (e >> s40) & tmask
            while tab[slot] != EMPTY: slot = (slot + np.uint64(1)) & tmask
            tab[slot] = e
        for j in range(b0, b1):
            e = B[j]; fp = e >> s40; slot = fp & tmask
            while tab[slot] != EMPTY:
                f = tab[slot]
                if (f >> s40) == fp:
                    nfp += 1
                    if sl:
                        x = (f >> ib) & im; y = f & im; c = (e >> ib) & im; d = e & im
                    else:
                        x = (e >> ib) & im; y = e & im; c = (f >> ib) & im; d = f & im
                    if hA[x] + hB[y] == hT - hC[c] - hD[d]:
                        if nh < maxhits:
                            hits[nh, 0] = np.int64(x); hits[nh, 1] = np.int64(y); hits[nh, 2] = np.int64(c); hits[nh, 3] = np.int64(d)
                        nh += 1
                slot = (slot + np.uint64(1)) & tmask
        for j in range(s0, s1):
            e = S[j]; slot = (e >> s40) & tmask
            while tab[slot] != EMPTY:
                tab[slot] = EMPTY; slot = (slot + np.uint64(1)) & tmask
    return nh, nfp

@njit(cache=True)
def _verify4(hits, n, SA, SB, SC, SD, T, out):
    L = T.shape[0]; m = 0
    for q in range(n):
        x = hits[q, 0]; y = hits[q, 1]; c = hits[q, 2]; d = hits[q, 3]
        good = True
        for j in range(L):
            if SA[x, j] + SB[y, j] + SC[c, j] + SD[d, j] != T[j]:
                good = False; break
        if good:
            out[m, 0] = x; out[m, 1] = y; out[m, 2] = c; out[m, 3] = d; m += 1
    return m

def join4p(S, T, order=None, P=4093, rbits=8, budget=1 << 22, seed=11, log=None, els=None, done_cb=None,
           maxhits=1 << 16, lists=None, stats_only=False, el_cb=None):
    """join4 with packed 8-byte entries (same contract; lists must have < 2^20 rows each).
    el_cb(EL, stats_of_EL) after each EL (pairs, key matches, solutions, seconds) for per-class accounting."""
    t0 = time.time()
    T = np.asarray(T, np.int64); E = int(T[0]); L = len(T)
    w = hash_w(L, seed); u = res_w(L, P, seed)
    hT = np.uint64((T[1:].astype(np.uint64) * w).sum(dtype=np.uint64)); rT = int((T[1:] * u).sum() % P)
    Ls = lists if lists is not None else [SL(s_, E, w, u, P) for s_ in S]
    assert max(len(l.idx) for l in Ls) <= IMASK + 1 and 1 <= rbits <= 16
    if order is None:
        best = None
        for o4 in PAIRINGS:
            a, b, c, d = o4
            cL = pair_counts(Ls[a], Ls[b], E, P).sum(1); cR = pair_counts(Ls[c], Ls[d], E, P).sum(1)
            both = (cL > 0) & (cR[::-1] > 0); wk = cL[both].sum() + cR[::-1][both].sum()
            if best is None or wk < best[0]: best = (wk, o4)
        order = best[1]
    a, b, c, d = order
    A, B, C, D = Ls[a], Ls[b], Ls[c], Ls[d]
    cL = pair_counts(A, B, E, P); cR = pair_counts(C, D, E, P)[:, (rT - np.arange(P)) % P]
    nparts = 1 << rbits
    tab = np.full(1 << 15, np.uint64(0xFFFFFFFFFFFFFFFF)); tbits = 15
    cap = int(budget * 1.25) + nparts * 64
    Lo = np.zeros(cap, np.uint64); Ro = np.zeros(cap, np.uint64)
    hits = np.zeros((maxhits, 4), np.int64)
    stats = dict(order=list(order), P=P, rbits=rbits, budget=budget, passes=0, left=0, right=0, fp_matches=0,
                 key_matches=0, overflow=0)
    wL = 0; wR = 0; per_el = {}
    for EL in range(E + 1):
        if cL[EL].sum() and cR[E - EL].sum(): wL += int(cL[EL].sum()); wR += int(cR[E - EL].sum())
    stats['work_log2'] = round(math.log2(wL + wR), 3) if wL + wR else None
    stats['left_total'] = wL; stats['right_total'] = wR
    if stats_only: return [], stats
    sols = []
    for EL in (range(E + 1) if els is None else els):
        ER = E - EL; te = time.time()
        rowL = cL[EL]; rowR = cR[ER]
        if rowL.sum() == 0 or rowR.sum() == 0:
            if done_cb: done_cb(EL, [])
            continue
        found = []; kmEL = 0
        p = 0
        mx = np.maximum(rowL, rowR); csum = np.concatenate([[0], np.cumsum(mx)])
        while p < P:
            g = int(np.searchsorted(csum, csum[p] + budget, 'right')) - 1 - p
            g = max(1, min(g, P - p))
            nL = int(rowL[p:p + g].sum()); nR = int(rowR[p:p + g].sum())
            if nL == 0 or nR == 0: p += g; continue
            need = max(nL, nR)
            if need * 1.25 + nparts * 64 > len(Lo):
                cap = int(need * 1.25) + nparts * 64; Lo = np.zeros(cap, np.uint64); Ro = np.zeros(cap, np.uint64)
            slL = int(nL / nparts * 1.15) + 64; slR = int(nR / nparts * 1.15) + 64
            posL = np.arange(nparts, dtype=np.int64) * slL; endL = posL + slL; stL = posL.copy()
            posR = np.arange(nparts, dtype=np.int64) * slR; endR = posR + slR; stR = posR.copy()
            okL = _emit8(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, rbits, posL, endL, Lo)
            okR = _emit8(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, rbits, posR, endR, Ro)
            if not (okL and okR):
                stats['overflow'] += 1
                sh = np.uint64(64 - rbits)
                cntA = np.zeros(nparts, np.int64); cntB = np.zeros(nparts, np.int64)
                _emit_count(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, sh, cntA)
                _emit_count(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, sh, cntB)
                stL = np.concatenate([[0], np.cumsum(cntA)[:-1]]).astype(np.int64); stR = np.concatenate([[0], np.cumsum(cntB)[:-1]]).astype(np.int64)
                posL = stL.copy(); endL = stL + cntA; posR = stR.copy(); endR = stR + cntB
                _emit8(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, rbits, posL, endL, Lo)
                _emit8(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, rbits, posR, endR, Ro)
            cntL_ = posL - stL; cntR_ = posR - stR
            assert int(cntL_.sum()) == nL and int(cntR_.sum()) == nR
            small = int(np.minimum(cntL_, cntR_).max())
            if 2 * small + 2 > len(tab):
                tbits = int(math.ceil(math.log2(2.6 * small + 4))); tab = np.full(1 << tbits, np.uint64(0xFFFFFFFFFFFFFFFF))
            while True:
                nh, nfp = _join8(Lo, stL, cntL_, Ro, stR, cntR_, nparts, tab, tbits, A.h, B.h, C.h, D.h, hT, hits, 0, len(hits))
                if nh <= len(hits): break
                hits = np.zeros((2 * nh, 4), np.int64)
            stats['passes'] += 1; stats['left'] += nL; stats['right'] += nR; stats['fp_matches'] += int(nfp); stats['key_matches'] += int(nh)
            kmEL += int(nh)
            if nh:
                out = np.zeros((nh, 4), np.int64)
                m = _verify4(hits, nh, A.S, B.S, C.S, D.S, T, out)
                for q in range(m):
                    t = [0] * 4
                    t[a] = int(A.idx[out[q, 0]]); t[b] = int(B.idx[out[q, 1]]); t[c] = int(C.idx[out[q, 2]]); t[d] = int(D.idx[out[q, 3]])
                    found.append(tuple(t))
            p += g
        sols.extend(found)
        if el_cb: el_cb(EL, dict(left=int(rowL.sum()), right=int(rowR.sum()), key_matches=kmEL, nsol=len(set(found)),
                                 secs=round(time.time() - te, 2)))
        if done_cb: done_cb(EL, sorted(set(found)))
        if log: log(json.dumps(dict(EL=EL, left=int(rowL.sum()), right=int(rowR.sum()), sols=len(sols), secs=round(time.time() - t0, 1))))
    stats['secs'] = round(time.time() - t0, 2); stats['solutions'] = len(set(sols))
    tot = stats['left'] + stats['right']
    stats['ns_per_pair'] = round(1e9 * stats['secs'] / tot, 2) if tot else None
    return sorted(set(sols)), stats

def selftest_p(trials=40):
    rng = np.random.default_rng(4); ok = 0
    for t in range(trials):
        L = int(rng.integers(2, 6)); N = [int(rng.integers(5, 150)) for _ in range(4)]
        S = []
        for n in N:
            s_ = rng.integers(-2, 3, (n, L)); s_[:, 0] = rng.integers(0, 15, n); S.append(s_)
        pick = [int(rng.integers(n)) for n in N]
        T = sum(S[k][pick[k]] for k in range(4))
        P = int(rng.choice([2, 3, 7, 31, 101]))
        sols, st = join4p(S, T, P=P, rbits=int(rng.integers(1, 6)), budget=int(rng.choice([16, 256, 4096])), seed=t)
        ok += (sols == _brute(S, T) and tuple(pick) in set(sols))
    return dict(selftest_packed=f"{ok}/{trials}")


# ---------------------------------------------------------------- optimised variant (join4q): branch-free emission,
# no per-x modulo, bitset-prefiltered partition join.  Same contract and same hits as join4p.
@njit(cache=True)
def _emitL(hX, rX, stX, enX, hY, csrY, EL, emaxX, emaxY, P, p0, G, rbits, pos, endp, out):
    ok = True
    sh = np.uint64(64 - rbits); ub = np.uint64(rbits); s40 = np.uint64(40)
    lo = EL - emaxY
    if lo < 0: lo = 0
    hi = EL if EL < emaxX else emaxX
    for ex in range(lo, hi + 1):
        ey = EL - ex
        if stX[ex] == enX[ex]: continue
        base = ey * P
        if csrY[base] == csrY[base + P]: continue
        for x in range(stX[ex], enX[ex]):
            s0 = p0 - rX[x]
            if s0 < 0: s0 += P
            kx = hX[x]; xi = np.uint64(x) << np.uint64(IB)
            s1 = s0 + G
            if s1 <= P:
                a0 = csrY[base + s0]; a1 = csrY[base + s1]; b0 = 0; b1 = 0
            else:
                a0 = csrY[base + s0]; a1 = csrY[base + P]; b0 = csrY[base]; b1 = csrY[base + s1 - P]
            for t in range(a0, a1):
                key = kx + hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    out[i] = ((key << ub) >> s40) << s40 | xi | np.uint64(t); pos[part] = i + 1
                else:
                    ok = False
            for t in range(b0, b1):
                key = kx + hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    out[i] = ((key << ub) >> s40) << s40 | xi | np.uint64(t); pos[part] = i + 1
                else:
                    ok = False
    return ok

@njit(cache=True)
def _emitR(hX, rX, stX, enX, hY, csrY, EL, emaxX, emaxY, P, p0, G, hT, rT, rbits, pos, endp, out):
    ok = True
    sh = np.uint64(64 - rbits); ub = np.uint64(rbits); s40 = np.uint64(40)
    lo = EL - emaxY
    if lo < 0: lo = 0
    hi = EL if EL < emaxX else emaxX
    c0 = rT - p0 - G + 1                     # s0 = (c0 - r_x) mod P, c0 in (-2P, P)
    while c0 < 0: c0 += P
    for ex in range(lo, hi + 1):
        ey = EL - ex
        if stX[ex] == enX[ex]: continue
        base = ey * P
        if csrY[base] == csrY[base + P]: continue
        for x in range(stX[ex], enX[ex]):
            s0 = c0 - rX[x]
            if s0 < 0: s0 += P
            kx = hT - hX[x]; xi = np.uint64(x) << np.uint64(IB)
            s1 = s0 + G
            if s1 <= P:
                a0 = csrY[base + s0]; a1 = csrY[base + s1]; b0 = 0; b1 = 0
            else:
                a0 = csrY[base + s0]; a1 = csrY[base + P]; b0 = csrY[base]; b1 = csrY[base + s1 - P]
            for t in range(a0, a1):
                key = kx - hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    out[i] = ((key << ub) >> s40) << s40 | xi | np.uint64(t); pos[part] = i + 1
                else:
                    ok = False
            for t in range(b0, b1):
                key = kx - hY[t]
                part = np.int64(key >> sh); i = pos[part]
                if i < endp[part]:
                    out[i] = ((key << ub) >> s40) << s40 | xi | np.uint64(t); pos[part] = i + 1
                else:
                    ok = False
    return ok

@njit(cache=True)
def _join8b(Lo, stL, cntL, Ro, stR, cntR, nparts, tab, tbits, bs, bbits, hA, hB, hC, hD, hT, hits, nh, maxhits):
    """as _join8, with an L1 bitset over the smaller side's fingerprints: the larger side probes the table only
    on a bitset hit"""
    tsize = np.int64(1) << tbits; tmask = np.uint64(tsize - 1)
    EMPTY = np.uint64(0xFFFFFFFFFFFFFFFF); s40 = np.uint64(40); im = np.uint64(IMASK); ib = np.uint64(IB)
    bmask = np.uint64((np.int64(1) << bbits) - 1); s6 = np.uint64(6); m63 = np.uint64(63); one = np.uint64(1)
    nfp = 0
    for i in range(nparts):
        nl = cntL[i]; nr = cntR[i]
        if nl == 0 or nr == 0: continue
        if nl <= nr:
            S = Lo; s0 = stL[i]; s1 = stL[i] + nl; B = Ro; b0 = stR[i]; b1 = stR[i] + nr; sl = True
        else:
            S = Ro; s0 = stR[i]; s1 = stR[i] + nr; B = Lo; b0 = stL[i]; b1 = stL[i] + nl; sl = False
        for j in range(s0, s1):
            e = S[j]; fp = e >> s40; slot = fp & tmask
            while tab[slot] != EMPTY: slot = (slot + one) & tmask
            tab[slot] = e
            h = fp & bmask; bs[h >> s6] |= one << (h & m63)
        for j in range(b0, b1):
            e = B[j]; fp = e >> s40
            h = fp & bmask
            if not ((bs[h >> s6] >> (h & m63)) & one): continue
            slot = fp & tmask
            while tab[slot] != EMPTY:
                f = tab[slot]
                if (f >> s40) == fp:
                    nfp += 1
                    if sl:
                        x = (f >> ib) & im; y = f & im; c = (e >> ib) & im; d = e & im
                    else:
                        x = (e >> ib) & im; y = e & im; c = (f >> ib) & im; d = f & im
                    if hA[x] + hB[y] == hT - hC[c] - hD[d]:
                        if nh < maxhits:
                            hits[nh, 0] = np.int64(x); hits[nh, 1] = np.int64(y); hits[nh, 2] = np.int64(c); hits[nh, 3] = np.int64(d)
                        nh += 1
                slot = (slot + one) & tmask
        for j in range(s0, s1):
            e = S[j]; fp = e >> s40; slot = fp & tmask
            while tab[slot] != EMPTY:
                tab[slot] = EMPTY; slot = (slot + one) & tmask
            h = fp & bmask; bs[h >> s6] = np.uint64(0)
    return nh, nfp

def join4q(S, T, order=None, P=4093, rbits=8, budget=1 << 22, seed=11, log=None, els=None, done_cb=None,
           maxhits=1 << 16, lists=None, stats_only=False, el_cb=None, bbits=18):
    """join4p with the optimised kernels (same contract, same hits)"""
    t0 = time.time()
    T = np.asarray(T, np.int64); E = int(T[0]); L = len(T)
    w = hash_w(L, seed); u = res_w(L, P, seed)
    hT = np.uint64((T[1:].astype(np.uint64) * w).sum(dtype=np.uint64)); rT = int((T[1:] * u).sum() % P)
    Ls = lists if lists is not None else [SL(s_, E, w, u, P) for s_ in S]
    assert max(len(l.idx) for l in Ls) <= IMASK + 1 and 1 <= rbits <= 16
    if order is None:
        best = None
        for o4 in PAIRINGS:
            a, b, c, d = o4
            cL = pair_counts(Ls[a], Ls[b], E, P).sum(1); cR = pair_counts(Ls[c], Ls[d], E, P).sum(1)
            both = (cL > 0) & (cR[::-1] > 0); wk = cL[both].sum() + cR[::-1][both].sum()
            if best is None or wk < best[0]: best = (wk, o4)
        order = best[1]
    a, b, c, d = order
    A, B, C, D = Ls[a], Ls[b], Ls[c], Ls[d]
    cL = pair_counts(A, B, E, P); cR = pair_counts(C, D, E, P)[:, (rT - np.arange(P)) % P]
    nparts = 1 << rbits
    tab = np.full(1 << 15, np.uint64(0xFFFFFFFFFFFFFFFF)); tbits = 15
    bs = np.zeros(1 << (bbits - 6), np.uint64)
    cap = int(budget * 1.25) + nparts * 64
    Lo = np.zeros(cap, np.uint64); Ro = np.zeros(cap, np.uint64)
    hits = np.zeros((maxhits, 4), np.int64)
    stats = dict(order=list(order), P=P, rbits=rbits, budget=budget, passes=0, left=0, right=0, fp_matches=0,
                 key_matches=0, overflow=0, kernel="q")
    wL = 0; wR = 0
    for EL in range(E + 1):
        if cL[EL].sum() and cR[E - EL].sum(): wL += int(cL[EL].sum()); wR += int(cR[E - EL].sum())
    stats['work_log2'] = round(math.log2(wL + wR), 3) if wL + wR else None
    stats['left_total'] = wL; stats['right_total'] = wR
    if stats_only: return [], stats
    sols = []
    for EL in (range(E + 1) if els is None else els):
        ER = E - EL; te = time.time()
        rowL = cL[EL]; rowR = cR[ER]
        if rowL.sum() == 0 or rowR.sum() == 0:
            if el_cb: el_cb(EL, dict(left=0, right=0, key_matches=0, nsol=0, secs=0.0))
            if done_cb: done_cb(EL, [])
            continue
        found = []; kmEL = 0
        p = 0
        mx = np.maximum(rowL, rowR); csum = np.concatenate([[0], np.cumsum(mx)])
        while p < P:
            g = int(np.searchsorted(csum, csum[p] + budget, 'right')) - 1 - p
            g = max(1, min(g, P - p))
            nL = int(rowL[p:p + g].sum()); nR = int(rowR[p:p + g].sum())
            if nL == 0 or nR == 0: p += g; continue
            need = max(nL, nR)
            if need * 1.25 + nparts * 64 > len(Lo):
                cap = int(need * 1.25) + nparts * 64; Lo = np.zeros(cap, np.uint64); Ro = np.zeros(cap, np.uint64)
            slL = int(nL / nparts * 1.15) + 64; slR = int(nR / nparts * 1.15) + 64
            posL = np.arange(nparts, dtype=np.int64) * slL; endL = posL + slL; stL = posL.copy()
            posR = np.arange(nparts, dtype=np.int64) * slR; endR = posR + slR; stR = posR.copy()
            okL = _emitL(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, rbits, posL, endL, Lo)
            okR = _emitR(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, hT, rT, rbits, posR, endR, Ro)
            if not (okL and okR):
                stats['overflow'] += 1
                sh = np.uint64(64 - rbits)
                cntA = np.zeros(nparts, np.int64); cntB = np.zeros(nparts, np.int64)
                _emit_count(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, False, hT, rT, sh, cntA)
                _emit_count(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, True, hT, rT, sh, cntB)
                stL = np.concatenate([[0], np.cumsum(cntA)[:-1]]).astype(np.int64); stR = np.concatenate([[0], np.cumsum(cntB)[:-1]]).astype(np.int64)
                posL = stL.copy(); endL = stL + cntA; posR = stR.copy(); endR = stR + cntB
                _emitL(A.h, A.r, A.st, A.en, B.h, B.csr, EL, A.emax, B.emax, P, p, g, rbits, posL, endL, Lo)
                _emitR(C.h, C.r, C.st, C.en, D.h, D.csr, ER, C.emax, D.emax, P, p, g, hT, rT, rbits, posR, endR, Ro)
            cntL_ = posL - stL; cntR_ = posR - stR
            assert int(cntL_.sum()) == nL and int(cntR_.sum()) == nR
            small = int(np.minimum(cntL_, cntR_).max())
            if 2 * small + 2 > len(tab):
                tbits = int(math.ceil(math.log2(2.6 * small + 4))); tab = np.full(1 << tbits, np.uint64(0xFFFFFFFFFFFFFFFF))
            while True:
                nh, nfp = _join8b(Lo, stL, cntL_, Ro, stR, cntR_, nparts, tab, tbits, bs, bbits, A.h, B.h, C.h, D.h, hT, hits, 0, len(hits))
                if nh <= len(hits): break
                hits = np.zeros((2 * nh, 4), np.int64)
            stats['passes'] += 1; stats['left'] += nL; stats['right'] += nR; stats['fp_matches'] += int(nfp); stats['key_matches'] += int(nh)
            kmEL += int(nh)
            if nh:
                out = np.zeros((nh, 4), np.int64)
                m = _verify4(hits, nh, A.S, B.S, C.S, D.S, T, out)
                for q in range(m):
                    t = [0] * 4
                    t[a] = int(A.idx[out[q, 0]]); t[b] = int(B.idx[out[q, 1]]); t[c] = int(C.idx[out[q, 2]]); t[d] = int(D.idx[out[q, 3]])
                    found.append(tuple(t))
            p += g
        sols.extend(found)
        if el_cb: el_cb(EL, dict(left=int(rowL.sum()), right=int(rowR.sum()), key_matches=kmEL, nsol=len(set(found)),
                                 secs=round(time.time() - te, 2)))
        if done_cb: done_cb(EL, sorted(set(found)))
        if log: log(json.dumps(dict(EL=EL, left=int(rowL.sum()), right=int(rowR.sum()), sols=len(sols), secs=round(time.time() - t0, 1))))
    stats['secs'] = round(time.time() - t0, 2); stats['solutions'] = len(set(sols))
    tot = stats['left'] + stats['right']
    stats['ns_per_pair'] = round(1e9 * stats['secs'] / tot, 2) if tot else None
    return sorted(set(sols)), stats

def selftest_q(trials=40):
    rng = np.random.default_rng(5); ok = 0
    for t in range(trials):
        L = int(rng.integers(2, 6)); N = [int(rng.integers(5, 150)) for _ in range(4)]
        S = []
        for n in N:
            s_ = rng.integers(-2, 3, (n, L)); s_[:, 0] = rng.integers(0, 15, n); S.append(s_)
        pick = [int(rng.integers(n)) for n in N]
        T = sum(S[k][pick[k]] for k in range(4))
        P = int(rng.choice([2, 3, 7, 31, 101]))
        sols, st = join4q(S, T, P=P, rbits=int(rng.integers(1, 6)), budget=int(rng.choice([16, 256, 4096])), seed=t,
                          bbits=int(rng.choice([6, 10, 18])))
        ok += (sols == _brute(S, T) and tuple(pick) in set(sols))
    return dict(selftest_q=f"{ok}/{trials}")

if __name__ == "__main__":
    if "--selftest" in sys.argv:
        print(json.dumps(dict(**selftest(), **selftest_p(), **selftest_q())), flush=True)
