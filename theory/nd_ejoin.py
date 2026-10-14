r"""Exact 4-list join with an energy split and numba hash tables [M tool].

Problem: four lists with integer signature matrices S_k (N_k x L); column 0 is a nonnegative 'energy'.  Find every
(i0, i1, i2, i3) with S_0[i0] + S_1[i1] + S_2[i2] + S_3[i3] = T (all L coordinates).
Method: the lists are sorted and bucketed by energy.  For every left energy EL = e_a + e_b, the left pairs go into an
open-addressing hash table keyed by a random 64-bit linear form of the remaining L - 1 coordinates (mod 2^64), and the
right pairs (energy T[0] - EL) are streamed against it; every key match is verified on the exact integer signature.
Large tables are split into P passes by the top bits of the key.  Optional energy-only restrictions (used for
symmetry reductions): e_b <= (T[0] - EL) * bmax_frac, and e_c <= e_d on the right.
Validated by plants (theory/nd_ejoin.py --selftest) against a brute force on small random lists."""
import sys, os, json, math, time
import numpy as np
from numba import njit

EMPTY = np.uint64(0)

def hashes(S, w):
    """64-bit linear form of the non-energy coordinates (wraps mod 2^64)"""
    return (S[:, 1:].astype(np.uint64) * w[None, :]).sum(1, dtype=np.uint64)

def buckets(e, emax):
    """e sorted ascending; start/end index of every energy value 0..emax"""
    st = np.searchsorted(e, np.arange(emax + 1), 'left').astype(np.int64)
    en = np.searchsorted(e, np.arange(emax + 1), 'right').astype(np.int64)
    return st, en

@njit(cache=True)
def _pass(ha, sa, ea, hb, sb, eb, hc, sc, ec, hd, sd, ed, EL, ER, emax, bmax, cled, hT, P, p, bits,
          out, maxout):
    size = np.int64(1) << bits; mask = np.uint64(size - 1)
    tk = np.zeros(size, np.uint64); tv = np.empty(size, np.int64)
    shift = np.uint64(58); Pu = np.uint64(P); pu = np.uint64(p); one = np.uint64(1)
    for e_a in range(0, EL + 1):
        e_b = EL - e_a
        if e_a > emax or e_b > emax or e_b > bmax: continue
        for i in range(sa[e_a], ea[e_a]):
            hai = ha[i]
            for j in range(sb[e_b], eb[e_b]):
                k = hai + hb[j]
                if (k >> shift) % Pu != pu: continue
                if k == EMPTY: k = one
                s = k & mask
                while tk[s] != EMPTY: s = (s + one) & mask
                tk[s] = k; tv[s] = i * np.int64(1 << 31) + j
    n = 0
    for e_c in range(0, ER + 1):
        e_d = ER - e_c
        if e_c > emax or e_d > emax: continue
        if cled and e_c > e_d: continue
        for i in range(sc[e_c], ec[e_c]):
            r0 = hT - hc[i]
            for j in range(sd[e_d], ed[e_d]):
                k = r0 - hd[j]
                if (k >> shift) % Pu != pu: continue
                if k == EMPTY: k = one
                s = k & mask
                while tk[s] != EMPTY:
                    if tk[s] == k:
                        if n < maxout:
                            out[n, 0] = tv[s] >> 31; out[n, 1] = tv[s] & ((1 << 31) - 1); out[n, 2] = i; out[n, 3] = j
                        n += 1
                    s = (s + one) & mask
    return n

def hash_weights(L, seed=11):
    return np.random.default_rng(seed).integers(1, 2 ** 62, L - 1, dtype=np.uint64) * np.uint64(2) + np.uint64(1)

def prepare(s_, E, w):
    """one list: kept rows (energy <= E) sorted by energy, their hashes and energy buckets"""
    keep = s_[:, 0] <= E
    idx = np.nonzero(keep)[0]; s2 = s_[idx]; o = np.argsort(s2[:, 0], kind='stable')
    idx = idx[o]; s2 = s2[o]
    st, en = buckets(s2[:, 0], E)
    return dict(idx=idx, S=s2, h=hashes(s2, w), st=st, en=en, cnt=np.bincount(s2[:, 0], minlength=E + 1)[:E + 1])

def ejoin4(S, T, order=None, bmax_frac=None, cled=False, max_table_bits=26, log=None, maxout=1 << 24, seed=11, els=None,
           done_cb=None, prepared=None):
    """S: list of 4 signature arrays (column 0 = energy >= 0); T: target.  order = (a, b, c, d) pairing (left a, b;
    right c, d); default: the pairing with the smallest estimated work.  Returns sorted list of (i0, i1, i2, i3) in the
    caller's list numbering, every one verified exactly, and a stats dict.  prepared: the four prepare() dicts (built
    with hash_weights(len(T), seed) and E = T[0]) to share them between calls or forked workers."""
    T = np.asarray(T, np.int64); E = int(T[0])
    w = hash_weights(len(T), seed)
    emax = E
    if prepared is None:
        for s_ in S: assert s_.shape[1] == len(T) and (s_[:, 0] >= 0).all()
        lists = [prepare(s_, E, w) for s_ in S]
    else:
        lists = prepared
    def work(a, b, c, d):
        ca, cb, cc, cd = [lists[x]['cnt'].astype(float) for x in (a, b, c, d)]
        L = np.convolve(ca, cb)[:E + 1]; R = np.convolve(cc, cd)[:E + 1]
        both = (L > 0) & (R[::-1] > 0)
        return float((L[both]).sum() + (R[::-1][both]).sum()), L, R
    if order is None:
        best = None
        for o4 in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2), (2, 3, 0, 1), (1, 3, 0, 2), (1, 2, 0, 3)):
            wk = work(*o4)[0]
            if best is None or wk < best[0]: best = (wk, o4)
        order = best[1]
    a, b, c, d = order
    wk, Lc, Rc = work(a, b, c, d)
    hT = np.uint64((T[1:].astype(np.uint64) * w).sum(dtype=np.uint64))
    A, B, C, D = [lists[x] for x in order]
    sols = set(); stats = dict(order=list(order), work_log2=round(math.log2(wk), 2) if wk > 0 else None, passes=0, key_matches=0)
    t0 = time.time()
    for EL in (range(E + 1) if els is None else els):
        ER = E - EL
        if Lc[EL] == 0 or Rc[ER] == 0:
            if done_cb: done_cb(EL, [])
            continue
        bmax = emax if bmax_frac is None else int(math.floor(ER * bmax_frac))
        # left size with the bmax restriction
        nL = sum(int(A['cnt'][ea]) * int(B['cnt'][EL - ea]) for ea in range(EL + 1) if EL - ea <= bmax)
        if nL == 0:
            if done_cb: done_cb(EL, [])
            continue
        P = 1
        while nL / P * 2 > (1 << max_table_bits) and P < 64: P *= 2
        bits = max(4, int(math.ceil(math.log2(nL / P * 2 + 1))))
        found = []
        for p in range(P):
            out = np.zeros((min(maxout, 1 << 20), 4), np.int64)
            n = _pass(A['h'], A['st'], A['en'], B['h'], B['st'], B['en'], C['h'], C['st'], C['en'], D['h'], D['st'], D['en'],
                      EL, ER, emax, bmax, cled, hT, P, p, bits, out, len(out))
            if n > len(out):
                if n > maxout: raise RuntimeError(f"too many key matches ({n}); raise maxout")
                out = np.zeros((n, 4), np.int64)
                n = _pass(A['h'], A['st'], A['en'], B['h'], B['st'], B['en'], C['h'], C['st'], C['en'], D['h'], D['st'], D['en'],
                          EL, ER, emax, bmax, cled, hT, P, p, bits, out, len(out))
            stats['passes'] += 1; stats['key_matches'] += int(n)
            for r in out[:n]:
                ia, ib, ic, id_ = int(r[0]), int(r[1]), int(r[2]), int(r[3])
                if (A['S'][ia] + B['S'][ib] + C['S'][ic] + D['S'][id_] == T).all():
                    t = [0] * 4
                    t[a], t[b], t[c], t[d] = int(A['idx'][ia]), int(B['idx'][ib]), int(C['idx'][ic]), int(D['idx'][id_])
                    sols.add(tuple(t)); found.append(tuple(t))
        if done_cb: done_cb(EL, sorted(set(found)))
        if log: log(json.dumps(dict(EL=EL, left=nL, passes=P, sols=len(sols), secs=round(time.time() - t0, 1))))
    stats['secs'] = round(time.time() - t0, 2); stats['solutions'] = len(sols)
    return sorted(sols), stats

def _selftest(trials=30):
    rng = np.random.default_rng(1); ok = 0
    for t in range(trials):
        L = int(rng.integers(2, 5)); N = [int(rng.integers(5, 40)) for _ in range(4)]
        S = []
        for n in N:
            s_ = rng.integers(-3, 4, (n, L)); s_[:, 0] = rng.integers(0, 12, n); S.append(s_)
        pick = [int(rng.integers(n)) for n in N]
        T = sum(S[k][pick[k]] for k in range(4))
        sols, st = ejoin4(S, T, max_table_bits=6)
        brute = sorted((i, j, k, l) for i in range(N[0]) for j in range(N[1]) for k in range(N[2]) for l in range(N[3])
                       if (S[0][i] + S[1][j] + S[2][k] + S[3][l] == T).all())
        ok += (sols == brute and tuple(pick) in set(sols))
    return ok, trials

if __name__ == "__main__":
    if "--selftest" in sys.argv:
        print(json.dumps(dict(selftest=_selftest())))
