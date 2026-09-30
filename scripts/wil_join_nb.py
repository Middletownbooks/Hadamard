"""Compiled stage-D join. Pairs failing PSD_1(f)+PSD_2(f) <= 4n at some f != 0 are skipped (exact: the four PSDs
are nonnegative and sum to 4n). Otherwise the (same exact MITM as wil_search.join4). Per lift: 64-bit linear hash of its half PAF
(h = P @ W mod 2^64) and its -1 mask on positions 1..(n-1)/2. A Williamson quadruple satisfies
(mA^mB)*K1 + (hA+hB) == (mC^mD^flip)*K1 - (hC+hD), so every true match collides; the caller re-verifies each hit."""
import numpy as np, numba as nb


@nb.njit(cache=True)
def _prep(F2):
    nf = F2.shape[1]; n2 = F2.shape[0]
    order = np.empty((nf, n2), np.int64); vals = np.empty((nf, n2), np.float64)
    for f in range(nf):
        o = np.argsort(F2[:, f]); order[f] = o
        for t in range(n2): vals[f, t] = F2[o[t], f]
    return order, vals

@nb.njit(cache=True)
def _partners(F1x, F2, order, vals, bound, buf):
    """indices y with F1x[f] + F2[y,f] <= bound for all f; returns count written to buf"""
    nf = F2.shape[1]; best = -1; bcnt = F2.shape[0] + 1
    for f in range(nf):
        c = np.searchsorted(vals[f], bound - F1x[f], side="right")
        if c < bcnt: bcnt = c; best = f
    k = 0
    for t in range(bcnt):
        y = order[best, t]; ok = True
        for f in range(nf):
            if F1x[f] + F2[y, f] > bound:
                ok = False; break
        if ok:
            buf[k] = y; k += 1
    return k

@nb.njit(cache=True)
def _join(h1, m1, F1, h2, m2, F2, f12, neg12, h3, m3, F3, h4, m4, F4, f34, neg34, K1, bound, kb, xb, yb, out, maxout):
    # pass 1: surviving (1,2) pairs -> key buffers
    o2, v2 = _prep(F2); buf = np.empty(h2.shape[0], np.int64); nk = 0
    for x in range(h1.shape[0]):
        kk = _partners(F1[x], F2, o2, v2, bound, buf)
        for t in range(kk):
            y = buf[t]
            hs = h1[x] + h2[y]
            if neg12: hs = np.uint64(0) - hs
            kb[nk] = ((m1[x] ^ m2[y] ^ f12) * K1) + hs; xb[nk] = x; yb[nk] = y; nk += 1
    # hash table sized to survivors
    size = 1
    while size < 2 * nk + 2: size <<= 1
    mask = size - 1
    tk = np.empty(size, np.uint64); ti = np.full(size, -1, np.int64)
    for r in range(nk):
        k = kb[r]
        s = np.int64((k ^ (k >> np.uint64(29))) & np.uint64(mask))
        while ti[s] != -1: s = (s + 1) & mask
        tk[s] = k; ti[s] = r
    # pass 2: probe with surviving (3,4) pairs
    o4, v4 = _prep(F4); buf4 = np.empty(h4.shape[0], np.int64); nout = 0
    for x in range(h3.shape[0]):
        kk = _partners(F3[x], F4, o4, v4, bound, buf4)
        for t in range(kk):
            y = buf4[t]
            hs = h3[x] + h4[y]
            if neg34: hs = np.uint64(0) - hs
            k = ((m3[x] ^ m4[y] ^ f34) * K1) + hs
            s = np.int64((k ^ (k >> np.uint64(29))) & np.uint64(mask))
            while ti[s] != -1:
                if tk[s] == k:
                    if nout < maxout:
                        r = ti[s]
                        out[nout, 0] = xb[r]; out[nout, 1] = yb[r]; out[nout, 2] = x; out[nout, 3] = y
                    nout += 1
                s = (s + 1) & mask
    return nout

class Joiner:
    def __init__(s): s.cap = 0
    def _ensure(s, need):
        if need > s.cap:
            s.cap = need; s.kb = np.empty(need, np.uint64); s.xb = np.empty(need, np.int32 if False else np.int64); s.yb = np.empty(need, np.int64)
    def run(s, hA, mA, FA, hB, mB, FB, hC, mC, FC, hD, mD, FD, flip, K1, bound, maxout=1 << 16):
        """returns candidate index quadruples (a,b,c,d)"""
        z = np.uint64(0); f = np.uint64(flip); K1 = np.uint64(K1)
        out = np.empty((maxout, 4), np.int64)
        if len(hA) * len(hB) <= len(hC) * len(hD):
            s._ensure(max(1, len(hA) * len(hB)))
            n = _join(hA, mA, FA, hB, mB, FB, z, False, hC, mC, FC, hD, mD, FD, f, True, K1, bound, s.kb, s.xb, s.yb, out, maxout)
            r = out[:min(n, maxout)]
        else:
            s._ensure(max(1, len(hC) * len(hD)))
            n = _join(hC, mC, FC, hD, mD, FD, f, True, hA, mA, FA, hB, mB, FB, z, False, K1, bound, s.kb, s.xb, s.yb, out, maxout)
            r = out[:min(n, maxout)][:, [2, 3, 0, 1]]
        assert n <= maxout, "candidate overflow"
        return r
