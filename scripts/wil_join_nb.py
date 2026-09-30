"""Compiled stage-D join (same exact MITM as wil_search.join4). Per lift: 64-bit linear hash of its half PAF
(h = P @ W mod 2^64) and its -1 mask on positions 1..(n-1)/2. A Williamson quadruple satisfies
(mA^mB)*K1 + (hA+hB) == (mC^mD^flip)*K1 - (hC+hD), so every true match collides; the caller re-verifies each hit."""
import numpy as np, numba as nb

@nb.njit(cache=True)
def _join(h1, m1, h2, m2, f12, neg12, h3, m3, h4, m4, f34, neg34, K1, tk, ti, tj, out, maxout):
    # build table on (1,2) side, probe with (3,4) side
    size = tk.shape[0]; mask = size - 1
    for s in range(size): ti[s] = -1
    for x in range(h1.shape[0]):
        for y in range(h2.shape[0]):
            hs = h1[x] + h2[y]
            if neg12: hs = np.uint64(0) - hs
            k = ((m1[x] ^ m2[y] ^ f12) * K1) + hs
            s = np.int64((k ^ (k >> np.uint64(29))) & np.uint64(mask))
            while ti[s] != -1: s = (s + 1) & mask
            tk[s] = k; ti[s] = x; tj[s] = y
    nout = 0
    for x in range(h3.shape[0]):
        for y in range(h4.shape[0]):
            hs = h3[x] + h4[y]
            if neg34: hs = np.uint64(0) - hs
            k = ((m3[x] ^ m4[y] ^ f34) * K1) + hs
            s = np.int64((k ^ (k >> np.uint64(29))) & np.uint64(mask))
            while ti[s] != -1:
                if tk[s] == k:
                    if nout < maxout:
                        out[nout, 0] = ti[s]; out[nout, 1] = tj[s]; out[nout, 2] = x; out[nout, 3] = y
                    nout += 1
                s = (s + 1) & mask
    return nout

class Joiner:
    def __init__(s): s.size = 0
    def _ensure(s, need):
        size = 1
        while size < 2 * need: size <<= 1
        if size > s.size:
            s.size = size; s.tk = np.empty(size, np.uint64); s.ti = np.empty(size, np.int32); s.tj = np.empty(size, np.int32)
    def run(s, hA, mA, hB, mB, hC, mC, hD, mD, flip, K1, maxout=1 << 16):
        """returns candidate index quadruples (a,b,c,d)"""
        z = np.uint64(0); f = np.uint64(flip); K1 = np.uint64(K1)
        out = np.empty((maxout, 4), np.int64)
        if len(hA) * len(hB) <= len(hC) * len(hD):
            s._ensure(max(1, len(hA) * len(hB)))
            n = _join(hA, mA, hB, mB, z, False, hC, mC, hD, mD, f, True, K1, s.tk[:s.size], s.ti[:s.size], s.tj[:s.size], out, maxout)
            r = out[:min(n, maxout)]
        else:
            s._ensure(max(1, len(hC) * len(hD)))
            n = _join(hC, mC, hD, mD, f, True, hA, mA, hB, mB, z, False, K1, s.tk[:s.size], s.ti[:s.size], s.tj[:s.size], out, maxout)
            r = out[:min(n, maxout)][:, [2, 3, 0, 1]]
        assert n <= maxout, "candidate overflow"
        return r
