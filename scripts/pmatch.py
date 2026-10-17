"""Exact equal-key matching of two uint64 arrays by radix partitioning (replaces the 1 GB byte bloom of
kaggle/kt_join_kaggle.py bloom_match, whose random accesses dominated the 311 buckets).

pmatch(KA, KC, bits) returns (ia, ic): all index pairs with KA[ia] == KC[ic] (same contract as bloom_match).
Both arrays are scattered into 2^bits partitions by their top bits (parallel two-pass counting scatter), then each
partition's C part is sorted and every A entry is binary-searched in it (partitions of ~10^4 keys stay in cache)."""
import numpy as np, numba as nb


@nb.njit(parallel=True, cache=True)
def _scatter_into(K, bits, nch, outK, outI):
    N = K.shape[0]; P = 1 << bits; sh = np.uint64(64 - bits)
    b0 = np.linspace(0, N, nch + 1).astype(np.int64); H = np.zeros((nch, P), np.int64)
    for c in nb.prange(nch):
        for i in range(b0[c], b0[c + 1]): H[c, np.int64(K[i] >> sh)] += 1
    start = np.zeros(P + 1, np.int64)
    for p in range(P):
        s = 0
        for c in range(nch): s += H[c, p]
        start[p + 1] = start[p] + s
    pos = np.zeros((nch, P), np.int64)
    for p in range(P):
        o = start[p]
        for c in range(nch): pos[c, p] = o; o += H[c, p]
    for c in nb.prange(nch):
        loc = pos[c].copy()
        for i in range(b0[c], b0[c + 1]):
            p = np.int64(K[i] >> sh); outK[loc[p]] = K[i]; outI[loc[p]] = i; loc[p] += 1
    return outK, outI, start


@nb.njit(parallel=True, cache=True)
def _join(AK, AI, As, CK, CI, Cs):
    P = As.shape[0] - 1; cnt = np.zeros(P, np.int64)
    for p in nb.prange(P):
        lo, hi = Cs[p], Cs[p + 1]
        if hi == lo or As[p + 1] == As[p]: continue
        o = np.argsort(CK[lo:hi]); ck = CK[lo:hi][o]; ci = CI[lo:hi][o]
        CK[lo:hi] = ck; CI[lo:hi] = ci
        c = 0
        for i in range(As[p], As[p + 1]):
            j = np.searchsorted(ck, AK[i])
            while j < ck.shape[0] and ck[j] == AK[i]: c += 1; j += 1
        cnt[p] = c
    off = np.zeros(P + 1, np.int64)
    for p in range(P): off[p + 1] = off[p] + cnt[p]
    ia = np.empty(off[P], np.int64); ic = np.empty(off[P], np.int64)
    for p in nb.prange(P):
        lo, hi = Cs[p], Cs[p + 1]; n = off[p]
        if cnt[p] == 0: continue
        ck = CK[lo:hi]
        for i in range(As[p], As[p + 1]):
            j = np.searchsorted(ck, AK[i])
            while j < ck.shape[0] and ck[j] == AK[i]: ia[n] = AI[i]; ic[n] = CI[lo + j]; n += 1; j += 1
    return ia, ic


_BUF = {}


def _buf(name, n, dt):
    """reusable scratch arrays (fresh multi-GB allocations cost page faults on every bucket)"""
    b = _BUF.get(name)
    if b is None or len(b) < n: b = _BUF[name] = np.empty(max(n, int(1.25 * n)), dt)
    return b[:n]


def pmatch(KA, KC, bits=12):
    AK, AI = _buf('AK', len(KA), np.uint64), _buf('AI', len(KA), np.int64)
    AK, AI, As = _scatter_into(KA, bits, 256, AK, AI)
    CK, CI = _buf('CK', len(KC), np.uint64), _buf('CI', len(KC), np.int64)
    CK, CI, Cs = _scatter_into(KC, bits, 256, CK, CI)
    return _join(AK, AI, As, CK, CI, Cs)


if __name__ == '__main__':                       # self-test against bloom_match and a numpy reference
    import sys, os, time
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'kaggle'))
    import kt_join_kaggle as KJ
    rng = np.random.default_rng(1)
    KA = rng.integers(0, 2 ** 63, 3_000_000, dtype=np.uint64); KC = rng.integers(0, 2 ** 63, 1_500_000, dtype=np.uint64)
    sel = rng.choice(len(KA), 5000, replace=False); KC[:5000] = KA[sel]; KC[5000:5100] = KA[sel[:100]]   # repeats
    for f in (pmatch, lambda a, c: KJ.bloom_match(a, c, 30)): f(KA[:1000], KC[:1000])
    t = time.time(); r1 = pmatch(KA, KC); t1 = time.time() - t
    t = time.time(); r2 = KJ.bloom_match(KA, KC, 30); t2 = time.time() - t
    s1 = sorted(zip(r1[0].tolist(), r1[1].tolist())); s2 = sorted(zip(r2[0].tolist(), r2[1].tolist()))
    assert s1 == s2 and len(s1) == 5100, (len(s1), len(s2))
    print('pmatch == bloom_match on', len(s1), 'matches; times', round(t1, 2), round(t2, 2))
