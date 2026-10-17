"""Flat pools with W-word orbit masks (r <= 64 W orbits), for b = 0 slices with more than 64 orbits (223: 75, 491: 71,
283: 95).  Same conventions as kaggle/kt_join_kaggle.py and scripts/flat_probe.py: bit q of the mask = orbit q is +1
(word q // 64, bit q % 64); bit 0 = orbit {0}; orbit 1 + j = g^j U; blocks stored with positive row sum.

Kernels: walk_w (tabu walk storing rho-flat blocks with an admissible |row sum|), paf_w (lag-class PAFs from masks),
rotate_w / canonical_w (multiplier closure: rotation of the e orbit bits), unique_w (row-wise unique), rowsum_w,
unpack_w (masks -> +-1 orbit labels).  Validated against the 64-bit path in scripts/flat_probe.py (selftest)."""
import os, sys
import numpy as np, numba as nb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kaggle'))
import kt_join_kaggle as KJ


def nwords(r): return (r + 63) // 64


@nb.njit(cache=True)
def _init(y, M):
    r = y.shape[0]; L = M.shape[0]; G = np.zeros((L, r), np.int64); P = np.zeros(L, np.int64)
    for c in range(L):
        for q in range(r):
            s = 0
            for p in range(r): s += M[c, q, p] * y[p]
            G[c, q] = s
        t = 0
        for q in range(r): t += y[q] * G[c, q]
        P[c] = t // 2
    return G, P


@nb.njit(cache=True)
def walk_w(y, M, F0, nsteps, tenure, seed, d, targets, out, nout):
    """tabu walk minimising F = sum_c PAF(c)^2; every state with F <= F0 and |row sum| in targets is stored (W words)"""
    np.random.seed(seed); r = y.shape[0]; L = M.shape[0]; W = out.shape[1]; G, P = _init(y, M)
    F = 0
    for c in range(L): F += P[c] * P[c]
    rs = y[0]
    for q in range(1, r): rs += d * y[q]
    tabu = np.zeros(r, np.int64); n = nout
    for step in range(nsteps):
        best = 1 << 62; bq = -1; nb_ = 0
        for q in range(r):
            if tabu[q] > step: continue
            yq = y[q]; dF = 0
            for c in range(L):
                dd = -2 * yq * G[c, q] + 2 * M[c, q, q]; dF += 2 * P[c] * dd + dd * dd
            if dF < best or (dF == best and np.random.random() < 1.0 / (nb_ + 1)):
                if dF < best: nb_ = 1
                else: nb_ += 1
                best = dF; bq = q
        yq = y[bq]
        for c in range(L):
            dd = -2 * yq * G[c, bq] + 2 * M[c, bq, bq]; F += 2 * P[c] * dd + dd * dd; P[c] += dd
        for c in range(L):
            for p in range(r): G[c, p] -= 2 * yq * M[c, p, bq]
        y[bq] = -yq; tabu[bq] = step + 1 + tenure + np.random.randint(tenure // 2 + 1)
        rs += -2 * yq * (1 if bq == 0 else d)
        if F <= F0:
            a = abs(rs); ok = False
            for t in targets:
                if a == t: ok = True
            if ok:
                sgn = 1 if rs > 0 else -1
                for wd in range(W): out[n, wd] = np.uint64(0)
                for q in range(r):
                    if y[q] * sgn > 0: out[n, q // 64] |= np.uint64(1) << np.uint64(q % 64)
                n += 1
                if n >= out.shape[0]: return n
    return n


@nb.njit(parallel=True, cache=True)
def paf_w(bits, cls, rows, cols, cnt, L):
    N = bits.shape[0]; out = np.zeros((N, L), np.int16)
    for n in nb.prange(N):
        for e in range(cls.shape[0]):
            a = rows[e]; b = cols[e]
            yq = 1 if (bits[n, a // 64] >> np.uint64(a % 64)) & np.uint64(1) else -1
            yp = 1 if (bits[n, b // 64] >> np.uint64(b % 64)) & np.uint64(1) else -1
            out[n, cls[e]] += cnt[e] * yq * yp
    return out


@nb.njit(parallel=True, cache=True)
def rotate_w(bits, s, e):
    """new orbit 1 + j = old orbit 1 + ((j + s) mod e); orbit 0 kept"""
    N, W = bits.shape; out = np.zeros_like(bits)
    for n in nb.prange(N):
        if bits[n, 0] & np.uint64(1): out[n, 0] |= np.uint64(1)
        for j in range(e):
            src = 1 + (j + s) % e; dst = 1 + j
            if (bits[n, src // 64] >> np.uint64(src % 64)) & np.uint64(1):
                out[n, dst // 64] |= np.uint64(1) << np.uint64(dst % 64)
    return out


@nb.njit(parallel=True, cache=True)
def canonical_w(bits, e):
    """least rotation (lexicographic on words from the top word down) of each mask"""
    N, W = bits.shape; out = bits.copy(); cur = np.zeros(W, np.uint64)
    for n in nb.prange(N):
        cur = np.zeros(W, np.uint64)
        for s in range(1, e):
            for wd in range(W): cur[wd] = np.uint64(0)
            if bits[n, 0] & np.uint64(1): cur[0] |= np.uint64(1)
            for j in range(e):
                src = 1 + (j + s) % e; dst = 1 + j
                if (bits[n, src // 64] >> np.uint64(src % 64)) & np.uint64(1):
                    cur[dst // 64] |= np.uint64(1) << np.uint64(dst % 64)
            less = False
            for wd in range(W - 1, -1, -1):
                if cur[wd] < out[n, wd]: less = True; break
                if cur[wd] > out[n, wd]: break
            if less:
                for wd in range(W): out[n, wd] = cur[wd]
    return out


def unique_w(bits):
    b = np.ascontiguousarray(bits)
    v = b.view(np.dtype((np.void, 8 * b.shape[1]))).ravel()
    _, idx = np.unique(v, return_index=True)
    return b[np.sort(idx)]


def isin_w(a, b):
    va = np.ascontiguousarray(a).view(np.dtype((np.void, 8 * a.shape[1]))).ravel()
    vb = np.ascontiguousarray(b).view(np.dtype((np.void, 8 * b.shape[1]))).ravel()
    return np.isin(va, vb)


def unpack_w(bits, r):
    q = np.arange(r)
    return np.where((bits[:, q // 64] >> (q % 64).astype(np.uint64)) & np.uint64(1), 1, -1).astype(np.int64)


def rowsum_w(bits, r, d):
    Y = unpack_w(bits, r); return Y[:, 0] + d * Y[:, 1:].sum(1)


def selftest():
    """W = 1 and W = 2 paths agree with the 64-bit path of flat_probe on random masks (v = 311, d = 5)"""
    sys.path.insert(0, HERE)
    from flat_probe import rotate, canonical
    R = KJ.Red(311, 5); e = (311 - 1) // 5; Nsp = KJ.sparse_N(R); rng = np.random.default_rng(1)
    w1 = rng.integers(0, 2 ** 63, 5000, dtype=np.uint64) & np.uint64((1 << (e + 1)) - 1)
    b1 = w1[:, None]; b2 = np.concatenate([b1, np.zeros((len(b1), 1), np.uint64)], 1)
    assert np.array_equal(paf_w(b1, *Nsp, R.L), KJ.paf_from_bits(w1, *Nsp, R.L))
    assert np.array_equal(paf_w(b2, *Nsp, R.L), KJ.paf_from_bits(w1, *Nsp, R.L))
    for s in (1, 5, e // 2): assert np.array_equal(rotate_w(b1, s, e)[:, 0], rotate(w1, s, e))
    assert np.array_equal(canonical_w(b1, e)[:, 0], canonical(w1, e))
    assert np.array_equal(canonical_w(b2, e)[:, 0], canonical(w1, e))
    # a 75-orbit slice (223, d = 3): rotation keeps F, canonical is rotation-invariant
    R2 = KJ.Red(223, 3); e2 = (223 - 1) // 3; N2 = KJ.sparse_N(R2)
    Y = rng.choice([-1, 1], size=(2000, R2.r)); B = np.zeros((2000, 2), np.uint64)
    for q in range(R2.r): B[Y[:, q] > 0, q // 64] |= np.uint64(1) << np.uint64(q % 64)
    assert np.array_equal(unpack_w(B, R2.r), Y)
    F = lambda bb: (paf_w(bb, *N2, R2.L).astype(np.int64) ** 2).sum(1)
    for s in (1, 7, e2 // 2): assert np.array_equal(F(B), F(rotate_w(B, s, e2)))
    assert np.array_equal(canonical_w(B, e2), canonical_w(rotate_w(B, 9, e2), e2))
    P = paf_w(B, *N2, R2.L); X = R2.expand(Y); Pd = R2.paf(X)[:, R2.lags]
    assert np.array_equal(P, Pd), 'paf_w disagrees with direct FFT PAF'
    print('flatw selftest: ok (W = 1, 2 against the 64-bit path at 311; rotation invariance and FFT PAF at 223)')


if __name__ == '__main__':
    selftest()
