"""Williamson helpers: symmetric +-1 sequences on Z_n, PAF, verification, row-sum cases.
A symmetric sequence is stored by its half h[0..(n-1)/2]; x[i] = x[n-i] = h[i]."""
import numpy as np

def expand(h, n):
    h = np.asarray(h); idx = np.minimum(np.arange(n), n - np.arange(n)); idx[0] = 0
    return h[..., idx]

def paf(x):
    """periodic autocorrelation along last axis (exact integers)"""
    n = x.shape[-1]
    return np.stack([np.sum(x * np.roll(x, -s, axis=-1), axis=-1) for s in range(n)], axis=-1)

def verify(Q):
    """independent check: four symmetric +-1 sequences with PAF sum 0 off zero"""
    Q = np.asarray(Q, dtype=np.int64); n = Q.shape[1]
    ok = Q.shape[0] == 4 and set(np.unique(Q)) <= {-1, 1}
    ok &= all(np.array_equal(x, x[(-np.arange(n)) % n]) for x in Q)
    R = sum(int(x[0] * 0) + np.array([np.dot(x, np.roll(x, -s)) for s in range(n)]) for x in Q)
    return bool(ok and R[0] == 4 * n and np.all(R[1:] == 0))

def rowsum_cases(n):
    """normalised row sums s_A >= s_B >= s_C >= s_D > 0, all odd (n odd), sum of squares 4n"""
    out = []; odd = range(1, int((4 * n) ** .5) + 1, 2)
    for a in odd:
        for b in odd:
            if b > a: break
            for c in odd:
                if c > b: break
                r = 4 * n - a * a - b * b - c * c
                if r <= 0: continue
                d = int(round(r ** .5))
                if d * d == r and d % 2 == 1 and d <= c: out.append((a, b, c, d))
    return out

def epsilon(n, sums):
    """product condition: a_j b_j c_j d_j = eps for all j != 0.
    From PAF_X(s) = n - 4 w_X + 4 sum_i c_i c_{i+s} and sum_i c_i c_{i+s} = c_{s/2} (mod 2) for symmetric c:
    the number of -1's among (a_j,b_j,c_j,d_j) has parity (W - n) mod 2, W = total number of -1's."""
    W = sum((n - s) // 2 for s in sums)
    return -1 if (W - n) % 2 else 1
