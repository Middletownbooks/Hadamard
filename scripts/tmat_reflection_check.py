"""Exact checks of the reflection lemma for T-matrices of odd prime order t.

Lemma: if T is invariant under x -> -x combined with a signed row permutation h (h an involution of B4), then h swaps
two rows i, j and fixes the other two; the fixed rows carry only position 0; T = {W, e*W^R, +-delta_0, 0} with W a
circulant weighing sequence of weight (t-1)/2 supported on a half-system.  Hence (t-1)/2 = s^2, and a prime q | s that
is self-conjugate mod t (-1 in <q>) rules it out.
Checks: (a) every T-matrix of order 3, 5, 7 (all 8^t labellings) tested against all 16 involutions h;
(b) T-matrices with all four rows symmetric (h = identity type) at t = 11, 13: count;
(c) t = 19: all W on half-systems with PAF(W) = 9*delta: count; (d) self-conjugacy table for t = 2s^2 + 1 prime."""
import itertools, math, json
import numpy as np

def paf_rows(T):
    F = np.fft.fft(T, axis=-1)
    return np.rint(np.fft.ifft(np.abs(F) ** 2, axis=-1).real).astype(np.int64)

def all_tmatrices(t):
    out = []
    idx = np.arange(8 ** t, dtype=np.int64)
    lab = np.stack([(idx >> (3 * x)) & 7 for x in range(t)], 1)               # N x t
    T = np.zeros((len(idx), 4, t), np.int8)
    for x in range(t):
        T[np.arange(len(idx)), lab[:, x] & 3, x] = np.where(lab[:, x] & 4, -1, 1)
    tot = paf_rows(T.astype(float)).sum(1)
    ok = (tot[:, 0] == t) & (tot[:, 1:] == 0).all(1)
    return T[ok]

INV = []                                                   # involutions of B4 (incl. identity type) as (perm, signs)
for p in itertools.permutations(range(4)):
    if any(p[p[i]] != i for i in range(4)): continue
    for e in itertools.product((1, -1), repeat=4):
        if all(e[i] * e[p[i]] == 1 for i in range(4)): INV.append((p, e))

def reflect_invariant(T, h):
    p, e = h; t = T.shape[1]; neg = [(-x) % t for x in range(t)]
    return all((T[p[i]][neg] == e[i] * T[i]).all() for i in range(4))

res = {}
for t in (3, 5, 7):
    Ts = all_tmatrices(t)
    inv = [T for T in Ts if any(reflect_invariant(T, h) for h in INV)]
    shapes = set()
    for T in inv:
        w = sorted(int(np.abs(r).sum()) for r in T); shapes.add(tuple(w))
    res[f"a_t{t}"] = {"T_matrices": len(Ts), "reflection_invariant": len(inv), "row_weights": sorted(shapes)}
print(json.dumps(res)); res = {}

for t in (11, 13):                                         # all rows symmetric: labels on 0 and on (t-1)/2 pairs
    h = (t - 1) // 2; n = 8 ** (h + 1); cnt = 0
    for s in range(0, n, 1 << 20):
        idx = np.arange(s, min(n, s + (1 << 20)), dtype=np.int64)
        lab = np.stack([(idx >> (3 * j)) & 7 for j in range(h + 1)], 1)
        T = np.zeros((len(idx), 4, t), np.float32); rr = np.arange(len(idx))
        for x in range(t):
            j = min(x, t - x)
            T[rr, lab[:, j] & 3, x] = np.where(lab[:, j] & 4, -1, 1)
        tot = paf_rows(T).sum(1); cnt += int(((tot[:, 0] == t) & (tot[:, 1:] == 0).all(1)).sum())
    res[f"b_symmetric_t{t}"] = cnt
print(json.dumps(res)); res = {}

t = 19; h = 9; cnt = 0                                     # W on a half-system: choose x or -x in each pair, and a sign
for side in range(2 ** h):
    pos = [(j + 1) if not (side >> j) & 1 else t - (j + 1) for j in range(h)]
    sg = np.array([[1 - 2 * ((s >> j) & 1) for j in range(h)] for s in range(2 ** h)], np.float64)
    W = np.zeros((2 ** h, t)); W[:, pos] = sg
    P = paf_rows(W); cnt += int((P[:, 1:] == 0).all(1).sum())
res["c_CW19_9_halfsystem"] = cnt

def selfconj(q, t):
    x = q % t
    for _ in range(t):
        if x == t - 1: return True
        x = x * q % t
    return False
def isprime(n): return n > 1 and all(n % p for p in range(2, math.isqrt(n) + 1))
tab = []
for s in range(1, 40):
    t = 2 * s * s + 1
    if not isprime(t): continue
    qs = [q for q in range(2, s + 1) if s % q == 0 and isprime(q)]
    tab.append({"t": t, "s": s, "self_conjugate_q": [q for q in qs if selfconj(q, t)], "other_q": [q for q in qs if not selfconj(q, t)]})
res["d_2s2p1"] = tab
print(json.dumps(res))
