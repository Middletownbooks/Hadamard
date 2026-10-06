"""Search for T-matrices of order t invariant under a multiplier subgroup K <= (Z/t)*, t prime.

T-matrices: four {0,+-1} sequences of length t with disjoint supports covering Z_t and sum of PAFs = t*delta.
Under K, each orbit gets a label in 0..7 (row = label & 3, sign = -1 if label & 4). PAF is constant on lag
classes s*K, so the conditions are f_c = 0 for each lag class c.
f_c = sum_{a,b} M[c,a,b] [row_a == row_b] sgn_a sgn_b,  M[c,a,b] = #{x in O_a : x + s_c in O_b}.
Local search: tabu on single-orbit relabels, minimising E = sum_c w_c f_c^2 (w_c = class size).
"""
import sys, json, time
import numpy as np
from numba import njit

def subgroup(t, k):
    for g in range(2, t):
        if pow(g, (t - 1) // k, t) != 1 or True:
            pass
    # generator of order exactly k
    for g in range(2, t):
        h = pow(g, (t - 1) // k, t)
        if all(pow(h, k // p, t) != 1 for p in prime_factors(k)) and pow(h, k, t) == 1:
            return sorted({pow(h, i, t) for i in range(k)})
    raise ValueError

def prime_factors(k):
    f, p = set(), 2
    while p * p <= k:
        while k % p == 0: f.add(p); k //= p
        p += 1
    if k > 1: f.add(k)
    return f

def setup(t, K):
    seen, orbits = set(), []
    for a in range(t):
        if a in seen: continue
        o = sorted({a * h % t for h in K}); seen.update(o); orbits.append(o)
    oid = np.zeros(t, dtype=np.int64)
    for i, o in enumerate(orbits): oid[o] = i
    reps = [o[0] for o in orbits[1:]]          # lag classes = nonzero orbits
    m = len(orbits)
    M = np.zeros((len(reps), m, m), dtype=np.int64)
    for c, s in enumerate(reps):
        for x in range(t):
            M[c, oid[x], oid[(x + s) % t]] += 1
    w = np.array([len(o) for o in orbits[1:]], dtype=np.int64)
    return orbits, oid, M, w

@njit(cache=True)
def fvals(L, M):
    C, m = M.shape[0], M.shape[1]
    f = np.zeros(C, dtype=np.int64)
    for c in range(C):
        s = 0
        for a in range(m):
            for b in range(m):
                if (L[a] & 3) == (L[b] & 3):
                    sa = -1 if L[a] & 4 else 1
                    sb = -1 if L[b] & 4 else 1
                    s += M[c, a, b] * sa * sb
        f[c] = s
    return f

@njit(cache=True)
def delta_f(L, M, a, new):
    # change in f when orbit a relabels L[a] -> new
    C, m = M.shape[0], M.shape[1]
    d = np.zeros(C, dtype=np.int64)
    old = L[a]
    so = -1 if old & 4 else 1
    sn = -1 if new & 4 else 1
    for c in range(C):
        s = 0
        for b in range(m):
            if b == a:
                continue
            sb = -1 if L[b] & 4 else 1
            mab = M[c, a, b] + M[c, b, a]
            if mab == 0:
                continue
            if (L[b] & 3) == (old & 3):
                s -= mab * so * sb
            if (L[b] & 3) == (new & 3):
                s += mab * sn * sb
        # a-a term: (row same, sign^2 = 1) unchanged
        d[c] = s
    return d

@njit(cache=True)
def tabu(L, M, w, iters, tenure, seed):
    np.random.seed(seed)
    C, m = M.shape[0], M.shape[1]
    f = fvals(L, M)
    E = 0
    for c in range(C): E += w[c] * f[c] * f[c]
    best_E = E
    tabu_until = np.zeros((m, 8), dtype=np.int64)
    for it in range(iters):
        if E == 0:
            return L, 0, it
        bE, ba, bn, ties = 1 << 62, -1, -1, 0
        bd = np.zeros(C, dtype=np.int64)
        for a in range(m):
            for nl in range(8):
                if nl == L[a]:
                    continue
                d = delta_f(L, M, a, nl)
                e = 0
                for c in range(C):
                    v = f[c] + d[c]; e += w[c] * v * v
                if tabu_until[a, nl] > it and e >= best_E:
                    continue
                if e < bE:
                    bE, ba, bn, ties, bd = e, a, nl, 1, d
                elif e == bE:
                    ties += 1
                    if np.random.randint(ties) == 0:
                        ba, bn, bd = a, nl, d
        if ba < 0:
            continue
        tabu_until[ba, L[ba]] = it + tenure + np.random.randint(tenure)
        L[ba] = bn
        for c in range(C): f[c] += bd[c]
        E = bE
        if E < best_E: best_E = E
    return L, best_E, iters

def verify(t, orbits, L):
    T = np.zeros((4, t), dtype=np.int64)
    for o, l in zip(orbits, L): T[l & 3, o] = -1 if l & 4 else 1
    tot = sum(np.array([int(np.dot(r, np.roll(r, -s))) for s in range(t)]) for r in T)
    return bool(tot[0] == t and (tot[1:] == 0).all()), T

if __name__ == "__main__":
    t, k = int(sys.argv[1]), int(sys.argv[2])
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
    restarts = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
    seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    K = subgroup(t, k)
    orbits, oid, M, w = setup(t, K)
    m = len(orbits)
    print(json.dumps({"t": t, "k": k, "K": K, "orbits": m, "classes": M.shape[0]}), flush=True)
    rng = np.random.default_rng(seed0)
    t0 = time.time()
    for r in range(restarts):
        L = rng.integers(0, 8, m).astype(np.int64)
        L, E, it = tabu(L, M, w, iters, max(3, m // 4), seed0 * 100003 + r)
        if E == 0:
            ok, T = verify(t, orbits, L)
            print(json.dumps({"FOUND": ok, "t": t, "k": k, "restart": r, "iters": int(it),
                              "secs": round(time.time() - t0, 1), "K": K, "labels": L.tolist()}), flush=True)
            if ok: break
        elif r % 10 == 0:
            print(json.dumps({"restart": r, "best_E": int(E), "secs": round(time.time() - t0, 1)}), flush=True)
