"""Turyn-family Williamson quadruple of order n=(q+1)/2 for q = p^k = 1 mod 4, built from GF(q^2) (W2-T63).
GF(q^2) = GF(p^(2k)) via a primitive polynomial; trace T(x) = x + x^q to GF(q); chi = quadratic character of GF(q).
s(j) = chi(T(g^j)), j in Z_{2(q+1)} (antiperiodic with period q+1). Candidate blocks are even/odd parts; the
independent verifier decides."""
import sys, itertools, json, numpy as np

def poly_mulmod(a, b, mod, p):
    d = len(mod) - 1; r = [0] * (2 * d - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
    for i in range(len(r) - 1, d - 1, -1):          # mod is monic, low-to-high coefficients
        c = r[i]
        if c:
            for t in range(d + 1): r[i - d + t] = (r[i - d + t] - c * mod[t]) % p
    return r[:d]

def field(p, m):
    """find a primitive monic polynomial of degree m; return list of powers g^0..g^(p^m-2) as tuples"""
    N = p ** m - 1
    for tail in itertools.product(range(p), repeat=m):
        if tail[0] == 0: continue
        mod = list(tail) + [1]
        x = [0, 1] + [0] * (m - 2); cur = [1] + [0] * (m - 1); pw = []
        seen1 = False
        for e in range(N):
            pw.append(tuple(cur)); cur = poly_mulmod(cur, x, mod, p)
            if cur == [1] + [0] * (m - 1) and e < N - 1: seen1 = True; break
        if not seen1 and cur == [1] + [0] * (m - 1): return pw
    raise RuntimeError

def turyn(p, k):
    q = p ** k; m = 2 * k; n = (q + 1) // 2
    pw = field(p, m); N = len(pw); log = {v: e for e, v in enumerate(pw)}
    zero = tuple([0] * m)
    def add(a, b): return tuple((x + y) % p for x, y in zip(a, b))
    def chi(v):
        if v == zero: return 0
        e = log[v]; assert e % (q + 1) == 0, "trace not in GF(q)"
        return 1 if (e // (q + 1)) % 2 == 0 else -1
    s = [chi(add(pw[j], pw[(j * q) % N])) for j in range(2 * (q + 1))]
    return n, np.array(s)

if __name__ == "__main__":
    p, k = int(sys.argv[1]), int(sys.argv[2])
    n, s = turyn(p, k)
    print("n", n, "zeros at", list(np.where(s == 0)[0]), "antiperiodic:", bool(np.all(s[n*2:] == -s[:n*2])))
    np.save(f"/tmp/claude-0/-home-user-Hadamard/e1244fe7-42e1-5e13-91cc-287618221334/scratchpad/turyn_s_{p}_{k}.npy", s)

def pool(n, s):
    """symmetric +-1 length-n sequences a_t(i) = (-1)^i s(2i+t), zero filled with +1 or -1"""
    out = {}
    L = len(s)
    for t in range(L):
        a = np.array([(-1) ** i * s[(2 * i + t) % L] for i in range(n)])
        zs = np.where(a == 0)[0]
        for fill in itertools.product([1, -1], repeat=len(zs)):
            b = a.copy(); b[zs] = fill
            if np.array_equal(b, b[(-np.arange(n)) % n]):
                out[b.tobytes()] = (t, fill, b)
    return list(out.values())

def find(n, s, limit=10):
    from wil_common import verify, paf
    P = pool(n, s)
    X = np.array([b for _, _, b in P]); A = paf(X)[:, 1:]
    tab = {}
    for i in range(len(X)):
        for j in range(i, len(X)): tab.setdefault(tuple(A[i] + A[j]), []).append((i, j))
    found = []
    for key, lst in tab.items():
        neg = tuple(-x for x in key)
        for (i, j) in lst:
            for (k, l) in tab.get(neg, []):
                Q = X[[i, j, k, l]]
                if verify(Q):
                    found.append(([P[x][0] for x in (i, j, k, l)], [list(P[x][1]) for x in (i, j, k, l)], Q))
                    if len(found) >= limit: return P, found
    return P, found
