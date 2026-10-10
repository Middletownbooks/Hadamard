r"""Stage-2 constraint budget by level for a generalised-Paley slice and a given stage-1 y [E] (+ exact parts).

Stage 2 is sum_k P_k(h) = 4v delta (theory/q523_lib.py).  Its 36 coordinates (F = 18, s = 4) split by
  * the order d | F of the Z_F-character (d = 1, 2, 3, 6, 9, 18); the S-level of the Z_2-split is d in {1, 3, 9}
    (= compression Z_18 -> Z_9), the D-level d in {2, 6, 18};
  * the Galois character r of Gal(K/Q) = Z_4 acting on eta-coordinates: r = 0 (rational: trace), r = 2 (quadratic
    subfield Q(sqrt q)), r = 1, 3 (the quartic part, K^-).
A level = (set of d) x (set of r); its coordinates are the integer projections sum_{d} c_d(h - h') (Ramanujan sums) and
likewise on Z_4.  For each level the constraint bits are -log2 P(level hits its target) under the lattice-coordinate
Gaussian (theory/q523_gauss.py), sampled from uniform configurations given y.  Conditional bits C(A | B) = C(A u B) - C(B).
Per-block entropies of compressions are exact (independent cosets; convolution of the per-coset alphabets).
usage: python3 theory/q523_levels.py F q s y.json [--N 40000] [--seed 1] [--block-only]"""
import sys, os, json, math, itertools as it
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, zcoords, pz, target_P
from q523_gauss import gauss_bits

def ramanujan(d, n, x):
    """c_d(x) for d | n, x integer array: sum over a of exact order d of exp(2 pi i a x / n), as integers"""
    a = np.array([t for t in range(n) if math.gcd(t, n) == n // d]) if d > 1 else np.array([0])
    v = np.cos(2 * np.pi * np.outer(np.asarray(x) % n, a) / n).sum(1)
    return np.rint(v).astype(np.int64)

def level_matrix(F, s, dset, rset):
    """integer projection matrix (F*s, F*s) for characters of order in dset (on Z_F) and Galois chars rset (on Z_s)"""
    h = np.arange(F); Mh = sum(ramanujan(d, F, (h[:, None] - h[None, :]).ravel()).reshape(F, F) for d in dset)
    l = np.arange(s); Ml = sum(ramanujan(r, s, (l[:, None] - l[None, :]).ravel()).reshape(s, s) for r in rset)
    return np.kron(Mh, Ml)

def sample_P(cy, Y, N, rng, chunk=4000):
    F = len(Y[0]); out = np.zeros((N, F, cy.s), np.int64)
    alph = [[zcoords(cy.alphabet(int(y))) for y in Yk] for Yk in Y]
    for c0 in range(0, N, chunk):
        n = min(chunk, N - c0)
        Z = np.zeros((n, 4, F, cy.s), np.int64)
        for k in range(4):
            for g in range(F):
                A = alph[k][g]; Z[:, k, g] = A[rng.integers(len(A), size=n)]
        out[c0:c0 + n] = pz(Z, cy.T).sum(1)
    return out

def dist_sum(dists):
    """distribution (Counter of tuples) of the sum of independent uniform choices from each list of tuples"""
    cur = Counter({tuple([0] * len(dists[0][0])): 1.0})
    for L in dists:
        w = 1.0 / len(L); nxt = Counter()
        for a, pa in cur.items():
            for b in L: nxt[tuple(x + y for x, y in zip(a, b))] += pa * w
        cur = nxt
    return cur

def H(c):
    return -sum(p * math.log2(p) for p in c.values() if p > 0)

def block_entropies(cy, Yk, F):
    """per block: entropy and log2 support of the Z_d-compressions (d | F) of Z; the S/D pair split for even F"""
    alph = [[tuple(r) for r in zcoords(cy.alphabet(int(y)))] for y in Yk]
    out = {}
    for d in [d for d in range(1, F + 1) if F % d == 0]:
        Ht = 0.0; sup = 0.0
        for r in range(d):
            c = dist_sum([alph[g] for g in range(r, F, d)])
            Ht += H(c); sup += math.log2(len(c))
        out[f"Z{d}"] = dict(H=round(Ht, 2), log2_support=round(sup, 2))
    out["full"] = dict(H=round(sum(math.log2(len(a)) for a in alph), 2))
    if F % 2 == 0:
        hD = 0.0
        for g in range(F // 2):
            A, B = alph[g], alph[g + F // 2]
            joint = Counter(); dS = Counter()
            for a in A:
                for b in B:
                    S = tuple(x + y for x, y in zip(a, b)); joint[(S, tuple(x - y for x, y in zip(a, b)))] += 1; dS[S] += 1
            n = len(A) * len(B)
            hD += sum((m / n) * math.log2(m) for m in dS.values())       # H(D | S) = E log2 mult(S)
        out["H(D|S)"] = round(hD, 2)
    return out

def main(F, q, s, Y, N=40000, seed=1, block_only=False):
    cy = Cyclo(q, s); m = F * q; v = m + 1
    res = dict(F=F, q=q, s=s, v=v, class_of_minus1=cy.cneg)
    res["blocks"] = [block_entropies(cy, Y[k], F) for k in range(4)]
    res["H2_total"] = round(sum(b["full"]["H"] for b in res["blocks"]), 2)
    if block_only: return res
    rng = np.random.default_rng(seed)
    P = sample_P(cy, Y, N, rng).reshape(N, -1)
    T = target_P(F, v, s).ravel()
    divs = [d for d in range(1, F + 1) if F % d == 0]
    rall = [r for r in range(s) if s % (r if r else s) == 0] if False else None
    rsets = {"all": [d for d in range(1, s + 1) if s % d == 0]}          # Galois characters by order: 1 (r=0), 2 (r=2), 4 (r=1,3)
    if s == 4: rsets.update({"rat": [1], "quad": [2], "quart": [4], "rat+quad": [1, 2]})
    elif s == 2: rsets.update({"rat": [1], "quad": [2]})
    levels = {}
    def lev(name, dset, rkey="all"):
        M = level_matrix(F, s, dset, rsets[rkey])
        V = P @ M.T; t = M @ T
        g = gauss_bits(V, t)
        levels[name] = dict(dset=dset, r=rkey, **g)
        return g
    lev("ALL", divs)
    for d in divs: lev(f"d={d}", [d])
    if F % 2 == 0:
        Sd = [d for d in divs if (F // 2) % d == 0]; Dd = [d for d in divs if d not in Sd]
        lev("S", Sd); lev("D", Dd)
    for e in divs:
        if e not in (1, F): lev(f"Z{e}", [d for d in divs if e % d == 0])
    if F % 6 == 0:
        lev("S+Z6", sorted(set([d for d in divs if (F // 2) % d == 0] + [d for d in divs if 6 % d == 0])))
    for rk in rsets:
        if rk != "all": lev(f"ALL[{rk}]", divs, rk)
    if F % 2 == 0:
        for rk in ("rat", "quad", "quart", "rat+quad"):
            if rk in rsets:
                lev(f"S[{rk}]", Sd, rk); lev(f"D[{rk}]", Dd, rk)
    res["levels"] = levels
    b = lambda k: levels[k]["bits"]
    if F % 2 == 0:
        res["C(D|S)"] = round(b("ALL") - b("S"), 2); res["C(S|D)"] = round(b("ALL") - b("D"), 2)
    res["yield_log2"] = round(res["H2_total"] - b("ALL"), 2)
    return res

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "40000", "--seed": "1"}; pos = []; flags = set(); i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        elif a[i].startswith("--"): flags.add(a[i]); i += 1
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y["y"] if isinstance(Y, dict) else Y
    r = main(F, q, s, Y, int(opt["--N"]), int(opt["--seed"]), "--block-only" in flags)
    print(json.dumps(r))
