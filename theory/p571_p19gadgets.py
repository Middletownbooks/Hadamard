r"""Gadget (constant-coset) configurations of the <61> = Paley(19) slice of v = 571 (report section 7), exactly.

Setting [P] (report section 7.1).  M = <61> fixes Z30 pointwise and acts on Z19 as the squares, so every block is, on
each Z19-coset g in Z30, either twisted (x = a at 0, +-u on QR / QNR) or constant (a gadget: x = a at 0, w elsewhere).
With A_k(g) = a, W_k(g) = w on gadgets (0 elsewhere), y_k = A_k + 18 W_k is the Z19-compression and S1 reads
    sum_k PAF(y_k)(h) = -76  (h != 0 in Z30),   sum_k sum_g y_k(g)^2 = 2208,   row sums r = (2, 0, 0, 0).
Proved there: (K, N-) in {(6, 1), (7, 6)} (K gadgets, N- of type J': a = -w); |t_k| <= 1 for t_k = sum_g W_k(g);
s(h) = sum_k PAF(W_k)(h) in {0, -1} (h != 0), s(15) = 0; the number of h in 1..14 with s(h) = -1 is (K - #odd)/2,
#odd = #{k : K_k odd}; Lemma 6 at the trivial character excludes #odd = 4; and at every lag the X-sum
    x_h = sum_k sum_{g in supp W_k} W_k(g) (A_k(g + h) + A_k(g - h))
lies in [8, 2K] if s(h) = -1 and in [-10, 2] if s(h) = 0 (2K terms, each +-1; a term is FORCED when g +- h is another
gadget of the same block: it equals eps(g') W(g) W(g'), eps = +1 for J, -1 for J').

This script enumerates all gadget configurations (W_0..W_3, J/J' labels) satisfying these conditions, up to per-block
translation and reflection, negation of blocks 1..3, permutations of blocks 1..3 and the global multipliers of Z30.
usage: python3 theory/p571_p19gadgets.py [--out file.pkl]"""
import sys, os, json, time, math, pickle, itertools as it
from collections import defaultdict
import numpy as np
import numba

n = 30; V = 571; K_OPTIONS = ((6, 1), (7, 6))         # (K, N-)
MULTS = [1, 7, 11, 13, 17, 19, 23, 29]

@numba.njit(cache=True)
def _paf_vec(pos, sg, n):
    out = np.zeros(16, np.int64)                        # PAF at h = 0..15
    K = pos.shape[0]
    for i in range(K):
        for j in range(K):
            if i == j: continue
            d = (pos[j] - pos[i]) % n
            if d <= 15: out[d] += sg[i] * sg[j]
    out[0] = K
    return out

@numba.njit(cache=True)
def _enum(K, n, tpar_ok, fix_first_sign, maxexcess, combs, out_pos, out_sg, out_P):
    """all signed K-subsets containing 0 (positions from combs rows), |t| <= 1, with excess <= maxexcess"""
    cnt = 0
    nsg = 1 << K
    for r in range(combs.shape[0]):
        pos = combs[r]
        for m in range(nsg):
            if fix_first_sign and (m & 1) == 0: continue          # gadget at 0 has sign +1 (negation symmetry)
            sg = np.empty(K, np.int64); t = 0
            for i in range(K):
                sg[i] = 1 if (m >> i) & 1 else -1; t += sg[i]
            if t > 1 or t < -1: continue
            P = _paf_vec(pos, sg, n)
            ex = 0
            for h in range(1, 15):
                if P[h] > 0: ex += P[h]
                elif P[h] < -1: ex += -1 - P[h]
            ex += abs(P[15]) // 2
            if ex > maxexcess: continue
            if cnt < out_pos.shape[0]:
                out_pos[cnt, :K] = pos; out_sg[cnt, :K] = sg; out_P[cnt] = P
            cnt += 1
    return cnt

def shapes(K, block0, maxexcess):
    """signed K-subsets of Z30 containing 0 with |t| <= 1 and excess <= maxexcess: (pos (N,K), sg (N,K), P (N,16))"""
    if K == 0:
        return np.zeros((1, 0), np.int64), np.zeros((1, 0), np.int64), np.zeros((1, 16), np.int64)
    combs = np.array([(0,) + c for c in it.combinations(range(1, n), K - 1)], np.int64)
    cap = 1 << 22
    while True:
        op = np.zeros((cap, K), np.int64); osg = np.zeros((cap, K), np.int64); oP = np.zeros((cap, 16), np.int64)
        c = _enum(K, n, True, not block0, maxexcess, combs, op, osg, oP)
        if c <= cap: return op[:c], osg[:c], oP[:c]
        cap = c

def to_vec(pos, sg):
    w = np.zeros(n, np.int64); w[list(pos)] = sg; return w

def _sparse_canon(items, allow_neg):
    """items: list of (g, w, e) for the gadgets of one block; canonical form under translation, reflection and (if
    allowed) negation w -> -w: lexicographically minimal sorted tuple of (g', w', e) with min position 0"""
    best = None
    for s_ in ((1, -1) if allow_neg else (1,)):
        for rf in (1, -1):
            pts = [((rf * g) % n, s_ * w, e) for g, w, e in items]
            for g0, _, _ in pts:
                c = tuple(sorted(((g - g0) % n, w, e) for g, w, e in pts))
                if best is None or c < best: best = c
    return best if best is not None else ()

def canon_block(w, allow_neg):
    """canonical form of a signed set (vector w) under translation, reflection (and negation)"""
    return _sparse_canon([(int(g), int(w[g]), 0) for g in np.nonzero(w)[0]], allow_neg)

def canon_config(Ws, eps):
    """canonical form of a configuration: blocks (W vectors) with J/J' labels (eps vectors, 0 off gadgets) under
    per-block translation/reflection, negation of blocks 1..3, permutation of blocks 1..3 and global multipliers"""
    best = None
    for mu in MULTS:
        enc = []
        for k in range(4):
            items = [((mu * int(g)) % n, int(Ws[k][g]), int(eps[k][g])) for g in np.nonzero(Ws[k])[0]]
            enc.append(_sparse_canon(items, k > 0))
        key = (enc[0],) + tuple(sorted(enc[1:]))
        if best is None or key < best: best = key
    return best

def lemma6_trivial(t):
    """sum_k (r_k - 19 t_k)^2 + 19 sum rho_k^2 = 4v with rho_k = 30 - K_k = t_k (mod 2), |rho_k| <= 30 - K_k (approx
    cap 30)"""
    r = (2, 0, 0, 0)
    num = 4 * V - sum((r[k] - 19 * t[k]) ** 2 for k in range(4))
    if num < 0 or num % 19: return False
    Q = num // 19; par = [abs(x) % 2 for x in t]
    for r0 in range(par[0], 31, 2):
        for r1 in range(par[1], 31, 2):
            for r2 in range(par[2], 31, 2):
                q = Q - r0 * r0 - r1 * r1 - r2 * r2
                if q < 0: break
                r3 = math.isqrt(q)
                if r3 * r3 == q and r3 % 2 == par[3]: return True
    return False

def x_feasible(Ws, eps, s, K):
    """per lag h = 1..15: forced X-terms and the reachable range of x_h must meet [8, 2K] (s = -1) or [-10, 2]"""
    for h in range(1, 16):
        F = 0; nfree = 0
        for k in range(4):
            w = Ws[k]; e = eps[k]
            for g in np.nonzero(w)[0]:
                for sgn in (1, -1):
                    g2 = (g + sgn * h) % n
                    if w[g2] != 0: F += e[g2] * w[g] * w[g2]
                    else: nfree += 1
        lo, hi = (8, 2 * K) if s[h] == -1 else (-10, 2)
        if F + nfree < lo or F - nfree > hi: return False
        # parity: x = F + (sum of nfree +-1) has parity of F + nfree
        ok = any((x - F - nfree) % 2 == 0 for x in range(max(lo, F - nfree), min(hi, F + nfree) + 1))
        if not ok: return False
    return True

def distributions(K):
    out = []
    for k0 in range(K + 1):
        for k1 in range(K + 1 - k0):
            for k2 in range(k1 + 1):
                k3 = K - k0 - k1 - k2
                if k3 < 0 or k3 > k2: continue
                odd = sum(x % 2 for x in (k0, k1, k2, k3))
                if K == 6 and odd not in (0, 2): continue
                if K == 7 and odd not in (1, 3): continue
                out.append((k0, k1, k2, k3))
    return out

def dedupe(lst, K, block0):
    """one representative per per-block class (translation, reflection, and negation for blocks 1..3)"""
    pos, sg, P = lst; seen = {}; keep = []
    for i in range(len(pos)):
        c = canon_block(to_vec(pos[i][:K], sg[i][:K]), not block0) if K else ()
        if c not in seen: seen[c] = i; keep.append(i)
    keep = np.array(keep, np.int64)
    return pos[keep], sg[keep], P[keep]

def join4(Ps, targets):
    """all (i0, i1, i2, i3) with sum_k Ps[k][i_k] in targets (2+2 split, dict on the larger side)"""
    sz = [len(x) for x in Ps]
    left, right = ((0, 1), (2, 3)) if sz[0] * sz[1] <= sz[2] * sz[3] else ((2, 3), (0, 1))
    R = defaultdict(list)
    for a in range(sz[right[0]]):
        for b in range(sz[right[1]]):
            R[tuple(Ps[right[0]][a] + Ps[right[1]][b])].append((a, b))
    out = []
    tg = [np.asarray(t) for t in targets]
    for a in range(sz[left[0]]):
        for b in range(sz[left[1]]):
            L = Ps[left[0]][a] + Ps[left[1]][b]
            for t in tg:
                for (c, d) in R.get(tuple(t - L), ()):
                    sel = [0] * 4; sel[left[0]], sel[left[1]], sel[right[0]], sel[right[1]] = a, b, c, d
                    out.append(tuple(sel))
    return out

def main(outpath=None):
    t0 = time.time(); classes = {}; stats = []
    shape_cache = {}
    for K, Nm in K_OPTIONS:
        for dist in distributions(K):
            pairs = [math.comb(x, 2) for x in dist]
            lists = []
            for k, Kk in enumerate(dist):
                maxex = sum(pairs) - pairs[k]
                key = (Kk, k == 0, maxex)
                if key not in shape_cache: shape_cache[key] = dedupe(shapes(Kk, k == 0, maxex), Kk, k == 0)
                lists.append(shape_cache[key])
            sizes = [len(l[0]) for l in lists]
            nspec = (K - sum(x % 2 for x in dist)) // 2
            targets = []
            for spec in it.combinations(range(1, 15), nspec):
                s = np.zeros(15, np.int64); s[np.array(spec) - 1] = -1; targets.append(s)
            hits = join4([l[2][:, 1:] for l in lists], targets)
            found = 0; n_l6 = 0; n_x = 0
            for sel in hits:
                Ws = []; ts = []
                for k in range(4):
                    pos, sg = lists[k][0][sel[k]], lists[k][1][sel[k]]
                    Ws.append(to_vec(pos[:dist[k]], sg[:dist[k]])); ts.append(int(sg[:dist[k]].sum()))
                if not lemma6_trivial(ts): n_l6 += 1; continue
                sfull = np.zeros(16, np.int64); sfull[1:] = sum(lists[k][2][sel[k]][1:] for k in range(4))
                gpos = [(k, g) for k in range(4) for g in np.nonzero(Ws[k])[0]]
                for jm in it.combinations(range(len(gpos)), Nm):
                    eps = [np.zeros(n, np.int64) for _ in range(4)]
                    for i_, (k, g) in enumerate(gpos): eps[k][g] = -1 if i_ in jm else 1
                    if not x_feasible(Ws, eps, sfull, K): n_x += 1; continue
                    c = canon_config(Ws, eps)
                    if c not in classes:
                        classes[c] = dict(K=K, dist=dist, t=ts, spec=[h for h in range(1, 15) if sfull[h] == -1],
                                          W=[w.tolist() for w in Ws], eps=[e.tolist() for e in eps])
                        found += 1
            stats.append(dict(K=K, dist=dist, shape_classes=sizes, gadget_hits=len(hits), lemma6_cut=n_l6,
                              x_cut_labelings=n_x, classes=found, secs=round(time.time() - t0, 1)))
            print(json.dumps(stats[-1]), flush=True)
    byK = defaultdict(int)
    for c, info in classes.items(): byK[(info['K'], info['dist'])] += 1
    print(json.dumps(dict(summary=True, classes=len(classes), by_K_dist={f'{k[0]}:{k[1]}': v for k, v in sorted(byK.items())},
                          secs=round(time.time() - t0, 1))), flush=True)
    if outpath:
        with open(outpath, 'wb') as fh: pickle.dump(classes, fh)
    return classes

if __name__ == "__main__":
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    main(out)
