r"""Exact stage-1 solver for b >= 3 Paley-coset cores with a fixed Q-compression and a fixed w-gadget.

Stage 1 (report/paley_design_answer.md section 2):  sum_k PAF(y_k) = T_J on G',  y = alpha + (p-1) w.
Input: per block the alpha coset sums over H' (alpha_Q; the border forces the compression, report section 5) and the
gadget w (positions, signs).  alpha ranges over all +-1 functions with those coset sums (alpha at gadget positions is
free: J' (alpha = -w) or J (alpha = +w); the norm identity is checked at the end).  The four blocks are joined by
meet-in-the-middle: the two smallest lists are summed and hashed (linear 64-bit hash of the PAF vector), the third
list is iterated and the fourth vectorised; every hash hit is re-checked exactly.

usage: python3 theory/paley_stage1.py check359     (published compression + gadget: lists all stage-1 solutions and
                                                     checks that the published y is among them)"""
import sys, os, json, time, itertools as it, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'theory'))
import paley_ls as LSm
from paley_algebra import load_core, components

def coset_index(P):
    Hs = sorted(P.H); rep = -np.ones(P.n, np.int64); q = 0
    for g in range(P.n):
        if rep[g] < 0:
            for h in Hs: rep[P.add[g, h]] = q
            q += 1
    return rep

def block_list(P, rep, alphaQ, w):
    """all alpha (+-1 on G') with coset sums alphaQ; returns array (N, n)"""
    per = []
    for q in range(rep.max() + 1):
        pos = np.nonzero(rep == q)[0]; k = len(pos); s = int(alphaQ[q])
        if (k - s) % 2 or abs(s) > k: return np.zeros((0, P.n), np.int64)
        nm = (k - s) // 2
        opts = []
        for minus in it.combinations(range(k), nm):
            a = np.ones(k, np.int64); a[list(minus)] = -1; opts.append(a)
        per.append((pos, np.array(opts)))
    N = math.prod(len(o) for _, o in per)
    out = np.zeros((N, P.n), np.int64)
    idx = np.indices([len(o) for _, o in per]).reshape(len(per), -1).T
    for j, (pos, opts) in enumerate(per):
        out[:, pos] = opts[idx[:, j]]
    return out

def paf_rows(P, Y, lags):
    return np.stack([(Y * Y[:, P.add[:, h]]).sum(1) for h in lags], 1)

def solve(P, alphaQ, W, max_print=10, verbose=True):
    p = P.p; rep = coset_index(P)
    lags = []                                                    # one representative per {h, -h}, h != 0
    seen = set()
    for h in range(1, P.n):
        if h in seen: continue
        seen.update({h, int(P.G.idx(-P.G.co[h]))}); lags.append(h)
    TJ = (P.T[:, 0] + (p - 1) * P.T[:, 1]).astype(np.int64)
    T = TJ[lags]
    t0 = time.time()
    L = []; Y = []
    for k in range(4):
        A = block_list(P, rep, alphaQ[k], W[k]); Yk = A + (p - 1) * W[k][None, :]
        L.append(A); Y.append(Yk)
    R = [paf_rows(P, Yk, lags) for Yk in Y]
    sizes = [len(r) for r in R]
    if verbose: print(json.dumps({"list_sizes": sizes, "log2_total": round(sum(math.log2(max(1, s)) for s in sizes), 1),
                                  "secs": round(time.time() - t0, 1)}), flush=True)
    order = np.argsort(sizes)                                   # smallest two hashed, third iterated, largest vectorised
    i1, i2, i3, i4 = [int(t) for t in order]
    wts = np.random.default_rng(7).integers(1, 2 ** 62, len(lags), dtype=np.int64).astype(np.uint64)
    key = lambda M: (M.astype(np.uint64) * wts[None, :]).sum(1)
    K12 = (key(R[i1])[:, None] + key(R[i2])[None, :]).reshape(-1)
    srt = np.argsort(K12, kind='stable'); Ks = K12[srt]
    tkey = (T.astype(np.uint64) * wts).sum()
    K4 = key(R[i4]); sols = []
    for j3 in range(sizes[i3]):
        need = tkey - key(R[i3][j3:j3 + 1])[0] - K4
        lo = np.searchsorted(Ks, need, 'left'); hi = np.searchsorted(Ks, need, 'right')
        for j4 in np.nonzero(hi > lo)[0]:
            for h in srt[lo[j4]:hi[j4]]:
                a1, a2 = divmod(int(h), sizes[i2])
                if not (R[i1][a1] + R[i2][a2] + R[i3][j3] + R[i4][j4] == T).all(): continue
                idx = {i1: a1, i2: a2, i3: j3, i4: int(j4)}
                sols.append([L[k][idx[k]] for k in range(4)])
    # full check (all lags incl. h' = 0) and norm identity
    good = []
    for s in sols:
        Ys = [s[k] + (p - 1) * W[k] for k in range(4)]
        tot = sum(np.array([(y * y[P.add[:, h]]).sum() for h in range(P.n)]) for y in Ys)
        if (tot == TJ).all(): good.append(s)
    if verbose: print(json.dumps({"stage1_solutions": len(good), "secs": round(time.time() - t0, 1)}), flush=True)
    return good

def check359():
    v, b, mods, X = load_core(os.path.join(ROOT, 'data/core_v359_b7.txt')); p = 11
    P = LSm.Prob(359, 7, 11, (2, 2, 8), [(0, 0, 2)])
    comps = [components(x, p) for x in X]
    A0 = [c[0].reshape(-1) for c in comps]; W0 = [c[2].reshape(-1) for c in comps]
    rep = coset_index(P)
    alphaQ = [np.array([int(A0[k][rep == q].sum()) for q in range(rep.max() + 1)]) for k in range(4)]
    sols = solve(P, alphaQ, W0)
    pub = any(all((s[k] == A0[k]).all() for k in range(4)) for s in sols)
    nminus = [sum(int(((W0[k] != 0) & (s[k] == -W0[k])).sum()) for k in range(4)) for s in sols]
    out = {"published_alpha_found": bool(pub), "solutions": len(sols), "N_minus_of_solutions": sorted(set(nminus))}
    print(json.dumps(out))
    np.save(os.path.join(ROOT, 'theory', 'logs', 'stage1_359_solutions.npy'), np.array(sols))
    return out

if __name__ == "__main__":
    {"check359": check359}[sys.argv[1]]()
