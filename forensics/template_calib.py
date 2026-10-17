"""Calibrate the Gaussian count of half-shift template halves against exact counts (small h).

For a support pattern sigma (S_k = {x_k(y) = x_k(y+h)}, D_k = complement, on Z_h, h odd), the b = 1 core splits into
    S-half: c_k = +-1 on S_k: sum_k PAF_{c_k}(t) = -2, t = 1..n;   D-half: d_k = +-1 on D_k: sum_k NPAF_{d_k}(t) = 0
(n = (h-1)/2; forensics/template_opt.py).  Here: parity-feasible templates (four arcs anchored at 0 with lengths near
h/2, |S_A| odd, |S_B|, |S_C|, |S_D| even, sum 2h - 1, plus 0-2 flipped positions, every V_t even), the exact number of
solutions of each half (all sign choices; 2+2 join with a 64-bit linear hash of the lag vector, so a collision could
only add), and the Gaussian of template_opt.evaluate.  Output: log2 exact vs log2 Gaussian per half.
usage: python3 forensics/template_calib.py [h1,h2,...] [--per-h K]"""
import os, sys, json, itertools, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from template_opt import lagterm_tables, evaluate, profile

W = np.random.default_rng(99).integers(1, 2 ** 62, 64).astype(np.uint64) | np.uint64(1)


def lagvecs(h, supp, nega):
    """all sign vectors on the support (first support element fixed to +1; caller multiplies by 2), lag vectors 1..n"""
    n = (h - 1) // 2; L = len(supp); N = 1 << (L - 1)
    bits = ((np.arange(N)[:, None] >> np.arange(L - 1)[None, :]) & 1).astype(np.int8)
    signs = np.concatenate([np.ones((N, 1), np.int8), 1 - 2 * bits], 1)
    X = np.zeros((N, h), np.int8); X[:, supp] = signs
    out = np.zeros((N, n), np.int32)
    for t in range(1, n + 1):
        sh = np.roll(X, -t, axis=1).astype(np.int32)
        if nega: sh[:, h - t:] *= -1
        out[:, t - 1] = (X.astype(np.int32) * sh).sum(1)
    return out


def hkey(V):
    return (V.astype(np.int64).astype(np.uint64) * W[:V.shape[1]][None, :]).sum(1, dtype=np.uint64)


def count_half(h, sets, nega, target):
    """exact number of sign choices on the four supports meeting sum_k lag_k(t) = target, t = 1..n"""
    vecs = [lagvecs(h, s, nega) if len(s) else np.zeros((1, (h - 1) // 2), np.int32) for s in sets]
    mult = np.prod([2 if len(s) else 1 for s in sets])
    def pairsums(a, b):
        ka, kb = hkey(a), hkey(b)
        tot = (ka[:, None] + kb[None, :]).ravel()
        u, c = np.unique(tot, return_counts=True); return u, c
    uA, cA = pairsums(vecs[0], vecs[1]); uB, cB = pairsums(vecs[2], vecs[3])
    tk = hkey(np.full((1, (h - 1) // 2), target, np.int32))[0]
    need = tk - uB                                                      # key(AB) must equal key(target) - key(CD)
    pos = np.searchsorted(uA, need); pos = np.minimum(pos, len(uA) - 1); hit = uA[pos] == need
    return float((cA[pos[hit]].astype(np.float64) * cB[hit]).sum()) * mult


def templates(h, per_h):
    """parity-feasible templates: arcs [0, L) per block with 0, 1 or 2 flipped positions"""
    n = (h - 1) // 2; Ls = range(n - 1, n + 3)                            # balanced arcs, as in the published cores
    def par(s): return tuple(int((s * np.roll(s, -t)).sum()) & 1 for t in range(1, n + 1))
    cands = []
    for L in Ls:
        base = np.zeros(h, np.int8); base[:L] = 1
        cands.append((base, 0))
        for y in range(h):
            s = base.copy(); s[y] ^= 1; cands.append((s, 1))
    want = [1, 0, 0, 0]
    byp = {}
    for s, f in cands:
        byp.setdefault((int(s.sum()) & 1), []).append((s, f, par(s)))
    found = []
    rng = np.random.default_rng(h)
    # search quadruples with small total flips: brute force over pairs (A,B) and (C,D) by parity XOR
    A = byp.get(1, []); B = byp.get(0, [])
    right = {}
    for i, (c2, f2, p2) in enumerate(B):
        for j, (c3, f3, p3) in enumerate(B):
            if f2 + f3 > 2: continue
            k = tuple(a ^ b for a, b in zip(p2, p3)); right.setdefault(k, []).append((c2, c3, f2 + f3))
    order = list(range(len(A))); rng.shuffle(order)
    for i in order:
        c0, f0, p0 = A[i]
        for c1, f1, p1 in B:
            if f0 + f1 > 2: continue
            k = tuple(a ^ b for a, b in zip(p0, p1))
            for c2, c3, f23 in right.get(k, []):
                S = np.array([c0, c1, c2, c3])
                if int(S.sum()) != 2 * h - 1 or f0 + f1 + f23 > 3: continue
                V = profile(S)[1:n + 1]
                if (V % 2).any() or (V < 2).any(): continue
                if any(np.array_equal(S, F_) for F_, _ in found): continue
                found.append((S, f0 + f1 + f23))
                if len(found) >= per_h: return found
    return found


def main():
    a = sys.argv[1:]
    hs = [int(x) for x in a[0].split(',')] if a and not a[0].startswith('--') else [11, 15, 19, 21, 23]
    per_h = int(a[a.index('--per-h') + 1]) if '--per-h' in a else 4
    out = []
    for h in hs:
        fS, fD = lagterm_tables(8 * h * h)
        T = templates(h, per_h)
        for S, flips in T:
            sets_S = [np.nonzero(S[k])[0] for k in range(4)]; sets_D = [np.nonzero(1 - S[k])[0] for k in range(4)]
            eS = count_half(h, sets_S, False, -2); eD = count_half(h, sets_D, True, 0)
            ev = evaluate(S, fS, fD)
            rec = dict(h=h, v=2 * h + 1, sizes=[int(x) for x in S.sum(1)], flips=flips,
                       exact_S=round(math.log2(eS), 2) if eS > 0 else None, gauss_S=round(ev['logES'], 2),
                       exact_D=round(math.log2(eD), 2) if eD > 0 else None, gauss_D=round(ev['logED'], 2))
            out.append(rec); print(json.dumps(rec), flush=True)
    json.dump(out, open(os.path.join(HERE, '..', 'logs', 'gap', 'template_calib.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
