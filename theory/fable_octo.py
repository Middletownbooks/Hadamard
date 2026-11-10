#!/usr/bin/env python3
"""All monomial (4,4,4) compositions over C with conjugates:
    z_j = sum_i  eps_ji * v_i^(e_ji) * w_{pi_j(i)}^(f_ji),   j = 1..4,
    sum_j |z_j|^2 = |v|^2 |w|^2   identically in v, w in C^4      (x^(+) = x, x^(-) = conj x).
If V = (V_1..V_4) and W = (W_1..W_4) are +-1 quadruples on groups G_t, G_v with PSD sums 4t and 4v, then
Z_j = (1/2) sum_i eps V_i(e tau) W_k(f g) has PSD sum 4tv, and Z_j is +-1 valued iff
    prod_i V_i(e_ji tau) = const  and  prod_i W_{pi_j(i)}(f_ji g) = const       for each j.
This script enumerates the compositions exhaustively (meet in the middle with 64-bit hashes, exact re-check)
and prints the distinct pairs of conjugation-pattern sets (for V and for W), up to the obvious symmetries.  [M]
"""
import itertools, sys, random
import numpy as np

random.seed(1)
# symbols: v_i^+ -> i, v_i^- -> 4+i ; same for w.  monomial = (sorted pair of v symbols, sorted pair of w symbols)
def vs(i, e): return i if e > 0 else 4 + i
MON = {}
def mon(a, b, c, d):
    key = (min(a, b), max(a, b), min(c, d), max(c, d))
    if key not in MON:
        MON[key] = random.getrandbits(62)
    return key

def row_terms(perm, eps, e, f):
    """off-diagonal terms of |z|^2 as dict monomial -> coeff ; z = sum_i eps_i v_i^{e_i} w_{perm[i]}^{f_i}"""
    out = {}
    for i in range(4):
        for j in range(4):
            if i == j: continue
            # term_i * conj(term_j)
            m = mon(vs(i, e[i]), vs(j, -e[j]), vs(perm[i], f[i]), vs(perm[j], -f[j]))
            out[m] = out.get(m, 0) + eps[i] * eps[j]
    return {k: c for k, c in out.items() if c}

def latin_squares():
    perms = list(itertools.permutations(range(4)))
    idp = (0, 1, 2, 3)
    res = []
    for p2 in perms:
        if any(p2[i] == idp[i] for i in range(4)): continue
        for p3 in perms:
            if any(p3[i] in (idp[i], p2[i]) for i in range(4)): continue
            for p4 in perms:
                if any(p4[i] in (idp[i], p2[i], p3[i]) for i in range(4)): continue
                if p2 < p3 < p4:
                    res.append((idp, p2, p3, p4))
    return res

def rows_for(perm):
    rows = []
    for eps in itertools.product((1, -1), repeat=3):
        eps = (1,) + eps
        for e in itertools.product((1, -1), repeat=3):
            e = (1,) + e                     # whole-row conjugation normalised: v_1 unconjugated
            for f in itertools.product((1, -1), repeat=4):
                t = row_terms(perm, eps, e, f)
                rows.append((eps, e, f, t))
    return rows

def hsh(t):
    h = sum(c * MON[m] for m, c in t.items()) % (1 << 64)
    return h - (1 << 64) if h >= (1 << 63) else h

def main():
    sols = []
    classes = {}
    for L in latin_squares():
        R = [rows_for(p) for p in L]
        H = [np.array([hsh(r[3]) for r in Rj], dtype=np.int64) for Rj in R]
        n = len(R[0])
        np.seterr(over='ignore')
        A = (H[0][:, None] + H[1][None, :]).reshape(-1)
        B = -(H[2][:, None] + H[3][None, :]).reshape(-1)
        oa = np.argsort(A, kind='stable'); ob = np.argsort(B, kind='stable')
        As = A[oa]; Bs = B[ob]
        common = np.intersect1d(As, Bs)
        cnt = 0
        for val in common:
            ia = oa[np.searchsorted(As, val, 'left'):np.searchsorted(As, val, 'right')]
            ib = ob[np.searchsorted(Bs, val, 'left'):np.searchsorted(Bs, val, 'right')]
            for x in ia:
                for y in ib:
                    idx = (x // n, x % n, y // n, y % n)
                    tot = {}
                    for j in range(4):
                        for m, c in R[j][idx[j]][3].items():
                            tot[m] = tot.get(m, 0) + c
                    if any(tot.values()): continue
                    cnt += 1
                    rows = [R[j][idx[j]] for j in range(4)]
                    # conjugation patterns: for V: e_j ; for W: g_j[k] = f at the slot i with perm[i] = k
                    pv = []; pw = []
                    for j in range(4):
                        eps, e, f, _ = rows[j]
                        g = [0] * 4
                        for i in range(4): g[L[j][i]] = f[i]
                        # patterns up to global sign (tau -> -tau)
                        e2 = tuple(e) if e[0] > 0 else tuple(-x for x in e)
                        g2 = tuple(g) if g[0] > 0 else tuple(-x for x in g)
                        pv.append(e2); pw.append(g2)
                    key = (tuple(sorted(set(pv))), tuple(sorted(set(pw))))
                    if key not in classes:
                        classes[key] = (L, [(r[0], r[1], r[2]) for r in rows])
        print('latin', L[1:], 'compositions', cnt, flush=True)
    print('distinct (V-pattern set, W-pattern set):', len(classes))
    for (pv, pw), ex in sorted(classes.items(), key=lambda kv: (len(kv[0][0]) + len(kv[0][1]), kv[0])):
        s = lambda P: ' '.join(''.join('+' if x > 0 else '-' for x in p) for p in P)
        print(f"  V: {s(pv):24s} | W: {s(pw)}")
    return classes

if __name__ == '__main__':
    main()
