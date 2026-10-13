r"""EXACT stage 1 of the 971 x81 slice (Z_10 fixed, index-8 classes on Z_97, f = 12) via the alpha reduction [P, M].

Letters y = alpha + 12 sigma, sigma even, so y = alpha + 24 j (j integer).  theory/gp_z2screen.py leaves one lag-0
composition, {23: 3, 25: 3, 1: 34}: j in {0, +-1}, six nonzero j, |y| = 23 iff alpha j = -1, |y| = 25 iff alpha j = +1
(three of each).  [P] PAF(alpha + 24 j)(h) = PAF(alpha)(h) + 24 C(h) + 576 PAF(j)(h), C(h) = sum_g (alpha_g j_{g+h} +
j_g alpha_{g+h}), and C(h) is even (it is = 2 sum j mod 2).  So mod 48 the target -388 forces sum_k PAF(alpha_k)(h) = -4
(|PAF| <= 40): the alpha_k form a b = 1 core on Z_10, with rows exactly (2, 0, 0, 0), so sum_g j_k(g) = 0 per block.
Then the remaining condition is sum_k v_k(h) = -384 with v_k = PAF(y_k) - PAF(alpha_k), computed directly per block.
Exhaustive: all alpha 4-tuples (2 + 2 join over all 2^10 sequences per block, block 0 row +2, the rest row 0), then
per alpha tuple all j (per block: an even number of nonzeros with sum 0), joined on (sum v, #(alpha j = +1), #nonzero).
alpha is taken canonical per block (shift, reversal, negation of row-0 blocks act on (alpha_k, j_k) jointly and
preserve every equation).  Every solution y is re-verified on Z_10.  `--selftest` plants random j on random alpha
cores (the target becomes the planted sum v) and checks recovery.
usage: python3 theory/q971_stage1.py [--max-alpha N] [--out file]"""
import sys, json, time, itertools as it
from collections import defaultdict
import numpy as np

F, q, f = 10, 97, 12
H = F // 2

def paf(a, h): return int(sum(a[g] * a[(g + h) % F] for g in range(F)))

def canon(a):
    """orbit representative of a block under cyclic shift, reversal and (row 0 only) negation; these act on
    (alpha_k, j_k) together and preserve every stage-1 equation [P], so alpha may be taken canonical per block"""
    c = []
    for b in (a, a[::-1]):
        for t in range(F):
            r = b[t:] + b[:t]; c.append(r)
            if sum(a) == 0: c.append(tuple(-x for x in r))
    return min(c)

def alpha_tuples():
    seqs = [tuple(1 if (n >> i) & 1 else -1 for i in range(F)) for n in range(1 << F)]
    by = defaultdict(list)
    for a in seqs:
        if a != canon(a): continue          # per-block symmetry: shift, reversal, negation of row-0 blocks
        by[(sum(a),) + tuple(paf(a, h) for h in range(1, H + 1))].append(a)
    k2 = [k for k in by if k[0] == 2]; k0 = [k for k in by if k[0] == 0]
    pair = defaultdict(list)
    for a, b in it.product(k0, k0): pair[tuple(x + y for x, y in zip(a[1:], b[1:]))].append((a, b))
    out = []
    for a in k2:
        for d in k0:
            need = tuple(-4 - x - y for x, y in zip(a[1:], d[1:]))
            for (b, c) in pair.get(need, []): out.append((a, b, c, d))
    return by, out      # key-level tuples (blocks' sequence lists in `by`)

def j_options(a):
    """per block: list of (vec C-part + 12 PAF(j) over h = 1..5, n_plus, n_nonzero, j)"""
    res = []
    for nz in (0, 2, 4, 6):
        for pos in it.combinations(range(F), nz):
            for sg in it.product((1, -1), repeat=nz):
                if sum(sg): continue
                j = [0] * F
                for p, e in zip(pos, sg): j[p] = e
                y = [a[g] + 24 * j[g] for g in range(F)]
                v = tuple(paf(y, h) - paf(a, h) for h in range(1, H + 1))
                npl = sum(1 for g in range(F) if a[g] * j[g] == 1)
                res.append((v, npl, nz, tuple(j)))
    return res

def main(max_alpha=None, out="", target=None, only=None, verify_full=True):
    t0 = time.time(); by, keyt = alpha_tuples()
    n_alpha = sum(len(by[a]) * len(by[b]) * len(by[c]) * len(by[d]) for a, b, c, d in keyt)
    jcache = {}
    def jo(a):
        if a not in jcache:
            d = defaultdict(list)
            for v, npl, nz, j in j_options(a): d[(v, npl, nz)].append(j)
            jcache[a] = d
        return jcache[a]
    sols = []; checked = 0
    target = tuple([-384] * H) if target is None else target
    for kt in (keyt if only is None else [only[0]]):
        for A in (it.product(*[by[k] for k in kt]) if only is None else [only[1]]):
            checked += 1
            if max_alpha and checked > max_alpha: break
            O = [jo(a) for a in A]
            L = defaultdict(list)
            for (v0, p0, n0) in O[0]:
                for (v1, p1, n1) in O[1]:
                    if n0 + n1 > 6 or p0 + p1 > 3: continue
                    L[(tuple(x + y for x, y in zip(v0, v1)), p0 + p1, n0 + n1)].append(((v0, p0, n0), (v1, p1, n1)))
            for (v2, p2, n2) in O[2]:
                for (v3, p3, n3) in O[3]:
                    nn = 6 - n2 - n3; pp = 3 - p2 - p3
                    if nn < 0 or pp < 0: continue
                    need = (tuple(t - x - y for t, x, y in zip(target, v2, v3)), pp, nn)
                    for (k0, k1) in L.get(need, []):
                        J = (O[0][k0][0], O[1][k1][0], O[2][(v2, p2, n2)][0], O[3][(v3, p3, n3)][0])
                        Y = [[A[k][g] + 24 * J[k][g] for g in range(F)] for k in range(4)]
                        ok = (not verify_full) or all(sum(paf(y, h) for y in Y) == -4 * q for h in range(1, H + 1)) and \
                             sum(paf(y, 0) for y in Y) == 4 * F * q - 4 * (q - 1) and [sum(y) for y in Y] == [2, 0, 0, 0]
                        sols.append(dict(y=Y, verified=ok))
                        if len(sols) >= 50: break
        if max_alpha and checked > max_alpha: break
    r = dict(v=971, F=F, q=q, s=8, alpha_key_tuples=len(keyt), alpha_tuples=n_alpha, alpha_checked=checked,
             solutions=len(sols), all_verified=all(s['verified'] for s in sols), complete=not (max_alpha and checked > max_alpha),
             secs=round(time.time() - t0, 1), examples=sols[:3])
    if out: json.dump(r, open(out, "w"))
    return r

def selftest(trials=20, seed=1):
    """plant: a canonical alpha tuple, random j (6 nonzeros, 3 + 3 by alpha j, sum 0 per block); target := its sum v"""
    rng = np.random.default_rng(seed); by, keyt = alpha_tuples(); ok = 0
    for t in range(trials):
        kt = keyt[rng.integers(len(keyt))]; A = tuple(by[k][0] for k in kt)
        while True:
            J = [[0] * F for _ in range(4)]
            for k, nz in zip(range(4), rng.permutation([2, 2, 2, 0])):
                pos = rng.choice(F, nz, replace=False)
                for i, p in enumerate(pos): J[k][p] = 1 if i % 2 == 0 else -1
            npl = sum(1 for k in range(4) for g in range(F) if A[k][g] * J[k][g] == 1)
            if npl == 3: break
        Y = [[A[k][g] + 24 * J[k][g] for g in range(F)] for k in range(4)]
        tgt = tuple(sum(paf(Y[k], h) - paf(A[k], h) for k in range(4)) for h in range(1, H + 1))
        r = main(target=tgt, only=(kt, A), verify_full=False)
        ok += r['solutions'] > 0
    return dict(selftest_trials=trials, recovered=ok)

if __name__ == "__main__" and "--selftest" in sys.argv:
    print(json.dumps(selftest())); sys.exit()
if __name__ == "__main__":
    a = sys.argv[1:]; mx = int(a[a.index('--max-alpha') + 1]) if '--max-alpha' in a else None
    out = a[a.index('--out') + 1] if '--out' in a else ""
    r = main(mx, out); r.pop('examples', None) if not r['solutions'] else None
    print(json.dumps(r))
