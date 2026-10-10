"""EXACT stage 1 for a generalised-Paley slice (b = 1): G = Z_F x Z_q, M fixes Z_F pointwise and acts on Z_q by the
index-s subgroup.  Stage 1 = coset sums y_k on Z_F with letters y = a + f * (sum of s signs), f = (q - 1)/s:
    sum_k PAF_{Z_F}(y_k)(h) = T_J(h)   (T_J(0) = 4m - 4(q - 1), T_J(h != 0) = -4q),   row sums (+-2, 0, 0, 0).
Exhaustive: lag 0 fixes the multiset of |y|-values (enumerated exactly); for each split of the large letters over the
four blocks, per-block lists of all y with that letter multiset are joined (2+2, exact on (row sum, PAF at lags
0..F/2)).  Returns every stage-1 solution up to the list order (blocks are labelled; no symmetry reduction).
Feasible for small F (F = 8 at 537).  usage: cyclo_stage1_exact.py F q s [--max N]"""
import sys, itertools as it, json, time
from collections import Counter, defaultdict
import numpy as np

def letters(q, s):
    f = (q - 1) // s
    return sorted({a + f * sum(c) for a in (1, -1) for c in it.product((1, -1), repeat=s)}), f

def paf(y):
    F = len(y); return tuple(int(sum(y[i] * y[(i + h) % F] for i in range(F))) for h in range(F // 2 + 1))

def main(F, q, s, maxsol=10 ** 6, only=None):
    m = F * q; T0 = 4 * m - 4 * (q - 1); Th = -4 * q
    L, f = letters(q, s); mags = sorted({abs(x) for x in L}); big = [a for a in mags if a > 1]
    t0 = time.time()
    # lag 0: total of 4F entries; counts n_a of each magnitude a with sum n_a a^2 = T0
    combos = []
    def rec(i, left, rem, cur):
        if i == len(big):
            if rem == left: combos.append({**cur, 1: left})            # rest are |y| = 1
            return
        a = big[i]; c = 0
        while c <= left and c * a * a + (left - c) <= rem + 0:
            rec(i + 1, left - c, rem - c * a * a, {**cur, a: c}); c += 1
    rec(0, 4 * F, T0, {})
    if only is not None: combos = [c for c in combos if all(c.get(a, 0) == only.get(a, 0) for a in big)]
    out = dict(F=F, q=q, s=s, letters=L, lag0_multisets=[{str(k): v for k, v in c.items() if v} for c in combos])
    sols = []
    target = tuple([T0] + [Th] * (F // 2))
    for c in combos:
        bigs = [a for a in big for _ in range(c[a])]
        if len(bigs) > 4 * F: continue
        # distribute the big magnitudes over blocks: multiset partitions into 4 labelled blocks (each <= F entries)
        dists = set()
        for assign in it.product(range(4), repeat=len(bigs)):
            per = [tuple(sorted(a for a, b in zip(bigs, assign) if b == k)) for k in range(4)]
            if all(len(p) <= F for p in per): dists.add(tuple(per))
        for per in dists:
            lists = []
            for k in range(4):
                bl = per[k]; L_k = defaultdict(list)
                nb = len(bl)
                for pos in it.combinations(range(F), nb):
                    for perm in set(it.permutations(bl)):
                        for sg in it.product((1, -1), repeat=nb):
                            for small in it.product((1, -1), repeat=F - nb):
                                y = [0] * F; si = 0
                                for i in range(F):
                                    if i in pos: y[i] = perm[pos.index(i)] * sg[pos.index(i)]
                                    else: y[i] = small[si]; si += 1
                                L_k[(sum(y),) + paf(y)].append(tuple(y))
                lists.append(L_k)
            # join: rows must be a permutation/sign pattern of (2,0,0,0): enumerate row-sum choices
            AB = defaultdict(list)
            for k0, v0 in lists[0].items():
                for k1, v1 in lists[1].items():
                    AB[tuple(a + b for a, b in zip(k0[1:], k1[1:]))].append((k0, k1))
            for k2, v2 in lists[2].items():
                for k3, v3 in lists[3].items():
                    need = tuple(t - a - b for t, a, b in zip(target, k2[1:], k3[1:]))
                    for (k0, k1) in AB.get(need, []):
                        rs = sorted(abs(x) for x in (k0[0], k1[0], k2[0], k3[0]))
                        if rs != [0, 0, 0, 2]: continue
                        for y0 in lists[0][k0][:2]:
                            for y1 in lists[1][k1][:2]:
                                for y2 in lists[2][k2][:2]:
                                    for y3 in lists[3][k3][:2]:
                                        sols.append((y0, y1, y2, y3))
                        if len(sols) >= maxsol: break
    # verify
    ok = all(tuple(sum(paf(y)[h] for y in Y) for h in range(F // 2 + 1)) == target for Y in sols)
    out.update(solutions_found=len(sols), all_verified=ok, examples=[list(map(list, Y)) for Y in sols[:3]],
               secs=round(time.time() - t0, 1))
    return out

if __name__ == "__main__":
    F, q, s = map(int, sys.argv[1:4])
    print(json.dumps(main(F, q, s)))
