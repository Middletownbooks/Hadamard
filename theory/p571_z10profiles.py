r"""Exact Z10-compression profiles of S1 (the y-level) for slices of Z_m = Z2 x Z3 x Z5 x Z_q' fixing Z6 pointwise and
mapping onto Z5^* (571 x27: q' = 19; v = 211 analogue: q' = 7).  [P]

y_k (a, t) in {+-1, +-3} (a in Z2, t in Z_q, q = 5q').  Its Z10-compression is determined by four odd integers per block
    Y(a, 0) = sum of y(a, .) over the c = 0 fibre (orbits o0, oC_i)        (q' terms)
    Y(a, *) = sum of y(a, .) over one c != 0 fibre (orbits o5, mixed)        (q' terms, equal for all c != 0)
The four characters of Z10 (the rational characters of Z_2q under M) give, for a b = 1 core with rows (2,0,0,0):
    (i)   r_k = Y(0,0) + Y(1,0) + 4 (Y(0,*) + Y(1,*))                    = (2, 0, 0, 0)_k        [trivial]
    (ii)  sum_k [ (Y(0,0) + Y(1,0)) - (Y(0,*) + Y(1,*)) ]^2               = 4v                    [theta_5, L1]
    (iii) sum_k [ (Y(0,0) - Y(1,0)) + 4 (Y(0,*) - Y(1,*)) ]^2            = 4v                    [(-1) x theta_0, L2]
    (iv)  sum_k [ (Y(0,0) - Y(1,0)) - (Y(0,*) - Y(1,*)) ]^2               = 4v                    [(-1) x theta_5, L2]
(the compressed S1: sum_k |Y^_k(chi)|^2 = T_J^(chi) = 4v at every nontrivial character chi of Z10).
Each Y is a weighted sum over orbits (c = 0 fibre: weights = #positions of the orbit in the fibre; same for c != 0),
so the reachable values factor over the two orbit groups.  This script enumerates ALL integer solutions with
reachable values, i.e. the complete list of compressed stage-1 profiles (an exact necessary condition).
usage: python3 theory/p571_z10profiles.py m gen   (e.g. 570 103, 210 97)"""
import sys, os, json, math, itertools as it
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst

def fibre_weights(inst):
    """per orbit of Z_q: (#positions with c = 0, #positions with c = c1 for a fixed c1 != 0)"""
    q = inst.q
    w0 = [sum(1 for t in o if t % 5 == 0) for o in inst.S.orbs]
    w1 = [sum(1 for t in o if t % 5 == 1) for o in inst.S.orbs]
    return w0, w1

def reach(weights):
    """value -> number of weighted y-configs (x-multiplicity 3 for |y| = 1) on the orbits with nonzero weight;
    orbits with weight 0 in this fibre are free (summed separately)"""
    D = {0: 1}
    for w in weights:
        if w == 0: continue
        E = defaultdict(int)
        for s, c in D.items():
            for yv, mult in ((1, 3), (-1, 3), (3, 1), (-3, 1)): E[s + w * yv] += c * mult
        D = E
    return dict(D)

def profiles(inst, rows=(2, 0, 0, 0)):
    v = inst.v; T = 4 * v
    w0, w1 = fibre_weights(inst)
    # c = 0 fibre and c != 0 fibre use disjoint orbit sets (an orbit lies in c = 0 or c != 0 part)
    assert all((a == 0) or (b == 0) for a, b in zip(w0, w1))
    R0 = reach(w0); R1 = reach(w1)               # values of Y(a, 0) and Y(a, *)
    V0 = sorted(R0); V1 = sorted(R1)
    # per block: all (Y00, Y0s, Y10, Y1s) with row sum r -> list of (q2, q3, q4) = squares of the three char values
    per = []
    for r in rows:
        lst = []
        for Y00 in V0:
            for Y10 in V0:
                rem = r - Y00 - Y10
                if rem % 4: continue
                s = rem // 4                         # = Y0s + Y1s
                for Y0s in V1:
                    Y1s = s - Y0s
                    if Y1s not in R1: continue
                    c2 = (Y00 + Y10) - (Y0s + Y1s); c3 = (Y00 - Y10) + 4 * (Y0s - Y1s); c4 = (Y00 - Y10) - (Y0s - Y1s)
                    if c2 * c2 > T or c3 * c3 > T or c4 * c4 > T: continue
                    wgt = R0[Y00] * R0[Y10] * R1[Y0s] * R1[Y1s]
                    lst.append((c2 * c2, c3 * c3, c4 * c4, (Y00, Y0s, Y10, Y1s), wgt))
        per.append(lst)
    idx = defaultdict(list)
    for e in per[3]: idx[(e[0], e[1], e[2])].append(e)
    sols = []
    for e0 in per[0]:
        for e1 in per[1]:
            a = (e0[0] + e1[0], e0[1] + e1[1], e0[2] + e1[2])
            if a[0] > T or a[1] > T or a[2] > T: continue
            for e2 in per[2]:
                b = (T - a[0] - e2[0], T - a[1] - e2[1], T - a[2] - e2[2])
                if min(b) < 0: continue
                for e3 in idx.get(b, []):
                    sols.append(((e0[3], e1[3], e2[3], e3[3]), e0[4] * e1[4] * e2[4] * e3[4]))
    return sols, (len(per[0]), len(per[1]))

if __name__ == "__main__":
    m, g = int(sys.argv[1]), [int(x) for x in sys.argv[2:]]
    inst = Inst(m, g)
    sols, sizes = profiles(inst)
    w0, w1 = fibre_weights(inst)
    tot = sum(w for _, w in sols)
    print(json.dumps(dict(v=inst.v, fibre_weights_c0=w0, fibre_weights_c1=w1, per_block_candidates=sizes,
                          profiles=len(sols), log2_weighted=(round(math.log2(tot), 1) if tot else None))), flush=True)
    sols.sort(key=lambda t: -t[1])
    for P, w in sols[:10]: print(json.dumps(dict(profile=P, log2w=round(math.log2(w), 1))))
