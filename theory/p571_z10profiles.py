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

# ---------------------------------------------------------------- lag-0 filter [P]
# Each block's reachable (J, rho) set (J = sum |O| [|y| = 3] over both halves, rho = sum |O| y0 y1) inside a profile
# is computed exactly as a 2-D boolean array; the c = 0 / c != 0 orbit groups are independent (OR-convolution), and
# within a group every a = 0 config pairs with every a = 1 config.  A profile survives iff sum_k J_k = J_tot (the
# J-coset count forced by S1 at lag 0) and sum_k rho_k = -6 (S1 at lag q, i.e. lag0(L1) - lag0(L2)) are reachable.
from numpy.fft import rfft2, irfft2

def _orconv(A, B):
    """boolean OR-convolution of 2-D arrays (index offsets add)"""
    s0 = A.shape[0] + B.shape[0] - 1; s1 = A.shape[1] + B.shape[1] - 1
    C = irfft2(rfft2(A.astype(float), (s0, s1)) * rfft2(B.astype(float), (s0, s1)), (s0, s1))
    return C > 0.5

class Lag0:
    def __init__(self, inst):
        from p571_hsolve import group_lists
        self.inst = inst; self.gl = group_lists
        self.w0, self.w1 = fibre_weights(inst)
        self.RMAX = 9 * inst.q                     # |rho| <= 9 q per block
        self.cache = {}; self.pair = {}
    def block(self, p):
        if p in self.cache: return self.cache[p]
        inst = self.inst; arr = None
        for g, w in ((0, self.w0), (1, self.w1)):
            orb = [i for i in range(inst.R) if w[i]]
            V0, _ = self.gl(w, p[0] if g == 0 else p[1]); V1, _ = self.gl(w, p[2] if g == 0 else p[3])
            if len(V0) == 0 or len(V1) == 0: self.cache[p] = None; return None
            sz = inst.sizes[orb]
            J0 = ((np.abs(V0) == 3) * sz).sum(1); J1 = ((np.abs(V1) == 3) * sz).sum(1)
            rho = (V0 * sz) @ V1.T
            Jg = (J0[:, None] + J1[None, :]).ravel(); Rg = rho.ravel()
            A = np.zeros((int(2 * sz.sum()) + 1, 2 * self.RMAX + 1), bool); A[Jg, Rg + self.RMAX] = True
            arr = A if arr is None else _orconv(arr, A)[:, self.RMAX:self.RMAX + 2 * self.RMAX + 1]
        self.cache[p] = arr; return arr
    def pairsum(self, p, q):
        key = (p, q)
        if key not in self.pair:
            A, B = self.block(p), self.block(q)
            self.pair[key] = None if A is None or B is None else _orconv(A, B)     # rho offset 2*RMAX
        return self.pair[key]
    def ok(self, P, Jtot, rhotot=-6):
        AB = self.pairsum(tuple(P[0]), tuple(P[1])); CD = self.pairsum(tuple(P[2]), tuple(P[3]))
        if AB is None or CD is None: return False
        R2 = 2 * self.RMAX
        # need (j, r) with AB[j, r + R2] and CD[Jtot - j, rhotot - r + R2]
        for j in np.nonzero(AB.any(1))[0]:
            jj = Jtot - j
            if jj < 0 or jj >= CD.shape[0]: continue
            ra = np.nonzero(AB[j])[0] - R2; rc = rhotot - ra + R2
            okm = (rc >= 0) & (rc < CD.shape[1])
            if CD[jj, rc[okm]].any(): return True
        return False

def lag0_filter(inst, sols):
    Jtot = (4 * inst.m - 8 - 8 * inst.q) // 8                    # sum y^2 = 2q*4 + 8 J = T_J(0) = 4m - 8
    L = Lag0(inst)
    keep = [(P, w) for P, w in sols if L.ok(P, Jtot)]
    return keep, Jtot
