r"""Independent check (separate code path from theory/p571_z10profiles.Lag0) of the stage-1 obstruction for the
571 x27 slice: no y on Z190 satisfies S1 at the four Z10 characters and at lags 0 and 95 simultaneously.

Setting [P]: G = Z570 = Z2 x Z3 x Z95, slice M = <103> (= x27 on Z95, identity on Z6).  For a b = 1 core the Z3-coset
sums y_k(a, t) (a in Z2, t in Z95) are in {+-1, +-3}, constant on M-orbits of Z95, and satisfy S1:
sum_k PAF(y_k)(h) = T_J(h) on Z190 (T_J(0) = 4m - 8 = 2272, T_J(h) = -12 otherwise).  Six of these conditions only see
per-block "summary" quantities:
  * the four Z10 characters (rational under M): with Y(a,0) / Y(a,*) the sums of y(a, .) over the c = 0 fibre / one
    c != 0 fibre,  rows (i), theta_5 (ii), (-1) x theta_0 (iii), (-1) x theta_5 (iv) as in p571_z10profiles.py;
  * lag 0:  sum_k J_k = 189   (J_k = #positions with |y| = 3:  sum y^2 = 760 + 8 J);
  * lag 95 (the Z2-shift): sum_k rho_k = -6   (rho_k = sum_t y_k(0, t) y_k(1, t)).
This script recomputes, from the orbit list alone, per block and per Z10-summary (Y00, Y0s, Y10, Y1s) the set of
reachable (J mod K, rho mod K) by a DP over orbits (each orbit: 16 choices of (y0, y1), contributing
(fibre weight * y0 to Y(0, .), fibre weight * y1 to Y(1, .), |O| [|y0| = 3] + |O| [|y1| = 3] to J, |O| y0 y1 to rho)),
enumerates all Z10 profiles satisfying (i)-(iv), and tests whether (189, -6) mod K is reachable.  It also reports the
exact (non-modular) result when K is large enough to make the residues exact (K > 2 max |value|).
usage: python3 theory/p571_empty_check.py [K]"""
import sys, os, json, time, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_slicelib import mgroup, orbits

def setup(m=570, gens=(103,)):
    q = m // 6; G = mgroup(m, list(gens)); Gq = sorted({t % q for t in G})
    orbs = orbits(q, Gq)
    w0 = [sum(1 for t in o if t % 5 == 0) for o in orbs]            # c = 0 fibre
    w1 = [sum(1 for t in o if t % 5 == 1) for o in orbs]            # c = 1 fibre (all c != 0 alike)
    return q, orbs, w0, w1

YV = (1, -1, 3, -3)

def group_table(orbs, w, K):
    """DP over the orbits with w > 0: dict (S0, S1) -> boolean array (K, K) of reachable (J mod K, rho mod K)"""
    idx = [i for i, x in enumerate(w) if x]
    T = {(0, 0): np.zeros((K, K), bool)}; T[(0, 0)][0, 0] = True
    for i in idx:
        sz = len(orbs[i]); wi = w[i]; N = {}
        for (s0, s1), A in T.items():
            for y0 in YV:
                for y1 in YV:
                    key = (s0 + wi * y0, s1 + wi * y1)
                    dj = (sz * ((abs(y0) == 3) + (abs(y1) == 3))) % K; dr = (sz * y0 * y1) % K
                    B = np.roll(np.roll(A, dj, 0), dr, 1)
                    if key in N: N[key] |= B
                    else: N[key] = B.copy()
        T = N
    return T

def cyc_or(A, B):
    """cyclic OR-convolution on Z_K x Z_K"""
    K = A.shape[0]; out = np.zeros_like(A)
    for (j, r) in zip(*np.nonzero(A)):
        out |= np.roll(np.roll(B, j, 0), r, 1)
    return out

def main(K=48, m=570, gens=(103,), v=571, rows=(2, 0, 0, 0)):
    t0 = time.time(); q, orbs, w0, w1 = setup(m, gens)
    T0 = group_table(orbs, w0, K); T1 = group_table(orbs, w1, K)
    V0 = sorted({s for s, _ in T0}); V1 = sorted({s for s, _ in T1})
    Jtot = (4 * m - 8 - 8 * q) // 8; target = (Jtot % K, (-6) % K)
    # per-block summaries satisfying the row sum, and their mod-K (J, rho) sets (combined c = 0 and c != 0 groups)
    cache = {}
    def bset(Y00, Y0s, Y10, Y1s):
        key = (Y00, Y0s, Y10, Y1s)
        if key not in cache:
            A = T0.get((Y00, Y10)); B = T1.get((Y0s, Y1s))
            cache[key] = None if A is None or B is None else cyc_or(A, B)
        return cache[key]
    Tq = 4 * v; per = []
    vals0 = sorted({a for a, b in T0} | {b for a, b in T0}); vals1 = sorted({a for a, b in T1} | {b for a, b in T1})
    for r in rows:
        lst = []
        for (Y00, Y10) in T0:
            rem = r - Y00 - Y10
            if rem % 4: continue
            s = rem // 4
            for (Y0s, Y1s) in T1:
                if Y0s + Y1s != s: continue
                c2 = (Y00 + Y10) - (Y0s + Y1s); c3 = (Y00 - Y10) + 4 * (Y0s - Y1s); c4 = (Y00 - Y10) - (Y0s - Y1s)
                if max(c2 * c2, c3 * c3, c4 * c4) > Tq: continue
                lst.append((c2 * c2, c3 * c3, c4 * c4, (Y00, Y0s, Y10, Y1s)))
        per.append(lst)
    from collections import defaultdict
    idx3 = defaultdict(list)
    for e in per[3]: idx3[e[:3]].append(e[3])
    n_prof = 0; n_ok = 0; examples = []
    for e0 in per[0]:
        for e1 in per[1]:
            a = (e0[0] + e1[0], e0[1] + e1[1], e0[2] + e1[2])
            if max(a) > Tq: continue
            AB = None
            for e2 in per[2]:
                b = (Tq - a[0] - e2[0], Tq - a[1] - e2[1], Tq - a[2] - e2[2])
                if min(b) < 0: continue
                for P3 in idx3.get(b, []):
                    n_prof += 1
                    if AB is None: AB = cyc_or(bset(*e0[3]), bset(*e1[3]))
                    CD = cyc_or(bset(*e2[3]), bset(*P3))
                    # reachable target: exists (j, r) in AB with (target - (j, r)) in CD
                    ok = False
                    for (j, rr) in zip(*np.nonzero(AB)):
                        if CD[(target[0] - j) % K, (target[1] - rr) % K]: ok = True; break
                    if ok:
                        n_ok += 1
                        if len(examples) < 3: examples.append([e0[3], e1[3], e2[3], P3])
    print(json.dumps(dict(K=K, fibre_c0=w0, fibre_c1=w1, Jtot=Jtot, profiles=n_prof, residue_feasible=n_ok,
                          examples=examples, secs=round(time.time() - t0, 1))), flush=True)

if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    main(K)
