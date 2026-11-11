#!/usr/bin/env python3
"""T-matrices of odd order t with the three 'half-amicability' relations that make the Williamson pattern a
multiplication for the type (A free, B, C, D symmetric):
    Z1 = T1 A + T2 B + T3 C + T4 D,  Z2 = T2 A - T1 B + T4 C - T3 D,
    Z3 = T3 A - T4 B - T1 C + T2 D,  Z4 = T4 A + T3 B - T2 C - T1 D
is a GS quadruple for every such (A,B,C,D) iff (T1..T4) are T-matrices and, with <X,Y>(s) the periodic
cross-correlation, the antisymmetric parts satisfy
    asym(T1 T2*) = asym(T3 T4*),  asym(T1 T3*) = -asym(T2 T4*),  asym(T1 T4*) = asym(T2 T3*).
Exhaustive over all 8^t (block, sign) assignments (cyclic T-matrices), numba.  [M]
Also checks every signed permutation variant of the block roles (the 3 relations up to relabelling T's).
"""
import sys, itertools, time
import numpy as np
from numba import njit

@njit(cache=True)
def search(t, maxout):
    out = np.zeros((maxout, t), dtype=np.int8)
    nout = 0; ntm = 0
    blk = np.zeros(t, dtype=np.int8); sg = np.zeros(t, dtype=np.int8)
    T = np.zeros((4, t), dtype=np.int8)
    total = 8 ** (t - 1)
    half = (t - 1) // 2
    for code in range(total):
        c = code
        # position 0: block 0 sign + (normalisation: translation not fixed, but block of position 0 := 0, sign +)
        blk[0] = 0; sg[0] = 1
        for i in range(1, t):
            d = c % 8; c //= 8
            blk[i] = d // 2; sg[i] = 1 - 2 * (d % 2)
        for j in range(4):
            for i in range(t): T[j, i] = 0
        for i in range(t): T[blk[i], i] = sg[i]
        ok = True
        for s in range(1, half + 1):
            acc = 0
            for j in range(4):
                for i in range(t):
                    acc += T[j, i] * T[j, (i + s) % t]
            if acc != 0:
                ok = False; break
        if not ok: continue
        ntm += 1
        # cross-correlation antisymmetric parts: c_xy(s) = sum_i X(i+s) Y(i);  asym = c(s) - c(-s)
        good = True
        for s in range(1, t):
            a12 = 0; a34 = 0; a13 = 0; a24 = 0; a14 = 0; a23 = 0
            for i in range(t):
                ip = (i + s) % t; im = (i - s) % t
                a12 += T[0, ip] * T[1, i] - T[0, im] * T[1, i]
                a34 += T[2, ip] * T[3, i] - T[2, im] * T[3, i]
                a13 += T[0, ip] * T[2, i] - T[0, im] * T[2, i]
                a24 += T[1, ip] * T[3, i] - T[1, im] * T[3, i]
                a14 += T[0, ip] * T[3, i] - T[0, im] * T[3, i]
                a23 += T[1, ip] * T[2, i] - T[1, im] * T[2, i]
            if a12 != a34 or a13 != -a24 or a14 != a23:
                good = False; break
        if good:
            if nout < maxout:
                for i in range(t): out[nout, i] = blk[i] * 2 + (1 - sg[i]) // 2
            nout += 1
    return ntm, nout, out

if __name__ == '__main__':
    for t in [int(x) for x in sys.argv[1:]]:
        t0 = time.time()
        ntm, nout, out = search(t, 20)
        print(f"t={t}: T-matrices with T(0) in block 1, sign +: {ntm}; with the three relations: {nout} ({time.time()-t0:.0f}s)", flush=True)
        for r in out[:min(nout, 5)]:
            print('   ', ' '.join(f"{'+-'[x % 2]}{x // 2 + 1}" for x in r))
