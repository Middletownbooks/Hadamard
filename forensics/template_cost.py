"""Gaussian cost of the half-shift arc template for b=1 (v = 167, 179), using the published cores' own supports.
With sigma_k fixed, x_k = [u_k | sigma_k u_k] and the b=1 conditions split exactly into
   periodic part   c_k = u_k on S_k = {sigma_k = +1}:   sum_k PAF_{c_k}(t)  = -2   (t = 1..n, n = (h-1)/2)
   negacyclic part d_k = u_k on S_k^c (twisted):        sum_k NPAF_{d_k}(t) =  0   (t = 1..n)
plus row sums (c: (+-1,0,0,0); d: sum of squares of D_k(-1) = 2h+1).  For fixed supports every lag sum has fixed
parity (lattice spacing 2), so each lag costs log2(sqrt(2 pi var)/2) bits.  Variances by Monte Carlo over random signs.
Compared with the unrestricted d=1 core (all 4*2h signs free, h lag classes over Z_2h, lattice 4)."""
import os, sys
from math import log2, pi, sqrt
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from halfshift_b1 import core_seqs
from load_originals import load_any
rng = np.random.default_rng(11)
def lag_bits(S, lattice):
    var = S.var(0); return float(np.sum(np.log2(np.sqrt(2 * pi * var) / lattice)))
def paf_mod(xs, h, nega):
    """sum over blocks of periodic (nega=False) or negacyclic (nega=True) autocorrelations, lags 1..n; xs: (N, 4, h)"""
    n = (h - 1) // 2; out = np.zeros((xs.shape[0], n))
    for t in range(1, n + 1):
        sh = np.roll(xs, -t, axis=2)
        if nega: sh[:, :, h - t:] *= -1
        out[:, t - 1] = (xs * sh).sum((1, 2))
    return out
if __name__ == '__main__':
    for name, v in [('matrix_H668_v167_h_b1', 167), ('matrix_H716_v179_h_b1', 179)]:
        H = load_any(name); m = v - 1; h = m // 2; n = (h - 1) // 2
        X = np.array([x for _, x in core_seqs(H, m)])
        sig = X[:, :h] * X[:, h:]
        Sc = (sig > 0).astype(int); Sd = 1 - Sc
        N = 20000
        signs = rng.choice([-1, 1], size=(N, 4, h))
        pc = paf_mod(signs * Sc[None], h, False); pdn = paf_mod(signs * Sd[None], h, True)
        bc, bd = lag_bits(pc, 2), lag_bits(pdn, 2)
        nc, nd = int(Sc.sum()), int(Sd.sum())
        # row-sum constraints: c: (sum c_A, sum c_B, ...) = (+-1, 0, 0, 0) ; d: sum_k D_k(-1)^2 = 2h+1 (use the same Gaussian bookkeeping)
        rc = sum(log2(sqrt(2 * pi * int(Sc[k].sum())) / 2) for k in range(4)) - 1          # 4 sums fixed up to the A sign
        # unrestricted d=1 core: 4 blocks of 2h signs, h lag classes over Z_2h, lattice 4 for b=1 (PAF(t) = m mod 4)
        xs = rng.choice([-1, 1], size=(4000, 4, m))
        full = paf_mod(xs, m, False)                      # lags 1..h-1 of Z_2h (t and 2h-t are the same equation)
        lagh = (xs * np.roll(xs, -h, axis=2)).sum((1, 2))  # lag h (self-paired)
        bfull = float(np.sum(np.log2(np.sqrt(2 * pi * full.var(0)) / 4))) + log2(sqrt(2 * pi * lagh.var()) / 4)
        print(f'v={v} (h={h}, n={n}):')
        print(f'   template, periodic part:   {nc} sign bits, {n} lag equations = {bc:.0f} bits (+ ~{rc:.0f} row-sum bits)  -> log2 E = {nc - bc - rc:+.0f}')
        print(f'   template, negacyclic part: {nd} sign bits, {n} lag equations = {bd:.0f} bits                       -> log2 E = {nd - bd:+.0f}')
        print(f'   unrestricted core:         {4*m} sign bits, {h} lag classes  = {bfull:.0f} bits                      -> log2 E = {4*m - bfull:+.0f}')
    # scaling: the same template with plain arcs of lengths (n, n+1, n+1, n) at other orders (exceptions ignored)
    print('scaling with plain arcs (expected solutions per half = 2^(sign bits - constraint bits)):')
    for v in [167, 179, 251, 307, 419, 443, 523, 739]:
        m = v - 1; h = m // 2; n = (h - 1) // 2
        Sc = np.zeros((4, h), int)
        for k, L in enumerate([n, n + 1, n + 1, n]): Sc[k, h - L:] = 1
        s = rng.choice([-1, 1], size=(3000 if h > 200 else 10000, 4, h))
        bc = lag_bits(paf_mod(s * Sc[None], h, False), 2); bd = lag_bits(paf_mod(s * (1 - Sc)[None], h, True), 2)
        print(f'   v={v}: h={h}: periodic half {int(Sc.sum())} - {bc:.0f} = {Sc.sum() - bc:+.0f};'
              f'  negacyclic half {int((1 - Sc).sum())} - {bd:.0f} = {(1 - Sc).sum() - bd:+.0f}')
