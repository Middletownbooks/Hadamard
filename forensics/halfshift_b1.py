"""Half-shift template of the b=1 cores.

For a core sequence x on Z_m, m = 2h (h odd), the half-shift product sigma(y) = x(y) x(y+h) is a function on Z_h.
Reported per block:
  * its cyclic sign-change count c and the null probability P(c' <= c) for a uniform random +-1 function on Z_h,
    P(c' = 2j) = 2 C(h, 2j) / 2^h  (sign changes are invariant under translation, negation and (-1)^y twists of x,
    not under multipliers, so they measure structure in the authors' own cyclic order);
  * in natural GS coordinates (blocks B, C, D un-reversed by y -> -1-y, undoing the back-identity R of the array),
    the arc [M, h-1] on which sigma = +1 and the exceptions, i.e. x = [p | q | -p | q] with |p| = M up to the exceptions;
  * the parity lemma check: sum_k PAF_{c_k}(1) must be even (it is -2), while four proper arcs would make it odd.
Inputs: the authors' original files (167, 179, 443) and, for 419, the repo core extracted from the normalised hoa64 copy."""
import os, sys
from math import comb, log10
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
def core_seqs(H, m):
    """4 border rows/cols first, then 4x4 blocks of size m; block (0,q) circulant f(j-i) or back-circulant f(i+j)"""
    M = H[4:, 4:]; I = np.arange(m); seqs = []
    for q in range(4):
        B = M[:m, q * m:(q + 1) * m]
        if np.array_equal(B, B[0][(I[None, :] - I[:, None]) % m]): seqs.append(('circ', B[0].copy()))
        else:
            assert np.array_equal(B, B[0][(I[None, :] + I[:, None]) % m]), 'unknown block type'
            seqs.append(('back', B[0].copy()))
    return seqs
def natural(seqs):
    m = len(seqs[0][1]); I = np.arange(m)
    return np.array([x if t == 'circ' else x[(-1 - I) % m] for t, x in seqs])
def sign_changes(s):
    return int((s != np.roll(s, -1)).sum())
def log10_null(h, c):
    return log10(2 * sum(comb(h, j) for j in range(0, c + 1, 2))) - h * log10(2)
def best_arc(s):
    """best fit of s by +1 on a cyclic arc [a, a+L) and -1 elsewhere: (#exceptions, a, L)"""
    h = len(s); best = (h + 1, 0, 0)
    for a in range(h):
        for L in range(1, h):
            pat = -np.ones(h, int); pat[(a + np.arange(L)) % h] = 1
            e = int((s != pat).sum())
            if e < best[0]: best = (e, a, L)
    return best
def analyse(label, X, verbose=True):
    m = X.shape[1]; h = m // 2; tot = 0.0
    print(f'{label}: m = {m}, h = {h}')
    lag1 = 0
    for k, x in enumerate(X):
        s = x[:h] * x[h:]
        c = sign_changes(s); lp = log10_null(h, c); tot += lp
        e, a, L = best_arc(s)
        arc = set(((a + np.arange(L)) % h).tolist())
        ex = [(int(y), 'hole in the +1 arc' if y in arc else 'isolated +1') for y in range(h) if (s[y] > 0) != (y in arc)]
        cpart = (x[:h] + x[h:]) // 2
        lag1 += int((cpart * np.roll(cpart, -1)).sum())
        line = f'   block {"ABCD"[k]}: sign changes {c:3d} (null mean {h / 2:.0f}), log10 P = {lp:6.1f}'
        if verbose and c < h // 4:
            line += f';  sigma = +1 on [{a}, {(a + L - 1) % h}] (length {L}), exceptions {ex}'
        print(line)
    print(f'   combined log10 P over the four blocks = {tot:.1f};  sum_k PAF_ck(1) = {lag1}')
    return tot
if __name__ == '__main__':
    from load_originals import load_any
    from hadlib import Core
    for name, v in [('matrix_H668_v167_h_b1', 167), ('matrix_H716_v179_h_b1', 179), ('matrix_H1772_v443_h_b1', 443)]:
        H = load_any(name)
        analyse(f'v={v} (authors\' original, natural coordinates)', natural(core_seqs(H, v - 1)))
    c = Core(os.path.join(HERE, '..', 'data', 'core_v419_b1.txt'))
    analyse('v=419 (repo core from the normalised hoa64 copy)', c.X)
