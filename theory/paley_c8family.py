r"""Inflated-border matrices for the c = 8, p = 11 family v = 8 (11 * 2^(j-1) + 1) - 1  (j = 2: 183, j = 3: 359, j = 4: 711)
over G = Z2 x Z2 x Z_{2^j} x Z11, H = <(0,0,2)> x Z11 (index c = 8, b = 7), reusing the PUBLISHED 359 border.

The sign-corrected published H(1436) (theory/paley_splice359.py) gives:  border rows K (28 x 32, constant on the 8 cosets of
each of the 4 core blocks), border columns Kc (28 x 32), and the core layout (sub-block (I,J) = rho(g) gam(h) s F_k(sigma)).
Q = G/H = Z2^3 has the same coordinates (a, b, e mod 2) for every j, so the layout and K, Kc transfer verbatim.
Given new core functions F'_k on G, the corner is FORCED by orthogonality of border rows and core rows:
      C' = -Kbig Core'^T Kcbig / (32 h),   h = |H| = 11 * 2^(j-1)   (Kcbig^T Kcbig = 32 h I),
and the matrix is Hadamard iff C' is +-1 and the core is a core (sum PAF(F'_k) = rho-twisted target) -- checked exactly.
usage: python3 theory/paley_c8family.py selftest        (j = 3 with the published functions reproduces H(1436))"""
import sys, os, json, itertools as it
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory'))
from paley_splice359 import Template359, OFF, M_ as M359

class C8Family:
    def __init__(s):
        T = Template359(); s.T = T; Mt = T.Mt
        co = T.g.coords; cos = co[:, 0] * 4 + co[:, 1] * 2 + (co[:, 2] % 2)
        # K[r, J*8+q], Kc[r, J*8+q] from the published (sign-corrected) matrix; check constancy on cosets
        s.K = np.zeros((OFF, 32), np.int64); s.Kc = np.zeros((OFF, 32), np.int64)
        for J in range(4):
            for q in range(8):
                blk = Mt[:OFF, OFF + J * M359:OFF + (J + 1) * M359][:, cos == q]
                assert (blk == blk[:, :1]).all(); s.K[:, J * 8 + q] = blk[:, 0]
                blkc = Mt[OFF + J * M359:OFF + (J + 1) * M359, :OFF][cos == q, :]
                assert (blkc == blkc[:1, :]).all(); s.Kc[:, J * 8 + q] = blkc[0, :]
        s.C359 = Mt[:OFF, :OFF].copy(); s.layout = T.layout

    def group(s, j):
        mods = (2, 2, 2 ** j, 11); n = int(np.prod(mods))
        co = np.array(np.unravel_index(np.arange(n), mods)).T
        idx = lambda c: np.ravel_multi_index(tuple(np.moveaxis(np.asarray(c) % np.array(mods), -1, 0)), mods)
        sig = {'h-g': idx(co[None, :, :] - co[:, None, :]), 'g+h': idx(co[:, None, :] + co[None, :, :]),
               'g-h': idx(co[:, None, :] - co[None, :, :]), '-g-h': idx(-co[:, None, :] - co[None, :, :])}
        chars = {''.join(map(str, e)): np.where(((co[:, 0] * e[0] + co[:, 1] * e[1] + co[:, 2] * e[2]) % 2) == 0, 1, -1)
                 for e in it.product((0, 1), repeat=3)}
        cos = co[:, 0] * 4 + co[:, 1] * 2 + (co[:, 2] % 2)
        return mods, n, co, sig, chars, cos

    def build(s, F, j):
        """F: (4, n) core functions on Z2 x Z2 x Z_{2^j} x Z11 in the published layout's coordinates"""
        mods, n, co, sig, chars, cos = s.group(j); h = n // 8; N = 4 * n + OFF
        Core = np.zeros((4 * n, 4 * n), np.int64)
        for key, L in s.layout.items():
            I_, J_ = int(key[0]), int(key[1])
            rho = L['rho_sign'] * chars[L['rho']]; gam = L['gam_sign'] * chars[L['gam']]
            Core[I_ * n:(I_ + 1) * n, J_ * n:(J_ + 1) * n] = np.outer(rho, gam) * np.asarray(F[L['k']])[sig[L['sigma']]]
        Kbig = np.concatenate([s.K[:, J * 8 + cos] for J in range(4)], axis=1)              # (28, 4n)
        Kcbig = np.concatenate([s.Kc[:, J * 8 + cos] for J in range(4)], axis=1).T          # (4n, 28)
        num = -(Kbig @ Core.T @ Kcbig)
        ok_div = bool((num % (32 * h) == 0).all())
        C = num // (32 * h) if ok_div else None
        if C is None or not (np.abs(C) == 1).all():
            return None, dict(corner_pm1=False, corner_divisible=ok_div)
        H = np.zeros((N, N), np.int64)
        H[:OFF, :OFF] = C; H[:OFF, OFF:] = Kbig; H[OFF:, :OFF] = Kcbig; H[OFF:, OFF:] = Core
        return H, dict(corner_pm1=True, hadamard=bool(np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64))), N=N)

def selftest():
    S = C8Family()
    H, info = S.build(S.T.F, 3)
    info['equals_published_corrected'] = bool(H is not None and np.array_equal(H, S.T.Mt))
    print(json.dumps(info))

if __name__ == "__main__":
    {"selftest": selftest}[sys.argv[1]]()
