r"""Exact algebra of the generalised-Paley (cyclotomic index s) slice of a b = 1 GS core, report/paley_523_answer.md.

G = Z_F x Z_q (q prime, gcd(F, q) = 1), the multiplier group fixes Z_F pointwise and acts on Z_q by the index-s subgroup
H of (Z/q)^*.  Classes C_j = g0^j H (g0 = least primitive root, j = 0..s-1), f = |C_j| = (q - 1)/s.  On every Z_q-coset g
of every block k:  x(g, 0) = alpha,  x(g, z) = c_j for z in C_j.

Gauss periods eta_j = sum_{z in C_j} zeta^z form a Z-basis of the ring of integers O_K of the degree-s subfield K of
Q(zeta_q), and 1 = -sum_j eta_j.  The Z_q-transform of a coset is therefore
    Z(g) = alpha + sum_j c_j eta_j = sum_j (c_j - alpha) eta_j,        eta-coordinates z_j = c_j - alpha in {0, +-2},
and the PAF of x, summed against zeta^z over Z_q, is P(h) = sum_g conj(Z(g)) Z(g + h)  (h in Z_F).  The b = 1 condition
at the characters nontrivial on Z_q is exactly
    sum_k P_k(h) = 4v delta_{h,0}   in O_K,  i.e. eta-coordinates (-4v, ..., -4v) at h = 0 and 0 elsewhere,
because PAF_x(h, 0) + sum_l R(h, l) eta_l = sum_l (R(h, l) - PAF_x(h, 0)) eta_l with R(h, l) = PAF_x(h, z in C_l).
Complex conjugation is eta_j -> eta_{j + c}, c = class of -1 (c = s/2 when (q - 1)/s is odd, e.g. q = 13, 29 at s = 4).

This module: classes, eta numerics, the multiplication table in the eta basis, per-letter coset alphabets, P(h) for
batches of configurations, and the cyclic-CRT expansion to the x-domain for FFT checks."""
import itertools as it
import numpy as np

def primroot(q):
    ps = [r for r in range(2, q) if (q - 1) % r == 0 and all(r % d for d in range(2, r))]
    return next(g for g in range(2, q) if all(pow(g, (q - 1) // r, q) != 1 for r in ps))

class Cyclo:
    """order-s cyclotomy on Z_q"""
    def __init__(self, q, s):
        assert (q - 1) % s == 0
        self.q, self.s, self.f = q, s, (q - 1) // s
        self.g0 = primroot(q)
        cls = np.full(q, -1, np.int64); x = 1
        for e in range(q - 1):
            cls[x] = e % s; x = x * self.g0 % q
        self.cls = cls
        self.C = [np.nonzero(cls == j)[0] for j in range(s)]
        self.cneg = int(cls[q - 1])                                   # class of -1
        zeta = np.exp(2j * np.pi * np.arange(q) / q)
        self.eta = np.array([zeta[self.C[j]].sum() for j in range(s)])  # embedding t = 0
        # eta_i eta_j = sum_l mult[i, j, l] eta_l   (constant term absorbed with 1 = -sum eta)
        mult = np.zeros((s, s, s), np.int64)
        for i in range(s):
            for j in range(s):
                cnt = np.zeros(q, np.int64)
                for a in self.C[i]:
                    for b in self.C[j]: cnt[(a + b) % q] += 1
                for l in range(s): mult[i, j, l] = cnt[self.C[l][0]] - cnt[0]
                assert all(cnt[z] == cnt[self.C[cls[z]][0]] for z in range(1, q))
        self.mult = mult
        # T[i, j, l]: coefficient of eta_l in eta_i * conj(eta_j) = eta_i eta_{j + cneg}
        self.T = np.stack([mult[:, (j + self.cneg) % s, :] for j in range(s)], 1)
        # cyclotomic numbers in the repo convention (scripts/cyclo_stage2_cpsat.py)
        reps = [int(self.C[t].min()) for t in range(s)]
        N = np.zeros((s, s, s), np.int64)
        for t, st in enumerate(reps):
            for z in range(1, q):
                w = (z + st) % q
                if w: N[t, cls[z], cls[w]] += 1
        self.N, self.reps = N, reps

    def num(self, z, t=0):
        """complex value of eta-coordinates z (..., s) in embedding t (eta_j -> eta_{j+t})"""
        e = np.roll(self.eta, -t)
        return np.asarray(z) @ e

    def letters(self):
        return sorted({a + self.f * sum(c) for a in (1, -1) for c in it.product((1, -1), repeat=self.s)})

    def alphabet(self, y):
        """all (alpha, c) with alpha + f * sum(c) = y, as an (n, 1 + s) int array"""
        out = [(a,) + c for a in (1, -1) for c in it.product((1, -1), repeat=self.s) if a + self.f * sum(c) == y]
        return np.array(out, np.int64).reshape(-1, 1 + self.s)

def zcoords(ac):
    """(..., 1 + s) rows (alpha, c_0..c_{s-1}) -> eta-coordinates c_j - alpha"""
    ac = np.asarray(ac)
    return ac[..., 1:] - ac[..., :1]

def pz(Z, T):
    """P(h) = sum_g conj(Z(g)) Z(g + h) for Z (..., F, s) eta-coordinates -> (..., F, s) eta-coordinates.
    CC[i, j](h) = sum_g Z_i(g + h) Z_j(g); P_l(h) = sum_ij T[i, j, l] CC[i, j](h)."""
    Z = np.asarray(Z, np.int64); F = Z.shape[-2]
    out = np.zeros(Z.shape, np.int64)
    for h in range(F):
        Zs = np.roll(Z, -h, axis=-2)                                   # Zs[g] = Z[g + h]
        CC = np.einsum('...gi,...gj->...ij', Zs, Z)
        out[..., h, :] = np.einsum('...ij,ijl->...l', CC, T)
    return out

def expand_x(AC, cy):
    """AC (4, F, 1 + s) -> x (4, F, q) +-1"""
    AC = np.asarray(AC); F = AC.shape[1]
    X = np.zeros((AC.shape[0], F, cy.q), np.int64)
    X[:, :, 0] = AC[:, :, 0]
    for z in range(1, cy.q): X[:, :, z] = AC[:, :, 1 + cy.cls[z]]
    return X

def paf2(X):
    """sum over blocks of the 2-D PAF on Z_F x Z_q"""
    return sum(np.rint(np.fft.ifft2(np.abs(np.fft.fft2(x)) ** 2).real).astype(np.int64) for x in X)

def to_cyclic(X):
    """(4, F, q) -> (4, Fq) on cyclic Z_{Fq} by CRT (index i <-> (i mod F, i mod q))"""
    k, F, q = X.shape; i = np.arange(F * q)
    return X[:, i % F, i % q]

def target_P(F, v, s):
    T = np.zeros((F, s), np.int64); T[0, :] = -4 * v
    return T
