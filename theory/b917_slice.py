r"""The index-s slice of b = 0 at v = 7q with U = <2> on Z_7 (s = 2s', f = (q-1)/s odd), exact algebra [P, M].

G = Z_7 x Z_q.  Rows: the U-orbits R0 = {0}, R1 = {1,2,4}, R2 = {3,5,6} of Z_7.  Columns: z = 0 and the s cosets
C_j = gamma^j M of the order-f subgroup M of F_q^* (gamma primitive).  f odd => -1 not in M and -C_j = C_{j+s/2}.
A block is a 3 x (s+1) sign array X[i, c]: X[i, 0] = x(R_i, 0), X[i, 1+j] = x(R_i, C_j).

Column (x0, x1, x2) -> w = x0 + 3 x1 + 3 x2 (trivial Z_7 character, weights = orbit sizes) and
u = x0 + eta7 x1 + conj(eta7) x2 (chi_1), eta7 = (-1 + i sqrt 7)/2.  [P] u = 2 u' with u' in {0, +-1, +-eta7, +-conj eta7}
and w^2 + 24 |u'|^2 = 49:
  sym  (x1 = x2 = e):  w = x0 + 6e in {+-5, +-7},  u' = (x0 - e)/2 in {0, +-1};
  anti (x1 = -x2 = e): w = x0 = +-1,               u' = (x0 + i sqrt7 e)/2.
With rho = Re u = x0 - p and m = Im u / sqrt 7 (p = e on sym columns, m = e on anti columns): w = rho + 7 p.

Characters: chi in {1, chi_1-orbit, chi_3-orbit}, psi = 1 or psi_t with t in C_a (a = 0..s-1).  For an M-invariant f on
Z_q (f0 at 0, f_j on C_j): fhat(psi_t) = f0 + sum_j f_j eta_{a+j} (Gauss periods eta_i = sum_{z in C_i} zeta_q^z).
  stage 1 (psi = 1):       row sums and sum_k |sum_z u_k|^2 = 4v;
  dual  (chi = 1):         D_a = A + sum_j w_j eta_{a+j},  sum_k |D_{k,a}|^2 = 4v  (a < s/2; D_{a+s/2} = conj D_a);
  cross (chi = chi_1):     X_a = uhat0 + sum_j u_j eta_{a+j}, sum_k |X_{k,a}|^2 = 4v  (a = 0..s-1, s real equations).
Stage 1 + dual + cross <=> sum_k PAF_G(x_k) = 4v delta (checked here by full_paf).
usage: python3 theory/b917_slice.py [--q 131 --s 26]   (self-test: identities on random blocks)"""
import sys, os, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import primitive_root

ETA7 = complex(-0.5, math.sqrt(7) / 2)
R_OF = np.array([0, 1, 1, 2, 1, 2, 2])            # row (U-orbit) of g in Z_7
WROW = np.array([1, 3, 3])                         # orbit sizes
# the 8 column states (x0, x1, x2), index = 4*(x0<0) + 2*(x1<0) + (x2<0)
STATES = np.array([[1 - 2 * ((t >> 2) & 1), 1 - 2 * ((t >> 1) & 1), 1 - 2 * (t & 1)] for t in range(8)], np.int64)
W_OF = STATES @ WROW                               # w of each state
U_OF = STATES[:, 0] + ETA7 * STATES[:, 1] + ETA7.conjugate() * STATES[:, 2]
SYM = STATES[:, 1] == STATES[:, 2]


class Slice917:
    def __init__(self, q=131, s=26):
        assert (q - 1) % s == 0 and s % 2 == 0
        self.q, self.s = q, s
        self.f = (q - 1) // s
        assert self.f % 2 == 1, "f must be odd (-1 not in M)"
        self.v = 7 * q; self.N4 = 4 * self.v
        g = primitive_root(q); self.gamma = g
        Mgen = pow(g, s, q)
        M = [1]
        while True:
            x = M[-1] * Mgen % q
            if x == 1: break
            M.append(x)
        assert len(M) == self.f
        self.M = M
        cls = -np.ones(q, np.int64); cls[0] = -1
        for j in range(s):
            r = pow(g, j, q)
            for m_ in M: cls[r * m_ % q] = j
        self.cls = cls                              # class of z (0 -> -1)
        h = s // 2
        assert all(cls[(q - z) % q] == (cls[z] + h) % s for z in range(1, q)), "-C_j != C_{j+s/2}"
        zeta = np.exp(2j * np.pi * np.arange(q) / q)
        self.eta = np.array([zeta[cls == i].sum() for i in range(s)])       # Gauss periods eta_i
        # period matrix P[a, j] = eta_{a+j}
        self.P = np.array([[self.eta[(a + j) % s] for j in range(s)] for a in range(s)])
        # column index c in 0..s of each z: 0 for z = 0, 1 + class otherwise
        self.col_of_z = np.where(cls < 0, 0, cls + 1)

    # ---------- blocks ----------
    def expand(self, X):
        """block(s) X (..., 3, s+1) -> x on G = Z_7 x Z_q, shape (..., 7, q)"""
        X = np.asarray(X)
        return X[..., R_OF, :][..., :, self.col_of_z]

    def full_paf(self, Xs):
        """sum_k PAF_G(x_k) over Z_7 x Z_q (exact integers), Xs (4, 3, s+1)"""
        x = self.expand(np.asarray(Xs)).astype(np.float64)
        F = np.fft.fft2(x)
        return np.rint(np.fft.ifft2((np.abs(F) ** 2).sum(0)).real).astype(np.int64)

    def is_core(self, Xs):
        P = self.full_paf(Xs)
        T = np.zeros_like(P); T[0, 0] = self.N4
        return bool((P == T).all())

    # ---------- levels ----------
    @staticmethod
    def wu(X):
        """per block: w (..., s+1) real, u (..., s+1) complex (column 0 included)"""
        X = np.asarray(X)
        w = X[..., 0, :] + 3 * X[..., 1, :] + 3 * X[..., 2, :]
        u = X[..., 0, :] + ETA7 * X[..., 1, :] + ETA7.conjugate() * X[..., 2, :]
        return w, u

    def dual_vals(self, w):
        """D_a = w0 + sum_j w_j eta_{a+j}, a = 0..s-1 (complex), w (..., s+1)"""
        w = np.asarray(w)
        return w[..., 0:1] + w[..., 1:] @ self.P.T

    def cross_vals(self, u):
        u = np.asarray(u)
        return u[..., 0:1] + u[..., 1:] @ self.P.T

    def ystage1(self, X):
        """the stage-1 compression y(R_i) = x(R_i, 0) + f sum_j x(R_i, C_j), per block (..., 3)"""
        X = np.asarray(X)
        return X[..., :, 0] + self.f * X[..., :, 1:].sum(-1)

    def residuals(self, Xs):
        """(row, stage-1 chi_1, dual a = 0..s/2-1, cross a = 0..s-1) residuals of sum_k |.|^2 - 4v"""
        w, u = self.wu(Xs)
        r = (w.sum(-1) * 0 + (w[:, 0] + self.f * w[:, 1:].sum(-1)))           # row sums r_k = sum_z w(z)
        st1 = u[:, 0] + self.f * u[:, 1:].sum(-1)
        D = self.dual_vals(w); Xc = self.cross_vals(u)
        h = self.s // 2
        return dict(row=float((r ** 2).sum() - self.N4), st1=float((np.abs(st1) ** 2).sum() - self.N4),
                    dual=(np.abs(D[:, :h]) ** 2).sum(0) - self.N4, cross=(np.abs(Xc) ** 2).sum(0) - self.N4)


def selftest(q=131, s=26, n=5, seed=1):
    S = Slice917(q, s)
    rng = np.random.default_rng(seed)
    h = s // 2
    out = []
    for t in range(n):
        Xs = rng.choice(np.array([-1, 1]), size=(4, 3, s + 1))
        P = S.full_paf(Xs)
        # character sums from the levels must equal the Fourier transform of the full PAF
        w, u = S.wu(Xs)
        D = S.dual_vals(w); Xc = S.cross_vals(u)
        zq = np.exp(-2j * np.pi * np.arange(q) / q); z7 = np.exp(-2j * np.pi * np.arange(7) / 7)
        FP = np.fft.fft2(P.astype(float))                       # = sum_k |x_k^(chi, psi)|^2
        errD = errX = 0.0
        for a in range(s):
            t_ = pow(S.gamma, a, q)                              # t in C_a
            errD = max(errD, abs(FP[0, t_] - (np.abs(D[:, a]) ** 2).sum()))
            errX = max(errX, abs(FP[1, t_] - (np.abs(Xc[:, a]) ** 2).sum()))
            # chi_1 has the same value on R1 (1, 2, 4): FP[2, t] = FP[4, t] = FP[1, t]
            errX = max(errX, abs(FP[2, t_] - FP[1, t_]), abs(FP[4, t_] - FP[1, t_]))
        # conjugate pairing of the dual: D_{a+h} = conj D_a
        errC = float(np.abs(D[:, h:] - D[:, :h].conj()).max())
        # half-degree identity: |D_a|^2 = (A + sum_j s_j Re eta_{a+j})^2 + (sum_j d_j Im eta_{a+j})^2, j < h
        sv = w[:, 1:1 + h] + w[:, 1 + h:]; dv = w[:, 1:1 + h] - w[:, 1 + h:]
        Re = np.array([[S.eta[(a + j) % s].real for j in range(h)] for a in range(h)])
        Im = np.array([[S.eta[(a + j) % s].imag for j in range(h)] for a in range(h)])
        lhs = np.abs(D[:, :h]) ** 2
        rhs = (w[:, :1] + sv @ Re.T) ** 2 + (dv @ Im.T) ** 2
        errH = float(np.abs(lhs - rhs).max())
        # u = 2u' with u' in {0, +-1, +-eta7, +-conj eta7}; w^2 + 24|u'|^2 = 49
        up = u / 2
        okU = bool(np.all(np.min(np.abs(up[..., None] - np.array([0, 1, -1, ETA7, -ETA7, ETA7.conjugate(), -ETA7.conjugate()])), -1) < 1e-12))
        okE = bool(np.allclose(w ** 2 + 24 * np.abs(up) ** 2, 49))
        out.append(dict(errD=float(errD), errX=float(errX), errConj=errC, errHalf=errH, u_alphabet=okU, energy49=okE))
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    q = int(a[a.index('--q') + 1]) if '--q' in a else 131
    s = int(a[a.index('--s') + 1]) if '--s' in a else 26
    for r in selftest(q, s): print(r)
