r"""Same-shape family of the 599 lead slice: c = 6, b = 5, G' = Z_2 x Z_3 x Z_n (n odd), -1 on Z_n, Paley(q) on Z_q
(q = 3 mod 4 prime, index 2, f = (q - 1)/2) [P, M].

  h = n q, m = 6 h, v = 6 (h + 1) - 1 = 6 n q + 5;  599 = (n, q) = (9, 11);  ladder (3, 11) -> 203, (5, 11) -> 335,
  (7, 11) -> 467, (9, 7) -> 383.
Block vectors: one value per G'-orbit, index o = J q_ + j with q_ = 3 e + t in Q = Z_2 x Z_3 and j in 0..J-1
(J = (n + 1)/2, j the representative of {j, -j} in Z_n); weight 1 (j = 0) or 2.  Letters alpha + f (p + n'), i.e.
{+-1, +-(2f - 1), +-(2f + 1)}.
Stage-1 target T1 on G' (report/alpoge_design_answer.md 1.2):  4m - 4b(q - 1) at 0, -4bq on H' \ 0, +4q off H'.
Lag classes of G' under negation and j -> -j: (e, t in {0, 1}, j in 0..J-1): 4 J classes, index 2 J e + J t + j.
"""
import os, sys, math, itertools, json
import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


@njit(cache=True)
def _sig_of(V, SH):
    n = V.shape[0]; nl = SH.shape[0]; ng = SH.shape[1]
    out = np.zeros((n, nl), np.int64)
    for i in range(n):
        for L in range(nl):
            acc = 0
            for g in range(ng): acc += V[i, g] * V[i, SH[L, g]]
            out[i, L] = acc
    return out


class Fam:
    def __init__(self, n=9, q=11):
        assert n % 2 == 1 and q % 4 == 3
        self.n, self.q = n, q
        self.f = (q - 1) // 2; f = self.f
        self.J = (n + 1) // 2; J = self.J
        self.b = 5; self.c = 6
        self.h = n * q; self.m = 6 * self.h; self.v = self.m + 5
        self.LET = np.array(sorted({a + f * (p + r) for a in (1, -1) for p in (1, -1) for r in (1, -1)}), np.int64)
        self.no = 6 * J; self.NG = 6 * n
        self.W = np.array(([1] + [2] * (J - 1)) * 6, np.int64)
        self.WJ = np.array([1] + [2] * (J - 1), np.int64)
        g = lambda e, t, j: 3 * n * e + n * (t % 3) + (j % n)
        self._g = g
        self.EXP = np.zeros(self.NG, np.int64)
        for e in range(2):
            for t in range(3):
                for j in range(n): self.EXP[g(e, t, j)] = J * (3 * e + t) + min(j, n - j)
        nl = 4 * J; self.nl = nl
        self.LAGREP = np.array([g(e, t, j) for e in range(2) for t in range(2) for j in range(J)], np.int64)
        self.SHIFT = np.zeros((nl, self.NG), np.int64)
        for L in range(nl):
            d = self.LAGREP[L]; de, dt, dj = d // (3 * n), (d // n) % 3, d % n
            for x in range(self.NG):
                e, t, j = x // (3 * n), (x // n) % 3, x % n
                self.SHIFT[L, x] = g((e + de) % 2, t + dt, j + dj)
        T = np.full(self.NG, 4 * q, np.int64)
        for j in range(1, n): T[g(0, 0, j)] = -4 * self.b * q
        T[0] = 4 * self.m - 4 * self.b * (q - 1)
        self.T1 = T; self.T1SIG = T[self.LAGREP].copy()
        # S/D on G'' = Z_3 x Z_n: cell = J t + j
        self.nc = 3 * J
        self.W15 = np.array(([1] + [2] * (J - 1)) * 3, np.int64)

    # --------------------------------------------------------- basic maps
    def expand(self, v):
        return np.asarray(v)[..., self.EXP]

    def sig(self, v):
        v = np.atleast_2d(np.asarray(v, np.int64))
        return _sig_of(np.ascontiguousarray(self.expand(v)), self.SHIFT)

    def qcomp(self, v):
        v = np.atleast_2d(np.asarray(v, np.int64))
        return (v * self.W).reshape(v.shape[0], 6, self.J).sum(2)

    def energy(self, v):
        v = np.atleast_2d(np.asarray(v, np.int64))
        return (v * v * self.W).sum(1)

    def paf_full(self, y):
        """exact PAF over G' (direct sums) of y (..., NG)"""
        y = np.asarray(y, np.int64); n = self.n; NG = self.NG
        out = np.zeros(y.shape, np.int64)
        for d in range(NG):
            de, dt, dj = d // (3 * n), (d // n) % 3, d % n
            idx = np.array([self._g((x // (3 * n) + de) % 2, (x // n) % 3 + dt, x % n + dj) for x in range(NG)])
            out[..., d] = (y * y[..., idx]).sum(-1)
        return out

    def check_stage1(self, Y):
        Y = np.asarray(Y, np.int64)
        if Y.shape != (4, self.no) or not np.isin(Y, self.LET).all(): return False
        return bool((self.paf_full(self.expand(Y)).sum(0) == self.T1).all())

    def sd_of(self, v):
        v = np.atleast_2d(np.asarray(v, np.int64)); h = 3 * self.J
        a, b = v[:, :h], v[:, h:]
        return (a + b) // 2, (a - b) // 2

    # --------------------------------------------------------- generic Slice config (theory/ab_sweep.py)
    def cfg(self):
        return dict(v=self.v, c=6, b=5, h=self.h, m=self.m, Q=[2, 3], H1=[self.n], q=self.q, s=2, f=self.f,
                    mods=[2, 3, self.n, self.q], nQ=2, U=[[2, 2, self.n - 1]])

    def selftest(self):
        import ab_sweep as Wp, ab_s1sd as SD
        sl = Wp.build(SD.cfg_from_json(json.dumps(self.cfg())))
        assert (sl.T1 == self.T1).all(), "T1"
        for o in range(self.no):
            els = np.nonzero(self.EXP == o)[0]
            assert len(set(sl.borb[els].tolist())) == 1 and len(els) == self.W[o]
        assert sorted({x for L in sl.letters for x in L}) == self.LET.tolist()
        rng = np.random.default_rng(1)
        V = self.LET[rng.integers(len(self.LET), size=(30, self.no))]
        P = sl.Gp.paf(self.expand(V))
        assert (P[:, self.LAGREP] == self.sig(V)).all() and (self.paf_full(self.expand(V)) == P).all()
        return True


if __name__ == "__main__":
    for (n, q) in ((3, 11), (5, 11), (7, 11), (9, 11), (9, 7), (3, 7), (5, 7)):
        F = Fam(n, q); F.selftest()
        print(json.dumps(dict(n=n, q=q, v=F.v, orbits=F.no, letters=F.LET.tolist(), T1=F.T1SIG[[0, 1, F.J]].tolist())))
