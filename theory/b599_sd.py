r"""Exact stage 1 of the c = 6 family (theory/b599_fam.py) for a fixed Q-compression quadruple z, by the S/D split [P, M].

S/D on the Z_2 axis of Q (report/alpoge_design_answer.md 1.2): per block, cell c = J t + j of G'' = Z_3 x Z_n (/ +-),
(sigma, delta)(c) = ((y(0, c) + y(1, c))/2, (y(0, c) - y(1, c))/2), with
  sum_k PAF_{G''}(sigma_k) = T_S,   sum_k PAF_{G''}(delta_k) = T_D   (T_S = (T1(0, .) + T1(1, .))/4, T_D likewise),
and the Q-compression fixes the Z_3-coset sums:  sum_j w_j sigma(t, j) = (z(0, t) + z(1, t))/2 =: zS(t),
sum_j w_j delta(t, j) = (z(0, t) - z(1, t))/2 =: zD(t).
Pipeline (every step exact):
  1. per block the distinct sigma vectors with the coset sums zS, S-energy <= its cap, and per coset a delta
     completion with sum zD (D-feasibility; the minimum D-energy is recorded for pruning);
  2. S-join: ordered 4-tuples with sum of S-signatures = T_S (64-bit hash 2 + 2 join, exact re-check);
  3. per S-solution: per block every delta compatible with sigma (cellwise), coset sums zD, D-energy cap; D-join
     (2 + 2 on D-signatures) = T_D; then y = (sigma + delta | sigma - delta) and the full stage-1 check.
usage: python3 theory/b599_sd.py --n 9 --q 11 --probe   (list sizes per z)   |   --run CLASS ROWIDX"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b599_fam import Fam, _sig_of


class SD:
    def __init__(self, F):
        self.F = F; J = F.J; n = F.n
        L = F.LET.tolist()
        # per-cell (sigma, delta) options
        self.pairs = sorted({((a + b) // 2, (a - b) // 2) for a in L for b in L})
        self.svals = sorted({s for s, _ in self.pairs})
        self.dopts = {s: sorted({d for s2, d in self.pairs if s2 == s}) for s in self.svals}
        # G'' = Z_3 x Z_n, cells J t + j; expansion and lag classes (t in {0,1}, j in 0..J-1)
        g = lambda t, j: n * (t % 3) + (j % n)
        self.EXP = np.zeros(3 * n, np.int64)
        for t in range(3):
            for j in range(n): self.EXP[g(t, j)] = J * t + min(j, n - j)
        self.LAG = np.array([g(t, j) for t in range(2) for j in range(J)], np.int64)
        self.SH = np.zeros((2 * J, 3 * n), np.int64)
        for Lc in range(2 * J):
            d = self.LAG[Lc]; dt, dj = d // n, d % n
            for x in range(3 * n):
                t, j = x // n, x % n
                self.SH[Lc, x] = g(t + dt, j + dj)
        T1 = F.T1; h = 3 * n
        self.TS = np.array([(T1[d] + T1[h + d]) // 4 for d in self.LAG], np.int64)
        self.TD = np.array([(T1[d] - T1[h + d]) // 4 for d in self.LAG], np.int64)
        self.Wc = np.array([1] + [2] * (J - 1), np.int64)

    def sig(self, C):
        C = np.atleast_2d(np.asarray(C, np.int64))
        return _sig_of(np.ascontiguousarray(C[:, self.EXP]), self.SH)

    # ------------------------------------------------------------- per-coset options
    def _coset_tables(self):
        """all sigma tuples of one coset (J cells) with their sums and S-energies (cached)"""
        if not hasattr(self, '_ct'):
            J = self.F.J; sv = np.array(self.svals, np.int64)
            idx = np.array(list(itertools.product(range(len(sv)), repeat=J)), np.int64)
            T = sv[idx]
            self._ct = (T, T @ self.Wc, (T * T) @ self.Wc)
            # delta options per sigma value as a padded table
            mx = max(len(v) for v in self.dopts.values())
            self._dtab = np.full((len(sv), mx), 10 ** 6, np.int64); self._dn = np.zeros(len(sv), np.int64)
            for i, s_ in enumerate(self.svals):
                self._dtab[i, :len(self.dopts[s_])] = self.dopts[s_]; self._dn[i] = len(self.dopts[s_])
            self._sidx = {s_: i for i, s_ in enumerate(self.svals)}
        return self._ct

    def coset_sigma(self, zS, zD, Scap):
        """all sigma tuples of one Z_3-coset (J cells) with sum zS, S-energy <= Scap, that admit a delta tuple with
        sum zD; returns (tuples (N x J), S-energies, min D-energies)"""
        T, sm, en = self._coset_tables()
        m = (sm == zS) & (en <= Scap)
        Ts = T[m]; es = en[m]
        sidx = np.searchsorted(np.array(self.svals), Ts)
        md = _min_denergy(sidx, self._dtab, self._dn, self.Wc, int(zD))
        ok = md < 10 ** 8
        return Ts[ok], es[ok], md[ok]

    def block_sigmas(self, zS, zD, Scap, Dcap):
        """distinct sigma vectors of a block (3 cosets) with S-energy <= Scap and min D-energy <= Dcap"""
        per = [self.coset_sigma(int(zS[t]), int(zD[t]), Scap) for t in range(3)]
        A, eA, dA = per[0]; B, eB, dB = per[1]; C, eC, dC = per[2]
        ia, ib = np.nonzero((eA[:, None] + eB[None, :] <= Scap) & (dA[:, None] + dB[None, :] <= Dcap))
        AB = np.concatenate([A[ia], B[ib]], 1); eAB = eA[ia] + eB[ib]; dAB = dA[ia] + dB[ib]
        out = []
        for c0 in range(0, len(AB), 4096):
            sl = slice(c0, c0 + 4096)
            i, j = np.nonzero((eAB[sl, None] + eC[None, :] <= Scap) & (dAB[sl, None] + dC[None, :] <= Dcap))
            out.append(np.concatenate([AB[sl][i], C[j]], 1))
        J = self.F.J
        return np.concatenate(out) if out else np.zeros((0, 3 * J), np.int64)


@njit(cache=True)
def _min_denergy(sidx, dtab, dn, W, zD):
    """per sigma tuple (indices into svals): minimum D-energy of a delta tuple with weighted sum zD (1e9 if none)"""
    N, J = sidx.shape
    OFF = 200; out = np.full(N, 10 ** 9, np.int64)
    cur = np.empty(2 * OFF + 1, np.int64); nxt = np.empty(2 * OFF + 1, np.int64)
    for i in range(N):
        for x in range(2 * OFF + 1): cur[x] = 10 ** 9
        cur[OFF] = 0
        for jj in range(J):
            for x in range(2 * OFF + 1): nxt[x] = 10 ** 9
            si = sidx[i, jj]
            for x in range(2 * OFF + 1):
                if cur[x] >= 10 ** 9: continue
                for t in range(dn[si]):
                    d = dtab[si, t]; y = x + W[jj] * d
                    if 0 <= y <= 2 * OFF:
                        e = cur[x] + W[jj] * d * d
                        if e < nxt[y]: nxt[y] = e
            for x in range(2 * OFF + 1): cur[x] = nxt[x]
        if 0 <= zD + OFF <= 2 * OFF: out[i] = cur[zD + OFF]
    return out


def zSD(z):
    z = np.asarray(z, np.int64)
    return (z[:, :3] + z[:, 3:]) // 2, (z[:, :3] - z[:, 3:]) // 2


def probe(F, z):
    S = SD(F)
    zs, zd = zSD(z)
    TS0, TD0 = int(S.TS[0]), int(S.TD[0])
    # minimum S- and D-energies per block for the caps
    mins = []; mind = []
    for k in range(4):
        cs = [S.coset_sigma(int(zs[k, t]), int(zd[k, t]), TS0) for t in range(3)]
        if any(len(c[0]) == 0 for c in cs): return None
        mins.append(sum(int(c[1].min()) for c in cs)); mind.append(sum(int(c[2].min()) for c in cs))
    out = []
    for k in range(4):
        Scap = TS0 - (sum(mins) - mins[k]); Dcap = TD0 - (sum(mind) - mind[k])
        Lk = S.block_sigmas(zs[k], zd[k], Scap, Dcap)
        out.append(len(Lk))
    return dict(sizes=out, minS=mins, minD=mind)


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    n, q = int(opt('--n', '9')), int(opt('--q', '11'))
    F = Fam(n, q)
    import b599_flevels as FL
    if '--probe' in a:
        quads = FL.class_quads(F)
        lim = int(opt('--lim', '20'))
        for ci, cl, c, z in quads[::max(1, len(quads) // lim)]:
            t0 = time.time(); r = probe(F, z)
            print(json.dumps(dict(cls=ci, c=list(c), rows=[int(x) for x in z.sum(1)], probe=r, secs=round(time.time() - t0, 1))), flush=True)
