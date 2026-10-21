r"""Exact gates and stratified counts for bordered-GS structures (theory/bt_engine.BStruct) [M gates, E counts].

Per structure:
1. coset-Parseval screen [P] (bt_engine.BStruct.coset_parseval);
2. Q-level gate [M]: every row pattern (sorted, nonnegative, sum of squares 4n), every Q-compression quadruple meeting
   sum_i PAF_{z_i} = F_K (exact join; blocks with a common S share their lists);
3. first compression levels [M]: for invariant subgroups N <= K with small G/N, the refined lists inside each Q-level
   solution and their exact join (count of solutions and their multiplicity = number of orbit-sign quadruples);
4. stratified count [E]: sum over the solutions of the finest exact level of (multiplicity) x (lattice Gaussian of the
   remaining class equations, read at the true target, from blocks sampled uniformly inside their fibres;
   theory/q523_gauss.gauss_bits).
"""
import math, time, json, sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bt_engine as E, bt_chain as C


def row_patterns(st):
    """sorted nonnegative row patterns (r_1 <= .. <= r_4) with sum r_i^2 = 4 n and r_i = |G| mod 2"""
    n4 = 4 * st.nH; par = st.G.n % 2
    lim = math.isqrt(n4); out = []
    for a in range(par, lim + 1, 2):
        for b in range(a, lim + 1, 2):
            for c in range(b, lim + 1, 2):
                d2 = n4 - a * a - b * b - c * c
                if d2 < c * c:
                    continue
                d = math.isqrt(d2)
                if d * d == d2 and d % 2 == par:
                    out.append((a, b, c, d))
    return out


def pattern_weight(pat):
    """number of ordered signed row vectors represented by a sorted nonnegative pattern"""
    from collections import Counter
    c = Counter(pat)
    perm = math.factorial(4)
    for v in c.values():
        perm //= math.factorial(v)
    return perm * 2 ** sum(1 for x in pat if x)


class Gates:
    def __init__(self, st, maxlen=4 * 10 ** 6, log=None, maxpairs=3e11):
        self.st = st
        self.maxlen = maxlen
        self.maxpairs = maxpairs
        self.qQ = st.quotient(st.Kel)
        self.chQ = C.char_matrix(self.qQ)
        self.log = log if log is not None else {}
        self._lq = {}

    def q_list(self, i, row):
        key = (i if self.st.split else 0, row)
        if key not in self._lq:
            self._lq[key] = C.q_level_list_two_step(self.st, self.qQ, i, row, self.chQ)
        return self._lq[key]

    def q_level(self, patterns=None, signs=False):
        """Q-level solutions per row pattern: list of (pattern-with-signs, [index quadruples], lists)"""
        st = self.st
        pats = row_patterns(st) if patterns is None else patterns
        out = []
        for pat in pats:
            sign_sets = [pat]
            if signs or st.split:
                sign_sets = sorted({tuple(s * p for s, p in zip(sg, pat))
                                    for sg in itertools.product((1, -1), repeat=4)})
            for rows in sign_sets:
                L = [self.q_list(i, rows[i]) for i in range(4)]
                if min(len(l[0]) for l in L) == 0:
                    continue
                P = [C.paf_reduced(self.qQ, i, L[i][0]) for i in range(4)]
                sol, no = C.join4(P, self.qQ.cls_target)
                if no:
                    out.append((rows, sol, L))
        return out

    def refine_level(self, qA, qsol, chA=None, cvA=None, count_only=False):
        """for one Q-level solution (rows, index quadruple, lists): the level-A lists and their join"""
        rows, s, L = qsol
        st = self.st
        chA = chA if chA is not None else C.char_matrix(qA)
        cvA = cvA if cvA is not None else [C.coord_values(qA, i) for i in range(4)]
        LA = [C.refine_list(self.qQ, L[i][0][s[i]], qA, i, chars=chA, cvals=cvA[i], maxlen=self.maxlen)
              for i in range(4)]
        if min(len(l[0]) for l in LA) == 0:
            return LA, np.zeros((0, 4), np.int64), 0, 0.0
        PA = [C.paf_reduced(qA, i, LA[i][0]) for i in range(4)]
        sz = sorted(len(l[0]) for l in LA)
        if float(sz[0]) * sz[3] + float(sz[1]) * sz[2] > self.maxpairs:
            raise MemoryError(f"level join pairs {float(sz[0]) * sz[3] + float(sz[1]) * sz[2]:.3g} > {self.maxpairs:.3g}")
        solA, noA = C.join4b(PA, qA.cls_target)
        mult = float(sum(np.prod([LA[i][1][a[i]] for i in range(4)]) for a in solA))
        return LA, solA, noA, mult


# ------------------------------------------------------------------------------------------------ sampling
class CellSampler:
    """uniform sampling of eps in {+-1}^k with sum_j w_j eps_j = v (subset-sum DP, backward sampling)"""

    def __init__(self, w, v):
        self.w = np.asarray(w, np.int64); self.v = int(v)
        k = len(self.w)
        tot = int(self.w.sum())
        self.off = tot
        # cnt[j][s + off] = number of ways for orbits j..k-1 to reach s (float64 to stay safe for big k)
        cnt = np.zeros((k + 1, 2 * tot + 1))
        cnt[k, tot] = 1.0
        for j in range(k - 1, -1, -1):
            wj = int(self.w[j])
            cnt[j, wj:] += cnt[j + 1, :2 * tot + 1 - wj]           # eps_j = +1: remainder s - w_j
            cnt[j, :2 * tot + 1 - wj] += cnt[j + 1, wj:]           # eps_j = -1: remainder s + w_j
        self.cnt = cnt
        self.total = cnt[0, self.v + tot] if abs(self.v) <= tot else 0.0

    def sample(self, N, rng):
        k = len(self.w); tot = self.off
        X = np.empty((N, k), np.int8)
        s = np.full(N, self.v, np.int64)
        for j in range(k):
            wj = int(self.w[j])
            c_plus = self.cnt[j + 1, s - wj + tot] if True else None
            c_all = self.cnt[j, s + tot]
            p = np.where(c_all > 0, c_plus / np.maximum(c_all, 1e-300), 0.0)
            e = np.where(rng.random(N) < p, 1, -1)
            X[:, j] = e
            s = s - wj * e
        assert (s == 0).all()
        return X


class FibreSampler:
    """uniform sampling of orbit sign vectors of block i inside the fibre of a reduced vector at level qb"""

    def __init__(self, st, i, qb, yb):
        self.st, self.i = st, i
        self.members = [np.flatnonzero(qb.omap[i] == U) for U in range(len(qb.bsize[i]))]
        self.samplers = [CellSampler(qb.w[i][m], yb[U]) for U, m in enumerate(self.members)]
        self.ok = all(cs.total > 0 for cs in self.samplers)
        self.log2_size = sum(math.log2(cs.total) for cs in self.samplers) if self.ok else -math.inf

    def sample(self, N, rng):
        r = self.st.r[self.i]
        X = np.empty((N, r), np.int8)
        for m, cs in zip(self.members, self.samplers):
            X[:, m] = cs.sample(N, rng)
        return X


def stratified_bits(st, qb, ys, N=4000, seed=1, gauss=None):
    """-log2 P(remaining class equations | compressions ys at qb) under the lattice Gaussian of the summed class
    PAF vectors (blocks sampled independently and uniformly in their fibres).  Returns dict from gauss_bits."""
    import q523_gauss as QG
    rng = np.random.default_rng(seed)
    V = np.zeros((N, st.ncls), np.int64)
    for i in range(4):
        fs = FibreSampler(st, i, qb, ys[i])
        if not fs.ok:
            return dict(bits=math.inf, note='empty fibre')
        X = fs.sample(N, rng)
        V += C.class_paf_orbits(st, i, X)
    return QG.gauss_bits(V, st.cls_target)
