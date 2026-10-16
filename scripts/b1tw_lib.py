"""Twisted b = 1 Goethals-Seidel slices on cyclic G = Z_m (m = v - 1): library.

Ansatz (report/prompt_b1_twist_design.md §2): u in (Z/m)* of order d, g a block permutation fixing A, and
    X_{pi(i)}(u x) = X_i(x).
Twists (signs are irrelevant up to block negation: 0 is a fixed point of u, so a g-fixed label at 0 forces the signs
around every cycle of g to multiply to +1, and a pi-fixed block cannot carry eps = -1):
    'id'  : all four blocks U-invariant (the usual multiplier slice);
    'BC'  : A, D U-invariant; C(u x) = B(x), B(u x) = C(x)  (B is U^2-invariant, C = B o u^-1; needs d even);
    'BCD' : A U-invariant; C(u x) = B(x), D(u x) = C(x), B(u x) = D(x)  (B is U^3-invariant; needs 3 | d).
An orbit of size s carries free bits: 'id' 4; 'BC' 2 + (2 if s even else 1, B = C there); 'BCD' 1 + (3 if 3 | s
else 1, B = C = D there).
b = 1 target: sum_k PAF_k(g) = -4 for g != 0, block sums (+-2, 0, 0, 0) with the +-2 block pi-fixed.
ΣPAF is constant on the classes of <U, -1> on Z_m minus {0}: one equation per class.  The class {m/2} has size 1.
"""
import math
import numpy as np


def units(m):
    return [a for a in range(1, m) if math.gcd(a, m) == 1]


def mult_order(u, m):
    k, x = 1, u % m
    while x != 1:
        x = x * u % m; k += 1
    return k


def cyclic_subgroups(m):
    """[(d, u, frozenset U)] one entry per cyclic subgroup of (Z/m)*, u = least generator"""
    seen = {}
    for a in units(m):
        S, x = [1], a % m
        while x != 1:
            S.append(x); x = x * a % m
        key = frozenset(S)
        if key not in seen: seen[key] = a
    return sorted([(len(k), u, k) for k, u in seen.items()], key=lambda t: (t[0], t[1]))


def orbits(u, m):
    seen, out = set(), []
    for x in range(m):
        if x in seen: continue
        o = [x]; y = x * u % m
        while y != x:
            o.append(y); y = y * u % m
        seen.update(o); out.append(o)
    return out


def classes(u, m):
    """reps and sizes of the <u, -1>-orbits on Z_m minus {0}"""
    seen, reps, sizes = set(), [], []
    for x in range(1, m):
        if x in seen: continue
        c, y = set(), x
        while True:
            c.add(y); c.add((-y) % m); y = y * u % m
            if y == x: break
        seen |= c; reps.append(x); sizes.append(len(c))
    return np.array(reps), np.array(sizes)


class Slice:
    """bit layout of a twisted slice: idx[k][y] = index of the free bit that gives block k at position y"""
    def __init__(self, m, u, tw):
        assert tw in ("id", "BC", "BCD")
        self.m, self.u, self.tw = m, u % m, tw
        self.d = mult_order(u, m)
        if tw == "BC": assert self.d % 2 == 0
        if tw == "BCD": assert self.d % 3 == 0
        self.orbits = orbits(u, m)
        self.reps, self.csize = classes(u, m)
        idx = -np.ones((4, m), np.int64)
        groups = {k: [] for k in ("A", "B", "C", "D")}       # bits of each row-sum group (group B carries C, D if twisted)
        nb = 0
        def new(g):
            nonlocal nb
            groups[g].append(nb); nb += 1; return nb - 1
        for o in self.orbits:
            s = len(o)
            a = new("A")
            for y in o: idx[0, y] = a
            if tw == "id":
                for k, g in ((1, "B"), (2, "C"), (3, "D")):
                    b = new(g)
                    for y in o: idx[k, y] = b
            elif tw == "BC":
                dd = new("D")
                for y in o: idx[3, y] = dd
                if s % 2 == 0:
                    b, c = new("B"), new("B")
                    for j, y in enumerate(o):
                        idx[1, y], idx[2, y] = (b, c) if j % 2 == 0 else (c, b)
                else:
                    b = new("B")
                    for y in o: idx[1, y] = idx[2, y] = b
            else:
                if s % 3 == 0:
                    b, c, dd = new("B"), new("B"), new("B")
                    for j, y in enumerate(o):
                        idx[1, y] = (b, dd, c)[j % 3]; idx[2, y] = (c, b, dd)[j % 3]; idx[3, y] = (dd, c, b)[j % 3]
                else:
                    b = new("B")
                    for y in o: idx[1, y] = idx[2, y] = idx[3, y] = b
        assert (idx >= 0).all()
        self.idx, self.nbits = idx, nb
        self.groups = {g: np.array(v, np.int64) for g, v in groups.items() if v}
        # row-sum weights: group -> (block whose sum it controls, weights per bit)
        blk = {"A": 0, "B": 1, "C": 2, "D": 3}
        self.gw = {}
        for g, bits in self.groups.items():
            k = blk[g]; w = np.bincount(idx[k], minlength=nb)[bits]
            self.gw[g] = (k, w.astype(np.int64))
        # +-2 block is A; symmetry factor for the other admissible placements of the +-2 block
        self.sym = {"id": 4, "BC": 2, "BCD": 1}[tw]
        self.r = len(self.orbits)
        self.L = len(self.reps)
        self.minus1_in_U = (m - 1) in {pow(self.u, j, m) for j in range(self.d)}
        h = m // 2
        self.half_class = int(np.flatnonzero(self.reps == h)[0]) if h in set(self.reps.tolist()) else None

    def expand(self, beta):
        """beta: (n, nbits) +-1 -> (n, 4, m)"""
        return beta[:, self.idx]

    def targets(self):
        """row-sum targets per group for the A = +-2 pattern (list of dicts)"""
        out = []
        for sA in (2, -2):
            t = {g: 0 for g in self.groups}; t["A"] = sA; out.append(t)
        return out


def count_signed(w, T):
    """number of eps in {+-1}^n with sum w_i eps_i = T (exact, Python ints)"""
    D = {0: 1}
    for x in w:
        E = {}
        for s, c in D.items():
            E[s + x] = E.get(s + x, 0) + c; E[s - x] = E.get(s - x, 0) + c
        D = E
    return D.get(T, 0)


class SignedSampler:
    """uniform sampler of eps in {+-1}^n with sum w_i eps_i = T (backward DP, float ratios)"""
    def __init__(self, w, T):
        self.w = np.asarray(w, np.int64); W = int(self.w.sum()); self.W = W; self.T = T
        n = len(self.w); C = np.zeros((n + 1, 2 * W + 1)); C[0, W] = 1.0
        for i in range(n):
            x = int(self.w[i]); prev = C[i]; cur = np.zeros_like(prev)
            cur[x:] += prev[:len(prev) - x] if x else prev
            if x: cur[:len(prev) - x] += prev[x:]
            else: cur += prev
            mx = cur.max(); C[i + 1] = cur / mx if mx > 0 else cur
        self.C = C
        self.ok = C[n, W + T] > 0

    def sample(self, N, rng):
        n = len(self.w); W = self.W; out = np.empty((N, n), np.int8); t = np.full(N, self.T, np.int64)
        for i in range(n - 1, -1, -1):
            x = int(self.w[i]); prev = self.C[i]
            pp = prev[np.clip(t - x + W, 0, 2 * W)] * ((t - x + W >= 0) & (t - x + W <= 2 * W))
            pm = prev[np.clip(t + x + W, 0, 2 * W)] * ((t + x + W >= 0) & (t + x + W <= 2 * W))
            p = pp / (pp + pm)
            e = np.where(rng.random(N) < p, 1, -1)
            out[:, i] = e; t = t - x * e
        assert (t == 0).all()
        return out


def class_values(sl, beta, batch=4000):
    """sum_k PAF_k at the class reps for each row of beta"""
    out = []
    for s in range(0, len(beta), batch):
        X = sl.expand(beta[s:s + batch].astype(np.float64))
        F = np.fft.rfft(X, axis=2); P = np.fft.irfft(F.real ** 2 + F.imag ** 2, n=sl.m, axis=2).sum(1)
        out.append(np.rint(P[:, sl.reps]).astype(np.int64))
    return np.vstack(out)


def sample_slice(sl, N, rng):
    """N uniform configurations of the slice with A = +-2 (sign uniform in proportion to counts), others 0"""
    nA = {sA: count_signed(sl.gw["A"][1], sA) for sA in (2, -2)}
    tot = nA[2] + nA[-2]
    beta = np.empty((N, sl.nbits), np.int8)
    k2 = int(rng.binomial(N, nA[2] / tot)) if tot else 0
    parts = [(2, k2), (-2, N - k2)]
    row = 0
    for sA, k in parts:
        if k == 0: continue
        for g, bits in sl.groups.items():
            T = sA if g == "A" else 0
            sm = SignedSampler(sl.gw[g][1], T)
            beta[row:row + k, bits] = sm.sample(k, rng)
        row += k
    return beta


def log2_rowsum_configs(sl):
    tot = 0
    for t in sl.targets():
        c = 1
        for g, (k, w) in sl.gw.items(): c *= count_signed(w, t[g])
        tot += c
    return math.log2(tot) if tot else float("-inf")


def gaussian(sl, N=40000, seed=1):
    """log2 E (all admissible placements of the +-2 block) and conditional constraint bits; Gaussian on the classes
    other than {m/2} (whose value is implied by the row sums), per-coordinate lattice steps from the samples"""
    rec = dict(m=sl.m, v=sl.m + 1, u=sl.u, d=sl.d, tw=sl.tw, orbits=sl.r, classes=sl.L, bits=sl.nbits,
               minus1_in_U=sl.minus1_in_U)
    lrs = log2_rowsum_configs(sl)
    rec["log2_rowsum_configs"] = round(lrs, 2)
    if lrs == float("-inf"):
        rec.update(feasible=False, why="row sums"); return rec
    rng = np.random.default_rng(seed)
    beta = sample_slice(sl, N, rng)
    F = class_values(sl, beta)
    keep = [c for c in range(sl.L) if c != sl.half_class]
    F = F[:, keep]; tgt = -4 * np.ones(len(keep))
    steps = np.array([max(1, int(np.gcd.reduce(np.abs(F[:, c] - F[0, c])))) for c in range(len(keep))])
    on = all(((tgt[c] - F[0, c]) % steps[c]) == 0 for c in range(len(keep)))
    if not on:
        rec.update(feasible=False, why="target off the per-class lattice"); return rec
    mu = F.mean(0); S = np.atleast_2d(np.cov(F.T))
    sign, logdet = np.linalg.slogdet(S)
    if sign <= 0:
        rec.update(feasible=None, why="degenerate covariance"); return rec
    dlt = tgt - mu; q = dlt @ np.linalg.solve(S, dlt)
    log2P = (np.log(steps).sum() - 0.5 * len(keep) * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q) / math.log(2)
    rec.update(feasible=True, n_cond=round(-log2P, 1), log2E=round(math.log2(sl.sym) + lrs + log2P, 1),
               steps=sorted(set(steps.tolist())), mahal=round(float(q), 2))
    return rec
