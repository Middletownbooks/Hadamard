"""Skew Goethals-Seidel multiplier slices over abelian groups of odd order v (prompt §2).

A slice = (G, H): G = Z_n1 x ... x Z_nk (|G| = v odd), H an odd-order group of automorphisms given by generators.
Blocks A, B, C, D are +-1 on G, constant on H-orbits; A is skew (A(0) = 1, A(-g) = -A(g)).
Sum of PAFs = 4v delta  <=>  skew-Hadamard matrix of order 4v (Goethals-Seidel array).

Contents
- Struct: orbits, negation pairing, lag classes (orbits of <H, -1> on G \\ 0), class sizes.
- rowsum_counts / sample_block: exact DP over orbit sizes (counts and uniform sampling with a given row sum).
- paf_classes: PAF of a batch of blocks at one representative per lag class (FFT, exact after rounding).
- gauss_estimate: lattice-Gaussian expected number of ordered (A, B, C, D) with fixed positive row sums [E].
- zp_gate: exact compression test to Z_p for a prime p | v along a cyclic factor (DP per coordinate, 2+2 join) [M].
Written for this repository (no code from Kumar's repository).
"""
import itertools, math
import numpy as np

# ---------------------------------------------------------------- groups and automorphisms
class Group:
    """G = prod Z_{n_i}; element index = mixed radix (first coordinate most significant)."""
    def __init__(self, mods):
        self.mods = tuple(mods); self.v = int(np.prod(mods))
        self.coords = np.array(list(itertools.product(*[range(n) for n in mods])), dtype=np.int64).reshape(self.v, len(mods))
        self.rad = np.array([int(np.prod(mods[i + 1:])) for i in range(len(mods))], dtype=np.int64)
    def index(self, C):
        C = np.asarray(C) % np.array(self.mods)
        return (C * self.rad).sum(-1)
    def neg(self):
        return self.index(-self.coords)
    def perm(self, aut):
        """aut: ('mat', block structure) given as function coords -> coords."""
        return self.index(aut(self.coords))

def aut_diag(mods, units):
    """x_i -> u_i x_i on each cyclic factor."""
    u = np.array(units, dtype=np.int64)
    return lambda C: C * u

def aut_block(mods, mats):
    """mats: list of (slice_start, slice_end, matrix) acting on coordinate blocks; coordinates outside are fixed."""
    def f(C):
        D = C.copy()
        for s, e, M in mats:
            D[:, s:e] = C[:, s:e] @ np.array(M, dtype=np.int64).T
        return D
    return f

class Struct:
    def __init__(self, G, gens, label=''):
        self.G, self.v, self.label = G, G.v, label
        v = self.v
        self.gens = [np.asarray(g, dtype=np.int64) for g in gens]
        parent = np.arange(v)
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for g in self.gens:
            for x in range(v):
                a, b = find(x), find(g[x])
                if a != b: parent[max(a, b)] = min(a, b)
        root = np.array([find(x) for x in range(v)])
        reps, self.orb = np.unique(root, return_inverse=True)
        self.r = len(reps)                                   # orbits per block
        self.osize = np.bincount(self.orb)
        neg = G.neg(); self.negidx = neg
        self.onegate = self.orb[neg[reps]]                   # orbit of -O
        assert self.osize[self.orb[0]] == 1
        self.zero = self.orb[0]
        if np.any((self.onegate == np.arange(self.r)) & (np.arange(self.r) != self.zero)):
            self.skew_ok = False                              # some orbit is closed under negation
        else:
            self.skew_ok = True
        # H order from the orbit of a generic element is not needed; record |H| as max orbit size
        self.Hord = int(self.osize.max())
        # lag classes: orbits of <H, -1> on G \ 0
        cl = np.minimum(np.arange(self.r), self.onegate)
        nz = [o for o in range(self.r) if o != self.zero]
        cids = sorted(set(int(cl[o]) for o in nz))
        self.C = len(cids)
        self.class_rep = np.array([reps[c] for c in cids])  # one lag per class
        cmap = {c: i for i, c in enumerate(cids)}
        self.class_of_orbit = np.array([cmap.get(int(cl[o]), -1) for o in range(self.r)])
        self.csize = np.array([self.osize[o] * (1 if self.onegate[o] == o else 2) for o in cids])
        # forced products: f(s) = #{g : g + s in the H-orbit of g} at each class representative s.
        # Every H-invariant block has PAF(s) >= 2 f(s) - v (skew A too: products with g + s in -Hg are forced -1 and
        # do not change the bound), so 2 f(s) > v at any class kills the slice [P] (gate L).
        Gc = G.coords; mods = np.array(G.mods)
        self.forced = np.array([int((self.orb[G.index(Gc + Gc[s])] == self.orb).sum()) for s in self.class_rep])
        self.gateL = bool(np.all(2 * self.forced <= v))
        # skew pairs: free sign per pair {O, -O} with O < -O
        self.pairs = [o for o in range(self.r) if o != self.zero and o < self.onegate[o]]

    # ------------------------------------------------------------ blocks
    def expand(self, S):
        """S: (n, r) orbit signs -> (n, v) blocks."""
        return S[:, self.orb]

    def skew_signs(self, P):
        """P: (n, len(pairs)) signs -> orbit signs of a skew block."""
        n = P.shape[0]; S = np.zeros((n, self.r), dtype=np.int8)
        S[:, self.zero] = 1
        S[:, self.pairs] = P
        S[:, self.onegate[self.pairs]] = -P
        return S

    def paf_classes(self, X):
        """X: (n, v) +-1 -> (n, C) PAF at the class representatives (exact integers)."""
        sh = self.G.mods
        Y = X.reshape((X.shape[0],) + sh).astype(np.float64)
        F = np.fft.rfftn(Y, axes=tuple(range(1, len(sh) + 1)))
        P = np.fft.irfftn(np.abs(F) ** 2, s=sh, axes=tuple(range(1, len(sh) + 1))).reshape(X.shape[0], self.v)
        return np.rint(P[:, self.class_rep]).astype(np.int64)

# ---------------------------------------------------------------- row sums
def rowsum_counts(sizes, fixed=None):
    """Counts of +-1 sign vectors on orbits of the given sizes by row sum. Returns dict sum -> count (python int)."""
    cnt = {0: 1}
    for s in sizes:
        nc = {}
        for k, c in cnt.items():
            nc[k + s] = nc.get(k + s, 0) + c
            nc[k - s] = nc.get(k - s, 0) + c
        cnt = nc
    return cnt

def block_rowsums(st):
    """Free block: dict row sum -> count."""
    return rowsum_counts(list(st.osize))

class RowSampler:
    """Uniform sampling of orbit-sign vectors of a free block with a given row sum (exact DP, float counts)."""
    def __init__(self, st, rsum):
        self.st, self.r = st, rsum
        sizes = list(st.osize); n = len(sizes); M = int(sum(sizes))
        self.off = M
        T = np.zeros((n + 1, 2 * M + 1)); T[n, M + 0] = 1.0     # T[i][s] = #ways orbits i.. reach remaining s
        for i in range(n - 1, -1, -1):
            s = sizes[i]
            T[i, s:] += T[i + 1, :2 * M + 1 - s]
            T[i, :2 * M + 1 - s] += T[i + 1, s:]
        self.T, self.sizes = T, sizes
        self.count = T[0, M + rsum]
    def sample(self, n, rng):
        sizes, T, M = self.sizes, self.T, self.off
        out = np.empty((n, len(sizes)), dtype=np.int8)
        rem = np.full(n, self.r)
        for i, s in enumerate(sizes):
            # choose +1 with prob T[i+1][rem - s] / T[i][rem]
            a = T[i + 1, np.clip(rem - s + M, 0, 2 * M)] * ((rem - s + M >= 0) & (rem - s + M <= 2 * M))
            b = T[i + 1, np.clip(rem + s + M, 0, 2 * M)] * ((rem + s + M >= 0) & (rem + s + M <= 2 * M))
            p = a / (a + b)
            up = rng.random(n) < p
            out[:, i] = np.where(up, 1, -1)
            rem = rem - np.where(up, s, -s)
        assert np.all(rem == 0)
        return out

def admissible_patterns(st, v):
    """Unordered positive (b, c, d) with b^2 + c^2 + d^2 = 4v - 1 and each attainable by a free block."""
    rs = block_rowsums(st)
    att = sorted(k for k in rs if k > 0 or (k == 0))
    att = [k for k in att if k >= 0]
    pats = []
    T = 4 * v - 1
    S = set(att)
    for b in att:
        for c in att:
            if c < b: continue
            d2 = T - b * b - c * c
            if d2 < c * c: continue
            d = int(math.isqrt(d2))
            if d * d == d2 and d in S: pats.append((b, c, d))
    return pats, rs

# ---------------------------------------------------------------- Gaussian estimate
def class_stats(st, X):
    P = st.paf_classes(X).astype(np.float64)
    return P.mean(0), np.cov(P, rowvar=False) if P.shape[0] > 1 else None

def gauss_estimate(st, pattern, nsamp=3000, rng=None, pools=None):
    """log10 expected number of ordered quadruples (A skew; B, C, D with row sums pattern) [E].
    Lattice: every class value of the sum lies in 4Z, on the hyperplane sum_c n_c P(c) = 0."""
    rng = rng or np.random.default_rng(1)
    v = st.v; C = st.C
    npairs = len(st.pairs)
    PA = rng.choice(np.array([-1, 1], dtype=np.int8), size=(nsamp, npairs))
    XA = st.expand(st.skew_signs(PA))
    mu, cov = class_stats(st, XA)
    logN = npairs * math.log10(2)
    for r in pattern:
        smp = RowSampler(st, r)
        if smp.count <= 0: return None
        logN += math.log10(smp.count)
        X = st.expand(smp.sample(nsamp, rng))
        m, c = class_stats(st, X)
        mu = mu + m; cov = cov + c
    n = st.csize.astype(np.float64)
    # orthonormal basis of the hyperplane n.y = 0
    Q, _ = np.linalg.qr(np.column_stack([n] + [np.eye(C)[:, i] for i in range(C)]))
    B = Q[:, 1:C]                                  # C x (C-1)
    S2 = B.T @ cov @ B
    z = B.T @ (0 - mu)
    sign, logdet = np.linalg.slogdet(S2)
    if sign <= 0: return None
    quad = float(z @ np.linalg.solve(S2, z))
    logdens = -0.5 * quad - 0.5 * (logdet + (C - 1) * math.log(2 * math.pi))
    g = np.gcd.reduce(st.csize.astype(np.int64))
    logcov = (C - 1) * math.log(4) + math.log(np.linalg.norm(n) / g)
    return logN + (logdens + logcov) / math.log(10), logN

# ---------------------------------------------------------------- exact Z_p compression gate (cyclic factor)
def zp_gate(st, p, axis, pattern, cap=300_000, joincap=3 * 10 ** 7):
    """p: any divisor d of the order of the cyclic factor `axis` (prime or not). Reduce coordinate `axis` mod p and
    sum over everything else. H must act on that factor by units (diagonal on it) -- checked by the caller.
    Returns (status, n_solutions_up_to_cap) with status 'pass', 'dead', 'undecided'."""
    G = st.G; v = st.v
    res = G.coords[:, axis] % p
    # H acts on Z_p through the induced action: x mod p -> (g x) mod p; get it from the generators
    hp = []
    for g in st.gens:
        img = G.coords[g][:, axis] % p
        # image of residue 1: take element with coords e_axis
        e = np.zeros(len(G.mods), dtype=np.int64); e[axis] = 1
        u = int(img[G.index(e)])
        if not np.all(img == (u * res) % p): return 'n/a', None, ['not equivariant']
        hp.append(u)
    # orbits of <hp> on Z_p
    seen = {}; porb = []
    for x in range(p):
        if x in seen: continue
        o = {x}; fr = [x]
        while fr:
            y = fr.pop()
            for u in hp:
                z = (u * y) % p
                if z not in o: o.add(z); fr.append(z)
        for y in o: seen[y] = len(porb)
        porb.append(sorted(o))
    k = len(porb); psz = np.array([len(o) for o in porb])
    if (k - 1) // 2 > 14: return 'undecided', None, ['too many coordinates %d' % k]
    pneg = [seen[(-o[0]) % p] for o in porb]
    # contribution of each G-orbit to coordinate o: |O| / |o|
    contrib = {}
    for O in range(st.r):
        x = int(np.flatnonzero(st.orb == O)[0])
        o = seen[int(res[x])]
        w = st.osize[O] // psz[o]
        contrib.setdefault(o, []).append((O, int(w)))
    def coord_sets(skew):
        sets = []
        for o in range(k):
            if skew:
                if o == seen[0]: sets.append({1}); continue
                if pneg[o] < o: sets.append(None); continue
                ws = [w for (O, w) in contrib.get(o, [])]
                s = {0}
                for w in ws: s = {a + w for a in s} | {a - w for a in s}
                sets.append(s)
            else:
                s = {0}
                for (O, w) in contrib.get(o, []): s = {a + w for a in s} | {a - w for a in s}
                sets.append(s)
        return sets
    E = 4 * v
    def enum(skew, rsum, emax):
        sets = coord_sets(skew)
        free = [o for o in range(k) if sets[o] is not None]
        out = []
        vals = [sorted(sets[o], key=abs) for o in free]
        # energy weight: coordinate o contributes |o| y_o^2 (and for skew its negative partner too)
        wt = [psz[o] * (2 if (skew and o != seen[0]) else 1) for o in free]
        rw = [psz[o] * (0 if skew else 1) for o in free]
        y = [0] * len(free)
        nf = len(free)
        smin = [0] * (nf + 1); smax = [0] * (nf + 1)
        for i in range(nf - 1, -1, -1):
            smin[i] = smin[i + 1] + min(rw[i] * a for a in vals[i]); smax[i] = smax[i + 1] + max(rw[i] * a for a in vals[i])
        def rec(i, en, rs):
            if len(out) > cap: return
            if not skew and not (smin[i] <= rsum - rs <= smax[i]): return
            if i == len(free):
                if (skew or rs == rsum) and en <= emax:
                    full = [0] * k
                    for j, o in enumerate(free):
                        full[o] = y[j]
                        if skew and o != seen[0]: full[pneg[o]] = -y[j]
                    out.append(tuple(full))
                return
            for a in vals[i]:
                e2 = en + wt[i] * a * a
                if e2 > emax: break
                y[i] = a
                rec(i + 1, e2, rs + rw[i] * a)
        rec(0, 0, 0)
        return out
    # PAF on Z_p of a vector constant on porb, at lag classes of <hp, -1>
    lagcls = sorted(set(min(o, pneg[o]) for o in range(k) if o != seen[0]))
    def pafs(L):
        Y = np.array(L, dtype=np.int64)
        full = Y[:, [seen[x] for x in range(p)]]
        P = np.stack([(full * np.roll(full, -porb[c][0], axis=1)).sum(1) for c in lagcls], 1)
        e0 = (full * full).sum(1)
        return np.column_stack([e0, P])
    lists = []
    emin_other = 0
    # energy lower bound for a block with row sum r on Z_p: r^2 / p
    rows = [1] + list(pattern)
    lbs = [ (r * r + p - 1) // p for r in rows]
    for i, r in enumerate(rows):
        emax = E - (sum(lbs) - lbs[i])
        L = enum(i == 0, r, emax)
        if len(L) > cap: return 'undecided', None, [len(x) for x in lists] + ['cap']
        if not L: return 'dead', 0, [len(x) for x in lists] + [0]
        lists.append(L)
    P = [np.unique(pafs(L), axis=0) for L in lists]
    tgt = np.zeros(P[0].shape[1], dtype=np.int64); tgt[0] = E
    # 2+2 join on (A,B) vs (C,D) by hashing
    nab, ncd = len(P[0]) * len(P[1]), len(P[2]) * len(P[3])
    if min(nab, ncd) > joincap or max(nab, ncd) > 100 * joincap:
        return 'undecided', None, [len(x) for x in P] + ['joincap']
    # chunked hash join: 64-bit random linear hash of the vectors, exact re-check of every hash match
    D = P[0].shape[1]
    R = np.random.default_rng(12345).integers(1, 2 ** 61, size=D, dtype=np.uint64)
    hsh = lambda X: (X.astype(np.uint64) * R).sum(1)
    if nab > ncd: P = [P[2], P[3], P[0], P[1]]
    small = (P[0][:, None, :] + P[1][None, :, :]).reshape(-1, D)
    h0 = hsh(small); order = np.argsort(h0); hs = h0[order]
    hD = hsh(P[3]); n = 0
    for i in range(len(P[2])):
        need = hsh((tgt - P[2][i])[None, :])[0] - hD          # hash of tgt - C_i - D_j (linear, mod 2^64)
        pos = np.searchsorted(hs, need); pos[pos >= len(hs)] = len(hs) - 1
        hit = np.flatnonzero(hs[pos] == need)
        for j in hit:                                        # exact re-check on the rows with an equal hash
            row = tgt - P[2][i] - P[3][j]; p0 = pos[j]
            while p0 < len(hs) and hs[p0] == need[j]:
                if np.array_equal(small[order[p0]], row): n += 1; break
                p0 += 1
    return ('pass' if n else 'dead'), n, [len(x) for x in P]

# ---------------------------------------------------------------- analytic in-list estimate (k-tree cost model)
def inlist_estimate(st, pattern, Ls, nsamp=4000, rng=None, samplers=None):
    """log10 expected solutions inside low-energy lists of sizes Ls = (L_A, L_B, L_C, L_D) [E].
    Model: block k's keys ~ N(mu_k, S_k) (uniform samples); a list keeps the members with |W_k (P - c_k)|^2 <= tau_k,
    c_k = mu_k + S_k S^+ (0 - sum mu) and tau_k the noncentral chi^2 quantile of L_k / N_k; the kept keys are taken as
    uniform in that ellipsoid (mean c_k, covariance S_k tau_k / (d + 2)); a complete list keeps N(mu_k, S_k).
    Returns (log10 E_in, log10 pair work L_A L_B + L_C L_D (best pairing), details)."""
    from scipy.stats import ncx2
    rng = rng or np.random.default_rng(1)
    C = st.C
    if samplers is not None:
        Ns = [s.count for s in samplers]; Xs = [st.expand(s.sample(nsamp, rng)) for s in samplers]
    else:
        PA = rng.choice(np.array([-1, 1], dtype=np.int8), size=(nsamp, len(st.pairs)))
        Xs = [st.expand(st.skew_signs(PA))]
        Ns = [2.0 ** len(st.pairs)]
        for r in pattern:
            smp = RowSampler(st, r)
            Ns.append(smp.count); Xs.append(st.expand(smp.sample(nsamp, rng)))
    stats = [class_stats(st, X) for X in Xs]
    mu = [s[0] for s in stats]; Sg = [s[1] for s in stats]
    Sp = np.linalg.pinv(sum(Sg), rcond=1e-9)
    cen = [mu[k] + Sg[k] @ Sp @ (0 - sum(mu)) for k in range(4)]
    means, covs, used = [], [], []
    for k in range(4):
        L = min(Ls[k], Ns[k]); used.append(L)
        if L >= Ns[k]:
            means.append(mu[k]); covs.append(Sg[k]); continue
        w, V = np.linalg.eigh(Sg[k]); keep = w > 1e-6 * w.max(); d = int(keep.sum())
        Wk = (V[:, keep] / np.sqrt(w[keep])) @ V[:, keep].T
        lam = float(np.sum((Wk @ (mu[k] - cen[k])) ** 2))
        tau = float(ncx2.ppf(L / Ns[k], d, lam))
        means.append(cen[k]); covs.append(Sg[k] * tau / (d + 2))
    n = st.csize.astype(np.float64)
    Q, _ = np.linalg.qr(np.column_stack([n] + [np.eye(C)[:, i] for i in range(C)])); B = Q[:, 1:C]
    S2 = B.T @ sum(covs) @ B; z = B.T @ (0 - sum(means))
    sign, logdet = np.linalg.slogdet(S2)
    logd = -0.5 * float(z @ np.linalg.solve(S2, z)) - 0.5 * (logdet + (C - 1) * math.log(2 * math.pi))
    g = np.gcd.reduce(st.csize.astype(np.int64))
    logE = (sum(math.log(L) for L in used) + logd + (C - 1) * math.log(4) + math.log(np.linalg.norm(n) / g)) / math.log(10)
    a, b, c, d_ = [math.log10(L) for L in used]
    W = min(math.log10(10 ** (a + b - 30) + 10 ** (c + d_ - 30)) + 30, math.log10(10 ** (a + c - 30) + 10 ** (b + d_ - 30)) + 30,
            math.log10(10 ** (a + d_ - 30) + 10 ** (b + c - 30)) + 30)
    return logE, W, dict(N=[math.log10(x) for x in Ns], L=[math.log10(x) for x in used])

# ---------------------------------------------------------------- stage-1 (Z_d compression) correction of the estimate
def stage1_ratio(st, d, axis, pattern, nsamp=20000, rng=None, cap=2_000_000):
    """R = (exact weighted number of stage-1 quadruples) / (Gaussian prediction of the same number) at the compression
    to Z_d along `axis` [M numerator / E denominator]. Exact: sum over compressed quadruples (y_A, y_B, y_C, y_D) with
    sum_X PAF_{y_X} = 4v delta on Z_d of prod_X N_X(y_X) (N_X(y) = number of H-invariant blocks with compression y and
    the row sum). Gaussian: the same lattice-Gaussian count applied to the compressed PAF vectors of uniform blocks.
    E_gauss * R is the staged estimate (stage 1 exact, the rest Gaussian), assuming the Gaussian is right given y."""
    from collections import defaultdict
    rng = rng or np.random.default_rng(1)
    G = st.G; v = st.v
    res = G.coords[:, axis] % d
    e = np.zeros(len(G.mods), dtype=np.int64); e[axis] = 1
    hp = []
    for g in st.gens:
        img = G.coords[g][:, axis] % d; u = int(img[G.index(e)])
        if not np.all(img == (u * res) % d): return None
        hp.append(u)
    seen = {}; porb = []
    for x in range(d):
        if x in seen: continue
        o = {x}; fr = [x]
        while fr:
            y = fr.pop()
            for u in hp:
                z = (u * y) % d
                if z not in o: o.add(z); fr.append(z)
        for y in o: seen[y] = len(porb)
        porb.append(sorted(o))
    k = len(porb); psz = np.array([len(o) for o in porb]); pneg = [seen[(-o[0]) % d] for o in porb]
    z0 = seen[0]
    contrib = defaultdict(list)
    for O in range(st.r):
        x = int(np.flatnonzero(st.orb == O)[0]); o = seen[int(res[x])]
        contrib[o].append(int(st.osize[O] // psz[o]))
    def counts(ws):
        c = {0: 1}
        for w in ws:
            n = defaultdict(int)
            for a, m in c.items(): n[a + w] += m; n[a - w] += m
            c = dict(n)
        return c
    # per-coordinate value distributions
    free_c = [counts(contrib[o]) for o in range(k)]
    def block_dist(skew, rsum):
        """dict y (tuple over coordinates) -> multiplicity, with the row sum (free) or skew structure; capped."""
        if skew:
            coords = [o for o in range(k) if o != z0 and o < pneg[o]]
            items = [(None, {(): 1})]
            dist = {(): 1}
            for o in coords:
                nd = {}
                for y, m in dist.items():
                    for a, ma in free_c[o].items():
                        nd[y + (a,)] = m * ma
                dist = nd
                if len(dist) > cap: return None
            out = {}
            for y, m in dist.items():
                full = [0] * k; full[z0] = 1
                for o, a in zip(coords, y): full[o] = a; full[pneg[o]] = -a
                out[tuple(full)] = m
            return out
        dist = {(): (0, 1)}
        dist = {((), 0): 1}
        for o in range(k):
            nd = defaultdict(int)
            for (y, rs), m in dist.items():
                for a, ma in free_c[o].items():
                    nd[(y + (a,), rs + psz[o] * a)] += m * ma
            dist = nd
            if len(dist) > cap: return None
        return {y: m for (y, rs), m in dist.items() if rs == rsum}
    lagcls = sorted(set(min(o, pneg[o]) for o in range(k) if o != z0))
    full_idx = [seen[x] for x in range(d)]
    def cpaf(Y):
        F = Y[:, full_idx]
        return np.column_stack([(F * F).sum(1)] + [(F * np.roll(F, -porb[c][0], axis=1)).sum(1) for c in lagcls])
    dists = [block_dist(True, 1)] + [block_dist(False, r) for r in pattern]
    if any(x is None for x in dists): return None
    P, M = [], []
    for D_ in dists:
        ys = np.array(list(D_.keys()), dtype=np.int64); ms = np.array([float(m) for m in D_.values()])
        pv = cpaf(ys)
        # merge equal PAF vectors
        u, inv = np.unique(pv, axis=0, return_inverse=True)
        w = np.zeros(len(u)); np.add.at(w, inv.ravel(), ms)
        P.append(u); M.append(w)
    tgt = np.zeros(P[0].shape[1], dtype=np.int64); tgt[0] = 4 * v
    # exact weighted join
    from collections import Counter
    AB = defaultdict(float)
    for i in range(len(P[0])):
        for j in range(len(P[1])): AB[tuple(P[0][i] + P[1][j])] += M[0][i] * M[1][j]
    if len(AB) > 5 * 10 ** 7: return None
    exact = 0.0
    for i in range(len(P[2])):
        for j in range(len(P[3])):
            key = tuple(tgt - P[2][i] - P[3][j])
            if key in AB: exact += AB[key] * M[2][i] * M[3][j]
    # Gaussian on the compressed PAF vectors of uniform blocks (same estimator as gauss_estimate, on Z_d)
    Xs = [st.expand(st.skew_signs(rng.choice(np.array([-1, 1], np.int8), size=(nsamp, len(st.pairs)))))]
    for r in pattern: Xs.append(st.expand(RowSampler(st, r).sample(nsamp, rng)))
    Ys = []
    for X in Xs:
        Y = np.zeros((len(X), k), dtype=np.int64)
        for o in range(k):
            x = porb[o][0]; Y[:, o] = X[:, res == x].sum(1)
        Ys.append(cpaf(Y).astype(float))
    mu = sum(y.mean(0) for y in Ys); cov = sum(np.cov(y, rowvar=False) for y in Ys)
    Ntot = sum(math.log(sum(m)) for m in M)
    # hyperplane: n . P = sum r^2 with n = (1, |class sizes| ...): lag 0 counts once, class c counts 2|o|
    n = np.array([1.0] + [float(2 * psz[c]) for c in lagcls])
    Cd = len(n)
    Q, _ = np.linalg.qr(np.column_stack([n] + [np.eye(Cd)[:, i] for i in range(Cd)])); B = Q[:, 1:Cd]
    S2 = B.T @ cov @ B; z = B.T @ (tgt - mu)
    sgn, ld = np.linalg.slogdet(S2)
    if sgn <= 0: return None
    # lattice of the summed compressed PAF vectors: gcd of differences per coordinate (empirical), projected
    tot = sum(y for y in Ys); diffs = (tot - tot[0]).astype(np.int64)
    gs = [int(np.gcd.reduce(np.abs(diffs[:, i]))) or 1 for i in range(Cd)]
    logd = -0.5 * float(z @ np.linalg.solve(S2, z)) - 0.5 * (ld + (Cd - 1) * math.log(2 * math.pi))
    logcov = sum(math.log(g) for g in gs) - math.log(gs[0]) + math.log(np.linalg.norm(n) / np.gcd.reduce(n.astype(np.int64)))
    gauss = Ntot + logd + logcov
    return (math.log10(exact) if exact > 0 else -np.inf) - gauss / math.log(10), dict(exact_log10=math.log10(exact) if exact else None,
                                                                              gauss_log10=gauss / math.log(10), k=k, lattice=gs)

# ---------------------------------------------------------------- stage-1 classes (exact) and class-conditional sampling
def compression_map(st, d, axis):
    """Coordinates of the compression to Z_d along `axis`: returns (k, coord_of_orbit[r], weight_of_orbit[r], z0, pneg)
    with y_o = sum_{O -> o} e_O w_O; None if H does not act on that coordinate."""
    G = st.G
    res = G.coords[:, axis] % d
    e = np.zeros(len(G.mods), dtype=np.int64); e[axis] = 1
    hp = []
    for g in st.gens:
        img = G.coords[g][:, axis] % d; u = int(img[G.index(e)])
        if not np.all(img == (u * res) % d): return None
        hp.append(u)
    seen = {}; porb = []
    for x in range(d):
        if x in seen: continue
        o = {x}; fr = [x]
        while fr:
            y = fr.pop()
            for u in hp:
                z = (u * y) % d
                if z not in o: o.add(z); fr.append(z)
        for y in o: seen[y] = len(porb)
        porb.append(sorted(o))
    k = len(porb); psz = np.array([len(o) for o in porb]); pneg = [seen[(-o[0]) % d] for o in porb]
    co = np.zeros(st.r, dtype=np.int64); wo = np.zeros(st.r, dtype=np.int64)
    for O in range(st.r):
        x = int(np.flatnonzero(st.orb == O)[0]); co[O] = seen[int(res[x])]; wo[O] = st.osize[O] // psz[co[O]]
    return dict(k=k, co=co, wo=wo, z0=seen[0], pneg=pneg, psz=psz, porb=porb, seen=seen, d=d)

def _sumdist(ws):
    from collections import defaultdict
    c = {0: 1}
    for w in ws:
        n = defaultdict(int)
        for a, m in c.items(): n[a + w] += m; n[a - w] += m
        c = dict(n)
    return c

def stage1_classes(st, d, axis, pattern, cap=2_000_000):
    """All stage-1 classes (y_A, y_B, y_C, y_D) at the compression to Z_d along `axis` with sum PAF_y = 4v delta,
    each with its exact weight prod_X N_X(y_X) (number of block quadruples in the class) [M]. Sorted by weight."""
    from collections import defaultdict
    cm = compression_map(st, d, axis)
    if cm is None: return None
    k, co, wo, z0, pneg, psz = cm['k'], cm['co'], cm['wo'], cm['z0'], cm['pneg'], cm['psz']
    fc = [_sumdist([int(wo[O]) for O in range(st.r) if co[O] == o]) for o in range(k)]
    def dist(skew, rsum):
        if skew:
            coords = [o for o in range(k) if o != z0 and o < pneg[o]]
            D = {(): 1}
            for o in coords:
                D = {y + (a,): m * ma for y, m in D.items() for a, ma in fc[o].items()}
                if len(D) > cap: return None
            out = {}
            for y, m in D.items():
                full = [0] * k; full[z0] = 1
                for o, a in zip(coords, y): full[o] = a; full[pneg[o]] = -a
                out[tuple(full)] = m
            return out
        D = {((), 0): 1}
        for o in range(k):
            nd = defaultdict(int)
            for (y, rs), m in D.items():
                for a, ma in fc[o].items(): nd[(y + (a,), rs + int(psz[o]) * a)] += m * ma
            D = nd
            if len(D) > cap: return None
        return {y: m for (y, rs), m in D.items() if rs == rsum}
    dists = [dist(True, 1)] + [dist(False, r) for r in pattern]
    if any(x is None for x in dists): return None
    lagcls = sorted(set(min(o, pneg[o]) for o in range(k) if o != z0))
    full_idx = [cm['seen'][x] for x in range(d)]
    def cpaf(y):
        F = np.array(y)[full_idx]
        return tuple([int((F * F).sum())] + [int((F * np.roll(F, -cm['porb'][c][0])).sum()) for c in lagcls])
    byp = []
    for D in dists:
        g = defaultdict(list)
        for y, m in D.items(): g[cpaf(y)].append((y, m))
        byp.append(g)
    tgt = np.zeros(len(lagcls) + 1, dtype=np.int64); tgt[0] = 4 * st.v
    AB = defaultdict(list)
    for pa in byp[0]:
        for pb in byp[1]: AB[tuple(np.array(pa) + np.array(pb))].append((pa, pb))
    out = []
    for pc in byp[2]:
        for pd in byp[3]:
            key = tuple(tgt - np.array(pc) - np.array(pd))
            for pa, pb in AB.get(key, []):
                for ya, ma in byp[0][pa]:
                    for yb, mb in byp[1][pb]:
                        for yc, mc in byp[2][pc]:
                            for yd, md in byp[3][pd]:
                                out.append((math.log10(ma) + math.log10(mb) + math.log10(mc) + math.log10(md), [ya, yb, yc, yd]))
    out.sort(key=lambda t: -t[0])
    return out, cm

class ClassSampler:
    """Uniform sampling of a block (orbit signs) with prescribed compression y (and hence row sum): per coordinate o,
    exact DP over the orbits mapping to o (skew block: over the pairs {O, -O} with O -> o, o < -o)."""
    def __init__(self, st, cm, y, skew):
        self.st, self.skew = st, skew
        k, co, wo, z0, pneg = cm['k'], cm['co'], cm['wo'], cm['z0'], cm['pneg']
        self.parts = []
        for o in range(k):
            if skew:
                if o == z0 or o > pneg[o]: continue
                units = [O for O in st.pairs if co[O] == o] + [int(st.onegate[O]) for O in st.pairs if co[st.onegate[O]] == o]
                # each pair contributes to o through exactly one of its two orbits
                sg = [1] * sum(1 for O in st.pairs if co[O] == o) + [-1] * sum(1 for O in st.pairs if co[st.onegate[O]] == o)
                pr = [O for O in st.pairs if co[O] == o] + [O for O in st.pairs if co[st.onegate[O]] == o]
                ws = [int(wo[u]) for u in units]
                self.parts.append((pr, sg, ws, y[o]))
            else:
                orbs = [O for O in range(st.r) if co[O] == o]
                self.parts.append((orbs, [1] * len(orbs), [int(wo[O]) for O in orbs], y[o]))
        self.tabs = []
        self.count = 1.0
        for items, sg, ws, t in self.parts:
            M = int(sum(ws)); n = len(ws)
            T = np.zeros((n + 1, 2 * M + 1)); T[n, M] = 1.0
            for i in range(n - 1, -1, -1):
                s = ws[i]; T[i, s:] += T[i + 1, :2 * M + 1 - s]; T[i, :2 * M + 1 - s] += T[i + 1, s:]
            self.tabs.append((T, M))
            self.count *= T[0, M + t] if abs(t) <= M else 0.0
    def sample(self, n, rng):
        st = self.st
        if self.skew:
            P = np.zeros((n, len(st.pairs)), dtype=np.int8); pidx = {O: i for i, O in enumerate(st.pairs)}
        else:
            S = np.zeros((n, st.r), dtype=np.int8)
        for (items, sg, ws, t), (T, M) in zip(self.parts, self.tabs):
            rem = np.full(n, t)
            for i, s in enumerate(ws):
                ia = rem - s + M; ib = rem + s + M
                a = np.where((ia >= 0) & (ia <= 2 * M), T[i + 1, np.clip(ia, 0, 2 * M)], 0.0)
                b = np.where((ib >= 0) & (ib <= 2 * M), T[i + 1, np.clip(ib, 0, 2 * M)], 0.0)
                up = rng.random(n) < a / (a + b)
                val = np.where(up, 1, -1).astype(np.int8)     # contribution sign of this orbit to coordinate o
                rem = rem - np.where(up, s, -s)
                if self.skew: P[:, pidx[items[i]]] = val * sg[i]
                else: S[:, items[i]] = val
            assert np.all(rem == 0)
        return st.skew_signs(P) if self.skew else S

# ---------------------------------------------------------------- exact compression gate onto a product of axes
def quotient_gate(st, axes, pattern, cap=2_000_000, joincap=3 * 10 ** 7):
    """Compression of every block onto Q = prod_{i in axes} Z_{n_i} (sum over the other coordinates), when H induces
    an action on Q. Exact join of the compressed PAFs on Q with lag-0 energy 4v [M]. Returns (status, n, sizes)."""
    G = st.G; v = st.v
    qmods = [G.mods[i] for i in axes]; Q = Group(qmods)
    proj = Q.index(G.coords[:, axes])
    # induced action
    qperm = []
    for g in st.gens:
        img = np.full(Q.v, -1)
        pg = proj[g]
        for x in range(v):
            q = proj[x]
            if img[q] < 0: img[q] = pg[x]
            elif img[q] != pg[x]: return 'n/a', None, ['no induced action']
        qperm.append(img)
    seen = {}; porb = []
    for x in range(Q.v):
        if x in seen: continue
        o = {x}; fr = [x]
        while fr:
            y = fr.pop()
            for pm in qperm:
                z = int(pm[y])
                if z not in o: o.add(z); fr.append(z)
        for y in o: seen[y] = len(porb)
        porb.append(sorted(o))
    k = len(porb); psz = np.array([len(o) for o in porb])
    qneg = Q.neg(); pneg = [seen[int(qneg[o[0]])] for o in porb]; z0 = seen[0]
    if (k - 1) // 2 > 16: return 'undecided', None, ['too many coordinates %d' % k]
    contrib = {}
    for O in range(st.r):
        x = int(np.flatnonzero(st.orb == O)[0]); o = seen[int(proj[x])]
        contrib.setdefault(o, []).append(int(st.osize[O] // psz[o]))
    def sets(skew):
        out = []
        for o in range(k):
            if skew and o == z0: out.append({1}); continue
            if skew and pneg[o] < o: out.append(None); continue
            s = {0}
            for w in contrib.get(o, []): s = {a + w for a in s} | {a - w for a in s}
            out.append(s)
        return out
    E4 = 4 * v
    def enum(skew, rsum, emax):
        S = sets(skew); free = [o for o in range(k) if S[o] is not None]
        vals = [sorted(S[o], key=abs) for o in free]
        wt = [int(psz[o]) * (2 if (skew and o != z0) else 1) for o in free]
        rw = [int(psz[o]) * (0 if skew else 1) for o in free]
        nf = len(free); smin = [0] * (nf + 1); smax = [0] * (nf + 1)
        for i in range(nf - 1, -1, -1):
            smin[i] = smin[i + 1] + min(rw[i] * a for a in vals[i]); smax[i] = smax[i + 1] + max(rw[i] * a for a in vals[i])
        out = []; y = [0] * nf
        def rec(i, en, rs):
            if len(out) > cap: return
            if not skew and not (smin[i] <= rsum - rs <= smax[i]): return
            if i == nf:
                if (skew or rs == rsum) and en <= emax:
                    full = [0] * k
                    for j, o in enumerate(free):
                        full[o] = y[j]
                        if skew and o != z0: full[pneg[o]] = -y[j]
                    out.append(tuple(full))
                return
            for a in vals[i]:
                e2 = en + wt[i] * a * a
                if e2 > emax: break
                y[i] = a; rec(i + 1, e2, rs + rw[i] * a)
        rec(0, 0, 0)
        return out
    lagcls = sorted(set(min(o, pneg[o]) for o in range(k) if o != z0))
    tshift = [Q.index(Q.coords + Q.coords[porb[c][0]]) for c in lagcls]   # q -> q + t
    full_idx = np.array([seen[x] for x in range(Q.v)])
    def pafs(L):
        Y = np.array(L, dtype=np.int64)[:, full_idx]
        return np.column_stack([(Y * Y).sum(1)] + [(Y * Y[:, ts]).sum(1) for ts in tshift])
    rows = [1] + list(pattern); lbs = [(r * r + Q.v - 1) // Q.v for r in rows]
    lists = []
    for i, r in enumerate(rows):
        L = enum(i == 0, r, E4 - (sum(lbs) - lbs[i]))
        if len(L) > cap: return 'undecided', None, [len(x) for x in lists] + ['cap']
        if not L: return 'dead', 0, [len(x) for x in lists] + [0]
        lists.append(L)
    P = [np.unique(pafs(L), axis=0) for L in lists]
    tgt = np.zeros(P[0].shape[1], dtype=np.int64); tgt[0] = E4
    nab, ncd = len(P[0]) * len(P[1]), len(P[2]) * len(P[3])
    if min(nab, ncd) > joincap or max(nab, ncd) > 100 * joincap: return 'undecided', None, [len(x) for x in P] + ['joincap']
    D = P[0].shape[1]; R = np.random.default_rng(12345).integers(1, 2 ** 61, size=D, dtype=np.uint64)
    hsh = lambda X: (X.astype(np.uint64) * R).sum(1)
    if nab > ncd: P = [P[2], P[3], P[0], P[1]]
    small = (P[0][:, None, :] + P[1][None, :, :]).reshape(-1, D)
    h0 = hsh(small); order = np.argsort(h0); hs = h0[order]; hD = hsh(P[3]); n = 0
    for i in range(len(P[2])):
        need = hsh((tgt - P[2][i])[None, :])[0] - hD
        pos = np.searchsorted(hs, need); pos[pos >= len(hs)] = len(hs) - 1
        for j in np.flatnonzero(hs[pos] == need):
            row = tgt - P[2][i] - P[3][j]; p0 = pos[j]
            while p0 < len(hs) and hs[p0] == need[j]:
                if np.array_equal(small[order[p0]], row): n += 1; break
                p0 += 1
    return ('pass' if n else 'dead'), n, [len(x) for x in P]
