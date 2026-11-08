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
    if len(P[0]) * len(P[1]) > joincap or len(P[2]) * len(P[3]) > joincap:
        return 'undecided', None, [len(x) for x in P] + ['joincap']
    AB = (P[0][:, None, :] + P[1][None, :, :]).reshape(-1, P[0].shape[1])
    CD = (P[2][:, None, :] + P[3][None, :, :]).reshape(-1, P[0].shape[1])
    sAB = set(map(tuple, AB)); n = 0
    for row in map(tuple, tgt - CD):
        if row in sAB: n += 1
    return ('pass' if n else 'dead'), n, [len(x) for x in P]
