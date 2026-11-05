r"""[E, M sizes] The lift chain H -> G/Z_m -> G for general hosts of the k4 bordered family (t = 3; Kumar / EKLO Theorem 7):
G = Z2^2 x Z_m x Z_p, K = Z_m x Z_p, multiplier t -> a t on Z_p (order f, e = (p-1)/f classes), identity on Z2^2 x Z_m.
Used for the analogues of the v = 863 S7 host (m = 5, p = 43, a = 4; set-up and certificate by OpenAI Codex) and the
transfer to the other open orders.  Independent of bt_engine and of the Codex scripts.

For a y (H-level compression; the Codex certificate at 863, or a synthetic y with the H-level lag-0 energy and random
positions elsewhere), with blocks sampled uniformly in the y-fibre:
  log2 fibre(y)          exact (class-sign rows with fixed sums);
  C(all | y)             lattice Gaussian of the full class PAF vector (read at the mean for synthetic y; at the true
                         target for a real y), rank;
  C(w | y)               the same for the G/Z_m-level part (lags summed over the Z_m component);
  H(w | y)               entropy of the G/Z_m compression w (table-weighted), exact by DP over rows when the state
                         space (m+1)^e is at most 2e6, otherwise a lattice-Gaussian entropy (flagged);
  derived [E]: log2 E[cores | y] = fibre - C(all); log2 E[#w | y] = H(w) - C(w); log2 fibre(y, w) = fibre - H(w);
               C(mixed | y, w) = C(all) - C(w);  D_def(Z_m split) = C(mixed) - fibre(y, w)  (cores per (y,w) = 2^-D_def);
  523-style Z2 S/D split along z = (0,1): h_D exact, C(D|S) = C(all) - C(S), D_def = C(D|S) - h_D.
usage: python3 theory/h863_analog.py m p a [--y cert.json] [--N samples] [--seed s]
"""
import sys, os, json, math, argparse, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q523_gauss as QG


def gauss_bits_mean(V, nlat=400):
    """-log2 of the lattice-Gaussian density at the sample mean (bits, rank): the target-free version of
    q523_gauss.gauss_bits (same lattice-coordinate construction)"""
    V = np.asarray(V, np.int64)
    D = (V[1:nlat + 1] - V[0]).tolist()
    while True:
        basis, piv = QG.echelon(D)
        X = (V - V[0]).astype(np.int64); bad = np.zeros(len(V), bool)
        for i, p in enumerate(piv):
            b = np.array(basis[i], np.int64); off = (X[:, p] % b[p]) != 0
            bad |= off; t = np.where(off, 0, X[:, p] // b[p]); X -= t[:, None] * b[None, :]
        bad |= X.any(1)
        if not bad.any(): break
        D = basis + (V[np.nonzero(bad)[0][:50]] - V[0]).tolist()
    r = len(basis)
    if r == 0: return dict(bits=0.0, rank=0)
    B = np.array(basis, dtype=np.int64)
    C = np.zeros((V.shape[0], r), np.int64); X = (V - V[0]).astype(np.int64)
    for i, p in enumerate(piv):
        t = X[:, p] // B[i, p]; C[:, i] = t; X -= t[:, None] * B[i][None, :]
    S = np.cov(C.astype(np.float64).T).reshape(r, r)
    w = np.linalg.eigvalsh(S)
    lp = -0.5 * r * math.log(2 * math.pi) - 0.5 * np.log(w).sum()
    return dict(bits=round(-lp / math.log(2), 2), rank=r)


class Host:
    def __init__(self, m, p, a):
        self.m, self.p, self.a = m, p, a
        sub = [1]; x = a % p
        while x != 1:
            sub.append(x); x = x * a % p
        self.sub = sub; self.f = len(sub); self.e = (p - 1) // self.f
        cls = np.full(p, -1, np.int64); reps = []
        for t in range(1, p):
            if cls[t] < 0:
                for h in sub: cls[t * h % p] = len(reps)
                reps.append(t)
        self.cls, self.reps = cls, reps
        self.neg_in_S = (p - 1) in sub
        self.nc = 4 * m                       # cells of H = Z2^2 x Z_m, index q*m + u
        self.G = 4 * m * p; self.n = self.G + 3
        rep = {}; lreps = []
        for q in range(4):
            for u in range(m):
                for t in range(p):
                    if (q, u, t) in rep: continue
                    orb = {(q, u, t * h % p) for h in sub} | {(q, (-u) % m, (-t * h) % p) for h in sub}
                    for g in orb: rep[g] = len(lreps)
                    lreps.append((q, u, t))
        self.lagcls, self.lagreps, self.L = rep, lreps, len(lreps)
        self.lagrep_idx = np.array([(q * m + u) * p + t for q, u, t in lreps], np.int64)
        self.target = np.array([4 * self.G if g == (0, 0, 0) else (-12 if g[0] == 0 else 4) for g in lreps], np.int64)
        # linear maps: w-level (sum over u0) and Z2-S level (sum over q0 and q0 ^ 1), as class -> class lists
        self.wmap = self._sum_map(lambda q, u, t: [(q, uu, t) for uu in range(m)], key=lambda q, u, t: (q, t))
        self.smap = self._sum_map(lambda q, u, t: [(q, u, t), (q ^ 1, u, t)], key=lambda q, u, t: (q >> 1, u, t))
        self.tfull = None

    def _sum_map(self, members, key):
        """matrix M (k x L): row = one coarse lag (key), entries = multiplicity of each fine class among the members of
        one representative (the coarse value is sum of the fine PAF values over the members)"""
        rows = {}
        for q in range(4):
            for u in range(self.m):
                for t in range(self.p):
                    k = key(q, u, t)
                    if k in rows: continue
                    v = np.zeros(self.L, np.int64)
                    for g in members(q, u, t): v[self.lagcls[g]] += 1
                    rows[k] = v
        M = np.unique(np.array(list(rows.values())), axis=0)
        return M

    def assemble(self, e, c):
        """e: (N, 4, nc) signs, c: (N, 4, nc, e) class signs -> (N, 4, G) int8"""
        N = e.shape[0]
        x = np.empty((N, 4, self.nc, self.p), np.int8)
        x[..., 0] = e
        x[..., 1:] = c[..., self.cls[1:]]
        return x.reshape(N, 4, self.G)

    def class_paf(self, x, chunk=2048):
        lead = x.shape[:-1]; xf = x.reshape(-1, self.G)
        out = np.empty((len(xf), self.L), np.int64)
        for a0 in range(0, len(xf), chunk):
            xc = xf[a0:a0 + chunk].astype(np.float64).reshape(-1, 2, 2, self.m, self.p)
            F = np.fft.rfftn(xc, axes=(-4, -3, -2, -1))
            pa = np.fft.irfftn(F.real ** 2 + F.imag ** 2, s=(2, 2, self.m, self.p), axes=(-4, -3, -2, -1))
            out[a0:a0 + chunk] = np.rint(pa.reshape(len(xc), self.G)[:, self.lagrep_idx]).astype(np.int64)
        return out.reshape(lead + (self.L,))

    # ---------------------------------------------------------------- y
    def h_target0(self):
        return 4 * self.G - 12 * (self.p - 1)

    def synthetic_y(self, rng):
        """(E, S) per block and cell: S = 0 except a set of cells with S = +-2 chosen so that sum y^2 equals the H-level
        lag-0 target (y = e + f S); positions and signs random"""
        f = self.f; tgt = self.h_target0() - 4 * self.nc
        lo, hi = 4 * f * f - 4 * f, 4 * f * f + 4 * f       # y^2 - 1 for S = +-2, e S = -2 / +2
        best = None
        for b in range(0, 4 * self.nc + 1):
            rest = tgt - b * hi
            if rest < 0: break
            if rest % lo == 0 and rest // lo + b <= 4 * self.nc:
                best = (rest // lo, b); break
        if best is None:
            raise ValueError('no composition with S in {0, +-2}')
        na, nb = best
        E = rng.choice([-1, 1], size=(4, self.nc)); S = np.zeros((4, self.nc), np.int64)
        cells = rng.permutation(4 * self.nc)[:na + nb]
        for k, idx in enumerate(cells):
            i, h = divmod(int(idx), self.nc)
            s = 2 * rng.choice([-1, 1])
            S[i, h] = s
            if k < na: E[i, h] = -np.sign(s)        # e S = -2  (|y| = 2f - 1)
            else: E[i, h] = np.sign(s)              # e S = +2  (|y| = 2f + 1)
        if self.e < 2: raise ValueError
        return E, S, dict(n_small=na, n_large=nb)

    def sample_y_fibre(self, S, N, rng):
        c = np.empty((N, 4, self.nc, self.e), np.int8)
        for i in range(4):
            for h in range(self.nc):
                n = (int(S[i, h]) + self.e) // 2
                subs = np.array([[1 if j in s_ else 0 for j in range(self.e)] for s_ in it.combinations(range(self.e), n)],
                                np.int8)
                c[:, i, h] = 2 * subs[rng.integers(len(subs), size=N)] - 1
        return c

    def log2_fibre_y(self, S):
        return sum(math.log2(math.comb(self.e, (int(s) + self.e) // 2)) for s in S.ravel())

    # ---------------------------------------------------------------- H(w | y)
    def H_w_cell(self, rowsums):
        """entropy (bits) of the column-sum vector of an m x e table with uniform rows of fixed sums: exact DP over
        rows with state = column plus-counts when (m+1)^e <= 2e6, else a Gaussian entropy (flag False)"""
        m, e = len(rowsums), self.e
        if (m + 1) ** e <= 2_000_000:
            base = m + 1; pw = base ** np.arange(e, dtype=np.int64)
            prob = {0: 1.0}
            for s in rowsums:
                n = (int(s) + e) // 2
                subs = [sum(int(pw[j]) for j in c_) for c_ in it.combinations(range(e), n)]
                w = 1.0 / len(subs); new = {}
                for st, pr in prob.items():
                    for d in subs:
                        k = st + d; new[k] = new.get(k, 0.0) + pr * w
                prob = new
            pv = np.array(list(prob.values()))
            return float(-(pv * np.log2(pv)).sum()), True
        # Gaussian entropy of the sum of m independent rows (d = e - 1 dims, lattice 2 Z^(e-1))
        var = sum(1 - (s / e) ** 2 for s in rowsums)
        Sig = var * (np.eye(e) - np.ones((e, e)) / e) * e / (e - 1)
        ev = np.linalg.eigvalsh(Sig)[1:]
        H = 0.5 * (e - 1) * math.log2(2 * math.pi * math.e) + 0.5 * np.log2(ev).sum() - (e - 1) * 1.0
        return float(H), False

    # ---------------------------------------------------------------- 523-style Z2 split h_D (pairs of cells q, q ^ 1)
    def hD_pair(self, s0, s1, _cache={}):
        key = (self.e, s0, s1)
        if key in _cache: return _cache[key]
        e = self.e
        n0, n1 = (s0 + e) // 2, (s1 + e) // 2
        R0 = list(it.combinations(range(e), n0)); R1 = list(it.combinations(range(e), n1))
        cnt = {}
        for A in R0:
            a = [1 if j in A else -1 for j in range(e)]
            for B in R1:
                k = tuple(a[j] + (1 if j in B else -1) for j in range(e)); cnt[k] = cnt.get(k, 0) + 1
        tot = len(R0) * len(R1)
        HS = -sum(v / tot * math.log2(v / tot) for v in cnt.values())
        _cache[key] = math.log2(tot) - HS
        return _cache[key]


def chain_numbers(host, E, S, N, rng, real_target=False):
    c = host.sample_y_fibre(S, N, rng)
    x = host.assemble(np.broadcast_to(E, (N,) + E.shape), c)
    V = host.class_paf(x).sum(1)
    if real_target:
        ga = QG.gauss_bits(V, host.target)
        gw = QG.gauss_bits(V @ host.wmap.T, host.wmap @ host.target)
        gs = QG.gauss_bits(V @ host.smap.T, host.smap @ host.target)
    else:
        ga = gauss_bits_mean(V); gw = gauss_bits_mean(V @ host.wmap.T); gs = gauss_bits_mean(V @ host.smap.T)
    lf = host.log2_fibre_y(S)
    Hw = 0.0; exact = True
    for i in range(4):
        for q in range(4):
            h, ex = host.H_w_cell([int(S[i, q * host.m + u]) for u in range(host.m)])
            Hw += h; exact &= ex
    hD = sum(host.hD_pair(int(S[i, q * host.m + u]), int(S[i, (q ^ 1) * host.m + u]))
             for i in range(4) for q in (0, 2) for u in range(host.m))
    Call, Cw, Cs = float(ga['bits']), float(gw['bits']), float(gs['bits'])
    out = dict(n=host.n, m=host.m, p=host.p, a=host.a, f=host.f, e=host.e, neg_in_S=host.neg_in_S,
               orbits_per_block=host.nc * (1 + host.e), lag_classes=host.L,
               log2_fibre_y=round(lf, 2), C_all=Call, rank_all=ga['rank'], C_w=Cw, rank_w=gw['rank'],
               H_w=round(Hw, 2), H_w_exact=exact,
               log2E_cores_per_y=round(lf - Call, 2), log2E_w_per_y=round(Hw - Cw, 2),
               log2_fibre_yw=round(lf - Hw, 2), C_mixed=round(Call - Cw, 2), rank_mixed=ga['rank'] - gw['rank'],
               D_def_Zm=round(Call - Cw - (lf - Hw), 2),
               Z2_C_S=Cs, Z2_h_D=round(hD, 2), Z2_D_def=round(Call - Cs - hD, 2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('m', type=int); ap.add_argument('p', type=int); ap.add_argument('a', type=int)
    ap.add_argument('--y', default=None, help='certificate (863 format) for a real y')
    ap.add_argument('--N', type=int, default=8000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--reps', type=int, default=1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    host = Host(a.m, a.p, a.a)
    rng = np.random.default_rng(a.seed)
    for r in range(a.reps):
        if a.y:
            d = json.load(open(a.y)); Y = np.array(d['y'], np.int64)
            E = np.where(np.mod(Y, 2 * host.f) == 1, 1, -1); S = (Y - E) // host.f
            info = dict(y=a.y); real = True
        else:
            E, S, info = host.synthetic_y(rng); real = False
        res = chain_numbers(host, E, S, a.N, rng, real_target=real)
        res.update(info); res.update(N=a.N, seed=a.seed, rep=r, read_at='target' if real else 'mean')
        s = json.dumps(res); print(s, flush=True)
        if a.out:
            with open(a.out, 'a') as fo: fo.write(s + '\n')


if __name__ == '__main__':
    main()
