r"""p = 3: the Paley(3)-coset ansatz is vacuous, and the right stage split is y | (Eisenstein Z), not (alpha, w) | u.

[P] For p = 3 the squares mod 3 are {1}, so EVERY +-1 triple (x0, x1, x2) on a Z_3-coset is alpha*delta + u*chi + w*eps
(x1 = x2: constant, w = x1; x1 = -x2: twisted, u = x1).  Write, on each coset g' of G = G' x Z_3,
    y(g') = x0 + x1 + x2 in {+-1, +-3},      Z(g') = x0 + x1 w + x2 w^2 in Z[w]   (w = exp(2 pi i / 3)).
|y| = 3 <=> the coset is constant <=> Z = 0;  |y| = 1 <=> exactly one minority position z* (x(z*) = -y) and
Z = -2 y w^{z*}.  So x <-> (y, z*) bijectively, and (Z_3-Fourier transform of the b = 1 target)
    S1:  sum_k PAF(y_k) = T_J   (4m - 8 at 0, -12 elsewhere)          [characters trivial on Z_3]
    S2:  sum_k PAF_C(Z'_k) = v delta,   Z'_k = y_k w^{z*_k} on |y_k| = 1, 0 on |y_k| = 3   [the other characters].
The repo's staged pipeline enumerated (gadget J/J', alpha) in stage 1 and u in stage 2.  For p = 3 a J' coset
(alpha, -alpha, -alpha) and a twisted coset (-alpha, u, -u) have the SAME y, so the bit "z* = 0 or not" was fixed in
stage 1, where S1 cannot see it: stage 1 was inflated by ~2^{n1} and stage 2 restricted to a binary choice.
Correct split: stage 1 = y only (distinct y's), stage 2 = z* in Z_3 (ternary) on every |y| = 1 orbit.

Slices: any multiplier group M <= (Z/m)^*, m = 3n', 3 !| n'.  Elements t = 1 mod 3 act on Z_n' only; t = -1 mod 3
also flip z -> -z, i.e. Z'(t g') = conj Z'(g') (z* -> -z*), and on M'-orbits fixed by such t, z* = 0 is forced.

This module: estimates [E] of stage 1 (# distinct y with S1, row sums (+-2,0,0,0)) and stage 2 (per y: 3^{n_free}
times the S2 Gaussian density), list sizes per block, plus the exact structure used by the join (p571_s2tern.py).
usage: python3 theory/p571_ternary.py m gen1 [gen2 ...] [--N 20000] [--ny 4] [--seed 1]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_slicelib import mgroup, orbits, lag_reps, gauss_logp

W = np.exp(2j * np.pi / 3)

class TernSlice:
    """structure of the M-slice of Z_m = Z_n x Z_3 (index i <-> (i mod n, i mod 3)) in (y, z*) coordinates"""
    def __init__(self, m, gens):
        assert m % 3 == 0 and (m // 3) % 3 != 0
        self.m = m; self.n = n = m // 3; self.v = m + 1
        self.G = mgroup(m, gens)
        self.Gp = sorted({t % n for t in self.G})                    # M' on Z_n
        self.G1p = sorted({t % n for t in self.G if t % 3 == 1})     # image of M_1 = {t = 1 mod 3}
        self.flip = any(t % 3 == 2 for t in self.G)
        self.yorbs = orbits(n, self.Gp)                              # y is constant on M'-orbits
        # stage-2 classes: an M'-orbit O' splits into M_1'-orbits; if two, they are swapped by the flips
        self.cls = []                                                # (rep orbit a, partner orbit b or None, free)
        for o in self.yorbs:
            sub = orbits_in(o, self.G1p, n)
            if len(sub) == 1:
                self.cls.append((sub[0], None, not self.flip))       # self-paired under a flip -> z* = 0 forced
            else:
                assert len(sub) == 2 and self.flip
                self.cls.append((sub[0], sub[1], True))
        self.reps1, _ = lag_reps(n, self.Gp)                         # S1 coordinates (orbits of <M', -1>)
        self.R = len(self.yorbs)

    def x_from(self, y, zs):
        """y: (4, R) orbit values in {+-1, +-3}; zs: (4, R) z* in {0,1,2} for the 'a' part (ignored where |y| = 3).
        Returns x (4, m) +-1 on Z_m (index i <-> (i mod n, i mod 3))."""
        n, m = self.n, self.m; X = np.zeros((4, n, 3), np.int64)
        for k in range(4):
            for j, (oa, ob, free) in enumerate(self.cls):
                yy = y[k, j]
                for part, sign in ((oa, 1), (ob, -1)):
                    if part is None: continue
                    if abs(yy) == 3: X[k, part, :] = yy // 3; continue
                    zz = zs[k, j] if free else 0                       # flip-fixed class: z* = 0 forced
                    z = zz % 3 if sign == 1 else (-zz) % 3
                    X[k, part, :] = yy; X[k, part, z] = -yy
        i = np.arange(m)
        return X[:, i % n, i % 3]

def orbits_in(o, G1p, n):
    seen, out = set(), []
    for g in o:
        if g in seen: continue
        s = sorted({g * t % n for t in G1p}); seen |= set(s); out.append(s)
    return out

# ---------------------------------------------------------------- stage 1: distinct y's
def y_count_tables(sizes, rows):
    """per block: dict (rowsum, J-count) -> # y-orbit labelings; y in {+-1, +-3} per orbit"""
    tabs = []
    for r in rows:
        D = {(0, 0): 1}
        for s in sizes:
            E = {}
            for (a, c), val in D.items():
                for yv in (1, -1, 3, -3):
                    key = (a + s * yv, c + (s if abs(yv) == 3 else 0))
                    E[key] = E.get(key, 0) + val
            D = E
        tabs.append({c: val for (a, c), val in D.items() if a == r})
    return tabs

class YSampler:
    """uniform y on orbits (values +-1, +-3) with sum s*y = r and J-count K: backward DP over (target, J-count)"""
    def __init__(self, sizes):
        self.sizes = list(sizes); R = len(sizes)
        cnt = [None] * (R + 1); cnt[R] = {(0, 0): 1.0}
        for i in range(R - 1, -1, -1):
            s = sizes[i]; E = {}
            for (t, c), val in cnt[i + 1].items():
                for yv in (1, -1, 3, -3):
                    key = (t + s * yv, c + (s if abs(yv) == 3 else 0))
                    E[key] = E.get(key, 0.0) + val
            cnt[i] = E
        self.cnt = cnt
    def sample(self, rng, r, K, N):
        if self.cnt[0].get((r, K), 0) == 0: return None
        out = np.zeros((N, len(self.sizes)), np.int64)
        for q in range(N):
            t, c = r, K
            for i, s in enumerate(self.sizes):
                opts = []; ws = []
                for yv in (1, -1, 3, -3):
                    w = self.cnt[i + 1].get((t - s * yv, c - (s if abs(yv) == 3 else 0)), 0.0)
                    if w > 0: opts.append(yv); ws.append(w)
                ws = np.array(ws); yv = opts[rng.choice(len(opts), p=ws / ws.sum())]
                out[q, i] = yv; t -= s * yv; c -= (s if abs(yv) == 3 else 0)
        return out

def stage1(S, N=20000, seed=1, rows=(2, 0, 0, 0)):
    """# distinct y with the norm count N+ = (n-1) J-cosets [exact], row sums [exact], S1 off 0 [Gaussian]"""
    n = S.n; sizes = [len(o) for o in S.yorbs]; Nplus = n - 1                 # 4n - K + 9 N+ ... : sum y^2 = T_J(0)
    assert 4 * n + 8 * Nplus == 4 * S.m - 8
    tabs = y_count_tables(sizes, rows)
    # distribution of (K_1..K_4) with sum = Nplus
    combos = {}
    for c0, v0 in tabs[0].items():
        for c1, v1 in tabs[1].items():
            if c0 + c1 > Nplus: continue
            for c2, v2 in tabs[2].items():
                c3 = Nplus - c0 - c1 - c2
                if c3 < 0 or c3 not in tabs[3]: continue
                combos[(c0, c1, c2, c3)] = v0 * v1 * v2 * tabs[3][c3]
    tot = sum(combos.values())
    if tot == 0: return dict(feasible=False)
    rng = np.random.default_rng(seed)
    keys = list(combos); pr = np.array([float(combos[k]) for k in keys]); pr /= pr.sum()
    pick = rng.choice(len(keys), size=N, p=pr)
    Y = np.zeros((N, 4, len(sizes)), np.int64); ys = YSampler(sizes)
    for idx in np.unique(pick):
        sel = np.nonzero(pick == idx)[0]; Ks = keys[idx]
        for k in range(4):
            Y[sel, k] = ys.sample(rng, rows[k], Ks[k], len(sel))
    Yf = np.zeros((N, 4, n), np.int64)
    for j, o in enumerate(S.yorbs): Yf[:, :, o] = Y[:, :, [j]]
    F = np.fft.rfft(Yf.astype(float), axis=2)
    P = np.rint(np.fft.irfft(np.abs(F) ** 2, n=n, axis=2)).astype(np.int64).sum(1)[:, S.reps1]
    cb, eqs = gauss_logp(P, np.full(len(S.reps1), -12))
    lE = math.log2(tot) + math.log2(8) - cb
    return dict(y_configs_log2=round(math.log2(tot) + 3, 1), C1=round(cb, 1), eqs1=eqs, stage1_log2E=round(lE, 1)), Y

# ---------------------------------------------------------------- stage 2: ternary z* given y
def gauss_cell(Smat, T, cell_log2):
    """-log2 Gaussian density at T plus the given lattice-cell log2 (no per-coordinate residue test: the residue class
    of the target is guaranteed by S1, see stage2)"""
    cols = [j for j in range(Smat.shape[1]) if Smat[:, j].std() > 0]
    S = Smat[:, cols].astype(float); Tt = np.asarray(T)[cols].astype(float)
    mu = S.mean(0); C = np.atleast_2d(np.cov(S.T)); w, V = np.linalg.eigh(C); kw = w > 1e-6 * w.max()
    pr = V.T @ (Tt - mu)
    lp = (cell_log2 * math.log(2) - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum()
          - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
    return -lp, int(kw.sum())

def stage2(S, y, N=20000, seed=2):
    """per-y stage-2 estimate [E]: free classes get uniform z* in Z_3 (z* = 0 on flip-fixed classes); S2 coordinates
    (a, b) of P(h) = sum_k PAF_C(Z'_k)(h) = a + b w at one lag per orbit of <M', -1>.
    Exact lattice [P]: Z' = y mod (1 - w), so P(h) = sum_k PAF(e_k)(h) mod (1 - w), e = y [|y| = 1]; and
    sum_k PAF(e_k) = sum_k PAF(y_k) = T_J = 0 mod 3 off 0 for every S1 solution y.  Hence (a, b) lies in the index-3
    lattice a + b = 0 mod 3 (b = 0 and a = 0 mod 3 where P(h) is forced real), the target 0 is in it for S1 solutions,
    and each lag contributes a cell of volume 3 (log2 3)."""
    rng = np.random.default_rng(seed); n = S.n
    free = np.zeros((4, S.R), bool)
    for k in range(4):
        for j, (oa, ob, fr) in enumerate(S.cls):
            free[k, j] = fr and abs(y[k, j]) == 1
    nfree = int(free.sum())
    Z = np.zeros((N, 4, n), np.complex128)
    for k in range(4):
        for j, (oa, ob, fr) in enumerate(S.cls):
            yy = y[k, j]
            if abs(yy) == 3: continue
            z = rng.integers(0, 3, N) if free[k, j] else np.zeros(N, np.int64)
            Z[:, k, oa] = (yy * W ** z)[:, None]
            if ob is not None: Z[:, k, ob] = (yy * W ** ((-z) % 3))[:, None]
    F = np.fft.fft(Z, axis=2)
    P = np.fft.ifft(np.abs(F) ** 2, axis=2).sum(1)
    reps = S.reps1
    b = np.rint(2 * P.imag / math.sqrt(3)); a = np.rint(P.real + b / 2)
    real_forced = [bool((b[:, h] == 0).all()) for h in reps]
    Smat = np.concatenate([a[:, reps], b[:, [h for h, rf in zip(reps, real_forced) if not rf]]], 1).astype(np.int64)
    nl = len(reps)
    cb, eqs = gauss_cell(Smat, np.zeros(Smat.shape[1]), nl * math.log2(3))
    per_block = [int(free[k].sum()) for k in range(4)]
    return dict(n_free=nfree, H2=round(nfree * math.log2(3), 1), C2=round(cb, 1), eqs2=eqs,
                lags=nl, complex_lags=nl - sum(real_forced),
                stage2_log2E_per_y=round(nfree * math.log2(3) - cb, 1),
                list_log2=[round(f * math.log2(3), 1) for f in per_block])

def main(m, gens, N=20000, ny=4, seed=1):
    S = TernSlice(m, gens)
    head = dict(m=m, v=S.v, gens=list(gens), order=len(S.G), flip=S.flip, y_orbits=S.R,
                free_classes=sum(1 for c in S.cls if c[2]), eqs1_lags=len(S.reps1))
    print(json.dumps(head), flush=True)
    r1 = stage1(S, N=N, seed=seed)
    if isinstance(r1, dict): print(json.dumps(r1)); return
    r1, Y = r1
    print(json.dumps(dict(stage='1', **r1)), flush=True)
    for q in range(ny):
        r2 = stage2(S, Y[q], N=N, seed=seed + 10 + q)
        r2['core_log2E_est'] = round(r1['stage1_log2E'] + r2['stage2_log2E_per_y'], 1)
        print(json.dumps(dict(stage='2', y=q, **r2)), flush=True)

if __name__ == "__main__" and sys.argv[1] != 'selftest':
    a = sys.argv[1:]; opt = {"--N": "20000", "--ny": "4", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    main(pos[0], pos[1:], N=int(opt["--N"]), ny=int(opt["--ny"]), seed=int(opt["--seed"]))

# ---------------------------------------------------------------- exact self-test
def decompose(S, x):
    """x (4, m) +-1 -> y (4, n) and Z' (4, n) complex (Z' = y w^{z*} on |y| = 1, 0 on |y| = 3)"""
    n, m = S.n, S.m; i = np.arange(m)
    X = np.zeros((4, n, 3), np.int64); X[:, i % n, i % 3] = x
    y = X.sum(2); Z = (X * W ** np.arange(3)).sum(2)              # Z = -2 y w^{z*}
    return y, Z / -2

def check_split(S, x, b=1):
    """exact: full PAF of x on Z_m vs S1 (PAF of y) and S2 (complex PAF of Z')"""
    n, m = S.n, S.m
    P = np.rint(np.fft.ifft(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=1).real).astype(np.int64).sum(0)
    full_ok = bool(P[0] == 4 * m and (P[1:] == -4 * b).all())
    y, Z = decompose(S, x)
    P1 = np.rint(np.fft.ifft(np.abs(np.fft.fft(y, axis=1)) ** 2, axis=1).real).astype(np.int64).sum(0)
    P2 = np.fft.ifft(np.abs(np.fft.fft(Z, axis=1)) ** 2, axis=1).sum(0)
    TJ = np.full(n, -12 * b); TJ[0] = 4 * m - 8 * b
    s1 = bool((P1 == TJ).all()); s2 = bool(np.allclose(P2, (m + b) * (np.arange(n) == 0)))
    units = bool(np.allclose(np.abs(Z)[np.abs(y) == 1], 1) and np.allclose(Z[np.abs(y) == 3], 0))
    return dict(full=full_ok, S1=s1, S2=s2, units=units, equiv=(full_ok == (s1 and s2)))

def load_core(path):
    L = open(path).read().split('\n'); hd = L[0]
    kv = dict(t.split('=') for t in hd.split(';')[0].split() if '=' in t)
    n0, p0, m = int(kv['n']), int(kv['p']), int(kv['m'])
    X = np.array([[1 if c == '+' else -1 for c in l.split('=', 1)[1].strip()] for l in L[1:] if l[:2] in ('A=', 'B=', 'C=', 'D=')])
    i = np.arange(m)
    return X[:, (i % n0) * p0 + (i % p0)], m                       # cyclic Z_m by CRT

def selftest():
    out = []
    for f in ('core_v43_b1_paley.txt', 'core_v85_b1_paley.txt'):
        x, m = load_core(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', f))
        S = TernSlice(m, []); r = check_split(S, x); r['core'] = f
        y, Z = decompose(S, x)
        zs = np.zeros((4, S.R), np.int64); yo = np.zeros((4, S.R), np.int64)
        for k in range(4):
            for j, (oa, ob, fr) in enumerate(S.cls):
                g = oa[0]; yo[k, j] = y[k, g]
                if abs(y[k, g]) == 1: zs[k, j] = int(np.argmin(np.abs(Z[k, g] - y[k, g] * W ** np.arange(3))))
        r['rebuild'] = bool((S.x_from(yo, zs) == x).all()); out.append(r)
    rng = np.random.default_rng(5); ok = True
    for (m, gens) in ((570, [217]), (570, [83]), (84, []), (570, [11, 49]), (390, [7])):
        S = TernSlice(m, gens)
        for _ in range(3):
            yo = rng.choice([1, -1, 3, -3], size=(4, S.R)); zs = rng.integers(0, 3, (4, S.R))
            x = S.x_from(yo, zs)
            inv = all(np.array_equal(x[:, (np.arange(m) * t) % m], x) for t in S.G)   # M-invariance
            y, Z = decompose(S, x)
            ok &= inv and bool(np.allclose(np.abs(Z)[np.abs(y) == 1], 1))
            P = np.fft.ifft(np.abs(np.fft.fft(x, axis=1)) ** 2, axis=1).real.sum(0)
            n = S.n; Pn = P.reshape(-1)
            # level split: sum over z of PAF(h', z) = PAF(y)(h'); sum_z w^{-z} PAF(h', z) = 4 * PAF_C(Z')(h')
            Ph = np.zeros((n, 3)); idx = np.arange(m); Ph[idx % n, idx % 3] = Pn
            P1 = np.fft.ifft(np.abs(np.fft.fft(y, axis=1)) ** 2, axis=1).real.sum(0)
            P2 = np.fft.ifft(np.abs(np.fft.fft(Z, axis=1)) ** 2, axis=1).sum(0)
            ok &= bool(np.allclose(Ph.sum(1), P1)) and bool(np.allclose((Ph * W ** (-np.arange(3))).sum(1), 4 * P2) or
                                                            np.allclose((Ph * W ** (np.arange(3))).sum(1), 4 * P2))
    out.append(dict(random_slices_invariant_and_level_split=bool(ok)))
    return out

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == 'selftest':
    for r in selftest(): print(json.dumps(r))
