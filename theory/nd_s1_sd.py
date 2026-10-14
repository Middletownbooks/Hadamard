r"""Exact compressed stage 1 of a multiplier slice by the alpha-congruence + S/D split (generalises
theory/nd_sym_stage1.py to any slice compressed onto Z_F with F = 2n, n odd) [P for the reduction, M for runs].

Slice: G = Z_m, x_k constant on the orbits of M = <mu>.  Compress along the Hall factor Z_K (n_F = F = m/K):
y_k(b) = sum_{g = b mod F} x_k(g) is M-invariant on Z_F, its value on an M-orbit ranges over the signed sums of the
sizes of the M_b-orbits that partition the fibre (theory/slice_compress_cpsat.py), and
    sum_k PAF(y_k)(0) = 4m - 4(K - 1),  sum_k PAF(y_k)(h) = -4K (h != 0),  rows (2, 0, 0, 0).
Reduction [P]:
 1. If every letter is = +-1 (mod 2 f0) with f0 >= 2, write y = alpha + 2 f0 t (alpha = +-1 unique).  Then
    PAF(y) = PAF(alpha) + 2 f0 (CC(alpha, t) + CC(t, alpha)) + 4 f0^2 PAF(t) and the middle sum is even, so
    sum_k PAF(alpha_k) = T (mod 4 f0) at every lag, and rows of alpha = R_k (mod 2 f0).
 2. Z_F = Z_2 x Z_n (CRT; every multiplier is odd so M acts trivially on Z_2).  With e0(b), e1(b) the lifts of b with
    even / odd parity, sigma(b) = (y(e0 b) + y(e1 b))/2 and delta(b) = (y(e0 b) - y(e1 b))/2 are M-invariant on Z_n and
        sum_k PAF_n(sigma_k) = ((1 + (n-1) v)/n, (1 - v)/n, ...),  rows sum_b sigma_k(b) = R_k / 2,
        sum_k PAF_n(delta_k) = (v, 0, ..., 0),                       v = m + 1.
 3. Given alpha the S-system involves only u = t(e0) + t(e1) and the D-system only w = t(e0) - t(e1).
Search (exhaustive): alpha-quadruples (block 0 has row +2; blocks 1..3 row 0, up to negation, in sorted order) ->
2+2 join of the S-system over u -> for each S-solution a 2+2 join of the D-system over w -> verify on Z_F.
No solution with "complete": true is an exact proof that the compressed stage 1, hence the slice, is empty.
usage: python3 theory/nd_s1_sd.py m mu K [--plant N] [--max-sol N] [--out FILE] [--max-quads N]
       python3 theory/nd_s1_sd.py --explicit C.json [...]   (compressed structure from theory/nd_fq_split.py --dump)"""
import sys, os, json, time, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def mgroup(m, mu):
    G = [1]; x = mu % m
    while x != 1: G.append(x); x = x * mu % m
    return G

def structure(m, mu, K):
    """orbits of M on Z_F (F = m/K) and the letter set of each (as in slice_compress_cpsat.structure)"""
    F = m // K; assert m % K == 0 and math.gcd(K, F) == 1
    Mg = mgroup(m, mu)
    orb = -np.ones(m, np.int64); sizes = []
    for g in range(m):
        if orb[g] >= 0: continue
        i = len(sizes); c = 0
        for a in Mg:
            x = g * a % m
            if orb[x] < 0: orb[x] = i; c += 1
        sizes.append(c)
    Fo = -np.ones(F, np.int64); reps = []
    for b in range(F):
        if Fo[b] >= 0: continue
        j = len(reps)
        for a in Mg: Fo[b * a % F] = j
        reps.append(b)
    letters = []
    for b in reps:
        parts = {}
        for g in range(b, m, F): parts[orb[g]] = parts.get(orb[g], 0) + 1
        L = {0}
        for p in parts.values(): L = {x + e * p for x in L for e in (1, -1)}
        letters.append(sorted(L))
    return F, Mg, Fo, reps, letters

def paf_lags(Y, lags):
    """PAF at the given lags (last axis is the group Z_N)"""
    return np.stack([(Y * np.roll(Y, -h, axis=-1)).sum(-1) for h in lags], -1)

def lag_reps(N, Mg):
    """representatives of the <M, -1>-orbits on Z_N (0 first)"""
    seen = set(); reps = []
    for h in range(N):
        if h in seen: continue
        reps.append(h)
        for a in Mg:
            seen.add(h * a % N); seen.add((-h * a) % N)
    return reps

def pair_join(L0, L1, L2, L3, T):
    """all (i, j, k, l) with L0[i] + L1[j] + L2[k] + L3[l] = T, exhaustive (sorted int64 keys on as many leading
    coordinates as fit in 62 bits, the rest checked exactly)"""
    if min(len(L0), len(L1), len(L2), len(L3)) == 0: return []
    d = L0.shape[1]
    lo_l = L0.min(0) + L1.min(0); hi_l = L0.max(0) + L1.max(0)
    lo_r = T - L2.max(0) - L3.max(0); hi_r = T - L2.min(0) - L3.min(0)
    lo = np.minimum(lo_l, lo_r); hi = np.maximum(hi_l, hi_r)
    bits = [int(hi[c] - lo[c]).bit_length() + 1 for c in range(d)]
    cols = []; tot = 0
    for c in range(d):
        if tot + bits[c] > 62: break
        cols.append(c); tot += bits[c]
    def key(X):
        k = np.zeros(X.shape[:-1], np.int64)
        for c in cols: k = (k << bits[c]) | (X[..., c] - lo[c])
        return k
    left = key(L0[:, None, :] + L1[None, :, :]).ravel()
    o = np.argsort(left, kind='stable'); ls = left[o]
    out = []; n1 = L1.shape[0]; n3 = L3.shape[0]; full = len(cols) == d
    step = max(1, (1 << 24) // max(1, n3))                      # right side in chunks of rows of L2
    for k0 in range(0, L2.shape[0], step):
        right = key(T[None, None, :] - L2[k0:k0 + step, None, :] - L3[None, :, :]).ravel()
        a = np.searchsorted(ls, right, 'left'); b = np.searchsorted(ls, right, 'right')
        for r in np.nonzero(b > a)[0]:
            k, l = divmod(int(r), n3); k += k0
            for p in o[a[r]:b[r]]:
                i, j = divmod(int(p), n1)
                if full or ((L0[i] + L1[j] + L2[k] + L3[l]) == T).all(): out.append((i, j, k, l))
    return out

class Problem:
    def __init__(self, m, mu, K, explicit=None):
        """explicit: dict(F, MF (multipliers on Z_F), Fo, freps, letters, T, m, K) for non-cyclic groups
        (theory/nd_fq_split.py --dump); then m = |G| and K = |compressed factor|"""
        self.m, self.mu, self.K = m, mu, K
        if explicit is None:
            F, Mg, Fo, reps, letters = structure(m, mu, K)
        else:
            F = explicit['F']; Mg = explicit['MF']; Fo = np.array(explicit['Fo']); reps = explicit['freps']
            letters = explicit['letters']
        assert F % 2 == 0 and (F // 2) % 2 == 1, "needs F = 2n with n odd"
        self.F, self.Mg, self.Fo, self.reps, self.letters = F, Mg, Fo, reps, letters
        n = F // 2; self.n = n; self.v = m + 1
        # f0: the largest f >= 2 with every letter = +-1 (mod 2f)
        allL = sorted({x for L in letters for x in L})
        assert all(x % 2 for x in allL), "letters must be odd"
        f0 = 1
        for f in range(2, max(abs(x) for x in allL) + 2):
            if all((x - 1) % (2 * f) == 0 or (x + 1) % (2 * f) == 0 for x in allL): f0 = f
        assert f0 >= 2, "no alpha-structure (f0 < 2)"
        self.f0 = f0
        # allowed t per (Z_F orbit j, alpha)
        self.tset = [{1: sorted((x - 1) // (2 * f0) for x in L if (x - 1) % (2 * f0) == 0),
                      -1: sorted((x + 1) // (2 * f0) for x in L if (x + 1) % (2 * f0) == 0)} for L in letters]
        # CRT lifts
        self.e0 = np.array([next(g for g in range(F) if g % 2 == 0 and g % n == b) for b in range(n)])
        self.e1 = np.array([next(g for g in range(F) if g % 2 == 1 and g % n == b) for b in range(n)])
        Mn = sorted({a % n for a in Mg}) if n > 1 else [0]
        self.Mn = Mn
        On = -np.ones(n, np.int64); self.nreps = []
        for b in range(n):
            if On[b] >= 0: continue
            j = len(self.nreps)
            for a in Mn: On[b * a % n] = j
            self.nreps.append(b)
        self.On = On; self.nw = np.bincount(On)
        self.nlags = lag_reps(n, Mn) if n > 1 else [0]
        self.Flags = lag_reps(F, [a % F for a in Mg])
        v = self.v
        self.T = np.array([4 * m - 4 * (K - 1)] + [-4 * K] * (len(self.Flags) - 1), np.int64)
        if explicit is not None: assert self.T[0] == explicit['T'][0] and self.T[1] == explicit['T'][1]
        assert (1 + (n - 1) * v) % n == 0 and (1 - v) % n == 0
        self.TS = np.array([(1 + (n - 1) * v) // n] + [(1 - v) // n] * (len(self.nlags) - 1), np.int64)
        self.TD = np.array([v] + [0] * (len(self.nlags) - 1), np.int64)

    # ---------------- alpha ----------------
    def alphas(self):
        """all M-invariant +-1 sequences on Z_F"""
        no = len(self.reps)
        if no > 22: raise RuntimeError(f"too many orbits for the alpha enumeration ({no})")
        A = np.empty((1 << no, self.F), np.int64)
        for msk in range(1 << no):
            vals = np.array([1 if (msk >> j) & 1 else -1 for j in range(no)])
            A[msk] = vals[self.Fo]
        return A

    def alpha_quads(self, T=None, rows=(2, 0, 0, 0), max_quads=10 ** 7):
        T = self.T if T is None else T
        f0 = self.f0; mod = 4 * f0; rmod = 2 * f0
        A = self.alphas(); P = paf_lags(A, self.Flags[1:]); R = A.sum(1)
        # only alphas that admit a t at every orbit (alpha value must have a letter)
        okA = np.array([all(len(self.tset[j][int(a[self.reps[j]])]) for j in range(len(self.reps))) for a in A])
        I2 = np.nonzero(((R - rows[0]) % rmod == 0) & okA)[0]
        I0 = np.nonzero((R % rmod == 0) & okA)[0]
        # canonical up to negation for row-0 blocks
        idx = {tuple(a): i for i, a in enumerate(A.tolist())}
        neg = np.array([idx[tuple((-a).tolist())] for a in A])
        I0c = I0[np.minimum(I0, neg[I0]) == I0]
        if len(I2) * len(I0c) > 3 * 10 ** 7 or len(I0c) ** 2 > 6 * 10 ** 7:
            raise RuntimeError(f"alpha lists too large ({len(I2)} x {len(I0c)})")
        nl = P.shape[1]
        base = np.array([mod ** i for i in range(nl)], dtype=object)
        if mod ** nl < 2 ** 62:
            base = np.array([mod ** i for i in range(nl)], np.int64); Kf = lambda M: (M % mod) @ base
        else:
            w = np.random.default_rng(7).integers(1, 2 ** 61, nl, dtype=np.int64)
            Kf = lambda M: ((M % mod) * w).sum(-1)          # hash; exact check below
        Tm = T[1:] % mod
        left = Kf(P[I2][:, None, :] + P[I0c][None, :, :]).ravel()
        pr = np.array([(a, b) for a in range(len(I0c)) for b in range(a, len(I0c))], np.int64).reshape(-1, 2)
        if len(pr) == 0 or len(left) == 0: return A, [], dict(row2=int(len(I2)), row0_classes=int(len(I0c)))
        rk = Kf(Tm[None, :] - P[I0c[pr[:, 0]]] - P[I0c[pr[:, 1]]])
        o = np.argsort(left); ls = left[o]
        lo = np.searchsorted(ls, rk, 'left'); hi = np.searchsorted(ls, rk, 'right')
        quads = set(); nI = len(I0c)
        for r in np.nonzero(hi > lo)[0]:
            a, b = pr[r]
            for p in o[lo[r]:hi[r]]:
                i0 = int(I2[p // nI]); i1 = int(I0c[p % nI])
                if ((P[i0] + P[i1] + P[I0c[a]] + P[I0c[b]] - T[1:]) % mod).any(): continue
                quads.add((i0,) + tuple(sorted((i1, int(I0c[a]), int(I0c[b])))))
                if len(quads) > max_quads: raise RuntimeError("too many alpha-quadruples")
        return A, sorted(quads), dict(row2=int(len(I2)), row0_classes=int(len(I0c)))

    # ---------------- per block ----------------
    def block(self, alpha, R, TS0, TD0):
        n = self.n; f0 = self.f0; nr = self.nreps
        aS = np.array([(alpha[self.e0[b]] + alpha[self.e1[b]]) // 2 for b in nr])
        aD = np.array([(alpha[self.e0[b]] - alpha[self.e1[b]]) // 2 for b in nr])
        T0 = [self.tset[self.Fo[self.e0[b]]][int(alpha[self.e0[b]])] for b in nr]
        T1 = [self.tset[self.Fo[self.e1[b]]][int(alpha[self.e1[b]])] for b in nr]
        pairs = [[(t0, t1) for t0 in a for t1 in c] for a, c in zip(T0, T1)]
        uopts = [sorted({t0 + t1 for t0, t1 in p}) for p in pairs]
        wts = self.nw
        # DFS over u with energy bound
        U = []
        def rec(j, e, cur):
            if j == len(nr):
                U.append(list(cur)); return
            for u in uopts[j]:
                sg = aS[j] + f0 * u; e2 = e + wts[j] * sg * sg
                if e2 <= TS0: cur.append(u); rec(j + 1, e2, cur); cur.pop()
        rec(0, 0, [])
        U = np.array(U, np.int64).reshape(-1, len(nr))
        sig = aS[None, :] + f0 * U
        ok = (wts[None, :] * sig).sum(1) * 2 == R
        U = U[ok]; sig = sig[ok]
        sigma = sig[:, self.On]
        Ssig = paf_lags(sigma, self.nlags) if len(U) else np.zeros((0, len(self.nlags)), np.int64)
        return dict(alpha=alpha, R=R, aS=aS, aD=aD, pairs=pairs, U=U, Ssig=Ssig, TD0=TD0)

    def dlist(self, B, u):
        f0 = self.f0; nr = self.nreps; wts = self.nw
        wopts = [sorted({t0 - t1 for t0, t1 in B['pairs'][j] if t0 + t1 == u[j]}) for j in range(len(nr))]
        W = []
        def rec(j, e, cur):
            if j == len(nr):
                W.append(list(cur)); return
            for w in wopts[j]:
                dl = B['aD'][j] + f0 * w; e2 = e + wts[j] * dl * dl
                if e2 <= B['TD0']: cur.append(w); rec(j + 1, e2, cur); cur.pop()
        rec(0, 0, [])
        W = np.array(W, np.int64).reshape(-1, len(nr))
        delta = (B['aD'][None, :] + f0 * W)[:, self.On]
        return W, (paf_lags(delta, self.nlags) if len(W) else np.zeros((0, len(self.nlags)), np.int64))

    def y_of(self, B, u, w):
        t = np.zeros(self.F, np.int64)
        for j, b in enumerate(self.nreps):
            t0 = (u[j] + w[j]) // 2; t1 = (u[j] - w[j]) // 2
            for a in self.Mn:
                bb = b * a % self.n
                t[self.e0[bb]] = t0; t[self.e1[bb]] = t1
        return B['alpha'] + 2 * self.f0 * t

    def verify(self, Y, T=None):
        T = self.T if T is None else T
        Y = np.array(Y)
        P = paf_lags(Y, self.Flags).sum(0)
        inv = all((Y[:, (np.arange(self.F) * a) % self.F] == Y).all() for a in self.Mg)
        if not inv: print("WARNING: not M-invariant", flush=True)
        lett = all(int(Y[k, self.reps[j]]) in self.letters[j] for k in range(4) for j in range(len(self.reps)))
        full = np.stack([(Y * np.roll(Y, -h, axis=1)).sum() for h in range(self.F)])
        fullT = np.array([T[0]] + [None] * 0)
        return bool((P == T).all() and inv and lett)

    def solve(self, max_sol=10 ** 9, log=print, T=None, TS=None, TD=None, rows=(2, 0, 0, 0), quads=None, max_quads=10 ** 7):
        T = self.T if T is None else T; TS = self.TS if TS is None else TS; TD = self.TD if TD is None else TD
        t0 = time.time()
        A, Q, info = self.alpha_quads(T, rows, max_quads) if quads is None else (self.alphas(), quads, {})
        log(json.dumps(dict(stage="alpha", m=self.m, mu=self.mu, K=self.K, F=self.F, f0=self.f0, alpha_quads=len(Q),
                            **info, secs=round(time.time() - t0, 1))))
        sols = []; nS = 0; cache = {}
        for quad in Q:
            blocks = []
            for k, ai in enumerate(quad):
                key = (ai, rows[k])
                if key not in cache: cache[key] = self.block(A[ai], rows[k], int(TS[0]), int(TD[0]))
                blocks.append(cache[key])
            Sh = pair_join(*[b['Ssig'] for b in blocks], TS)
            nS += len(Sh)
            for (i, j, k, l) in Sh:
                us = [blocks[0]['U'][i], blocks[1]['U'][j], blocks[2]['U'][k], blocks[3]['U'][l]]
                Ds = [self.dlist(b, u) for b, u in zip(blocks, us)]
                for (a, b_, c, d) in pair_join(*[x[1] for x in Ds], TD):
                    ws = [Ds[0][0][a], Ds[1][0][b_], Ds[2][0][c], Ds[3][0][d]]
                    Y = [self.y_of(B, u, w).tolist() for B, u, w in zip(blocks, us, ws)]
                    sols.append(dict(alpha_quad=[int(x) for x in quad], y=Y, verified=self.verify(Y, T)))
                    if len(sols) >= max_sol: break
                if len(sols) >= max_sol: break
            if len(sols) >= max_sol: break
        out = dict(stage="done", m=self.m, v=self.m + 1, mu=self.mu, K=self.K, F=self.F, f0=self.f0, alpha_quads=len(Q),
                   S_solutions=nS, stage1_solutions=len(sols), all_verified=all(x['verified'] for x in sols),
                   complete=len(sols) < max_sol, secs=round(time.time() - t0, 1))
        log(json.dumps(out))
        return out, sols

    def planted_real(self, seed, log=print):
        """plant drawn uniformly at the slice's exact lag-0 energy and rows (theory/nd_compress_est.quad_sampler);
        targets = its own PAF sums (lag 0 = the real one); the exhaustive search must return it"""
        import nd_compress_est as CE
        C = CE.compressed_cyclic(self.m, self.mu, self.K) if not hasattr(self, 'C_explicit') else self.C_explicit
        Y = CE.quad_sampler(C, 500 + seed)()
        T = paf_lags(Y, self.Flags).sum(0)
        sig = (Y[:, self.e0] + Y[:, self.e1]) // 2; dl = (Y[:, self.e0] - Y[:, self.e1]) // 2
        TS = paf_lags(sig, self.nlags).sum(0); TD = paf_lags(dl, self.nlags).sum(0)
        out, sols = self.solve(log=lambda x: None, T=T, TS=TS, TD=TD)
        V = [{tuple(Y[0])}] + [{tuple(Y[k]), tuple(-Y[k])} for k in (1, 2, 3)]
        hit = any(tuple(x['y'][0]) in V[0] and any(all(tuple(x['y'][1 + i]) in V[1 + p[i]] for i in range(3))
                                                   for p in itertools.permutations(range(3))) for x in sols)
        res = dict(plant_seed=seed, alpha_quads=out['alpha_quads'], S_solutions=out['S_solutions'],
                   solutions=len(sols), planted_found=hit, secs=out['secs'])
        log(json.dumps(res)); return res

    def planted(self, seed, log=print):
        """random M-invariant quadruple with letters from the slice (small letters favoured), rows (2, 0, 0, 0) by
        rejection; targets = its own sums; the exhaustive search must return it (block order / negation aside)"""
        rng = np.random.default_rng(seed)
        def rnd(R):
            for _ in range(200000):
                vals = []
                for j, L in enumerate(self.letters):
                    small = [x for x in L if abs(x) <= sorted(abs(z) for z in L)[min(3, len(L) - 1)]]
                    vals.append(int(rng.choice(small if rng.random() < 0.85 else L)))
                y = np.array(vals)[self.Fo]
                if y.sum() == R: return y
            return None
        Y = [rnd(2)] + [rnd(0) for _ in range(3)]
        if any(y is None for y in Y): log(json.dumps(dict(plant_seed=seed, note="no plant"))); return None
        Y = np.array(Y)
        T = paf_lags(Y, self.Flags).sum(0)
        sig = (Y[:, self.e0] + Y[:, self.e1]) // 2; dl = (Y[:, self.e0] - Y[:, self.e1]) // 2
        TS = paf_lags(sig, self.nlags).sum(0); TD = paf_lags(dl, self.nlags).sum(0)
        out, sols = self.solve(log=lambda x: None, T=T, TS=TS, TD=TD)
        V = [{tuple(Y[0])}] + [{tuple(Y[k]), tuple(-Y[k])} for k in (1, 2, 3)]
        hit = any(tuple(x['y'][0]) in V[0] and any(all(tuple(x['y'][1 + i]) in V[1 + p[i]] for i in range(3))
                                                   for p in itertools.permutations(range(3))) for x in sols)
        res = dict(plant_seed=seed, alpha_quads=out['alpha_quads'], S_solutions=out['S_solutions'],
                   solutions=len(sols), planted_found=hit, secs=out['secs'])
        log(json.dumps(res)); return res

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--plant": None, "--plant-real": None, "--max-sol": str(10 ** 9), "--out": "", "--max-quads": str(10 ** 7), "--explicit": ""}
    pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    if opt["--explicit"]:
        C = json.load(open(opt["--explicit"]))
        UF = sorted({u % C['F'] for (u, e) in C['gens']})
        MF = [1]
        while True:
            new = sorted({x * u % C['F'] for x in MF for u in UF} | set(MF))
            if new == MF: break
            MF = new
        C['MF'] = MF; m, mu, K = C['m'], 0, C['q']
        P = Problem(m, mu, K, explicit=C); P.C_explicit = C
    else:
        m, mu, K = pos
        P = Problem(m, mu, K)
    print(json.dumps(dict(m=m, mu=mu, K=K, F=P.F, n=P.n, f0=P.f0, Z_F_orbits=len(P.reps), Z_n_orbits=len(P.nreps),
                          letters=[P.letters[j] for j in range(len(P.reps))][:3], lag_classes=len(P.Flags) - 1)), flush=True)
    if opt["--plant-real"] is not None:
        r = [P.planted_real(sd, log=lambda x: print(x, flush=True)) for sd in range(int(opt["--plant-real"]))]
        print(json.dumps(dict(planted_real=len(r), found=sum(x['planted_found'] for x in r))), flush=True)
        sys.exit()
    if opt["--plant"] is not None:
        r = [P.planted(sd, log=lambda x: print(x, flush=True)) for sd in range(int(opt["--plant"]))]
        r = [x for x in r if x]
        print(json.dumps(dict(planted=len(r), found=sum(x['planted_found'] for x in r))), flush=True)
        sys.exit()
    out, sols = P.solve(int(opt["--max-sol"]), log=lambda x: print(x, flush=True), max_quads=int(opt["--max-quads"]))
    if opt["--out"]:
        json.dump(dict(result=out, solutions=sols[:1000]), open(opt["--out"], "w"))
