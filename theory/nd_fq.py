r"""Generalised-Paley slices over a finite field F_q, q = p^2 (non-cyclic G = Z_F x F_q, additive group Z_p^2) [P, M, E].

b = 1 bordered GS core on the abelian group G = Z_F x (F_q, +) of order m = F q (v = m + 1).  The slice: Z_F fixed
pointwise, the index-s subgroup H of F_q^* acting by field multiplication; per Z_F-point g and block k the values are
alpha(g) at z = 0 and c_i(g) on the class C_i = gamma^i H (i = 0..s-1).  Stage 1 (the F_q-compression) is the same
problem as for a prime q: y = alpha + f sum_i c_i on Z_F, sum_k PAF(y_k) = (4m - 4(q-1), -4q, ...), rows (2,0,0,0).

Uniform cyclotomy [P, standard; checked numerically by `check`]: if -1 is a power of p modulo s (here s | p + 1), the
Gauss periods eta_i(a) = sum_{z in C_i} zeta_p^{Tr(a z)} take the value A on one class and B on the others, with
A - B = p (for the sign that makes them integers), B = (-1 - p)/s.  Hence for the additive character psi_a
    xhat(g; psi_a) = alpha(g) + B sigma(g) + p c_{i(a)}(g),     sigma = sum_i c_i = (y - alpha)/f,
and stage 2 is equivalent to s integer systems on Z_F, one per class i:
    sum_k PAF_{Z_F}(w_{k,i}) = 4v delta,   w_{k,i}(g) = alpha_k(g) + B sigma_k(g) + p c_{k,i}(g),
coupled only by sum_i c_{k,i}(g) = sigma_k(g).
This file: field arithmetic, the period check, a full-group PAF (3-D FFT on Z_F x Z_p x Z_p, no use of the
decomposition) to verify cores and the decomposition, and the stage-2 Gaussian given y (lattice-exact, as in
theory/nd_stage2_est.py).
usage: python3 theory/nd_fq.py check p s
       python3 theory/nd_fq.py s2est F p s Yfile [--idx i] [--N 20000]"""
import sys, os, json, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits

class Fq:
    """F_{p^2} = F_p[t]/(t^2 - r), r a non-residue; element a + b t <-> integer a + p b"""
    def __init__(self, p):
        self.p = p; q = p * p; self.q = q
        r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        self.r = r
        a = np.arange(q) % p; b = np.arange(q) // p
        A1, A2 = np.meshgrid(a, a, indexing='ij'); B1, B2 = np.meshgrid(b, b, indexing='ij')
        self.mul = (((A1 * A2 + r * B1 * B2) % p) + p * ((A1 * B2 + A2 * B1) % p)).astype(np.int64)
        self.add = (((A1 + A2) % p) + p * ((B1 + B2) % p)).astype(np.int64)
        self.neg = ((-a) % p) + p * ((-b) % p)
        # generator of F_q^*
        for g in range(2, q):
            x = 1; seen = 0
            for e in range(1, q - 1):
                x = self.mul[x, g]
                if x == 1: break
            else:
                self.gen = g; break
        self.log = -np.ones(q, np.int64); x = 1
        for e in range(q - 1): self.log[x] = e; x = self.mul[x, self.gen]

    def classes(self, s):
        """class index (log mod s) of every nonzero element; -1 at 0"""
        c = np.where(self.log >= 0, self.log % s, -1); return c

    def trace(self, z):
        """absolute trace F_{p^2} -> F_p: Tr(a + b t) = 2a"""
        return (2 * (z % self.p)) % self.p

    def periods(self, s):
        """eta_i(a) for a = 1: sum_{z in C_i} zeta_p^{Tr(z)} (complex)"""
        cl = self.classes(s); z = np.arange(self.q)
        w = np.exp(2j * np.pi * self.trace(z) / self.p)
        return np.array([w[cl == i].sum() for i in range(s)])

def check(p, s):
    K = Fq(p); eta = K.periods(s)
    vals = sorted(set(np.round(eta.real, 9)))
    imag = float(np.abs(eta.imag).max())
    out = dict(p=p, q=p * p, s=s, periods=[round(float(x), 6) for x in eta.real], distinct=vals, max_imag=imag,
               integral=bool(imag < 1e-9 and all(abs(x - round(x)) < 1e-9 for x in vals)))
    print(json.dumps(out)); return out

class Slice:
    """the GP slice on G = Z_F x F_q: configurations (alpha, c) -> x on Z_F x Z_p x Z_p"""
    def __init__(self, F, p, s):
        self.F, self.p, self.s = F, p, s; self.K = Fq(p); q = p * p; self.q = q
        self.f = (q - 1) // s; self.m = F * q; self.v = self.m + 1
        self.cl = self.K.classes(s)
        # lag orbits of <H, -1> on G: (h, z) ~ (h, u z) (u in H) ~ (-h, -z)
        H = [x for x in range(1, q) if self.K.log[x] % s == 0]
        seen = set(); reps = []
        for h in range(F):
            for z in range(q):
                if (h, z) in seen: continue
                reps.append((h, z))
                for u in H:
                    uz = int(self.K.mul[u, z]); seen.add((h, uz)); seen.add(((-h) % F, int(self.K.neg[uz])))
        self.lags = reps
        self.T = np.array([4 * self.m if (h, z) == (0, 0) else -4 for (h, z) in reps], np.int64)

    def x_of(self, alpha, C):
        """alpha: (F,), C: (F, s) -> x on (F, p, p) [index z = a + p b -> (a, b)]"""
        F, p = self.F, self.p
        X = np.empty((F, self.q), np.int64)
        X[:, 0] = alpha
        nz = np.arange(1, self.q); X[:, nz] = C[:, self.cl[nz]]
        return X.reshape(F, p, p, order='F') if False else X.reshape(F, p, p).transpose(0, 2, 1)

    def paf_full(self, Xs):
        """sum_k PAF on G of the blocks Xs (list of (F, p, p) arrays, axis 1 = a, axis 2 = b)"""
        acc = 0
        for X in Xs:
            fx = np.fft.fftn(X.astype(np.float64)); acc = acc + np.fft.ifftn(fx * fx.conj()).real
        return np.rint(acc).astype(np.int64)

    def lagvec(self, P):
        return np.array([P[h, z % self.p, z // self.p] for (h, z) in self.lags], np.int64)

    def patterns(self, sig):
        """sign vectors c in {+-1}^s with sum = sig"""
        s = self.s; npos = (sig + s) // 2
        out = []
        for pos in itertools.combinations(range(s), npos):
            c = -np.ones(s, np.int64); c[list(pos)] = 1; out.append(c)
        return np.array(out, np.int64).reshape(-1, s)

    def s2est(self, Y, N=20000, seed=1, chunk=500):
        rng = np.random.default_rng(seed); Y = np.array(Y, np.int64); f = self.f
        alpha = np.where((Y - 1) % (2 * f) == 0, 1, -1) if (2 * f) > 2 else None
        sig = (Y - alpha) // f
        pats = [[self.patterns(int(sig[k, g])) for g in range(self.F)] for k in range(4)]
        if any(len(P) == 0 for pk in pats for P in pk): return dict(note="coset sum not realisable")
        H2 = sum(math.log2(len(P)) for pk in pats for P in pk)
        V = np.zeros((N, len(self.lags)), np.int64)
        for n in range(N):
            Xs = []
            for k in range(4):
                C = np.stack([pats[k][g][rng.integers(len(pats[k][g]))] for g in range(self.F)])
                Xs.append(self.x_of(alpha[k], C))
            V[n] = self.lagvec(self.paf_full(Xs))
        gb = gauss_bits(V, self.T)
        return dict(H2=round(H2, 2), lag_classes=len(self.lags), **gb,
                    stage2_log2E=round(H2 - gb['bits'], 2) if gb.get('in_lattice') else -math.inf)

if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "check":
        check(int(a[1]), int(a[2])); sys.exit()
    if a[0] == "s2est":
        opt = {"--idx": "-1", "--N": "20000"}; pos = []; i = 1
        while i < len(a):
            if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
            else: pos.append(a[i]); i += 1
        F, p, s = map(int, pos[:3]); d = json.load(open(pos[3]))
        sols = d['solutions'] if 'solutions' in d else [d]
        S = Slice(F, p, s)
        idx = range(len(sols)) if int(opt["--idx"]) < 0 else [int(opt["--idx"])]
        for ix in idx:
            print(json.dumps(dict(F=F, p=p, s=s, idx=ix, **S.s2est(sols[ix]['y'], int(opt["--N"])))), flush=True)

def check_decomposition(F, p, s, trials=5, seed=0):
    """[M] for random slice configurations: the full-group DFT at (chi, psi_a) equals the Z_F-DFT of
    w_i = alpha + B sigma + p c_i with i = i*(a) (the class whose period is A); B = (-1 - p)/s"""
    S = Slice(F, p, s); K = S.K; q = S.q; rng = np.random.default_rng(seed)
    B = (-1 - p) // s; assert (-1 - p) % s == 0
    eta = K.periods(s); istar1 = int(np.argmax(eta.real))           # class with period A for a = 1
    worst = 0.0
    for t in range(trials):
        alpha = rng.choice([-1, 1], F); C = rng.choice([-1, 1], (F, s))
        X = S.x_of(alpha, C).astype(float)                           # (F, a, b)
        sigma = C.sum(1)
        for a in range(1, q):
            # psi_a(z) = zeta_p^{Tr(a z)}
            az = K.mul[a, np.arange(q)]; ph = np.exp(2j * np.pi * K.trace(az) / p)
            Xz = np.empty((F, q))
            zi = np.arange(q); Xz[:, zi] = X[:, zi % p, zi // p]
            xhat_g = Xz @ ph                                             # sum_z x(g, z) psi_a(z), per g
            i = (istar1 - int(K.log[a])) % s
            w = alpha + B * sigma + p * C[:, i]
            worst = max(worst, float(np.abs(xhat_g - w).max()))
    out = dict(F=F, p=p, s=s, B=B, A=B + p, trials=trials, max_abs_diff=worst, ok=worst < 1e-8)
    print(json.dumps(out)); return out

def class_system(F, p, s, Y):
    """the single class system Phi(c) = R of a uniform-cyclotomy slice given stage-1 y (theory/nd_fq.py docstring):
    w = beta + p c with beta = alpha + B sigma; Phi_k(c)(h) = p (CC(beta,c) + CC(c,beta))(h) + p^2 PAF(c)(h), h = 0..F/2;
    R = 4 v delta - sum_k PAF(beta_k).  Returns (alpha, sigma, beta, R)"""
    q = p * p; f = (q - 1) // s; B = (-1 - p) // s; v = F * q + 1
    Y = np.array(Y, np.int64)
    alpha = np.where((Y - 1) % (2 * f) == 0, 1, -1); sigma = (Y - alpha) // f
    assert ((Y - alpha) % f == 0).all() and (np.abs(sigma) <= s).all() and ((sigma + s) % 2 == 0).all()
    beta = alpha + B * sigma
    H = F // 2
    Pb = np.array([[int((beta[k] * np.roll(beta[k], -h)).sum()) for h in range(H + 1)] for k in range(4)])
    R = np.array([4 * v if h == 0 else 0 for h in range(H + 1)], np.int64) - Pb.sum(0)
    return alpha, sigma, beta, R

def phi_block(beta_k, C, p):
    """Phi_k for many c (rows of C, shape (n, F)) and one block's beta; returns (n, F/2 + 1)"""
    F = C.shape[1]; H = F // 2; out = np.zeros((C.shape[0], H + 1), np.int64)
    for h in range(H + 1):
        bc = (beta_k[None, :] * np.roll(C, -h, axis=1)).sum(1) + (C * np.roll(beta_k, -h)[None, :]).sum(1)
        out[:, h] = p * bc + p * p * (C * np.roll(C, -h, axis=1)).sum(1)
    return out

def class_bits(F, p, s, Y, N=40000, seed=1):
    """[E] constraint bits X of one class system and log2 |S_Phi| = 4F - X (c uniform in {+-1}^(4F)); coupling bits
    and the implied stage-2 yield H2 - s X (exact H2 and coupling)"""
    rng = np.random.default_rng(seed)
    alpha, sigma, beta, R = class_system(F, p, s, Y)
    V = np.zeros((N, F // 2 + 1), np.int64)
    for k in range(4):
        C = rng.choice([-1, 1], (N, F)); V += phi_block(beta[k], C, p)
    gb = gauss_bits(V, R)
    H2 = sum(math.log2(math.comb(s, (s + int(x)) // 2)) for x in sigma.ravel())
    coupling = s * 4 * F - H2
    X = gb['bits']
    return dict(X=round(X, 2), log2_S_Phi=round(4 * F - X, 2), rank=gb['rank'], in_lattice=gb['in_lattice'],
                H2=round(H2, 2), coupling_bits=round(coupling, 2), stage2_log2E_indep=round(H2 - s * X, 2),
                sigma_hist={int(a): int(b) for a, b in zip(*np.unique(np.abs(sigma), return_counts=True))})

def _ramanujan(d, F):
    out = np.zeros(F)
    for j in range(F):
        if F // math.gcd(j, F) == d: out += np.cos(2 * np.pi * j * np.arange(F) / F)
    return np.rint(out).astype(np.int64)

def levels_classes(F, p, s, Y, N=20000, seed=1, levels=None):
    """[E] stage-2 constraint bits by Z_F-character level for a uniform-cyclotomy slice, from the class systems:
    signature = for each class i the summed Z_F PAF of w_{.,i} (lags 0..F-1), projected by Ramanujan sums (all h);
    patterns drawn uniformly per (block, position) given y.  Also h_D = H(D | S) for the Z_2 split (F = 2n):
    S = (c(g) + c(g + n)), D = (c(g) - c(g + n)) per position pair, exact from the pattern lists."""
    rng = np.random.default_rng(seed); q = p * p; f = (q - 1) // s; B = (-1 - p) // s; v = F * q + 1
    alpha, sigma, beta, R = class_system(F, p, s, Y)
    pats = {}
    def P(sg):
        if sg not in pats:
            out = []
            for pos in itertools.combinations(range(s), (s + sg) // 2):
                c = -np.ones(s, np.int64); c[list(pos)] = 1; out.append(c)
            pats[sg] = np.array(out)
        return pats[sg]
    H2 = sum(math.log2(len(P(int(x)))) for x in sigma.ravel())
    # samples: C[n, k, g, i]
    sig = np.zeros((N, s, F), np.int64)
    for k in range(4):
        C = np.zeros((N, F, s), np.int64)
        for g in range(F):
            L = P(int(sigma[k, g])); C[:, g, :] = L[rng.integers(len(L), size=N)]
        W = beta[k][None, :, None] + p * C                       # (N, F, s)
        for h in range(F):
            sig[:, :, h] += (W * np.roll(W, -h, axis=1)).sum(1)
    T = np.zeros(F, np.int64); T[0] = 4 * v
    ds = sorted({F // math.gcd(j, F) for j in range(F)})
    hh = np.arange(F)
    Cm = {d: _ramanujan(d, F)[(hh[:, None] - hh[None, :]) % F] for d in ds}
    def bits(level):
        Vs = []; Ts = []
        for d in level:
            Vs.append(np.einsum('hj,nij->nih', Cm[d], sig).reshape(N, -1)); Ts.append(np.tile(Cm[d] @ T, s))
        V = np.concatenate(Vs, 1); TT = np.concatenate(Ts)
        _, keep = np.unique(V[:300].T, axis=0, return_index=True); keep = np.sort(keep)
        return gauss_bits(V[:, keep], TT[keep])
    out = dict(F=F, p=p, s=s, H2=round(H2, 2), orders=ds, all=bits(ds))
    for lv in (levels or []): out['level_' + ','.join(map(str, lv))] = bits(lv)
    if F % 2 == 0:
        n = F // 2; hD = 0.0
        from collections import Counter
        for k in range(4):
            for g in range(n):
                A, Bb = P(int(sigma[k, g])), P(int(sigma[k, g + n]))
                cnt = Counter()
                for a in A:
                    for b in Bb: cnt[tuple(a + b)] += 1
                tot = len(A) * len(Bb)
                hD += sum((m_ / tot) * math.log2(m_) for m_ in cnt.values())
        out['h_D'] = round(hD, 2)
    return out

def real_char_screen(F, p, s, Y):
    """[P, M] exact necessary condition at the real characters chi (chi^2 = 1) of Z_F for a uniform-cyclotomy slice:
    for every class i, sum_k (chihat(beta_k) + p * chihat(c_{k,i}))^2 = 4v, where chihat(c_{k,i}) ranges over the
    values of a +-1 sum with the right parity and |.| <= F, and sum_i chihat(c_{k,i}) = chihat(sigma_k) per block.
    Checked jointly for chi = 1 and chi = (-1)^g (F even) via the per-position pattern multisets (exact DP).
    Returns the number of admissible (t, u) assignments (0 => y has no stage-2 completion)."""
    alpha, sigma, beta, R = class_system(F, p, s, Y)
    v = F * p * p + 1
    chis = [np.ones(F, np.int64)] + ([np.where(np.arange(F) % 2 == 0, 1, -1)] if F % 2 == 0 else [])
    # per block: the set of achievable vectors (t_i, u_i)_{i} (t = sum_g c_i(g), u = sum_g (-1)^g c_i(g)) from the
    # per-position patterns, by DP over positions
    def pats(sg):
        out = []
        for pos in itertools.combinations(range(s), (s + sg) // 2):
            c = -np.ones(s, np.int64); c[list(pos)] = 1; out.append(tuple(c))
        return out
    blocks = []
    for k in range(4):
        cur = {tuple([0] * (len(chis) * s)): 1}
        for g in range(F):
            nxt = {}
            for key, cnt in cur.items():
                for c in pats(int(sigma[k, g])):
                    nk = list(key)
                    for ci, ch in enumerate(chis):
                        for i in range(s): nk[ci * s + i] += int(ch[g]) * c[i]
                    nk = tuple(nk); nxt[nk] = nxt.get(nk, 0) + cnt
            cur = nxt
        bh = [int((ch * beta[k]).sum()) for ch in chis]
        blocks.append((cur, bh))
    # need: for every chi and class i: sum_k (bh_k(chi) + p * T_k(chi, i))^2 = 4v
    # join blocks 0,1 and 2,3 on the vector of partial sums of squares
    def sq(key, bh):
        return tuple((bh[ci] + p * key[ci * s + i]) ** 2 for ci in range(len(chis)) for i in range(s))
    left = {}
    for k0, c0 in blocks[0][0].items():
        a = sq(k0, blocks[0][1])
        for k1, c1 in blocks[1][0].items():
            b = sq(k1, blocks[1][1]); key = tuple(x + y for x, y in zip(a, b))
            left[key] = left.get(key, 0) + c0 * c1
    target = tuple([4 * v] * (len(chis) * s))
    tot = 0
    rq = {}
    for k2, c2 in blocks[2][0].items():
        a = sq(k2, blocks[2][1])
        for k3, c3 in blocks[3][0].items():
            b = sq(k3, blocks[3][1]); need = tuple(t - x - y for t, x, y in zip(target, a, b))
            tot += left.get(need, 0) * c2 * c3
    return tot

def dual_row_screen(F, p, s, Y, chis=None):
    """[P, M] fast exact screen (necessary): for each real character chi of Z_F (chi = 1, and (-1)^g if F is even), the
    class values T_{k,i} = chihat(c_{k,i}) must satisfy, for every class i,
        sum_k (chihat(beta_k) + p T_{k,i})^2 = 4v,
    with T_{k,i} = F (mod 2), |T_{k,i}| <= F, and sum_i T_{k,i} = chihat(sigma_k) for every block k.
    The 4-tuples (T_{0,i}, .., T_{3,i}) solving the per-class equation are listed; then an exact DP over the s classes
    checks whether s of them (with repetition) add up to (chihat(sigma_k))_k.  Characters are checked separately
    (necessary, weaker than the joint screen).  Returns {chi: number of admissible class assignments (ordered)}."""
    alpha, sigma, beta, R = class_system(F, p, s, Y)
    v = F * p * p + 1
    allchis = {'one': np.ones(F, np.int64)}
    if F % 2 == 0: allchis['alt'] = np.where(np.arange(F) % 2 == 0, 1, -1)
    out = {}
    for name, ch in allchis.items():
        bh = [int((ch * beta[k]).sum()) for k in range(4)]
        sh = tuple(int((ch * sigma[k]).sum()) for k in range(4))
        vals = range(-F, F + 1, 2)
        quads = [t for t in itertools.product(vals, repeat=4) if sum((bh[k] + p * t[k]) ** 2 for k in range(4)) == 4 * v]
        # DP over classes: count ordered s-tuples of quads summing to sh
        cur = {(0, 0, 0, 0): 1}
        for i in range(s):
            nxt = {}
            for key, c in cur.items():
                for t in quads:
                    nk = (key[0] + t[0], key[1] + t[1], key[2] + t[2], key[3] + t[3])
                    if all(abs(nk[k] - sh[k]) <= F * (s - i - 1) for k in range(4)):
                        nxt[nk] = nxt.get(nk, 0) + c
            cur = nxt
        out[name] = dict(quads=len(quads), assignments=cur.get(sh, 0))
    return out
