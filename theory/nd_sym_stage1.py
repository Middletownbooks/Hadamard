r"""Exact (exhaustive) stage 1 of a symmetric generalised-Paley slice: decides 739 x215 [P for the reduction, M].

Slice (739 x215: F = 18, q = 41, s = 8, f = 5).  G = Z_F x Z_q with F = 2n, n odd; the multiplier group is
<-1 on Z_F> x <index-s subgroup of Z_q^*>, s even.  Stage 1 = the Z_q-coset sums y_k on Z_F:
    y_k symmetric (y(-g) = y(g)),  y = alpha + 2f t  (alpha = +-1, |t| <= s/2),
    sum_k PAF(y_k)(0) = 4Fq - 4(q - 1),  sum_k PAF(y_k)(h) = -4q (h != 0),  rows (+-2, 0, 0, 0).

Reduction [P] (proved in report/next_design_answer.md, section 2):
 1. alpha-congruence: PAF(y) = PAF(alpha) + 4f X + 4f^2 PAF(t) with X(h) = sum_g alpha(g) t(g + h) (the two cross
    sums agree for symmetric alpha, t).  Hence sum_k PAF(alpha_k)(h) = T(h) (mod 4f) and rows of alpha = R_k (mod 2f).
    The alpha_k are symmetric +-1 sequences (2^(n+1) each): a finite list of alpha-quadruples.
 2. S/D split on Z_F = Z_2 x Z_n: sigma(b) = (y(b) + y(b+n))/2, delta(b) = (-1)^b (y(b) - y(b+n))/2 (b = 0..n-1).
    Both are symmetric on Z_n, and stage 1 is equivalent to
        sum_k PAF_n(sigma_k) = (1 + 2q(n-1), -2q, ..., -2q),   rows sum_b sigma_k(b) = R_k / 2,
        sum_k PAF_n(delta_k) = (v, 0, ..., 0)                   (v = Fq + 1).
 3. Given alpha: sigma(b) = (alpha(b) + alpha(b+n))/2 + f u(b), delta(b) = (-1)^b ((alpha(b) - alpha(b+n))/2 + f w(b))
    with u = t(b) + t(b+n), w = t(b) - t(b+n): u = w (mod 2), |u| + |w| <= s.  The S-system involves only u, the
    D-system only w.
Search: every alpha-quadruple (each block up to its own shift by n and, for row-0 blocks, negation; row-0 blocks in
sorted order) -> exhaustive 2+2 join of the S-system over u -> for each S-solution an exhaustive 2+2 join of the
D-system over w (w restricted by u) -> every hit verified on Z_F directly.  Exhaustive: no solution printed and
"complete": true means stage 1 is empty.
usage: python3 theory/nd_sym_stage1.py F q s [--plant N | --plant-real N] [--max-sol N] [--out FILE]"""
import sys, os, json, time, math, itertools
import numpy as np

def sym_seqs(F):
    """all symmetric +-1 sequences on Z_F (orbit values at the negation-orbit representatives 0..F/2)"""
    reps = list(range(F // 2 + 1)); out = []
    for msk in range(1 << len(reps)):
        a = np.empty(F, np.int64)
        for i, g in enumerate(reps):
            v = 1 if (msk >> i) & 1 else -1; a[g] = v; a[(-g) % F] = v
        out.append(a)
    return np.array(out)

def paf(a, H=None):
    F = a.shape[-1]; H = F // 2 if H is None else H
    return np.stack([(a * np.roll(a, -h, axis=-1)).sum(-1) for h in range(0, H + 1)], -1)

def targets(F, q):
    n = F // 2; v = F * q + 1
    T = np.array([4 * F * q - 4 * (q - 1)] + [-4 * q] * (F // 2), np.int64)
    TS = np.array([1 + 2 * q * (n - 1)] + [-2 * q] * ((n - 1) // 2), np.int64)
    TD = np.array([v] + [0] * ((n - 1) // 2), np.int64)
    return T, TS, TD

def shift_class(F, A):
    """canonical index of each symmetric sequence under shift by n and (optionally) negation"""
    n = F // 2; idx = {tuple(a): i for i, a in enumerate(A.tolist())}
    A = [np.array(a) for a in A.tolist()]
    sh = np.array([idx[tuple(np.roll(a, -n))] for a in A])
    ng = np.array([idx[tuple(-a)] for a in A])
    return sh, ng

ENC = 4096

def enc(S):
    """int64 key of signature rows (components offset by ENC/2, 12 bits each)"""
    k = np.zeros(S.shape[:-1], np.int64)
    for c in range(S.shape[-1]): k = (k << 12) | (S[..., c] + ENC // 2)
    return k

def pair_join(L0, L1, L2, L3, T):
    """all (i, j, k, l) with L0[i] + L1[j] + L2[k] + L3[l] = T (signature arrays), exhaustive"""
    if min(len(L0), len(L1), len(L2), len(L3)) == 0: return []
    left = enc(L0[:, None, :] + L1[None, :, :]).ravel()
    right = enc(T[None, None, :] - L2[:, None, :] - L3[None, :, :]).ravel()
    o = np.argsort(left, kind='stable'); ls = left[o]
    lo = np.searchsorted(ls, right, 'left'); hi = np.searchsorted(ls, right, 'right')
    out = []; n1 = L1.shape[0]; n3 = L3.shape[0]
    for r in np.nonzero(hi > lo)[0]:
        for p in o[lo[r]:hi[r]]:
            out.append((int(p // n1), int(p % n1), int(r // n3), int(r % n3)))
    return out

def zn_expand(vals, n):
    """symmetric sequence on Z_n from orbit values at 0..(n-1)/2 (last axis)"""
    h = (n - 1) // 2; idx = [0] + list(range(1, h + 1)) + list(range(h, 0, -1))
    return vals[..., idx]

class Block:
    """per-block data for a fixed symmetric alpha on Z_F (F = 2n)"""
    def __init__(self, alpha, F, f, s, R, TS0, TD0):
        n = F // 2; h = (n - 1) // 2; self.n = n
        b = np.arange(h + 1)
        self.aS = (alpha[b] + alpha[b + n]) // 2                  # sigma = aS + f u
        self.aD = (alpha[b] - alpha[b + n]) // 2                  # delta = (-1)^b (aD + f w)
        self.sgn = np.where(b % 2 == 0, 1, -1)
        self.f, self.s, self.R, self.alpha = f, s, R, alpha
        wts = np.array([1] + [2] * h)
        # S-list: u in [-s, s]^(h+1), sigma energy <= TS0, row = R/2
        grids = []
        for i in range(h + 1):
            u = np.arange(-s, s + 1); sg = self.aS[i] + f * u
            grids.append(u[wts[i] * sg * sg <= TS0])
        U = np.array(list(itertools.product(*grids)), np.int64).reshape(-1, h + 1)
        sig = self.aS[None, :] + f * U
        ok = ((wts[None, :] * sig * sig).sum(1) <= TS0) & ((wts[None, :] * sig).sum(1) * 2 == R)
        self.U = U[ok]; self.sigma = zn_expand(sig[ok], n)
        self.Ssig = paf(self.sigma, h) if len(self.U) else np.zeros((0, h + 1), np.int64)
        self.TD0 = TD0; self.wts = wts

    def dlist(self, u):
        """D-list given u: w = u (mod 2), |u| + |w| <= s, delta energy <= TD0"""
        h = len(u) - 1; grids = []
        for i in range(h + 1):
            lim = self.s - abs(int(u[i])); w = np.arange(-lim, lim + 1); w = w[(w - u[i]) % 2 == 0]
            dl = self.aD[i] + self.f * w
            grids.append(w[self.wts[i] * dl * dl <= self.TD0])
        W = np.array(list(itertools.product(*grids)), np.int64).reshape(-1, h + 1)
        dl = self.sgn[None, :] * (self.aD[None, :] + self.f * W)
        ok = (self.wts[None, :] * dl * dl).sum(1) <= self.TD0
        W = W[ok]; delta = zn_expand(dl[ok], self.n)
        return W, (paf(delta, h) if len(W) else np.zeros((0, h + 1), np.int64))

    def y(self, u, w):
        n = self.n; h = (n - 1) // 2; F = 2 * n
        t = np.zeros(F, np.int64)
        for b in range(h + 1):
            tb = (u[b] + w[b]) // 2; tbn = (u[b] - w[b]) // 2
            for g, val in ((b, tb), (b + n, tbn)):
                t[g % F] = val; t[(-g) % F] = val
        return self.alpha + 2 * self.f * t

def verify(Y, F, q):
    T, _, _ = targets(F, q)
    Y = np.array(Y); P = paf(Y).sum(0)
    rows = sorted(np.abs(Y.sum(1)).tolist())
    sym = all((y == np.roll(y[::-1], 1)).all() for y in Y)
    return bool((P == T).all() and rows == [0, 0, 0, 2] and sym)

def alpha_quads(F, q, s, T, plant=None):
    f = (q - 1) // s; mod = 4 * f; rmod = 2 * f
    A = sym_seqs(F); P = paf(A)[:, 1:]; R = A.sum(1)
    sh, ng = shift_class(F, A)
    canon2 = np.minimum(np.arange(len(A)), sh)                                  # block 0: shift only
    canon0 = np.minimum.reduce([np.arange(len(A)), sh, ng, sh[ng]])             # row-0 blocks: shift and negation
    I2 = np.nonzero((R - 2) % rmod == 0)[0]; I0 = np.nonzero(R % rmod == 0)[0]
    I2 = I2[canon2[I2] == I2]; I0c = I0[canon0[I0] == I0]
    Tm = T[1:] % mod
    base = np.array([mod ** i for i in range(P.shape[1])], np.int64)     # mod^9 < 2^63 for mod <= 100
    K = lambda M: (M % mod) @ base
    # left: (block0 in I2, block1 in I0c); right: (block2, block3 in I0c, block2 <= block3)
    left = K(P[I2][:, None, :] + P[I0c][None, :, :]).ravel()
    pr = [(a, b) for a in range(len(I0c)) for b in range(a, len(I0c))]
    rk = K(Tm[None, :] - P[I0c[[a for a, b in pr]]] - P[I0c[[b for a, b in pr]]])
    o = np.argsort(left); ls = left[o]
    lo = np.searchsorted(ls, rk, 'left'); hi = np.searchsorted(ls, rk, 'right')
    quads = set()
    for r in np.nonzero(hi > lo)[0]:
        a, b = pr[r]
        for p in o[lo[r]:hi[r]]:
            i0 = int(I2[p // len(I0c)]); i1 = int(I0c[p % len(I0c)])
            trip = tuple(sorted((i1, int(I0c[a]), int(I0c[b]))))
            quads.add((i0,) + trip)
    return A, sorted(quads), dict(row2=int(len(I2)), row0_classes=int(len(I0c)))

def solve(F, q, s, max_sol=10 ** 9, log=print, T=None, TS=None, TD=None, rows=(2, 0, 0, 0), quads=None):
    f = (q - 1) // s; n = F // 2; h = (n - 1) // 2
    T0, TS0, TD0 = targets(F, q)
    T = T0 if T is None else T; TS = TS0 if TS is None else TS; TD = TD0 if TD is None else TD
    t0 = time.time()
    A, Q, info = alpha_quads(F, q, s, T)
    if quads is not None: Q = quads
    log(json.dumps(dict(stage="alpha", F=F, q=q, s=s, f=f, alpha_quads=len(Q), **info, secs=round(time.time() - t0, 1))))
    sols = []; nS = 0; cache = {}
    for qi, quad in enumerate(Q):
        blocks = []
        for k, ai in enumerate(quad):
            key = (ai, rows[k])
            if key not in cache: cache[key] = Block(A[ai], F, f, s, rows[k], int(TS[0]), int(TD[0]))
            blocks.append(cache[key])
        Sh = pair_join(*[b.Ssig for b in blocks], TS)
        nS += len(Sh)
        for (i, j, k, l) in Sh:
            us = [blocks[0].U[i], blocks[1].U[j], blocks[2].U[k], blocks[3].U[l]]
            Ds = [b.dlist(u) for b, u in zip(blocks, us)]
            Dh = pair_join(*[d[1] for d in Ds], TD)
            for (a, b_, c, d) in Dh:
                ws = [Ds[0][0][a], Ds[1][0][b_], Ds[2][0][c], Ds[3][0][d]]
                Y = [blk.y(u, w).tolist() for blk, u, w in zip(blocks, us, ws)]
                ok = verify(Y, F, q) if T is T0 else bool((paf(np.array(Y)).sum(0) == T).all())
                sols.append(dict(alpha_quad=[int(x) for x in quad], y=Y, verified=ok))
                if len(sols) >= max_sol: break
            if len(sols) >= max_sol: break
        if len(sols) >= max_sol: break
    out = dict(stage="done", F=F, q=q, s=s, alpha_quads=len(Q), S_solutions=nS, stage1_solutions=len(sols),
               all_verified=all(x['verified'] for x in sols), complete=len(sols) < max_sol, secs=round(time.time() - t0, 1))
    log(json.dumps(out))
    return out, sols

def planted(F, q, s, seed, log=print):
    """plant a random symmetric quadruple (block 0 row 2, others row 0 by rejection), set the targets to its own
    PAF sums, and check that the exhaustive search returns it"""
    rng = np.random.default_rng(seed); f = (q - 1) // s; n = F // 2; h = (n - 1) // 2
    reps = list(range(F // 2 + 1))
    def rnd(R):
        while True:
            y = np.empty(F, np.int64)
            for g in reps:
                a = rng.choice([-1, 1]); t = int(rng.choice([-1, 0, 0, 0, 0, 1])) if rng.random() < 0.6 else 0
                y[g] = a + 2 * f * t; y[(-g) % F] = y[g]
            if y.sum() == R: return y
    Y = np.array([rnd(2), rnd(0), rnd(0), rnd(0)])
    T = paf(Y).sum(0)
    sig = (Y[:, :n] + Y[:, n:]) // 2; dl = np.where(np.arange(n) % 2 == 0, 1, -1)[None, :] * (Y[:, :n] - Y[:, n:]) // 2
    TS = paf(sig, h).sum(0); TD = paf(dl, h).sum(0)
    # alpha-quadruples for this target, then search all of them
    A, Q, info = alpha_quads(F, q, s, T)
    out, sols = solve(F, q, s, log=lambda x: None, T=T, TS=TS, TD=TD, quads=Q)
    canon = lambda Z: sorted(tuple(z) for z in Z)
    # the planted y is found up to the per-block symmetries (shift by n, negation of row-0 blocks) and block order
    def variants(y, R):
        vs = [y, np.roll(y, -n)]
        if R == 0: vs += [-y, -np.roll(y, -n)]
        return {tuple(v) for v in vs}
    V = [variants(Y[k], [2, 0, 0, 0][k]) for k in range(4)]
    hit = False
    for x in sols:
        Z = [tuple(z) for z in x['y']]
        if Z[0] in V[0] and any(all(Z[1 + i] in V[1 + p[i]] for i in range(3)) for p in itertools.permutations(range(3))):
            hit = True; break
    res = dict(plant_seed=seed, alpha_quads=len(Q), solutions=len(sols), planted_found=hit, secs=out['secs'])
    log(json.dumps(res)); return res

def planted_real(F, q, s, seed, log=print):
    """planted test with realistic plants: block configurations drawn uniformly (exact DP sampler of
    theory/nd_sym_s1est.py) from the real lag-0 energy and rows (2, 0, 0, 0); targets = the plant's own PAF sums"""
    import nd_sym_s1est as E
    rng = np.random.default_rng(1000 + seed); f = (q - 1) // s; n = F // 2; h = (n - 1) // 2
    L, _ = E.letters(q, s); g = 0
    for y in L: g = math.gcd(g, y * y - 1)
    unit = {y: (y * y - 1) // g for y in L}
    T0, _, _ = targets(F, q); Etot = (int(T0[0]) - 4 * F) // g
    reps, tabs = E.block_dp(F, L, Etot, unit)
    c = {}
    for (e, r), cnt in tabs[0].items(): c[(e, r)] = cnt
    while True:
        cut = np.sort(rng.integers(0, Etot + 1, 3)); es = [cut[0], cut[1] - cut[0], cut[2] - cut[1], Etot - cut[2]]
        if all(c.get((int(e), r), 0) for e, r in zip(es, (2, 0, 0, 0))): break
    Y = np.array([E.expand(F, reps, E.sample_block(reps, tabs, L, unit, int(e), r, rng)) for e, r in zip(es, (2, 0, 0, 0))])
    T = paf(Y).sum(0)
    sig = (Y[:, :n] + Y[:, n:]) // 2; dl = np.where(np.arange(n) % 2 == 0, 1, -1)[None, :] * (Y[:, :n] - Y[:, n:]) // 2
    TS = paf(sig, h).sum(0); TD = paf(dl, h).sum(0)
    A, Q, info = alpha_quads(F, q, s, T)
    out, sols = solve(F, q, s, log=lambda x: None, T=T, TS=TS, TD=TD, quads=Q)
    def variants(y, R):
        vs = [y, np.roll(y, -n)]
        if R == 0: vs += [-y, -np.roll(y, -n)]
        return {tuple(v) for v in vs}
    V = [variants(Y[k], [2, 0, 0, 0][k]) for k in range(4)]
    hit = False
    for x in sols:
        Z = [tuple(z) for z in x['y']]
        if Z[0] in V[0] and any(all(Z[1 + i] in V[1 + p[i]] for i in range(3)) for p in itertools.permutations(range(3))):
            hit = True; break
    res = dict(plant_seed=seed, energies=[int(e) for e in es], letters_used=sorted({abs(int(v)) for v in Y.ravel()}),
               alpha_quads=len(Q), S_solutions=out['S_solutions'], solutions=len(sols), planted_found=hit, secs=out['secs'])
    log(json.dumps(res)); return res

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--plant": None, "--plant-real": None, "--max-sol": str(10 ** 9), "--out": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    F, q, s = pos
    if opt["--plant"] is not None:
        for sd in range(int(opt["--plant"])): planted(F, q, s, sd, log=lambda x: print(x, flush=True))
        sys.exit()
    if opt["--plant-real"] is not None:
        r = [planted_real(F, q, s, sd, log=lambda x: print(x, flush=True)) for sd in range(int(opt["--plant-real"]))]
        print(json.dumps(dict(planted_real=len(r), found=sum(x['planted_found'] for x in r))), flush=True)
        sys.exit()
    out, sols = solve(F, q, s, int(opt["--max-sol"]), log=lambda x: print(x, flush=True))
    if opt["--out"]:
        json.dump(dict(result=out, solutions=sols[:1000]), open(opt["--out"], "w"))
