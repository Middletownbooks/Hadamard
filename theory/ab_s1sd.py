r"""Exact stage 1 of a slice on G' = Z_2 x G'' (S/D split on axis 0), any target T1 (b = 0, 1, >= 3) [P, M].

Generalises theory/nd_s1_sd.py (cyclic Z_F, b = 1) to any G' with a Z_2 direct factor on axis 0 (U acting on G'' only),
and to the two-valued b >= 3 targets.  [P]
 1. alpha-congruence: every letter is = +-1 (mod 2 f0); y = alpha + 2 f0 t (alpha = +-1); the cross term
    2 f0 (CC(alpha, t) + CC(t, alpha)) is 2 f0 times an even number, so sum_k PAF(alpha_k) = T1 (mod 4 f0) at every
    lag, and rows of alpha = rows of y (mod 2 f0).
 2. S/D: sigma(b) = (y(0, b) + y(1, b))/2, delta(b) = (y(0, b) - y(1, b))/2 on G'' satisfy
        sum_k PAF_{G''}(sigma_k) = T_S = (T1(0, .) + T1(1, .))/4,   sum_k PAF_{G''}(delta_k) = T_D = (T1(0, .) - T1(1, .))/4,
    and given alpha, sigma depends only on u = t(0, .) + t(1, .), delta only on w = t(0, .) - t(1, .) (u = w mod 2).
Search (exhaustive): alpha-quadruples (each block up to sign, blocks sorted) with rows from an admissible row quadruple
(sum r^2 = R4) -> 2+2 join of the S-system over u -> per S-solution a 2+2 join of the D-system over w -> verify.
usage: python3 theory/ab_s1sd.py CFGJSON [--out FILE] [--max-sol N] [--plant N]"""
import sys, os, json, time, itertools, math
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_lib import AG, alpha_congruence_mod
from ab_gauss import lag_reps
from nd_s1_sd import pair_join


def cfg_from_json(cfgs):
    cfg = json.load(open(cfgs[1:])) if cfgs.startswith('@') else json.loads(cfgs)
    cfg['mods'] = tuple(cfg['mods']); cfg['Q'] = tuple(cfg['Q']); cfg['H1'] = tuple(cfg['H1'])
    cfg['U'] = [tuple(u) if u[0] != 'neg' else ('neg', tuple(u[1])) for u in cfg['U']]
    return cfg


def ejoin4(L0, L1, L2, L3, T, max_pairs=3 * 10 ** 8, max_out=10 ** 7):
    """all (i, j, k, l) with L0[i] + L1[j] + L2[k] + L3[l] = T (exact), coordinate 0 = an energy (non-negative,
    additive): pairs are formed only within matching energy splits; keys are 64-bit linear hashes, every match is
    re-checked exactly"""
    if min(len(L0), len(L1), len(L2), len(L3)) == 0: return []
    d = L0.shape[1]
    rng = np.random.default_rng(99)
    wv = rng.integers(1, 2 ** 62, size=d, dtype=np.int64)
    H = [L @ wv for L in (L0, L1, L2, L3)]
    E = [L[:, 0] for L in (L0, L1, L2, L3)]
    T0 = int(T[0]); th = np.int64(int((T.astype(np.int64) @ wv)))
    def groups(e):
        o = np.argsort(e, kind='stable'); es = e[o]
        u, st = np.unique(es, return_index=True)
        en = np.r_[st[1:], len(es)]
        return {int(a): o[b:c] for a, b, c in zip(u, st, en)}
    G = [groups(e) for e in E]
    out = []; work = 0
    # left energies EL = e0 + e1, right = T0 - EL
    lsum = defaultdict(list)
    for e0, i0 in G[0].items():
        for e1, i1 in G[1].items():
            if e0 + e1 <= T0: lsum[e0 + e1].append((i0, i1))
    rsum = defaultdict(list)
    for e2, i2 in G[2].items():
        for e3, i3 in G[3].items():
            if e2 + e3 <= T0: rsum[e2 + e3].append((i2, i3))
    for EL, lefts in lsum.items():
        rights = rsum.get(T0 - EL)
        if not rights: continue
        nl = sum(len(a) * len(b) for a, b in lefts); nr = sum(len(a) * len(b) for a, b in rights)
        work += nl + nr
        if work > max_pairs: raise RuntimeError(f"ejoin4 work cap ({work})")
        LH = np.concatenate([(H[0][a][:, None] + H[1][b][None, :]).ravel() for a, b in lefts])
        LI = np.concatenate([np.stack(np.broadcast_arrays(a[:, None], b[None, :]), -1).reshape(-1, 2) for a, b in lefts])
        RH = np.concatenate([(th - (H[2][a][:, None] + H[3][b][None, :])).ravel() for a, b in rights])
        RI = np.concatenate([np.stack(np.broadcast_arrays(a[:, None], b[None, :]), -1).reshape(-1, 2) for a, b in rights])
        o = np.argsort(LH); LS = LH[o]
        lo = np.searchsorted(LS, RH, 'left'); hi = np.searchsorted(LS, RH, 'right')
        for r in np.nonzero(hi > lo)[0]:
            k, l = RI[r]
            for p in o[lo[r]:hi[r]]:
                i, j = LI[p]
                if ((L0[i] + L1[j] + L2[k] + L3[l]) == T).all():
                    out.append((int(i), int(j), int(k), int(l)))
                    if len(out) > max_out: raise RuntimeError("ejoin4 output cap")
    return out


class S1SD:
    def __init__(self, cfg, T1=None, sl=None, zallow=None, nQ=None):
        """zallow: optional list of allowed Q-compressions (length-c vectors, Q index = ravel of the first nQ axes);
        each block's compression must be +-one of them (prunes alphas, S-lists and D-lists by coset sums)"""
        self.cfg = cfg
        sl = W.build(cfg) if sl is None else sl; self.sl = sl
        Gp = sl.Gp; self.Gp = Gp
        assert Gp.mods[0] == 2, "axis 0 of G' must be Z_2"
        for u in sl.Sp: assert u[0] % 2 == 1
        self.T1 = sl.T1.copy() if T1 is None else np.array(T1, np.int64)
        self.R4 = int(self.T1.sum())
        self.Gpp = AG(Gp.mods[1:])
        # G''-orbits under U and the pairing with G'-orbits
        Spp = sorted({tuple(u[1:]) for u in sl.Sp})
        self.Upp = [self.Gpp.act(u) for u in Spp] if len(self.Gpp.mods) else [np.zeros(1, np.int64)]
        n2 = self.Gpp.n
        Oo = -np.ones(n2, np.int64); reps = []
        for b in range(n2):
            if Oo[b] >= 0: continue
            for P in self.Upp: Oo[P[b]] = len(reps)
            reps.append(b)
        self.Oo, self.oreps = Oo, reps; self.ow = np.bincount(Oo)
        # G' element (e, b) has index e * n2 + b (axis 0 slowest)
        self.e0 = np.array([b for b in reps]); self.e1 = np.array([n2 + b for b in reps])
        self.o0 = sl.borb[self.e0]; self.o1 = sl.borb[self.e1]
        self.f0 = alpha_congruence_mod(sl.letters)
        assert self.f0 >= 2
        f0 = self.f0
        self.tset = [{a: sorted((y - a) // (2 * f0) for y in L if (y - a) % (2 * f0) == 0) for a in (1, -1)}
                     for L in sl.letters]
        self.reps, _ = lag_reps(Gp, [Gp.act(u) for u in sl.Sp])
        self.preps, _ = lag_reps(self.Gpp, self.Upp)
        T1 = self.T1
        self.TS = np.array([(T1[b] + T1[n2 + b]) // 4 for b in range(n2)], np.int64)
        self.zallow = None
        if zallow is not None:
            nQ = cfg['nQ'] if nQ is None else nQ
            assert nQ >= 1
            Za = np.array(zallow, np.int64); Za = np.concatenate([Za, -Za])
            self.zallow = np.unique(Za, axis=0)
            Qmods = Gp.mods[:nQ]; cq = int(np.prod(Qmods)); cpp = cq // 2
            # Q'' = Q / Z_2 (axis 0): coset of b'' in G'' by its first nQ-1 coordinates
            self.qpp = AG(Qmods[1:]).idx(self.Gpp.co[:, :nQ - 1]) if nQ > 1 else np.zeros(n2, np.int64)
            self.cpp = cpp
            z0 = self.zallow[:, :cpp]; z1 = self.zallow[:, cpp:]
            self.zS = {tuple(r) for r in ((z0 + z1) // 2).tolist()}
            self.zD = {tuple(r) for r in ((z0 - z1) // 2).tolist()}
            self.zSD = {(tuple(((a + b) // 2).tolist()), tuple(((a - b) // 2).tolist())) for a, b in zip(z0, z1)}
            self.qfull = AG(Qmods).idx(Gp.co[:, :nQ])
        self.TD = np.array([(T1[b] - T1[n2 + b]) // 4 for b in range(n2)], np.int64)
        assert ((T1[:n2] + T1[n2:]) % 4 == 0).all() and ((T1[:n2] - T1[n2:]) % 4 == 0).all()

    def alphas(self, keep=None, chunk=1 << 20):
        """all alpha vectors (one sign per G'-orbit, restricted by tset); keep(vals) -> bool mask filters each chunk,
        so that L1 up to ~30 is possible when a strong filter (the border class) is given"""
        L1 = len(self.sl.breps)
        allowed = [[a for a in (1, -1) if self.tset[j][a]] for j in range(L1)]
        free = [j for j in range(L1) if len(allowed[j]) == 2]
        if len(free) > 24 and keep is None: raise RuntimeError("too many G'-orbits for alpha enumeration")
        n = 1 << len(free)
        out = []
        for c0 in range(0, n, chunk):
            idx = np.arange(c0, min(n, c0 + chunk), dtype=np.int64)
            vals = np.empty((len(idx), L1), np.int64)
            for j in range(L1):
                if len(allowed[j]) == 1: vals[:, j] = allowed[j][0]
            for t, j in enumerate(free): vals[:, j] = 2 * ((idx >> t) & 1) - 1
            out.append(vals[keep(vals)] if keep is not None else vals)
        return np.concatenate(out) if out else np.zeros((0, L1), np.int64)

    def alpha_quads(self, rowquads):
        """canonical alpha-quadruples: alpha_k up to sign (alpha at orbit of 0 = +1 if allowed), sorted"""
        f0 = self.f0; mod, rmod = 4 * f0, 2 * f0; sl = self.sl
        if self.zallow is not None:
            cq = self.zallow.shape[1]; md = 2 * f0
            Wq = np.zeros((len(sl.breps), cq), np.int64)          # orbit j -> number of its elements in coset q
            for g_ in range(len(sl.borb)): Wq[sl.borb[g_], self.qfull[g_]] += 1
            base = md ** np.arange(cq, dtype=np.int64)
            zcode = np.unique(((self.zallow % md) @ base))
            def keepf(A):
                return np.isin(((A @ Wq) % md) @ base, zcode)
            A = self.alphas(keep=keepf)
        else:
            A = self.alphas()
        # canonical sign: keep alpha and its negative once (if -alpha is also admissible)
        key = {tuple(a): i for i, a in enumerate(A.tolist())}
        keep = []
        for i, a in enumerate(A.tolist()):
            j = key.get(tuple(-x for x in a))
            if j is None or i <= j: keep.append(i)
        A = A[keep]
        full = A[:, sl.borb]
        P = self.Gp.paf(full)[:, self.reps] % mod
        R = full.sum(1)
        T = self.T1[self.reps] % mod
        # admissible residues of rows: each block row = +-r for r in a row quadruple (block order free)
        rq = set(tuple(sorted(q)) for q in rowquads)
        okrows = set()
        for q in rq:
            for perm in set(itertools.permutations(q)):
                for sg in itertools.product((1, -1), repeat=4):
                    okrows.add(tuple(sorted(((s_ * r) % rmod) for s_, r in zip(sg, perm))))
        # alpha row residues can be +-; since alpha sign is canonical we allow both signs per block
        # pair join on P
        nA = len(A)
        keyP = [tuple(p) for p in P.tolist()]
        byP = defaultdict(list)
        for i in range(nA): byP[keyP[i]].append(i)
        pk = list(byP)
        Pk = np.array(pk, np.int64)
        pairs = defaultdict(list)
        for i in range(len(pk)):
            S = (Pk[i] + Pk[i:]) % mod
            for jj, srow in enumerate(S): pairs[tuple(srow.tolist())].append((i, i + jj))
        quads = set()
        for s1, lst in pairs.items():
            need = tuple(((T - np.array(s1)) % mod).tolist())
            if need not in pairs: continue
            for (a, b) in lst:
                for (c, d) in pairs[need]:
                    if (c, d) < (a, b): continue
                    for ia in byP[pk[a]]:
                        for ib in byP[pk[b]]:
                            for ic in byP[pk[c]]:
                                for idd in byP[pk[d]]:
                                    qd = tuple(sorted((ia, ib, ic, idd)))
                                    if qd in quads: continue
                                    rr = tuple(sorted(int(R[x]) % rmod for x in qd))
                                    # rows: each block row = +-r_k (mod rmod); test all sign choices
                                    ok = False
                                    for sg in itertools.product((1, -1), repeat=4):
                                        if tuple(sorted((s_ * int(R[x])) % rmod for s_, x in zip(sg, qd))) in okrows:
                                            ok = True; break
                                    if ok: quads.add(qd)
        return A, sorted(quads)

    def block(self, alpha, TS0, TD0):
        f0 = self.f0
        a0 = alpha[self.o0]; a1 = alpha[self.o1]
        aS = (a0 + a1) // 2; aD = (a0 - a1) // 2
        T0 = [self.tset[self.o0[j]][int(a0[j])] for j in range(len(self.oreps))]
        T1_ = [self.tset[self.o1[j]][int(a1[j])] for j in range(len(self.oreps))]
        pairs = [[(t0, t1) for t0 in x for t1 in y] for x, y in zip(T0, T1_)]
        uopts = [sorted({t0 + t1 for t0, t1 in p}) for p in pairs]
        w = self.ow
        nr = len(self.oreps)
        minrest = np.zeros(nr + 1, np.int64)
        for j in range(nr - 1, -1, -1):
            minrest[j] = minrest[j + 1] + w[j] * min((aS[j] + f0 * u) ** 2 for u in uopts[j])
        U = []; cur = [0] * nr
        def rec(j, e):
            if j == nr:
                U.append(list(cur)); return
            for u in uopts[j]:
                sg = aS[j] + f0 * u; e2 = e + w[j] * sg * sg
                if e2 + minrest[j + 1] <= TS0: cur[j] = u; rec(j + 1, e2)
        rec(0, 0)
        U = np.array(U, np.int64).reshape(-1, nr)
        sig = aS[None, :] + f0 * U
        sigma = sig[:, self.Oo]
        if self.zallow is not None and len(U):
            cs = np.zeros((len(U), self.cpp), np.int64)
            for t in range(self.cpp): cs[:, t] = sigma[:, self.qpp == t].sum(1)
            keep = np.array([tuple(r) in self.zS for r in cs.tolist()])
            U = U[keep]; sig = sig[keep]; sigma = sigma[keep]
        Ssig = self.Gpp.paf(sigma)[:, self.preps] if len(U) else np.zeros((0, len(self.preps)), np.int64)
        return dict(alpha=alpha, aS=aS, aD=aD, pairs=pairs, U=U, Ssig=Ssig, TD0=TD0)

    def dlist(self, B, u):
        f0 = self.f0; w = self.ow; nr = len(self.oreps)
        wopts = [sorted({t0 - t1 for t0, t1 in B['pairs'][j] if t0 + t1 == u[j]}) for j in range(nr)]
        minrest = np.zeros(nr + 1, np.int64)
        for j in range(nr - 1, -1, -1):
            minrest[j] = minrest[j + 1] + w[j] * min((B['aD'][j] + f0 * x) ** 2 for x in wopts[j])
        Wl = []; cur = [0] * nr
        def rec(j, e):
            if j == nr:
                Wl.append(list(cur)); return
            for x in wopts[j]:
                dl = B['aD'][j] + f0 * x; e2 = e + w[j] * dl * dl
                if e2 + minrest[j + 1] <= B['TD0']: cur[j] = x; rec(j + 1, e2)
        rec(0, 0)
        Wl = np.array(Wl, np.int64).reshape(-1, nr)
        delta = (B['aD'][None, :] + f0 * Wl)[:, self.Oo]
        if self.zallow is not None and len(Wl):
            sg = (B['aS'] + f0 * np.asarray(u))[self.Oo]
            sS = tuple(int(sg[self.qpp == t].sum()) for t in range(self.cpp))
            okD = {d for (a, d) in self.zSD if a == sS}
            cs = np.zeros((len(Wl), self.cpp), np.int64)
            for t in range(self.cpp): cs[:, t] = delta[:, self.qpp == t].sum(1)
            keep = np.array([tuple(r) in okD for r in cs.tolist()])
            Wl = Wl[keep]; delta = delta[keep]
        return Wl, (self.Gpp.paf(delta)[:, self.preps] if len(Wl) else np.zeros((0, len(self.preps)), np.int64))

    def y_of(self, B, u, w):
        n2 = self.Gpp.n; y = np.zeros(2 * n2, np.int64)
        t0 = (u + w) // 2; t1 = (u - w) // 2
        alpha_full = B['alpha'][self.sl.borb]
        for j, b in enumerate(self.oreps):
            for bb in np.nonzero(self.Oo == j)[0]:
                y[bb] = alpha_full[bb] + 2 * self.f0 * t0[j]
                y[n2 + bb] = alpha_full[n2 + bb] + 2 * self.f0 * t1[j]
        return y

    def verify(self, Y):
        Y = np.array(Y)
        ok = (self.Gp.paf(Y).sum(0) == self.T1).all()
        inv = all((Y[:, self.Gp.act(u)] == Y).all() for u in self.sl.Sp)
        lett = all(int(Y[k, b]) in self.sl.letters[self.sl.borb[b]] for k in range(4) for b in range(self.Gp.n))
        return bool(ok and inv and lett)

    def solve(self, rowquads, max_sol=10 ** 9, log=print, shuffle=0, progress=0, only=None):
        t0 = time.time()
        A, Q = self.alpha_quads(rowquads)
        if only is not None: Q = [q for q in Q if q in set(only)]
        log(json.dumps(dict(stage="alpha", f0=self.f0, alphas=len(A), alpha_quads=len(Q), secs=round(time.time() - t0, 1))))
        TS0, TD0 = int(self.TS[0]), int(self.TD[0])
        cache = {}; sols = []; nS = 0
        TS = self.TS[self.preps]; TD = self.TD[self.preps]
        if shuffle:
            rr = np.random.default_rng(shuffle); Q = [Q[i] for i in rr.permutation(len(Q))]
        for qi, quad in enumerate(Q):
            if progress and qi % progress == 0:
                log(json.dumps(dict(stage="progress", quads_done=qi, of=len(Q), S_solutions=nS, solutions=len(sols),
                                    secs=round(time.time() - t0, 1))))
            if len(cache) > 3000: cache.clear()
            blocks = []
            for ai in quad:
                if ai not in cache: cache[ai] = self.block(A[ai], TS0, TD0)
                blocks.append(cache[ai])
            Sh = ejoin4(*[b['Ssig'] for b in blocks], TS)
            nS += len(Sh)
            for (i, j, k, l) in Sh:
                us = [blocks[0]['U'][i], blocks[1]['U'][j], blocks[2]['U'][k], blocks[3]['U'][l]]
                Ds = [self.dlist(b, u) for b, u in zip(blocks, us)]
                for (a, b_, c, d) in ejoin4(*[x[1] for x in Ds], TD):
                    ws = [Ds[0][0][a], Ds[1][0][b_], Ds[2][0][c], Ds[3][0][d]]
                    Y = [self.y_of(B, u, w_).tolist() for B, u, w_ in zip(blocks, us, ws)]
                    sols.append(dict(alpha_quad=[int(x) for x in quad], y=Y, verified=self.verify(Y)))
                    if len(sols) >= max_sol: break
                if len(sols) >= max_sol: break
            if len(sols) >= max_sol: break
        out = dict(stage="done", alpha_quads=len(Q), S_solutions=nS, stage1_solutions=len(sols),
                   all_verified=all(x['verified'] for x in sols), complete=len(sols) < max_sol, secs=round(time.time() - t0, 1))
        log(json.dumps(out))
        return out, sols


def row_quads_of(cfg, sl=None):
    sl = W.build(cfg) if sl is None else sl
    A = W.screen_A(sl, cfg)
    return A.get('row_quads', []), A


def planted(cfg, seed=0, log=print, sl=None, bias=0.85):
    """random U-invariant quadruple with the slice's letters (small letters favoured); targets = its own PAF;
    the exhaustive solver must return it (up to block order and sign)"""
    P0 = S1SD(cfg, sl=sl)
    sl = P0.sl; rng = np.random.default_rng(seed)
    def rnd():
        vals = []
        for L in sl.letters:
            Ls = sorted(L, key=abs); small = [x for x in Ls if abs(x) <= abs(Ls[min(3, len(Ls) - 1)])]
            vals.append(int(rng.choice(small if rng.random() < bias else Ls)))
        return np.array(vals)[sl.borb]
    Y = np.array([rnd() for _ in range(4)])
    T1 = P0.Gp.paf(Y).sum(0)
    P = S1SD(cfg, T1=T1, sl=sl)
    rq = [tuple(sorted(abs(int(r)) for r in Y.sum(1)))]
    out, sols = P.solve(rq, log=lambda s: None)
    canon = lambda Qd: sorted(tuple(int(x) * (1 if r[0] > 0 else -1) for x in r) for r in np.asarray(Qd).tolist())
    found = any(canon(s_['y']) == canon(Y) for s_ in sols)
    res = dict(seed=seed, found=bool(found), solutions=len(sols), alpha_quads=out['alpha_quads'], secs=out['secs'])
    log(json.dumps(res)); return res



if __name__ == "__main__":
    a = sys.argv[1:]
    cfg = cfg_from_json(a[0])
    max_sol = int(a[a.index('--max-sol') + 1]) if '--max-sol' in a else 10 ** 9
    if '--plant' in a:
        n = int(a[a.index('--plant') + 1])
        r = [planted(cfg, sd) for sd in range(n)]
        print(json.dumps(dict(planted=n, found=sum(x['found'] for x in r)))); sys.exit()
    P = S1SD(cfg)
    rq, A = row_quads_of(cfg)
    print(json.dumps(dict(cfg=dict(v=cfg['v'], Q=cfg['Q'], H1=cfg['H1'], q=cfg['q'], s=cfg['s'], U=[list(u) for u in cfg['U']]),
                          f0=P.f0, Gp=P.Gp.mods, Gpp_orbits=len(P.oreps), TS=P.TS[P.preps].tolist(), TD=P.TD[P.preps].tolist(),
                          row_quads=len(rq), A=A.get('A'))), flush=True)
    shuffle = int(a[a.index('--shuffle') + 1]) if '--shuffle' in a else 0
    out, sols = P.solve(rq, max_sol, shuffle=shuffle, progress=200)
    if '--out' in a:
        json.dump(dict(cfg=dict(cfg, U=[list(u) for u in cfg['U']]), result=out, solutions=sols), open(a[a.index('--out') + 1], 'w'))

