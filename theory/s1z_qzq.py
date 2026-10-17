r"""Exact compression screens of stage 2 at the level of block TYPES, for any c (character-free) [P, M].

Setting: G = Q x H' x Z_q (Z_q the staging axis), a core x_k (k = 0..3), invariant under S = U x M.  For a quotient
Qbar of Q (keeping some of Q's axes; all of them = the Q x Z_q level, none = the Z_q level with the trivial character),
compress along the kernel and H':
    xbar_k(p, t) = sum of x_k over the points of G above (p, t) in Qbar x Z_q,    sum_k PAF_{Qbar x Z_q}(xbar_k) = T',
    T'(d, tau) = sum of T over the points of G above (d, tau)                                                     [P]
(sum the core identity over the fibres of G -> Qbar x Z_q; report/paley6_answer.md 3.1, 3.5, 3.6 for c = 4).
For a stage-1 solution y (letters y_kj on the G'-orbits j) the lifts are independent fibre patterns per orbit; S acts
trivially on the staging coordinate of G', so every point of orbit j carries the same pattern, and orbit j adds
n_{j,p} times its pattern to xbar_k(p, .) (n_{j,p} = number of points of orbit j above p).  So the distribution of
xbar_k over the lifts depends only on the multiset of (orbit class, letter) pairs, orbit class = (profile n_{j,.},
fibre parts): the block's TYPE.  The exact number of lift quadruples meeting the condition is a function of the type
quadruple; zero kills every stage-1 solution with that type quadruple [M].
Per type: an exact DP over its orbits (state = xbar at the representatives of the orbits of the S-action on
Qbar x Z_q, exact integer counts); PAF signatures at the lag-class representatives (lag 0 first); per type quadruple a
2+2 join keyed by (lag-0 energy, 64-bit linear hash of the signature), counts summed in float64: a hash collision can
only add to a count, so a zero count is exact.
usage (library): Screen(cfg, keep=[axes of Q kept]).count([y_0, y_1, y_2, y_3])"""
import sys, os, json, math, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_lib import AG


class Screen:
    def __init__(self, cfg, keep=None, sl=None):
        sl = W.build(cfg) if sl is None else sl; self.sl = sl; self.cfg = cfg
        G = sl.G; Gp = sl.Gp; nQ = cfg['nQ']; qn = sl.q
        assert sl.qaxis == len(G.mods) - 1
        keep = list(range(nQ)) if keep is None else list(keep); self.keep = keep
        self.Pmods = tuple(Gp.mods[a] for a in keep); self.P = AG(self.Pmods)
        self.PZ = AG(self.Pmods + (qn,)); PZ = self.PZ
        proj = PZ.idx(np.concatenate([G.co[:, keep], G.co[:, -1:]], 1)) if keep else G.co[:, -1]
        self.T_full = np.zeros(PZ.n, np.int64); np.add.at(self.T_full, proj, np.asarray(sl.T, np.int64))
        perms = {}
        for u in sl.S:
            uu = tuple(int(u[a]) for a in keep) + (int(u[-1]),)
            if uu not in perms: perms[uu] = PZ.act(uu)
        self.perms = list(perms.values())
        orb = -np.ones(PZ.n, np.int64); reps = []
        for x in range(PZ.n):
            if orb[x] >= 0: continue
            for Pm in self.perms: orb[Pm[x]] = len(reps)
            reps.append(x)
        self.orb = orb; self.oreps = np.array(reps)
        self.lreps, _ = AGs.lag_reps(PZ, self.perms)
        assert self.lreps[0] == 0
        self.T = self.T_full[self.lreps]
        L1 = len(sl.breps); self.L1 = L1
        pfull = self.P.idx(Gp.co[:, keep]) if keep else np.zeros(Gp.n, np.int64)
        prof = np.zeros((L1, self.P.n), np.int64); np.add.at(prof, (np.asarray(sl.borb), pfull), 1)
        self.M = []; self.ocls = []
        for j, b in enumerate(sl.breps):
            parts = [o for o, _ in sl.parts[j]]; pidx = {o: i for i, o in enumerate(parts)}
            Mj = np.zeros((len(parts), PZ.n), np.int64)
            for t in range(qn):
                pi = pidx[int(sl.orb[int(b) * qn + t])]
                for pp in np.nonzero(prof[j])[0]:
                    Mj[pi, PZ.idx(np.array(list(self.P.co[pp]) + [t]))] += prof[j, pp]
            self.M.append(Mj[:, self.oreps])
            self.ocls.append((tuple(prof[j].tolist()), tuple(s for _, s in sl.parts[j])))
        self._pat = {}; self._tcache = {}
        self.wv = np.random.default_rng(12345).integers(1, 2 ** 62, size=len(self.lreps), dtype=np.int64)
        self.th = int((self.T.astype(np.int64) @ self.wv))

    def set_target(self, TG):
        """replace the core target by TG (on G), e.g. a planted one"""
        sl = self.sl; G = sl.G; keep = self.keep
        proj = self.PZ.idx(np.concatenate([G.co[:, keep], G.co[:, -1:]], 1)) if keep else G.co[:, -1]
        self.T_full = np.zeros(self.PZ.n, np.int64); np.add.at(self.T_full, proj, np.asarray(TG, np.int64))
        self.T = self.T_full[self.lreps]; self.th = int((self.T.astype(np.int64) @ self.wv))

    def patterns(self, j, y):
        key = (j, int(y))
        if key not in self._pat:
            Pt, _ = AGs.fibre_patterns(self.sl, j, int(y))
            V = Pt.astype(np.int64) @ self.M[j]
            cnt = defaultdict(int)
            for v in map(tuple, V.tolist()): cnt[v] += 1
            self._pat[key] = dict(cnt)
        return self._pat[key]

    def type_key(self, yb):
        return tuple(sorted((self.ocls[j], int(yb[j])) for j in range(self.L1)))

    def block_dist(self, yb):
        """exact distribution of xbar (at the orbit reps) over the lifts of the block letters yb (per G'-orbit)"""
        D = {tuple([0] * len(self.oreps)): 1}
        for j in range(self.L1):
            cnt = self.patterns(j, yb[j]); nd = defaultdict(int)
            for s, n in D.items():
                for v, m in cnt.items(): nd[tuple(a + b for a, b in zip(s, v))] += n * m
            D = dict(nd)
        return D

    def sigs(self, yb):
        """(lag-0 energies, hashes, counts) of the block's lift distribution, aggregated by (energy, hash)"""
        key = self.type_key(yb)
        if key in self._tcache: return self._tcache[key]
        D = self.block_dist(yb)
        X = np.array(list(D.keys()), np.int64); C = np.array([float(n) for n in D.values()])
        Pf = self.PZ.paf(X[:, self.orb])[:, self.lreps]
        E = Pf[:, 0]; Hh = Pf @ self.wv
        agg = defaultdict(float)
        for e, h, c in zip(E.tolist(), Hh.tolist(), C.tolist()): agg[(e, h)] += c
        k = np.array(list(agg.keys()), np.int64).reshape(-1, 2); cn = np.array(list(agg.values()))
        self._tcache[key] = (k[:, 0], k[:, 1], cn)
        return self._tcache[key]

    def count(self, Y4):
        """count of lift quadruples meeting the condition (zero is exact); Y4: four blocks of letters per G'-orbit"""
        S = [self.sigs(yb) for yb in Y4]
        T0 = int(self.T[0])
        def pairs(a, b):
            Ea, Ha, Ca = a; Eb, Hb, Cb = b
            E = (Ea[:, None] + Eb[None, :]).ravel(); H = (Ha[:, None] + Hb[None, :]).ravel()
            C = (Ca[:, None] * Cb[None, :]).ravel(); ok = E <= T0
            return E[ok], H[ok], C[ok]
        EL, HL, CL = pairs(S[0], S[1]); ER, HR, CR = pairs(S[2], S[3])
        if not len(EL) or not len(ER): return 0.0
        kl = HL * np.int64(1000003) + EL
        kr = (np.int64(self.th) - HR) * np.int64(1000003) + (np.int64(T0) - ER)
        o = np.argsort(kl, kind='stable'); kls = kl[o]; cs = np.concatenate([[0.0], np.cumsum(CL[o])])
        lo = np.searchsorted(kls, kr, 'left'); hi = np.searchsorted(kls, kr, 'right')
        m = hi > lo
        if not m.any(): return 0.0
        return float(((cs[hi[m]] - cs[lo[m]]) * CR[m]).sum())

    def count_exact_small(self, Y4):
        """independent check for small cases: every lift of every block enumerated, PAF on Qbar x Z_q computed directly
        on the full vector, exact integer pair join (no hashing)"""
        PZ = self.PZ; T = self.T_full
        per = []
        for yb in Y4:
            lists = []
            for j in range(self.L1):
                Pt, _ = AGs.fibre_patterns(self.sl, j, int(yb[j]))
                lists.append(Pt.astype(np.int64) @ self.M[j])
            acc = defaultdict(int)
            for combo in itertools.product(*[range(len(l)) for l in lists]):
                v = sum(lists[j][i] for j, i in enumerate(combo))
                full = np.asarray(v)[self.orb]
                acc[tuple(int(x) for x in full)] += 1
            per.append(acc)
        def paf(full):
            F = np.array(full).reshape(PZ.mods)
            out = np.zeros(PZ.n, np.int64)
            for d in range(PZ.n):
                sh = tuple(int(x) for x in PZ.co[d])
                out[d] = int((F * np.roll(F, tuple(-s for s in sh), axis=tuple(range(len(PZ.mods))))).sum())
            return out
        sig = [{v: tuple(paf(v).tolist()) for v in acc} for acc in per]
        L = defaultdict(int)
        for v0, n0 in per[0].items():
            for v1, n1 in per[1].items():
                L[tuple(a + b for a, b in zip(sig[0][v0], sig[1][v1]))] += n0 * n1
        tot = 0
        for v2, n2 in per[2].items():
            for v3, n3 in per[3].items():
                need = tuple(int(t) - a - b for t, a, b in zip(T, sig[2][v2], sig[3][v3]))
                tot += n2 * n3 * L.get(need, 0)
        return tot
