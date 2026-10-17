r"""Stage 1 of a b >= 3 GP slice counted conditioned on the Q-compressions, over ALL compression quadruples up to
symmetry [counts exact; remaining lags E].

Setting (theory/ab_lib.py, report/alpoge_design_answer.md 1): G' = Q x H' (Q = the first nQ axes), y_k constant on the
G'-orbits of U with letters L_j, Q-compression z_k(q) = sum of y_k over the Q-coset q.  Exact facts [P]:
  * sum_k PAF_Q(z_k) = T_Q (the stage-1 identity summed over the lags with a given Q-part), hence sum_k r_k^2 = 4v;
  * per block |zhat_k(psi)|^2 <= That_Q(psi) for psi != 1 (a sum of four non-negative terms), so
    sum_q (z_k(q) - zbar_k)^2 <= sum_{psi != 1} That_Q(psi) / c (Parseval);
  * the stage-1 solution set is closed under, per block, y -> -y, y -> y(. - t) for t in G' fixed by U, and
    y -> y(-.); and globally (all blocks at once) under the automorphisms of G' that preserve H' and commute with U
    (here: those induced by Aut(Q), acting as the identity on H'); and under permutations of the blocks.
    All of them map compression quadruples to compression quadruples and keep every count below.
Method:
  1. per coset component (cosets joined by U-orbits), the exact number of letter assignments with a given
     (component compression, energy) by DP over its orbits;
  2. rows: every achievable row sum r with r^2 <= 4v occurring in some quadruple with sum r_k^2 = 4v (no truncation);
     compressions: every z with these rows, the Parseval bound and |zhat(psi)|^2 <= That_Q(psi), each with its exact
     energy distribution;
  3. classes of compressions under the per-block group; class quadruples (multisets) with sum PAF_Q = T_Q by a 2+2
     join; their exact lag-0 count N0 = #{y-quadruples with these compressions, sum_k E_k = T1(0)} (exact integers);
     the number of ordered compression quadruples each multiset stands for; orbits of multisets under Aut(Q);
     N0 = 0 for every multiset (or no multiset) => stage 1 empty [M];
  4. per orbit representative: the lattice Gaussian (theory/q523_gauss.gauss_bits) of the remaining lag classes (those
     with a nonzero H'-part) under y-quadruples drawn uniformly from the N0 set:
     E(stage 1 | multiset) = N0 * 2^-bits; optional finite-sample correction of the log-det bias;
  5. exact counts per multiset for calibration (every block completion with the compression enumerated, energy-split
     hashed 2+2 join on the full PAF vectors (theory/ab_s1sd.ejoin4), every solution re-checked).
usage: see theory/s1z_run.py"""
import sys, os, json, math, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W
from ab_lib import AG
from ab_gauss import lag_reps
from q523_gauss import gauss_bits


def conv_dict(A, B, cap=None):
    out = defaultdict(int)
    for a, x in A.items():
        for b, y in B.items():
            if cap is None or a + b <= cap: out[a + b] += x * y
    return dict(out)


class ZOrb:
    def __init__(self, cfg, sl=None):
        self.cfg = cfg; sl = W.build(cfg) if sl is None else sl; self.sl = sl
        Gp = sl.Gp; self.Gp = Gp; nQ = cfg['nQ']; assert nQ >= 1
        self.v, self.b = cfg['v'], cfg['b']
        self.Q = AG(Gp.mods[:nQ]); self.c = self.Q.n
        self.qfull = self.Q.idx(Gp.co[:, :nQ])
        self.L1 = len(sl.breps); self.borb = np.asarray(sl.borb)
        self.letters = [sorted(int(x) for x in L) for L in sl.letters]
        self.w = np.asarray(sl.bsz, np.int64)
        self.prof = np.zeros((self.L1, self.c), np.int64); np.add.at(self.prof, (self.borb, self.qfull), 1)
        self.T1 = np.asarray(sl.T1, np.int64); self.T0 = int(self.T1[0])
        self.Up = [Gp.act(u) for u in sl.Sp]
        self.reps, self.lcls = lag_reps(Gp, self.Up)
        self.TQ = np.zeros(self.c, np.int64); np.add.at(self.TQ, self.qfull, self.T1)
        Hpart = Gp.co[:, nQ:]
        self.rest = [g for g in self.reps if g != 0 and (Hpart[g] != 0).any()]
        self.That = np.rint(np.real(np.fft.fftn(self.TQ.reshape(self.Q.mods)))).ravel().astype(np.int64)
        # coset components
        par = list(range(self.c))
        def find(a):
            while par[a] != a: par[a] = par[par[a]]; a = par[a]
            return a
        for j in range(self.L1):
            qs = np.nonzero(self.prof[j])[0]
            for q in qs[1:]: par[find(int(q))] = find(int(qs[0]))
        comp = defaultdict(list)
        for q in range(self.c): comp[find(q)].append(q)
        self.comps = [sorted(x) for x in comp.values()]
        self.corbs = [[j for j in range(self.L1) if self.prof[j, C].sum() > 0] for C in self.comps]
        self.Emin = int(sum(int(self.w[j]) * min(y * y for y in self.letters[j]) for j in range(self.L1)))
        self.Ecap = self.T0 - 3 * self.Emin                     # largest energy a block can have
        self._fwd = [self._forward(i) for i in range(len(self.comps))]
        # per-block symmetry on Q: translations by Q-projections of U-fixed elements, reversal, negation
        fixed = [g for g in range(Gp.n) if all(P[g] == g for P in self.Up)]
        self.tQ = sorted({int(self.qfull[g]) for g in fixed})
        self.qneg = self.Q.idx(-self.Q.co)
        self.qshift = {t: self.Q.idx(self.Q.co - self.Q.co[t]) for t in self.tQ}   # z(. - t) = z[qshift[t]]
        self.autQ = self._autQ()

    # ------------------------------------------------------------------ exact per-component tables
    def _forward(self, i):
        """forward DP over the orbits of component i: list of dicts (zc tuple, E) -> count after t orbits"""
        C, orbs = self.comps[i], self.corbs[i]
        tabs = [{(tuple(0 for _ in C), 0): 1}]
        for j in orbs:
            pj = [int(x) for x in self.prof[j, C]]; wj = int(self.w[j]); nd = defaultdict(int)
            for (zc, e), n in tabs[-1].items():
                for y in self.letters[j]:
                    e2 = e + wj * y * y
                    if e2 > self.Ecap: continue
                    nd[(tuple(a + y * b for a, b in zip(zc, pj)), e2)] += n
            tabs.append(dict(nd))
        return tabs

    def comp_values(self, i):
        """{zc: {E: count}} of component i"""
        out = defaultdict(dict)
        for (zc, e), n in self._fwd[i][-1].items(): out[zc][e] = n
        return dict(out)

    def _autQ(self):
        Q = self.Q; mods = Q.mods; perms = set()
        if all(m == 2 for m in mods):
            r = len(mods)
            for M in itertools.product((0, 1), repeat=r * r):
                A = np.array(M).reshape(r, r)
                if round(abs(np.linalg.det(A))) % 2 == 0: continue
                perms.add(tuple(Q.idx((Q.co @ A.T) % 2).tolist()))
        else:
            opts = [[u for u in range(1, m) if math.gcd(u, m) == 1] if m > 1 else [0] for m in mods]
            for us in itertools.product(*opts): perms.add(tuple(Q.act(us).tolist()))
        # keep the ones that map U-invariant compressions to U-invariant compressions and permute the components
        compset = {tuple(C) for C in self.comps}
        ok = []
        for p in perms:
            p = np.array(p)
            if all(tuple(sorted(int(p[q]) for q in C)) in compset for C in self.comps): ok.append(p)
        return ok

    # ------------------------------------------------------------------ rows and compressions
    def rows(self):
        """achievable block row sums r (r^2 <= 4v) that occur in a quadruple with sum r_k^2 = 4v"""
        R4 = 4 * self.v
        sums = {0}
        for i in range(len(self.comps)):
            vals = {sum(zc) for zc in self.comp_values(i)}
            sums = {a + b for a in sums for b in vals if abs(a + b) <= math.isqrt(R4) + 200}
        ach = sorted(r for r in sums if r * r <= R4)
        sq = sorted({r * r for r in ach})
        ok2 = set()
        pair = defaultdict(list)
        for a in sq:
            for b in sq:
                if a + b <= R4: pair[a + b].append((a, b))
        good = set()
        for s, lst in pair.items():
            if R4 - s in pair:
                for (a, b) in lst: good.add(a); good.add(b)
        return [r for r in ach if r * r in good]

    def compressions(self, rows=None):
        """{z tuple: {E: count}} over the admissible rows (Parseval + per-character bound)"""
        rows = self.rows() if rows is None else rows
        c = self.c
        dev = float(sum(self.That[1:])) / c
        rad = math.sqrt(dev) + 1e-9
        cvals = [self.comp_values(i) for i in range(len(self.comps))]
        out = {}
        bound = self.That.astype(float) + 1e-6
        for r in rows:
            zb = r / c
            cand = [[(zc, ed) for zc, ed in cv.items() if all(abs(x - zb) <= rad for x in zc)] for cv in cvals]
            if any(not x for x in cand): continue
            # recursive enumeration with partial-deviation pruning
            order = sorted(range(len(cand)), key=lambda i: len(cand[i]))
            def rec(t, parts, s, dv):
                if t == len(order):
                    if s != r: return
                    z = np.zeros(c, np.int64)
                    for i, (zc, ed) in parts: z[self.comps[i]] = zc
                    P = np.abs(np.fft.fftn(z.reshape(self.Q.mods).astype(float))).ravel() ** 2
                    if (P[1:] <= bound[1:]).all():
                        Ed = {0: 1}
                        for i, (zc, ed) in parts: Ed = conv_dict(Ed, ed, self.Ecap)
                        if Ed: out[tuple(int(x) for x in z)] = Ed
                    return
                i = order[t]
                for zc, ed in cand[i]:
                    d2 = dv + sum((x - zb) ** 2 for x in zc)
                    if d2 > dev + 1e-9: continue
                    rec(t + 1, parts + [(i, (zc, ed))], s + sum(zc), d2)
            rec(0, [], 0, 0.0)
        return out

    def paf_q(self, z):
        z = np.asarray(z, np.int64)
        return tuple(int((z * z[self.Q.add(np.arange(self.c), np.full(self.c, d))]).sum()) for d in range(self.c))

    # ------------------------------------------------------------------ classes and quadruples
    def images(self, z):
        z = np.asarray(z, np.int64); out = set()
        for zz in (z, z[self.qneg]):
            for t in self.tQ:
                y = zz[self.qshift[t]]
                out.add(tuple(int(x) for x in y)); out.add(tuple(int(-x) for x in y))
        return out

    def classes(self, Z):
        """list of dicts: rep, size, sig, Edist; and a map z -> class index"""
        cls = []; idx = {}
        for z in Z:
            if z in idx: continue
            im = self.images(z); im = [x for x in im if x in Z]
            rep = min(im)
            k = len(cls)
            for x in im: idx[x] = k
            cls.append(dict(rep=rep, size=len(im), sig=self.paf_q(rep), E=Z[rep]))
        return cls, idx

    def class_quads(self, cls, max_multisets=None):
        """multisets (i <= j <= k <= l) of classes with sum PAF_Q = T_Q, with the number of ordered compression
        quadruples they stand for"""
        sig = [np.array(c_['sig'], np.int64) for c_ in cls]
        pairs = defaultdict(list)
        n = len(cls)
        for i in range(n):
            for j in range(i, n):
                pairs[tuple((sig[i] + sig[j]).tolist())].append((i, j))
        TQ = self.TQ
        ms = set()
        for key, lst in pairs.items():
            other = tuple((TQ - np.array(key)).tolist())
            if other not in pairs or other < key: continue
            for (i, j) in lst:
                for (k, l) in pairs[other]:
                    ms.add(tuple(sorted((i, j, k, l))))
                    if max_multisets and len(ms) > max_multisets: raise RuntimeError('too many multisets')
        out = []
        for m in sorted(ms):
            cnt = defaultdict(int)
            for x in m: cnt[x] += 1
            perm = math.factorial(4)
            for x in cnt.values(): perm //= math.factorial(x)
            out.append(dict(m=m, ordered=perm * math.prod(cls[x]['size'] for x in m)))
        return out

    def n0(self, cls, m, cache=None):
        """exact number of y-quadruples with these compressions (one per class, any order) and sum E_k = T1(0)"""
        cache = {} if cache is None else cache
        a, b, c_, d = m
        def pc(i, j):
            key = (i, j)
            if key not in cache: cache[key] = conv_dict(cls[i]['E'], cls[j]['E'], self.T0)
            return cache[key]
        L = pc(a, b); R = pc(c_, d)
        return sum(n * R.get(self.T0 - e, 0) for e, n in L.items())

    def aut_orbits(self, cls, idx, quads):
        """group the multisets into orbits of Aut(Q) (acting on all blocks at once)"""
        cmap = []
        for p in self.autQ:
            mp = []
            for c_ in cls:
                im = tuple(int(x) for x in np.asarray(c_['rep'])[p])
                mp.append(idx.get(im, -1))
            cmap.append(mp)
        mkey = {}
        for t, q in enumerate(quads):
            if t in mkey: continue
            imgs = set()
            for mp in cmap:
                im = tuple(sorted(mp[x] for x in q['m']))
                if -1 not in im: imgs.add(im)
            for im in imgs: mkey[im] = t
            mkey[q['m']] = t
        orb = defaultdict(list)
        for t, q in enumerate(quads): orb[mkey.get(q['m'], t)].append(t)
        return orb

    # ------------------------------------------------------------------ sampling
    def _comp_backward(self, i, zc, e, rng):
        """uniform letter assignment of component i's orbits with compression zc and energy e"""
        C, orbs = self.comps[i], self.corbs[i]; tabs = self._fwd[i]
        out = {}; rz = list(zc); re_ = e
        for t in range(len(orbs) - 1, -1, -1):
            j = orbs[t]; pj = [int(x) for x in self.prof[j, C]]; wj = int(self.w[j])
            ys = self.letters[j]; ws = []
            for y in ys:
                key = (tuple(a - y * b for a, b in zip(rz, pj)), re_ - wj * y * y)
                ws.append(tabs[t].get(key, 0))
            tot = sum(ws); r = rng.random() * tot; acc = 0
            for y, wgt in zip(ys, ws):
                acc += wgt
                if acc > r and wgt > 0: break
            out[j] = y; rz = [a - y * b for a, b in zip(rz, pj)]; re_ -= wj * y * y
        assert re_ == 0 and not any(rz)
        return out

    def block_sampler(self, z):
        """returns (Edist, draw(E, rng) -> y per orbit) for compression z"""
        z = np.asarray(z, np.int64)
        zcs = [tuple(int(z[q]) for q in C) for C in self.comps]
        per = []
        for i, zc in enumerate(zcs):
            ed = {e: n for (zz, e), n in self._fwd[i][-1].items() if zz == zc}
            per.append(ed)
        suff = [None] * (len(per) + 1); suff[-1] = {0: 1}
        for i in range(len(per) - 1, -1, -1): suff[i] = conv_dict(per[i], suff[i + 1], self.Ecap)
        def draw(E, rng):
            y = np.zeros(self.L1, np.int64); rem = E
            for i in range(len(per)):
                opts = [(e1, n1 * suff[i + 1].get(rem - e1, 0)) for e1, n1 in per[i].items()]
                tot = sum(wt for _, wt in opts); r = rng.random() * tot; acc = 0
                for e1, wt in opts:
                    acc += wt
                    if acc > r and wt > 0: break
                for j, yv in self._comp_backward(i, zcs[i], e1, rng).items(): y[j] = yv
                rem -= e1
            assert rem == 0
            return y
        return suff[0], draw

    def sample_quads(self, zq, N, rng):
        """N y-quadruples drawn uniformly from those with compressions zq and sum_k E_k = T1(0)"""
        samp = [self.block_sampler(z) for z in zq]
        Ns = [s[0] for s in samp]
        P01 = defaultdict(list)
        for e0, c0 in Ns[0].items():
            for e1, c1 in Ns[1].items():
                if e0 + e1 <= self.T0: P01[e0 + e1].append((e0, e1, c0 * c1))
        splits = []; wts = []
        for e2, c2 in Ns[2].items():
            for e3, c3 in Ns[3].items():
                for (e0, e1, c01) in P01.get(self.T0 - e2 - e3, []):
                    splits.append((e0, e1, e2, e3)); wts.append(c01 * c2 * c3)
        tot = sum(wts)
        if not tot: return None, 0
        p = np.array([x / tot for x in wts], float); p /= p.sum()
        pick = rng.choice(len(splits), size=N, p=p)
        Y = np.zeros((N, 4, self.L1), np.int64)
        for t, si in enumerate(pick):
            for k in range(4): Y[t, k] = samp[k][1](splits[si][k], rng)
        return Y, tot

    def estimate(self, zq, N=2000, seed=1, bias_correct=True):
        rng = np.random.default_rng(seed)
        Y, n0 = self.sample_quads(zq, N, rng)
        if not n0: return dict(log2N0=None, log2E=None)
        Pf = self.Gp.paf(Y.reshape(4 * N, self.L1)[:, self.borb]).reshape(N, 4, -1).sum(1)
        V = Pf[:, self.rest]
        gb = gauss_bits(V, self.T1[self.rest])
        out = dict(log2N0=round(math.log2(n0), 3), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'),
                   mahal=gb.get('mahal'))
        b = gb.get('bits')
        if b is None or not np.isfinite(b):
            out['bits'] = None if b is None else float(b); out['log2E'] = -math.inf if gb.get('in_lattice') is False else None
            return out
        out['bits_raw'] = float(b)
        if bias_correct and gb.get('rank'):
            out['bits'] = round(float(b) + logdet_bias_bits(gb['rank'], N), 3)
        else: out['bits'] = float(b)
        out['log2E'] = round(math.log2(n0) - out['bits'], 3)
        return out

    # ------------------------------------------------------------------ exact counts (calibration)
    def enumerate_block(self, z, Ecap, cap=1 << 23):
        """all y (per orbit) with compression z and energy <= Ecap"""
        z = np.asarray(z, np.int64)
        per = []
        for i, C in enumerate(self.comps):
            zc = tuple(int(z[q]) for q in C); orbs = self.corbs[i]; tabs = self._fwd[i]
            emin_pre = [min(e for (_, e) in tb) for tb in tabs]
            zsets = [{k[0] for k in tb} for tb in tabs]
            res = []
            def dfs(t, rz, cur, e):
                if t == 0:
                    if not any(rz): res.append((list(reversed(cur)), e))
                    return
                j = orbs[t - 1]; pj = [int(x) for x in self.prof[j, C]]; wj = int(self.w[j])
                for y in self.letters[j]:
                    rz2 = tuple(a - y * b for a, b in zip(rz, pj)); e2 = e + wj * y * y
                    if e2 + emin_pre[t - 1] > Ecap: continue
                    if rz2 not in zsets[t - 1]: continue
                    cur.append(y); dfs(t - 1, rz2, cur, e2); cur.pop()
            dfs(len(orbs), zc, [], 0)
            per.append((orbs, res))
        tot = math.prod(len(r) for _, r in per)
        if tot > cap: raise RuntimeError(f'block enumeration cap ({tot})')
        Ys = [(np.zeros(self.L1, np.int64), 0)]
        for orbs, res in per:
            new = []
            for y, e in Ys:
                for vals, e1 in res:
                    if e + e1 > Ecap: continue
                    y2 = y.copy(); y2[orbs] = vals; new.append((y2, e + e1))
            Ys = new
        return np.array([y for y, _ in Ys], np.int64).reshape(-1, self.L1)

    def exact_count(self, zq, max_pairs=6 * 10 ** 8):
        """exact number of ordered y-quadruples with compressions zq (in this order) solving stage 1"""
        from ab_s1sd import ejoin4
        mins = []
        lists = []
        for z in zq:
            Ed = self.block_sampler(z)[0]; mins.append(min(Ed) if Ed else None)
        if None in mins: return 0, [0, 0, 0, 0], []
        for k, z in enumerate(zq):
            cap = self.T0 - (sum(mins) - mins[k])
            Yb = self.enumerate_block(z, cap)
            P = self.Gp.paf(Yb[:, self.borb])[:, self.reps].astype(np.int64)
            lists.append((Yb, P))
        T = self.T1[self.reps]
        sols = ejoin4(lists[0][1], lists[1][1], lists[2][1], lists[3][1], T, max_pairs=max_pairs)
        Ys = []
        for (a, b, c_, d) in sols:
            Y = np.array([lists[0][0][a], lists[1][0][b], lists[2][0][c_], lists[3][0][d]])
            assert (self.Gp.paf(Y[:, self.borb]).sum(0) == self.T1).all()
            Ys.append(Y)
        return len(sols), [len(l[0]) for l in lists], Ys


def logdet_bias_bits(d, N):
    """expected underestimate (in bits of -log2 density) of the log-det term from a sample covariance of N draws in d
    dimensions (Wishart: E ln det S = ln det Sigma + sum_i psi((N - i)/2) + d ln(2/(N - 1)))"""
    from math import log
    try:
        from scipy.special import digamma
    except Exception:
        return 0.0
    s = sum(digamma((N - i) / 2.0) for i in range(1, d + 1)) + d * log(2.0 / (N - 1))
    return float(-0.5 * s / log(2.0))


class SDCount:
    """exact stage-1 count for a prescribed ordered compression quadruple by the S/D split on axis 0 of G' (= Z_2, U
    trivial on it) [P, M]:  sigma = (y(0, .) + y(1, .))/2, delta = (y(0, .) - y(1, .))/2 on G'' satisfy
    sum_k PAF_{G''}(sigma_k) = T_S, sum_k PAF_{G''}(delta_k) = T_D (theory/ab_s1sd.py), and the compressions split as
    zs = (z(0, .) + z(1, .))/2, zd = (z(0, .) - z(1, .))/2 on Q'' = Q / Z_2.  Every block's sigma-list (energy-capped,
    compression zs) is enumerated; S-solutions by the energy-split 2+2 join; per S-solution the delta-lists
    (compatible with sigma and zd) and the D-join; every solution is re-checked on G'."""
    def __init__(self, Z):
        self.Z = Z; Gp = Z.Gp; assert Gp.mods[0] == 2
        self.n2 = Gp.n // 2; n2 = self.n2
        self.Gpp = AG(Gp.mods[1:])
        Spp = sorted({tuple(u[1:]) for u in Z.sl.Sp})
        self.Upp = [self.Gpp.act(u) for u in Spp]
        Oo = -np.ones(n2, np.int64); reps = []
        for b in range(n2):
            if Oo[b] >= 0: continue
            for P in self.Upp: Oo[P[b]] = len(reps)
            reps.append(b)
        self.Oo = Oo; self.oreps = reps; self.ow = np.bincount(Oo)
        self.j0 = [int(Z.borb[b]) for b in reps]; self.j1 = [int(Z.borb[n2 + b]) for b in reps]
        nQ = Z.cfg['nQ']
        self.Qpp = AG(Gp.mods[1:nQ]); cpp = self.Qpp.n; self.cpp = cpp
        qpp = self.Qpp.idx(self.Gpp.co[:, :nQ - 1]) if nQ > 1 else np.zeros(n2, np.int64)
        self.pro = np.zeros((len(reps), cpp), np.int64); np.add.at(self.pro, (Oo, qpp), 1)
        self.preps, _ = lag_reps(self.Gpp, self.Upp)
        T1 = Z.T1
        self.TS = np.array([(T1[b] + T1[n2 + b]) // 4 for b in range(n2)], np.int64)
        self.TD = np.array([(T1[b] - T1[n2 + b]) // 4 for b in range(n2)], np.int64)
        # per G''-orbit: sigma -> list of delta
        self.sd = []
        for o in range(len(reps)):
            m = defaultdict(set)
            for a in Z.letters[self.j0[o]]:
                for b in Z.letters[self.j1[o]]:
                    if (a + b) % 2 == 0: m[(a + b) // 2].add((a - b) // 2)
            self.sd.append({s: sorted(d) for s, d in m.items()})

    def _vectors(self, opts, target, Ecap):
        """all vectors (per G''-orbit value from opts[o]) with compression `target` on Q'' and energy <= Ecap"""
        nr = len(self.oreps); w = self.ow; pro = self.pro
        order = list(range(nr))
        minrest = np.zeros(nr + 1, np.int64)
        for t in range(nr - 1, -1, -1): minrest[t] = minrest[t + 1] + w[order[t]] * min(x * x for x in opts[order[t]])
        # remaining achievable compression ranges for pruning
        lo = np.zeros((nr + 1, self.cpp), np.int64); hi = np.zeros((nr + 1, self.cpp), np.int64)
        for t in range(nr - 1, -1, -1):
            o = order[t]; mn = min(opts[o]); mx = max(opts[o])
            lo[t] = lo[t + 1] + pro[o] * mn; hi[t] = hi[t + 1] + pro[o] * mx
        out = []; cur = [0] * nr; tgt = np.asarray(target, np.int64)
        def rec(t, acc, e):
            if t == nr:
                if (acc == tgt).all(): out.append(list(cur))
                return
            o = order[t]
            for x in opts[o]:
                e2 = e + int(w[o]) * x * x
                if e2 + minrest[t + 1] > Ecap: continue
                a2 = acc + pro[o] * x
                rem = tgt - a2
                if (rem < lo[t + 1]).any() or (rem > hi[t + 1]).any(): continue
                cur[o] = x; rec(t + 1, a2, e2)
        rec(0, np.zeros(self.cpp, np.int64), 0)
        return np.array(out, np.int64).reshape(-1, nr)

    def count(self, zq, max_pairs=4 * 10 ** 8, keep=False):
        from ab_s1sd import ejoin4
        Z = self.Z; cq = Z.c; cpp = self.cpp
        zq = [np.asarray(z, np.int64) for z in zq]
        zs = []; zd = []
        for z in zq:
            z0 = z[:cpp]; z1 = z[cpp:]                      # Q index = e * |Q''| + q''
            if ((z0 + z1) % 2).any(): return 0, []
            zs.append((z0 + z1) // 2); zd.append((z0 - z1) // 2)
        TS0 = int(self.TS[0]); TD0 = int(self.TD[0])
        sopts = [sorted(self.sd[o]) for o in range(len(self.oreps))]
        Sl = [self._vectors(sopts, zs[k], TS0) for k in range(4)]
        if min(len(s) for s in Sl) == 0: return 0, []
        P = [self.Gpp.paf(S[:, self.Oo])[:, self.preps].astype(np.int64) for S in Sl]
        Ssol = ejoin4(P[0], P[1], P[2], P[3], self.TS[self.preps], max_pairs=max_pairs)
        total = 0; sols = []
        dcache = {}
        for (a, b, c_, d) in Ssol:
            sig = [Sl[0][a], Sl[1][b], Sl[2][c_], Sl[3][d]]
            Dl = []; DP = []
            for k in range(4):
                key = (k, tuple(sig[k].tolist()))
                if key not in dcache:
                    dopts = [self.sd[o][int(sig[k][o])] for o in range(len(self.oreps))]
                    Dv = self._vectors(dopts, zd[k], TD0)
                    dcache[key] = (Dv, self.Gpp.paf(Dv[:, self.Oo])[:, self.preps].astype(np.int64) if len(Dv) else None)
                Dl.append(dcache[key][0]); DP.append(dcache[key][1])
            if min(len(x) for x in Dl) == 0: continue
            Dsol = ejoin4(DP[0], DP[1], DP[2], DP[3], self.TD[self.preps], max_pairs=max_pairs)
            for (e, f, g, h) in Dsol:
                dl = [Dl[0][e], Dl[1][f], Dl[2][g], Dl[3][h]]
                Y = np.zeros((4, Z.Gp.n), np.int64)
                for k in range(4):
                    s_ = sig[k][self.Oo]; d_ = dl[k][self.Oo]
                    Y[k, :self.n2] = s_ + d_; Y[k, self.n2:] = s_ - d_
                assert (Z.Gp.paf(Y).sum(0) == Z.T1).all()
                total += 1
                if keep: sols.append(Y)
        return total, sols


def log2_symmetry(cfg, sl=None):
    """log2 of the order of the symmetry group acting on ordered stage-1 solutions (generic orbit size) [P]:
    per block, translations by the U-fixed points of G', the sign, and the reversal g -> -g when it is not already a
    multiplier of the slice; the 24 block orders; the per-axis unit scalars of G' modulo <U, -1> and, for Q = Z_2^r,
    GL(r, 2) on Q (they preserve H', commute with U and keep T1).  Also returns s (the multipliers of Z_q outside M,
    which fix y and permute its lifts)."""
    import ab_sweep as W_
    sl = W_.build(cfg) if sl is None else sl
    Gp = sl.Gp; mods = Gp.mods
    Up = [Gp.act(u) for u in sl.Sp]
    nfix = sum(1 for g in range(Gp.n) if all(P[g] == g for P in Up))
    neg = tuple((m - 1) % m if m > 1 else 0 for m in mods)
    Sred = {tuple(int(a) % m if m > 1 else 0 for a, m in zip(u, mods)) for u in sl.Sp}
    rev_trivial = neg in Sred
    # <U, -1> as per-axis scalars
    grp = set(Sred) | {tuple((a * b) % m if m > 1 else 0 for a, b, m in zip(u, neg, mods)) for u in Sred}
    units = 1
    for m in mods: units *= sum(1 for x in range(1, m) if math.gcd(x, m) == 1) if m > 1 else 1
    gl = 1
    nQ = cfg['nQ']
    if all(m == 2 for m in mods[:nQ]):
        r = nQ; gl = 1
        for i in range(r): gl *= (2 ** r - 2 ** i)
    per_block = nfix * 2 * (1 if rev_trivial else 2)
    lg = 4 * math.log2(per_block) + math.log2(24) + math.log2(units / len(grp)) + math.log2(gl)
    return dict(log2_sym=round(lg, 3), fix=nfix, reversal=not rev_trivial, global_auts=units // len(grp) * gl, s=cfg['s'])
