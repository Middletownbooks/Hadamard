r"""Stage-1 count of a GP-type slice conditioned on the Q-compressions [counts exact; remaining lags E].

Why: whole-slice and typical-y estimates of stage 1 sample y with the lag-0 energy only and can be wildly off: 599 at
c = 6 looked like +48 bits and holds about 2^4.9 (report/b599_answer.md, theory/b599_flevels.py).  This generalises
the z-conditioned count of theory/b599_flevels.py to any configuration of theory/ab_sweep.py with a Q part (nQ >= 1).

Setting: G' = Q x H' (Q = the first nQ axes, |Q| = c), y_k constant on G'-orbits with letters L_j, Q-compression
z_k(q) = sum of y_k over the Q-coset q.  Exact facts [P]:
  * summing the stage-1 identity over the lags with Q-part delta gives sum_k PAF_Q(z_k)(delta) = T_Q(delta), where
    T_Q(delta) = sum of T1 over the G'-lags with Q-part delta;
  * by Parseval on Q, sum_{psi != 1} |zhat_k(psi)|^2 = c sum_q (z_k(q) - zbar_k)^2, and every |zhat_k(psi)|^2 is at
    most That_Q(psi) (a sum of four non-negative terms), so sum_q (z_k(q) - zbar_k)^2 <= sum_{psi != 1} That_Q(psi) / c.
Method:
  1. per block, the exact number of letter assignments with a given (z, E) (E = lag-0 energy of the block): DP over
     the orbits of each connected component of Q-cosets (orbits can join cosets when U acts on Q);
  2. compressions z enumerated around each admissible row-sum mean with the Parseval bound, then the exact per-character
     bound; z-quadruples by a 2+2 join on (PAF_Q(z), |z|^2) with sum = T_Q, energies by exact convolution
     (sum_k E_k = T1(0)): N(z-quad) = exact number of y-quadruples meeting the Q-level and lag-0 equations;
  3. the remaining stage-1 equations (G'-lag classes with a nonzero H'-part) by the lattice Gaussian
     (theory/q523_gauss.gauss_bits) of uniform samples of y given (z_k, E_k);
  log2 E(stage 1 | z-quad) = log2 N - bits; the slice total is the sum over z-quadruples.
Calibration: --calibrate K picks K random z-quadruples and counts their stage-1 solutions exactly (every block
completion enumerated, then a hashed 2+2 join on the full PAF vectors), against the estimate.
usage: python3 theory/ab_s1_zcount.py CFG.json [--N 2000] [--max-quads 200] [--calibrate K] [--out FILE]"""
import sys, os, json, math, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs
from ab_lib import AG
from q523_gauss import gauss_bits


class ZCount:
    def __init__(self, cfg, sl=None):
        self.cfg = cfg; self.sl = sl = W.build(cfg) if sl is None else sl
        Gp = sl.Gp; self.Gp = Gp; nQ = cfg['nQ']; assert nQ >= 1
        self.Q = AG(Gp.mods[:nQ]); self.c = self.Q.n
        self.qfull = self.Q.idx(Gp.co[:, :nQ]); self.L1 = len(sl.breps); borb = np.asarray(sl.borb); self.borb = borb
        self.letters = [np.array(sorted(int(x) for x in L), np.int64) for L in sl.letters]
        self.prof = np.zeros((self.L1, self.c), np.int64); np.add.at(self.prof, (borb, self.qfull), 1)
        self.bsz = np.asarray(sl.bsz, np.int64)
        self.T1 = np.asarray(sl.T1, np.int64)
        Up = [Gp.act(u) for u in sl.Sp]
        self.reps, _ = AGs.lag_reps(Gp, Up)
        self.TQ = np.zeros(self.c, np.int64); np.add.at(self.TQ, self.qfull, self.T1)      # T_Q(delta)
        # remaining lag representatives: those with a nonzero H'-part (the Q-only lags are implied by z)
        Hpart = Gp.co[:, nQ:]
        self.rest = [c for c in self.reps if c != 0 and (Hpart[c] != 0).any()]
        # coset components (orbits joining cosets)
        par = list(range(self.c))
        def find(a):
            while par[a] != a: par[a] = par[par[a]]; a = par[a]
            return a
        for j in range(self.L1):
            qs = np.nonzero(self.prof[j])[0]
            for q in qs[1:]: par[find(int(q))] = find(int(qs[0]))
        comp = defaultdict(list)
        for q in range(self.c): comp[find(q)].append(q)
        self.comps = [sorted(v) for v in comp.values()]
        self.corbs = [[j for j in range(self.L1) if self.prof[j, C].sum() > 0] for C in self.comps]
        self.compdist = [self._comp_dp(C, O) for C, O in zip(self.comps, self.corbs)]

    def _comp_dp(self, C, orbs):
        """{(z restricted to C, E): count} over the letter assignments of the orbits of component C"""
        D = {(tuple(0 for _ in C), 0): 1}
        for j in orbs:
            pj = self.prof[j, C]; nd = defaultdict(int)
            for (zc, e), n in D.items():
                for y in self.letters[j]:
                    nd[(tuple(a + int(y) * int(b) for a, b in zip(zc, pj)), e + int(self.bsz[j]) * int(y) * int(y))] += n
            D = dict(nd)
        return D

    def block_z(self, dev_bound, rowsums):
        """compressions z with sum_q z(q) = r for an admissible block row sum r and sum_q (z(q) - zbar)^2 <= dev_bound:
        {z tuple: {E: count}} (each component's values are restricted to |z(q) - zbar| <= sqrt(dev_bound) first)"""
        per = []
        for D in self.compdist:
            m = defaultdict(lambda: defaultdict(int))
            for (zc, e), n in D.items(): m[zc][e] += n
            per.append({zc: dict(v) for zc, v in m.items()})
        out = {}; rad = math.sqrt(dev_bound) + 1e-9
        for r in sorted(set(rowsums)):
            zb = r / self.c
            cand = [{zc: ed for zc, ed in P.items() if all(abs(x - zb) <= rad for x in zc)} for P in per]
            if any(not c_ for c_ in cand): continue
            if np.prod([float(len(c_)) for c_ in cand]) > 2e7: raise RuntimeError('too many compression combinations')
            def rec(i, zparts, Edist, ssum):
                if i == len(cand):
                    if ssum != r: return
                    z = np.zeros(self.c, np.int64)
                    for C, zc in zip(self.comps, zparts): z[C] = zc
                    if ((z - zb) ** 2).sum() <= dev_bound + 1e-9: out[tuple(int(x) for x in z)] = Edist
                    return
                for zc, ed in cand[i].items():
                    nd = defaultdict(int)
                    for e1, n1 in Edist.items():
                        for e2, n2 in ed.items(): nd[e1 + e2] += n1 * n2
                    rec(i + 1, zparts + [zc], dict(nd), ssum + sum(zc))
            rec(0, [], {0: 1}, 0)
        return out

    def paf_q(self, z):
        z = np.asarray(z, np.int64); return np.array([int((z * z[self.Q.add(np.arange(self.c), np.full(self.c, d))]).sum()) for d in range(self.c)])

    def zquads(self, Z, max_quads=None):
        """ordered quadruples (z0..z3) of compressions with sum_k PAF_Q(z_k) = T_Q (all delta, incl. delta = 0)"""
        keys = list(Z); P = np.array([self.paf_q(z) for z in keys], np.int64)
        idx = {}
        for i in range(len(keys)):
            for j in range(len(keys)):
                idx.setdefault(tuple((P[i] + P[j]).tolist()), []).append((i, j))
        out = []
        for s, pairs in idx.items():
            other = tuple((self.TQ - np.array(s)).tolist())
            if other not in idx: continue
            for (i, j) in pairs:
                for (k, l) in idx[other]:
                    out.append((keys[i], keys[j], keys[k], keys[l]))
                    if max_quads and len(out) >= max_quads: return out
        return out

    def exact_q_count(self, zq, Z):
        """number of y-quadruples with these compressions and sum_k E_k = T1(0) (exact)"""
        T0 = int(self.T1[0]); acc = {0: 1}
        for z in zq:
            nd = defaultdict(int)
            for e1, n1 in acc.items():
                for e2, n2 in Z[z].items():
                    if e1 + e2 <= T0: nd[e1 + e2] += n1 * n2
            acc = dict(nd)
        return acc.get(T0, 0)

    def _suffix(self, i, zc):
        """suffix DP tables of component i for target compression zc (cached)"""
        key = (i, zc)
        if not hasattr(self, '_sufcache'): self._sufcache = {}
        if key in self._sufcache: return self._sufcache[key]
        C, orbs = self.comps[i], self.corbs[i]
        suf = [None] * (len(orbs) + 1); suf[-1] = {(tuple(0 for _ in C), 0): 1}
        for t in range(len(orbs) - 1, -1, -1):
            j = orbs[t]; pj = self.prof[j, C]; nd = defaultdict(int)
            for (zz, e), n in suf[t + 1].items():
                for y in self.letters[j]:
                    nd[(tuple(a + int(y) * int(b) for a, b in zip(zz, pj)), e + int(self.bsz[j]) * int(y) * int(y))] += n
            suf[t] = dict(nd)
        compE = {e: cnt for (zz, e), cnt in suf[0].items() if zz == zc}
        self._sufcache[key] = (suf, compE); return self._sufcache[key]

    def sample_block(self, z, E, rng):
        """one uniform y (per orbit) with compression z and energy E"""
        y = np.zeros(self.L1, np.int64); zcs = [tuple(int(z[q]) for q in C) for C in self.comps]
        tabs = [self._suffix(i, zc) for i, zc in enumerate(zcs)]
        # energy split across components (suffix convolution of the components' energy counts)
        suff = [None] * (len(tabs) + 1); suff[-1] = {0: 1}
        for i in range(len(tabs) - 1, -1, -1):
            nd = defaultdict(int)
            for e1, n1 in tabs[i][1].items():
                for e2, n2 in suff[i + 1].items(): nd[e1 + e2] += n1 * n2
            suff[i] = dict(nd)
        rem = E
        for i, (C, orbs) in enumerate(zip(self.comps, self.corbs)):
            opts = [(e1, n1 * suff[i + 1].get(rem - e1, 0)) for e1, n1 in tabs[i][1].items()]
            tot = sum(w for _, w in opts); r = rng.random() * tot; acc = 0; ei = None
            for e1, w in opts:
                acc += w
                if acc >= r and w > 0: ei = e1; break
            rem -= ei
            suf = tabs[i][0]; rz, re_ = list(zcs[i]), ei
            for t, j in enumerate(orbs):
                pj = self.prof[j, C]; ch = []
                for yv in self.letters[j]:
                    z2 = tuple(a - int(yv) * int(b) for a, b in zip(rz, pj)); e2 = re_ - int(self.bsz[j]) * int(yv) * int(yv)
                    ch.append((int(yv), suf[t + 1].get((z2, e2), 0)))
                tot = sum(w for _, w in ch); r = rng.random() * tot; acc = 0
                for yv, w in ch:
                    acc += w
                    if acc >= r and w > 0:
                        y[j] = yv; rz = [a - yv * int(b) for a, b in zip(rz, pj)]; re_ -= int(self.bsz[j]) * yv * yv; break
        return y

    def enumerate_block(self, z, cap=1 << 22):
        """all y (per orbit) with compression z (any energy), by depth-first search on the suffix tables"""
        zcs = [tuple(int(z[q]) for q in C) for C in self.comps]
        tabs = [self._suffix(i, zc) for i, zc in enumerate(zcs)]
        per = []
        for i, (C, orbs) in enumerate(zip(self.comps, self.corbs)):
            suf = tabs[i][0]; out = []
            def dfs(t, rz, re_, cur):
                if t == len(orbs):
                    out.append(list(cur)); return
                j = orbs[t]; pj = self.prof[j, C]
                for yv in self.letters[j]:
                    z2 = tuple(a - int(yv) * int(b) for a, b in zip(rz, pj))
                    # any energy: the suffix table must reach z2 with some energy
                    if any(zz == z2 for (zz, e) in suf[t + 1]):
                        cur.append(int(yv)); dfs(t + 1, z2, 0, cur); cur.pop()
            dfs(0, zcs[i], 0, [])
            per.append((orbs, out))
            if np.prod([float(len(o)) for _, o in per]) > cap: raise RuntimeError('block enumeration cap')
        Ys = [np.zeros(self.L1, np.int64)]
        for orbs, opts in per:
            new = []
            for y in Ys:
                for o in opts:
                    y2 = y.copy(); y2[orbs] = o; new.append(y2)
            Ys = new
        return np.array(Ys)

    def exact_count(self, zq):
        """exact number of y-quadruples with these compressions solving stage 1 (all lag classes; 2+2 join)"""
        lists = []
        for z in zq:
            Yb = self.enumerate_block(np.array(z))
            P = self.Gp.paf(Yb[:, self.borb])[:, self.reps]
            lists.append(P.astype(np.int64))
        T = self.T1[self.reps]
        assert self.reps[0] == 0                                  # coordinate 0 = lag 0 = block energy (ejoin4 splits on it)
        sols = SD.ejoin4(lists[0], lists[1], lists[2], lists[3], T, max_pairs=4 * 10 ** 8)
        return len(sols), [len(l) for l in lists]

    def estimate(self, zq, Z, N=2000, seed=1):
        """log2 of the expected number of y-quadruples with these compressions solving stage 1"""
        n = self.exact_q_count(zq, Z)
        if n == 0: return dict(log2N=None, log2E=None)
        rng = np.random.default_rng(seed); T0 = int(self.T1[0])
        # energy splits by exact weights, then blocks
        splits = []; W_ = []
        Edists = [Z[z] for z in zq]
        for e0, n0 in Edists[0].items():
            for e1, n1 in Edists[1].items():
                for e2, n2 in Edists[2].items():
                    e3 = T0 - e0 - e1 - e2
                    if e3 in Edists[3]: splits.append((e0, e1, e2, e3)); W_.append(n0 * n1 * n2 * Edists[3][e3])
        W_ = np.array(W_, float); W_ /= W_.sum()
        pick = rng.choice(len(splits), size=N, p=W_)
        V = np.zeros((N, len(self.rest)), np.int64)
        for t in range(N):
            es = splits[pick[t]]; acc = np.zeros(self.Gp.n, np.int64)
            for k in range(4):
                y = self.sample_block(np.array(zq[k]), es[k], rng)
                acc += self.Gp.paf(y[self.borb][None, :])[0]
            V[t] = acc[self.rest]
        gb = gauss_bits(V, self.T1[self.rest])
        return dict(log2N=round(math.log2(n), 2), bits=round(gb['bits'], 2) if math.isfinite(gb['bits']) else None,
                    rank=gb['rank'], log2E=round(math.log2(n) - gb['bits'], 2) if math.isfinite(gb['bits']) else None)


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read())
    N = int(a[a.index('--N') + 1]) if '--N' in a else 2000
    mq = int(a[a.index('--max-quads') + 1]) if '--max-quads' in a else 200
    t0 = time.time(); log = lambda **k: print(json.dumps(dict(k, secs=round(time.time() - t0, 1))), flush=True)
    Zc = ZCount(cfg)
    That = np.real(np.fft.fft(Zc.TQ)) if len(Zc.Q.mods) == 1 else None
    dev = float(sum(abs(x) for x in np.real(np.fft.fftn(Zc.TQ.reshape(Zc.Q.mods))).ravel()[1:])) / Zc.c
    rq, A = SD.row_quads_of(cfg, sl=Zc.sl)
    rows = sorted({int(s_) * int(x) for q in rq for x in q for s_ in (1, -1)})
    Z = Zc.block_z(dev, rows)
    log(stage='blocks', L1=Zc.L1, c=Zc.c, components=[len(C) for C in Zc.comps], rest_lags=len(Zc.rest), compressions=len(Z),
        dev_bound=round(dev, 2), TQ=Zc.TQ.tolist())
    Q4 = Zc.zquads(Z, max_quads=None)
    log(stage='zquads', ordered=len(Q4))
    tot = 0.0; recs = []
    for i, zq in enumerate(Q4[:mq]):
        r = Zc.estimate(zq, Z, N=N, seed=i + 1)
        if r['log2E'] is not None: tot += 2.0 ** r['log2E']
        recs.append(dict(z=[list(z) for z in zq], **r))
        log(stage='zquad', i=i, **r)
    log(stage='total_sampled_quads', n=min(mq, len(Q4)), log2_sum=round(math.log2(tot), 2) if tot > 0 else None,
        scaled_to_all=round(math.log2(tot * len(Q4) / min(mq, len(Q4))), 2) if tot > 0 and Q4 else None)
    if '--calibrate' in a:                       # random z-quadruples: exact count (2+2 join) against the estimate
        ncal = int(a[a.index('--calibrate') + 1]) if len(a) > a.index('--calibrate') + 1 and a[a.index('--calibrate') + 1].isdigit() else 20
        rng = np.random.default_rng(7)
        for i in rng.choice(len(Q4), size=min(ncal, len(Q4)), replace=False):
            zq = Q4[int(i)]
            try: ex, sizes = Zc.exact_count(zq)
            except RuntimeError as e: log(stage='calib', i=int(i), skipped=str(e)); continue
            r = Zc.estimate(zq, Z, N=N, seed=int(i) + 1)
            log(stage='calib', i=int(i), exact=ex, log2_exact=round(math.log2(ex), 2) if ex else None, block_lists=sizes, **r)
    if '--out' in a: json.dump(dict(cfg=a[0], quads=recs), open(a[a.index('--out') + 1], 'w'))


if __name__ == '__main__':
    main()
