r"""Hierarchical stage-2 solver for GP-type slices (any b), the 523 method generalised [P for the identities, M for runs].

Setting: G = G' x Z_q, G' = Z_2 x G'' (axis 0 of G' is the S/D axis), slice group S = U x M (U on G'', M <= F_q^* of
index s on Z_q).  Given a stage-1 solution y, every G'-orbit j of every block k carries a fibre pattern
p = (alpha, c_0, ..., c_{s-1}) (signs on the parts of the fibre, signed size-sum y_k(b_j)).
Identities [P] (x_S = x(0, .) + x(1, .), x_D = x(0, .) - x(1, .) on G'' x Z_q):
    PAF_{G''xZq}(x_S)(g) = PAF_G(x)(0, g) + PAF_G(x)(1, g),   PAF_{G''xZq}(x_D)(g) = PAF_G(x)(0, g) - PAF_G(x)(1, g),
so the core condition sum_k PAF_G(x_k) = T is equivalent to
    (S) sum_k PAF(x_{S,k}) = T_S := T(0, .) + T(1, .)     and     (D) sum_k PAF(x_{D,k}) = T_D := T(0, .) - T(1, .).
Given x_S on a fibre pair (orbit o of G''), the pairs (p0, p1) with p0 + p1 = s_o are the D-options (d = p0 - p1).
Middle level: L = an axis of G'' (U-stable); x_M = compression of x_S along L on (G''/L) x Z_q; its PAF is the sum of
the S-level PAF over L-cosets of lags (characters trivial on Z_2 x L), so (M) sum_k PAF(x_{M,k}) = T_M is implied by (S).
Search: parents = per-block x_M values; for a middle value c, left = parent pairs of blocks (0, 1) with
PAF-signature sum c, right = parent pairs of blocks (2, 3) with sum T_M - c; fibre products are joined on the full
S-signature; every S-solution gets an exact D-check (2+2 join over its D-options), and every D-solution is assembled
and verified on G.
usage: python3 theory/ab_stage2.py CFGJSON YSOURCE --stats      (sizes of rows, parents, fibres)"""
import sys, os, json, math, itertools, time
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_lib import AG
from ab_s1sd import cfg_from_json


class Stage2:
    def __init__(self, cfg, Y, Laxis):
        """Laxis: index of the middle-level axis among the axes of G' (>= 1)"""
        self.cfg = cfg; sl = W.build(cfg); self.sl = sl
        self.Y = np.asarray(Y)
        G, Gp = sl.G, sl.Gp
        assert (Gp.paf(self.Y).sum(0) == sl.T1).all(), "y is not a stage-1 solution"
        assert Gp.mods[0] == 2
        self.q = sl.q
        # half group H2 = G'' x Z_q (axes 1.. of G' then Z_q)
        self.Gpp_mods = Gp.mods[1:]
        self.Hh = AG(self.Gpp_mods + (sl.q,))
        n2 = Gp.n // 2; self.n2 = n2
        # element of G with (e, b'', z): G' index e*n2 + b'' ; G index = bidx/zco mapping
        # build map: G element -> (e, half index)
        e = sl.bco[:, 0]
        bpp = sl.Gp.idx(np.concatenate([np.zeros((G.n, 1), np.int64), sl.bco[:, 1:]], 1)) if len(Gp.mods) > 1 else np.zeros(G.n, np.int64)
        half = self.Hh.idx(np.concatenate([sl.bco[:, 1:], sl.zco[:, None]], 1))
        self.g_e, self.g_half = e, half
        # targets on the half group
        T = sl.T
        TS = np.zeros(self.Hh.n, np.int64); TD = np.zeros(self.Hh.n, np.int64)
        np.add.at(TS, half, T); np.add.at(TD, half, np.where(e == 0, T, -T))
        self.TS, self.TD = TS, TD
        # S-invariant structure on the half group: orbits of S restricted (U on G'', M on Z_q)
        self.Hperms = []
        for u in sl.S:
            uu = tuple(u[1:])   # drop axis 0 of G' (identity there)
            self.Hperms.append(self.Hh.act(uu))
        self.hreps, self.hcls = AGs.lag_reps(self.Hh, self.Hperms)
        # G''-orbits (pairs of G'-orbits) and fibre parts
        self.pairs = []        # list of (j0, j1, G''-orbit elements)
        seen = set()
        for j, bb in enumerate(sl.breps):
            b0 = int(bb) % n2
            j0 = int(sl.borb[b0]); j1 = int(sl.borb[b0 + n2])
            if (j0, j1) in seen: continue
            seen.add((j0, j1))
            self.pairs.append((j0, j1))
        # fibre parts: the S-orbits of each G'-orbit's fibre, as sets of G elements; for the half group we need, for
        # pair (j0, j1), the half-group elements of each part (parts of j0 and j1 correspond by z-class)
        self.part_half = []
        for (j0, j1) in self.pairs:
            parts0 = [np.nonzero(sl.orb == o)[0] for o, _ in sl.parts[j0]]
            parts1 = [np.nonzero(sl.orb == o)[0] for o, _ in sl.parts[j1]]
            h0 = [sorted(set(self.g_half[p].tolist())) for p in parts0]
            h1 = [sorted(set(self.g_half[p].tolist())) for p in parts1]
            # match parts by their half-group support
            idx1 = {tuple(h): i for i, h in enumerate(h1)}
            perm = [idx1[tuple(h)] for h in h0]
            self.part_half.append((h0, perm))
        # middle axis
        self.Laxis = Laxis                    # axis index in G' (>= 1) -> axis Laxis - 1 in G''
        La = Laxis - 1
        self.Lmask_half = np.ones(self.Hh.n, bool)
        for a in range(len(self.Hh.mods)):
            if a != La: self.Lmask_half &= (self.Hh.co[:, a] == 0)
        self.Mcos, self.nM = AGs.subgroup_compress_index(self.Hh, self.Lmask_half)
        TM = np.zeros(self.nM, np.int64); np.add.at(TM, self.Mcos, TS); self.TM = TM

    def block_pairs(self, k):
        """per G''-orbit: list of S-values (vectors over the half-group elements of the pair's fibre) with the
        D-option lists"""
        sl = self.sl; out = []
        for pi, (j0, j1) in enumerate(self.pairs):
            P0, _ = AGs.fibre_patterns(sl, j0, int(self.Y[k, sl.breps[j0]]))
            P1, _ = AGs.fibre_patterns(sl, j1, int(self.Y[k, sl.breps[j1]]))
            h0, perm = self.part_half[pi]
            P1p = P1[:, perm]
            byS = defaultdict(list)
            for a in range(len(P0)):
                for b in range(len(P1p)):
                    s_ = tuple((P0[a].astype(np.int64) + P1p[b]).tolist())
                    byS[s_].append((a, b))
            out.append(dict(P0=P0, P1=P1p, byS=byS, h0=h0))
        return out

    def stats(self):
        res = []
        for k in range(4):
            bp = self.block_pairs(k)
            nS = [len(d['byS']) for d in bp]
            res.append(dict(block=k, S_values_per_orbit=nS, log2_S_configs=round(sum(math.log2(x) for x in nS), 2),
                            log2_pairs=round(sum(math.log2(len(d['P0']) * len(d['P1'])) for d in bp), 2)))
        return res


if __name__ == "__main__":
    a = sys.argv[1:]
    cfg = cfg_from_json(a[0])
    src = a[1]
    sl = W.build(cfg)
    if src.startswith('core:'):
        from ab_verify_published import load_core, recoordinate
        hdr, X0 = load_core(src[5:])
        mo = tuple(int(t) for t in hdr.split('group=')[1].split()[0].split('x'))
        amap = json.loads(a[a.index('--amap') + 1])
        X = recoordinate(X0, mo, cfg['mods'], amap)
        Y = sl.compress(X)
    else:
        d = json.load(open(src)); Y = np.array(d['solutions'][0]['y'] if 'solutions' in d else d['y'])
    Lax = int(a[a.index('--Laxis') + 1]) if '--Laxis' in a else 2
    S2 = Stage2(cfg, Y, Lax)
    print(json.dumps(dict(half=S2.Hh.mods, lag_classes_half=len(S2.hreps), pairs=len(S2.pairs), middle_cosets=S2.nM)))
    for r in S2.stats(): print(json.dumps(r))


def row_structure(S2):
    """rows = U-orbits on G''/L (L = middle axis); for each row, the G''-orbit pairs it contains and, for each, the
    number of its elements lying over the row representative coset"""
    Hh = S2.Hh; sl = S2.sl
    La = S2.Laxis - 1
    # G'' coordinates of the pair orbits' representatives
    Gpp = AG(S2.Gpp_mods)
    nax = len(Gpp.mods)
    # U action on G''/L: drop axis La
    keep = [a for a in range(nax) if a != La]
    Qm = AG(tuple(Gpp.mods[a] for a in keep))
    Uq = sorted({tuple(u[1:][a] for a in keep) for u in sl.Sp})
    Qperm = [Qm.act(u) for u in Uq]
    rorb = -np.ones(Qm.n, np.int64); rreps = []
    for t in range(Qm.n):
        if rorb[t] >= 0: continue
        for P in Qperm: rorb[P[t]] = len(rreps)
        rreps.append(t)
    rows = [dict(rep=t, orbits=[]) for t in rreps]
    for pi, (j0, j1) in enumerate(S2.pairs):
        el = np.nonzero(sl.borb == j0)[0]                   # G' elements with axis-0 coordinate 0
        cop = sl.Gp.co[el][:, 1:]
        tq = Qm.idx(cop[:, keep])
        r = {int(rorb[t]) for t in tq}
        assert len(r) == 1
        r = r.pop()
        cnt = int((tq == rreps[r]).sum())
        rows[r]['orbits'].append((pi, cnt))
    return rows, Qm, rorb


if __name__ == "__main__" and '--rows' in sys.argv:
    rows, Qm, rorb = row_structure(S2)
    print(json.dumps(dict(rows=[dict(rep=r['rep'], orbits=r['orbits']) for r in rows])))
    for k in range(4):
        bp = S2.block_pairs(k)
        for ri, r in enumerate(rows):
            # row values: sum over orbits of cnt * S-value (vector over parts)
            vals = Counter(); wts = Counter()
            lists = [list(bp[pi]['byS'].items()) for pi, c in r['orbits']]
            for combo in itertools.product(*lists):
                v_ = np.zeros(len(combo[0][0]), np.int64); w_ = 1
                for (sv, opts), (pi, c) in zip(combo, r['orbits']):
                    v_ += c * np.array(sv); w_ *= len(opts)
                key = tuple(v_.tolist()); vals[key] += 1; wts[key] += w_
            nconf = sum(vals.values())
            top = sorted(wts.items(), key=lambda x: -x[1])
            totw = sum(wts.values())
            cum = np.cumsum([w for _, w in top]) / totw
            print(json.dumps(dict(block=k, row=ri, configs=nconf, distinct_values=len(vals),
                                  log2_max_fibre=round(math.log2(max(vals.values())), 2),
                                  values_for_50pct_weight=int(np.searchsorted(cum, 0.5)) + 1,
                                  values_for_90pct_weight=int(np.searchsorted(cum, 0.9)) + 1)))


# ---------------------------------------------------------------------------------------------------------------
# bilinear PAF tensors on the half group, per pair of G''-orbit pairs and fibre parts
def paf_tensor(S2, lag_reps):
    """C[o, o2, p, p2, g] = #{(b, z): b in orbit o (G''-part), z in part p; (b, z) + lag_g in (orbit o2, part p2)}
    where 'orbit o' = the half-group elements of G''-orbit-pair o (all its G'' elements), part p = the fibre part
    index (z-class) common to both halves.  Then PAF(x)(g) = sum C[o,o2,p,p2,g] s_o[p] s_o2[p2] for x constant
    s_o[p] on (orbit o, part p)."""
    Hh = S2.Hh
    lab = -np.ones(Hh.n, np.int64)          # half element -> (orbit pair index, part index)
    no = len(S2.pairs); npart = None
    for pi in range(no):
        h0, perm = S2.part_half[pi]
        npart = len(h0)
        for p, els in enumerate(h0):
            for e in els: lab[e] = pi * npart + p
    assert (lab >= 0).all()
    nlab = no * npart
    C = np.zeros((nlab, nlab, len(lag_reps)), np.int64)
    for gi, g in enumerate(lag_reps):
        sh = Hh.idx(Hh.co + Hh.co[g])
        np.add.at(C[:, :, gi], (lab, lab[sh]), 1)
    return C.reshape(no, npart, no, npart, len(lag_reps)).transpose(0, 2, 1, 3, 4), lab, npart


class Solver:
    def __init__(self, cfg, Y, Laxis=2, log=print):
        self.S2 = Stage2(cfg, Y, Laxis); S2 = self.S2; self.log = log
        self.rows, self.Qm, self.rorb = row_structure(S2)
        self.lr = S2.hreps
        self.Ct, self.hlab, self.npart = paf_tensor(S2, self.lr)
        self.TS = S2.TS[self.lr]; self.TD = S2.TD[self.lr]
        self.blocks = [S2.block_pairs(k) for k in range(4)]
        no = len(S2.pairs)
        # per block per orbit: arrays of S-values (n_o x npart) and the pair lists
        self.Sv = [[np.array(list(bp[o]['byS'].keys()), np.int64) for o in range(no)] for bp in self.blocks]
        self.Smult = [[np.array([len(v) for v in bp[o]['byS'].values()], np.int64) for o in range(no)] for bp in self.blocks]
        # per-block bilinear tables M[o][o2] (n_o x n_o2 x L)
        self.Mt = []
        for k in range(4):
            Mk = {}
            for o in range(no):
                for o2 in range(o, no):
                    A = self.Sv[k][o]; B = self.Sv[k][o2]
                    T_ = np.einsum('ip,pql,jq->ijl', A, self.Ct[o, o2], B)
                    if o2 != o:
                        T_ = T_ + np.einsum('jp,pql,iq->ijl', B, self.Ct[o2, o], A)
                    Mk[(o, o2)] = T_
            self.Mt.append(Mk)
        # middle level: x_M on (G''/L) x Z_q -- parts per row
        log(json.dumps(dict(stage="setup", lag_classes=len(self.lr), orbits=no, parts=self.npart,
                            rows=[r['orbits'] for r in self.rows])))

    def sig_configs(self, k, cfgs):
        """S-signatures of configurations (array n x no of S-value indices) of block k"""
        no = cfgs.shape[1]; out = np.zeros((len(cfgs), len(self.lr)), np.int64)
        for o in range(no):
            for o2 in range(o, no):
                out += self.Mt[k][(o, o2)][cfgs[:, o], cfgs[:, o2]]
        return out

    def row_tables(self, k, ri, topR, K, rng):
        """row ri of block k: the topR row values by weight, each with up to K configurations (weighted sample)"""
        r = self.rows[ri]
        orbs = [pi for pi, c in r['orbits']]; cnts = [c for pi, c in r['orbits']]
        lists = [range(len(self.Sv[k][pi])) for pi in orbs]
        vals = defaultdict(list); wts = defaultdict(float)
        for combo in itertools.product(*lists):
            v_ = sum(c * self.Sv[k][pi][i] for pi, c, i in zip(orbs, cnts, combo))
            w_ = float(np.prod([self.Smult[k][pi][i] for pi, i in zip(orbs, combo)]))
            key = tuple(v_.tolist()); vals[key].append(combo); wts[key] += w_
        keys = list(wts); w = np.array([wts[x] for x in keys], float)
        if getattr(self, 'row_select', 'sample') == 'top':
            top = sorted(wts, key=lambda x: -wts[x])[:topR]
        else:                                  # weighted sampling without replacement (typical energies)
            pick = rng.choice(len(keys), size=min(topR, len(keys)), replace=False, p=w / w.sum())
            top = [keys[i] for i in pick]
        out = []
        for key in top:
            cl = vals[key]
            if len(cl) > K:
                w = np.array([np.prod([self.Smult[k][pi][i] for pi, i in zip(orbs, cb)]) for cb in cl], float)
                pick = rng.choice(len(cl), size=K, replace=False, p=w / w.sum())
                cl = [cl[i] for i in pick]
            out.append((np.array(key), np.array(cl, np.int64), wts[key]))
        return orbs, out

    # ------------------------------------------------------------------ middle level
    def middle_setup(self):
        S2 = self.S2; Hh = S2.Hh
        La = S2.Laxis - 1
        keep = [a for a in range(len(S2.Gpp_mods)) if a != La]
        QH = AG(tuple(S2.Gpp_mods[a] for a in keep) + (S2.q,))
        hq = QH.idx(np.concatenate([Hh.co[:, keep], Hh.co[:, -1:]], 1))
        TM = np.zeros(QH.n, np.int64); np.add.at(TM, hq, S2.TS)
        perms = []
        for u in S2.sl.S:
            uu = tuple(u[1:]); uq = tuple(uu[a] for a in keep) + (uu[-1],)
            perms.append(QH.act(uq))
        mreps, _ = AGs.lag_reps(QH, perms)
        # z -> part (from pair 0)
        h0, _ = S2.part_half[0]
        zpart = -np.ones(S2.q, np.int64)
        for p, els in enumerate(h0):
            for e in els: zpart[Hh.co[e, -1]] = p
        assert (zpart >= 0).all()
        tq = Hh.co[:, keep]
        Qm = self.Qm
        lab = np.array([self.rorb[Qm.idx(QH.co[i, :-1])] * self.npart + zpart[QH.co[i, -1]] for i in range(QH.n)])
        nr = len(self.rows)
        CM = np.zeros((nr * self.npart, nr * self.npart, len(mreps)), np.int64)
        for gi, g in enumerate(mreps):
            sh = QH.idx(QH.co + QH.co[g])
            np.add.at(CM[:, :, gi], (lab, lab[sh]), 1)
        self.CM = CM.reshape(nr, self.npart, nr, self.npart, len(mreps)).transpose(0, 2, 1, 3, 4)
        self.mreps = mreps; self.TM = TM[mreps]; self.QH = QH

    def parent_sigs(self, rowvals):
        """rowvals: list over rows of arrays (n_r x npart); returns the middle signature of every combination
        (cartesian product, row 0 slowest)"""
        nr = len(rowvals)
        sizes = [len(v) for v in rowvals]
        idx = np.array(list(itertools.product(*[range(n) for n in sizes])), np.int64).reshape(-1, nr)
        out = np.zeros((len(idx), len(self.mreps)), np.int64)
        for r in range(nr):
            for r2 in range(nr):
                A = rowvals[r][idx[:, r]]; B = rowvals[r2][idx[:, r2]]
                out += np.einsum('np,pql,nq->nl', A, self.CM[r, r2], B)
        return idx, out

    # ------------------------------------------------------------------ the S-join for one middle value
    def build_side(self, kA, kB, partsA, partsB, pairlist):
        """S-signature sums and config index pairs for the fibre products of the given parent pairs"""
        sigs = []; idx = []
        for (a, b) in pairlist:
            cA, sA = partsA[a]; cB, sB = partsB[b]
            S = (sA[:, None, :] + sB[None, :, :]).reshape(-1, sA.shape[1])
            ii = np.stack(np.broadcast_arrays(np.arange(len(cA))[:, None], np.arange(len(cB))[None, :]), -1).reshape(-1, 2)
            sigs.append(S); idx.append(np.concatenate([np.full((len(ii), 1), a), np.full((len(ii), 1), b), ii], 1))
        if not sigs: return np.zeros((0, len(self.lr)), np.int64), np.zeros((0, 4), np.int64)
        return np.concatenate(sigs), np.concatenate(idx)

    def prepare_parents(self, topR, K, seed=0):
        """per block: parents (row-value combinations) with their configurations and S-signatures"""
        rng = np.random.default_rng(seed)
        self.middle_setup()
        self.P = []
        for k in range(4):
            rowinfo = [self.row_tables(k, ri, topR[ri], K[ri], rng) for ri in range(len(self.rows))]
            rowvals = [np.array([v for v, c, w in ri_[1]]) for ri_ in rowinfo]
            pidx, psig = self.parent_sigs(rowvals)
            no = len(self.S2.pairs)
            parts = []
            for comb in pidx:
                # configurations: product of the rows' config lists
                lists = [rowinfo[ri][1][comb[ri]][1] for ri in range(len(self.rows))]
                orbs = [rowinfo[ri][0] for ri in range(len(self.rows))]
                prod = itertools.product(*[range(len(L)) for L in lists])
                cf = []
                for t in prod:
                    full = np.zeros(no, np.int64)
                    for ri in range(len(self.rows)):
                        full[orbs[ri]] = lists[ri][t[ri]]
                    cf.append(full)
                cf = np.array(cf, np.int64)
                parts.append((cf, self.sig_configs(k, cf)))
            self.P.append(dict(pidx=pidx, psig=psig, parts=parts))
        self.log(json.dumps(dict(stage="parents", parents=[len(p['parts']) for p in self.P],
                                 configs=[int(sum(len(c) for c, s in p['parts'])) for p in self.P])))

    def middle_candidates(self, top=20, seed=1, max_pairs=8 * 10 ** 7):
        """exact: all left pair sums psig0 + psig1 and right complements TM - psig2 - psig3 (64-bit hashes);
        middle values present on both sides, ranked by (left fibre mass) x (right fibre mass)"""
        P = self.P
        rng = np.random.default_rng(seed)
        wv = rng.integers(1, 2 ** 62, size=len(self.mreps), dtype=np.int64)
        h = [p['psig'] @ wv for p in P]
        fm = [np.array([len(c) for c, s_ in p['parts']], np.float64) for p in P]
        if len(h[0]) * len(h[1]) > max_pairs or len(h[2]) * len(h[3]) > max_pairs: raise RuntimeError("too many parent pairs")
        L = (h[0][:, None] + h[1][None, :]).ravel(); LW = (fm[0][:, None] * fm[1][None, :]).ravel()
        R = (np.int64(self.TM @ wv) - (h[2][:, None] + h[3][None, :])).ravel(); RW = (fm[2][:, None] * fm[3][None, :]).ravel()
        ul, il = np.unique(L, return_inverse=True); wl = np.bincount(il.ravel(), weights=LW)
        ur, ir = np.unique(R, return_inverse=True); wr = np.bincount(ir.ravel(), weights=RW)
        common, a_, b_ = np.intersect1d(ul, ur, return_indices=True)
        score = wl[a_] * wr[b_]
        order = np.argsort(-score)[:top]
        # recover an explicit c vector for each chosen hash
        out = []
        for oi in order:
            hv = common[oi]
            k = int(np.nonzero(L == hv)[0][0]); a, b = divmod(k, len(h[1]))
            out.append((P[0]['psig'][a] + P[1]['psig'][b], float(np.log2(score[oi]))))
        return out, len(common)

    def pairs_for(self, kA, kB, c):
        """all parent pairs (a, b) of blocks kA, kB with psig_A[a] + psig_B[b] = c"""
        PA, PB = self.P[kA]['psig'], self.P[kB]['psig']
        d = defaultdict(list)
        for b, row in enumerate(PB): d[row.tobytes()].append(b)
        out = []
        for a, row in enumerate(PA):
            for b in d.get((c - row).tobytes(), []): out.append((a, b))
        return out

    def s_join(self, c, max_pairs=6 * 10 ** 7):
        left = self.pairs_for(0, 1, c)
        right = self.pairs_for(2, 3, self.TM - c)
        nl = sum(len(self.P[0]['parts'][a][0]) * len(self.P[1]['parts'][b][0]) for a, b in left)
        nr = sum(len(self.P[2]['parts'][a][0]) * len(self.P[3]['parts'][b][0]) for a, b in right)
        if nl + nr > max_pairs: raise RuntimeError(f"side too large ({nl}, {nr})")
        SL, IL = self.build_side(0, 1, self.P[0]['parts'], self.P[1]['parts'], left)
        SR, IR = self.build_side(2, 3, self.P[2]['parts'], self.P[3]['parts'], right)
        rng = np.random.default_rng(5)
        wv = rng.integers(1, 2 ** 62, size=SL.shape[1] if len(SL) else len(self.lr), dtype=np.int64)
        hl = SL @ wv; hr = (self.TS @ wv) - SR @ wv
        o = np.argsort(hl); hs = hl[o]
        lo = np.searchsorted(hs, hr, 'left'); hi = np.searchsorted(hs, hr, 'right')
        sols = []
        for r in np.nonzero(hi > lo)[0]:
            for p in o[lo[r]:hi[r]]:
                if ((SL[p] + SR[r]) == self.TS).all(): sols.append((IL[p], IR[r]))
        return sols, dict(left_parent_pairs=len(left), right_parent_pairs=len(right), left=int(nl), right=int(nr))

    def configs_of(self, sol):
        """the four blocks' S-configurations (orbit S-value indices) of an S-solution"""
        IL, IR = sol
        a, b, i, j = IL; a2, b2, i2, j2 = IR
        return [self.P[0]['parts'][a][0][i], self.P[1]['parts'][b][0][j], self.P[2]['parts'][a2][0][i2],
                self.P[3]['parts'][b2][0][j2]]

    # ------------------------------------------------------------------ D-check and assembly
    def d_check(self, cfgs, max_list=200000):
        """exact D-level for an S-solution: all D-configurations, 2+2 join on the D-signature"""
        no = len(self.S2.pairs)
        lists = []
        for k in range(4):
            opts = []
            for o in range(no):
                sv = tuple(self.Sv[k][o][cfgs[k][o]].tolist())
                ab = self.blocks[k][o]['byS'][sv]
                P0 = self.blocks[k][o]['P0']; P1 = self.blocks[k][o]['P1']
                opts.append([(P0[a].astype(np.int64) - P1[b], a, b) for a, b in ab])
            n = int(np.prod([len(x) for x in opts]))
            if n > max_list: raise RuntimeError("D-list too large")
            D = []; meta = []
            for combo in itertools.product(*opts):
                D.append(np.stack([d for d, a, b in combo])); meta.append([(a, b) for d, a, b in combo])
            D = np.array(D)                                  # n x no x npart
            sig = np.einsum('nop,oqprl,nqr->nl', D, self.Ct, D)
            lists.append((sig, meta))
        L0, L1, L2, L3 = [x[0] for x in lists]
        rng = np.random.default_rng(9); wv = rng.integers(1, 2 ** 62, size=L0.shape[1], dtype=np.int64)
        h01 = (L0 @ wv)[:, None] + (L1 @ wv)[None, :]
        h23 = (self.TD @ wv) - ((L2 @ wv)[:, None] + (L3 @ wv)[None, :])
        o = np.argsort(h01.ravel()); hs = h01.ravel()[o]
        lo = np.searchsorted(hs, h23.ravel(), 'left'); hi = np.searchsorted(hs, h23.ravel(), 'right')
        sols = []
        n1, n3 = L1.shape[0], L3.shape[0]
        for r in np.nonzero(hi > lo)[0]:
            c_, d_ = divmod(int(r), n3)
            for p in o[lo[r]:hi[r]]:
                a_, b_ = divmod(int(p), n1)
                if ((L0[a_] + L1[b_] + L2[c_] + L3[d_]) == self.TD).all():
                    sols.append([lists[0][1][a_], lists[1][1][b_], lists[2][1][c_], lists[3][1][d_]])
        return sols, [len(x[0]) for x in lists]

    def assemble(self, dsol):
        """x on G for the four blocks from the per-orbit (a, b) pattern-pair choices; verified"""
        S2 = self.S2; sl = S2.sl
        X = np.zeros((4, sl.G.n), np.int64)
        for k in range(4):
            lab = np.zeros(sl.r, np.int64)
            for o, (j0, j1) in enumerate(S2.pairs):
                a, b = dsol[k][o]
                p0 = self.blocks[k][o]['P0'][a]
                h0, perm = S2.part_half[o]
                p1p = self.blocks[k][o]['P1'][b]                 # permuted to j0's part order
                p1 = np.empty_like(p1p); p1[perm] = p1p          # back to j1's own part order
                for i, (orb, sz) in enumerate(sl.parts[j0]): lab[orb] = p0[i]
                for i, (orb, sz) in enumerate(sl.parts[j1]): lab[orb] = p1[i]
            X[k] = lab[sl.orb]
        ok = bool((sl.G.paf(X).sum(0) == sl.T).all())
        return X, ok


class FastSolver(Solver):
    """parents as (row value index per row); configurations and S-signatures built lazily; vectorised products"""
    def prepare_rows(self, R, K, seed=0):
        rng = np.random.default_rng(seed)
        self.middle_setup()
        self.rowdata = []         # per block: list over rows of (orbs, values (n x npart), [config arrays])
        for k in range(4):
            rd = []
            for ri in range(len(self.rows)):
                orbs, out = self.row_tables(k, ri, R[ri], K[ri], rng)
                vals = np.array([v for v, c, w in out]); cfgs = [c for v, c, w in out]
                rd.append((orbs, vals, cfgs))
            self.rowdata.append(rd)
        self.pidx = []; self.psig = []
        for k in range(4):
            idx, sig = self.parent_sigs([rd[1] for rd in self.rowdata[k]])
            self.pidx.append(idx); self.psig.append(sig)
        self.fib = [np.array([np.prod([len(self.rowdata[k][r][2][i]) for r, i in enumerate(comb)]) for comb in self.pidx[k]], np.int64)
                    for k in range(4)]
        self._sigcache = [dict() for _ in range(4)]
        self.log(json.dumps(dict(stage="rows", parents=[len(p) for p in self.pidx],
                                 mean_fibre=[round(float(f.mean()), 1) for f in self.fib])))

    def parent_configs(self, k, a):
        if a in self._sigcache[k]: return self._sigcache[k][a]
        comb = self.pidx[k][a]; no = len(self.S2.pairs)
        lists = [self.rowdata[k][r][2][i] for r, i in enumerate(comb)]
        orbs = [self.rowdata[k][r][0] for r in range(len(comb))]
        grids = np.meshgrid(*[np.arange(len(L)) for L in lists], indexing='ij')
        n = grids[0].size
        cf = np.zeros((n, no), np.int64)
        for r, L in enumerate(lists):
            cf[:, orbs[r]] = L[grids[r].ravel()]
        res = (cf, self.sig_configs(k, cf))
        if len(self._sigcache[k]) > 3000: self._sigcache[k].clear()
        self._sigcache[k][a] = res
        return res

    def choose_middles(self, nsample=4 * 10 ** 6, top=30, seed=2):
        rng = np.random.default_rng(seed)
        wv = rng.integers(1, 2 ** 62, size=len(self.mreps), dtype=np.int64)
        h = [p @ wv for p in self.psig]
        pw = [f / f.sum() for f in self.fib]
        i = rng.choice(len(h[0]), size=nsample, p=pw[0]); j = rng.choice(len(h[1]), size=nsample, p=pw[1])
        L = h[0][i] + h[1][j]
        i2 = rng.choice(len(h[2]), size=nsample, p=pw[2]); j2 = rng.choice(len(h[3]), size=nsample, p=pw[3])
        R = np.int64(self.TM @ wv) - (h[2][i2] + h[3][j2])
        ul, cl = np.unique(L, return_counts=True); ur, cr = np.unique(R, return_counts=True)
        common, a_, b_ = np.intersect1d(ul, ur, return_indices=True)
        score = cl[a_].astype(float) * cr[b_]
        order = np.argsort(-score)[:top]
        out = []
        for oi in order:
            k_ = int(np.nonzero(L == common[oi])[0][0])
            out.append((self.psig[0][i[k_]] + self.psig[1][j[k_]], float(score[oi])))
        self._hv = (wv, h)
        return out, len(common)

    def pairs_for(self, kA, kB, c):
        wv, h = self._hv
        if not hasattr(self, '_hd'): self._hd = {}
        if kB not in self._hd:
            d = defaultdict(list)
            for b, x in enumerate(h[kB].tolist()): d[x].append(b)
            self._hd[kB] = d
        d = self._hd[kB]
        hc = int(c @ wv)
        out = []
        for a, x in enumerate(h[kA].tolist()):
            for b in d.get((hc - x + 2 ** 63) % 2 ** 64 - 2 ** 63, []):
                if (self.psig[kA][a] + self.psig[kB][b] == c).all(): out.append((a, b))
        return out

    def s_join(self, c, max_side=1.5 * 10 ** 8):
        left = self.pairs_for(0, 1, c)
        right = self.pairs_for(2, 3, self.TM - c)
        nl = int(sum(self.fib[0][a] * self.fib[1][b] for a, b in left))
        nr = int(sum(self.fib[2][a] * self.fib[3][b] for a, b in right))
        info = dict(left_parent_pairs=len(left), right_parent_pairs=len(right), left=nl, right=nr)
        if nl == 0 or nr == 0: return [], info
        if nl + nr > max_side: raise RuntimeError(f"side too large {info}")
        rng = np.random.default_rng(5)
        wv = rng.integers(1, 2 ** 62, size=len(self.lr), dtype=np.int64)
        def side(kA, kB, pl):
            H = []; I = []
            for (a, b) in pl:
                cA, sA = self.parent_configs(kA, a); cB, sB = self.parent_configs(kB, b)
                hA = sA @ wv; hB = sB @ wv
                H.append((hA[:, None] + hB[None, :]).ravel())
                ii = np.stack(np.broadcast_arrays(np.arange(len(cA))[:, None], np.arange(len(cB))[None, :]), -1).reshape(-1, 2)
                I.append(np.concatenate([np.full((len(ii), 1), a), np.full((len(ii), 1), b), ii], 1))
            return np.concatenate(H), np.concatenate(I)
        HL, IL = side(0, 1, left); HR, IR = side(2, 3, right)
        HR = np.int64(self.TS @ wv) - HR
        o = np.argsort(HL); hs = HL[o]
        lo = np.searchsorted(hs, HR, 'left'); hi = np.searchsorted(hs, HR, 'right')
        sols = []
        for r in np.nonzero(hi > lo)[0]:
            for p in o[lo[r]:hi[r]]:
                a, b, i, j = IL[p]; a2, b2, i2, j2 = IR[r]
                cf = [self.parent_configs(0, a)[0][i], self.parent_configs(1, b)[0][j],
                      self.parent_configs(2, a2)[0][i2], self.parent_configs(3, b2)[0][j2]]
                sig = sum(self.sig_configs(k, cf[k][None, :])[0] for k in range(4))
                if (sig == self.TS).all(): sols.append(cf)
        return sols, info


class FastSolver2(FastSolver):
    """memory-lean S-join: packed indices (parent a: 22 bits, parent b: 22 bits, config i, j: 10 bits each),
    right side processed in chunks; configuration S-signatures cached per parent"""
    def s_join(self, c, max_side=2.5 * 10 ** 8, chunk=4 * 10 ** 6):
        left = self.pairs_for(0, 1, c)
        right = self.pairs_for(2, 3, self.TM - c)
        nl = int(sum(self.fib[0][a] * self.fib[1][b] for a, b in left))
        nr = int(sum(self.fib[2][a] * self.fib[3][b] for a, b in right))
        info = dict(left_parent_pairs=len(left), right_parent_pairs=len(right), left=nl, right=nr)
        if nl == 0 or nr == 0: return [], info
        if nl > max_side or nr > max_side: raise RuntimeError(f"side too large {info}")
        rng = np.random.default_rng(5)
        wv = rng.integers(1, 2 ** 62, size=len(self.lr), dtype=np.int64)
        def pack(a, b, i, j): return (a << 42) | (b << 20) | (i << 10) | j
        def side_iter(kA, kB, pl):
            for (a, b) in pl:
                cA, sA = self.parent_configs(kA, a); cB, sB = self.parent_configs(kB, b)
                hA = sA @ wv; hB = sB @ wv
                H = (hA[:, None] + hB[None, :]).ravel()
                ii, jj = np.meshgrid(np.arange(len(cA), dtype=np.int64), np.arange(len(cB), dtype=np.int64), indexing='ij')
                yield H, pack(np.int64(a), np.int64(b), ii.ravel(), jj.ravel())
        HL = []; PL = []
        for H, P in side_iter(0, 1, left): HL.append(H); PL.append(P)
        HL = np.concatenate(HL); PL = np.concatenate(PL)
        o = np.argsort(HL); HL = HL[o]; PL = PL[o]
        tgt = np.int64(self.TS @ wv)
        sols = []
        bufH = []; bufP = []; nb = 0
        def flush():
            if not bufH: return
            HR = tgt - np.concatenate(bufH); PR = np.concatenate(bufP)
            lo = np.searchsorted(HL, HR, 'left'); hi = np.searchsorted(HL, HR, 'right')
            for r in np.nonzero(hi > lo)[0]:
                for p in range(lo[r], hi[r]):
                    pl, pr = int(PL[p]), int(PR[r])
                    a, b, i, j = pl >> 42, (pl >> 20) & (2 ** 22 - 1), (pl >> 10) & 1023, pl & 1023
                    a2, b2, i2, j2 = pr >> 42, (pr >> 20) & (2 ** 22 - 1), (pr >> 10) & 1023, pr & 1023
                    cf = [self.parent_configs(0, a)[0][i], self.parent_configs(1, b)[0][j],
                          self.parent_configs(2, a2)[0][i2], self.parent_configs(3, b2)[0][j2]]
                    sig = sum(self.sig_configs(k, cf[k][None, :])[0] for k in range(4))
                    if (sig == self.TS).all(): sols.append(cf)
            bufH.clear(); bufP.clear()
        for H, P in side_iter(2, 3, right):
            bufH.append(H); bufP.append(P); nb += len(H)
            if nb >= chunk: flush(); nb = 0
        flush()
        return sols, info
