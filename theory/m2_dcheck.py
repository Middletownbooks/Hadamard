r"""Compiled exact D-check for the b >= 3 hierarchical stage-2 solver of theory/ab_stage2.py [M].

Same algorithm and output as ab_stage2._dcheck_fast (the numpy reference): given an S-solution (per block, per
G''-orbit pair, the index of its S-value), every block's D-configurations are the products of its per-orbit D-option
lists (D = P0[a] - P1[b] over the (a, b) with P0[a] + P1[b] = S); the D-level signature of a configuration is the
bilinear form sum_{o <= o2} d_o^T Ct[o, o2] d_o2 at the half-group lag classes; a 2+2 join on a 64-bit linear hash of
the signatures (left pair keys sorted, right needs binary-searched) finds candidate quadruples and every candidate is
re-checked exactly on the full integer signature.  Solutions are returned as [[(a, b) per orbit] per block], the
format of ab_stage2.Solver.assemble.
Speed: per block the per-orbit-pair tables (n_o x n_o2 x L) are built once per S-solution and each configuration costs
no (no + 1)/2 table rows; the join is O(n0 n1 log(n0 n1) + n2 n3 log(n0 n1)).
usage: python3 theory/m2_dcheck.py test347 [--n 200]    (compares with _dcheck_fast on real 347 S-solutions)"""
import sys, os, json, time
import numpy as np
import numba
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


@numba.njit(cache=True)
def _block_sigs(dval, st, cnt, Ct, w):
    """dval (n_tot, npart); st, cnt (no,): this block's option ranges; Ct (no, no, npart, npart, L); w (L,) int64.
    returns sig (n_cfg, L) and key (n_cfg,) for all configurations (mixed radix, last orbit fastest)"""
    no = st.shape[0]; npart = dval.shape[1]; L = Ct.shape[4]
    n = 1
    for o in range(no): n *= cnt[o]
    # tables T[o, o2] (o <= o2): flat storage with offsets
    off = np.zeros((no, no), np.int64); tot = 0
    for o in range(no):
        for o2 in range(o, no):
            off[o, o2] = tot; tot += cnt[o] * cnt[o2] * L
    T = np.zeros(tot, np.int64)
    for o in range(no):
        for o2 in range(o, no):
            base = off[o, o2]
            for i in range(cnt[o]):
                di = st[o] + i
                for j in range(cnt[o2]):
                    dj = st[o2] + j
                    for l in range(L):
                        acc = 0
                        for p in range(npart):
                            a = dval[di, p]
                            if a == 0: continue
                            for q in range(npart):
                                b = dval[dj, q]
                                if b != 0: acc += a * Ct[o, o2, p, q, l] * b
                        if o2 != o:
                            for p in range(npart):
                                a = dval[dj, p]
                                if a == 0: continue
                                for q in range(npart):
                                    b = dval[di, q]
                                    if b != 0: acc += a * Ct[o2, o, p, q, l] * b
                        T[base + (i * cnt[o2] + j) * L + l] = acc
    sig = np.zeros((n, L), np.int64); key = np.zeros(n, np.int64)
    dig = np.zeros(no, np.int64)
    for c in range(n):
        rem = c
        for o in range(no - 1, -1, -1):
            dig[o] = rem % cnt[o]; rem //= cnt[o]
        for o in range(no):
            for o2 in range(o, no):
                base = off[o, o2] + (dig[o] * cnt[o2] + dig[o2]) * L
                for l in range(L): sig[c, l] += T[base + l]
        k = 0
        for l in range(L): k += sig[c, l] * w[l]
        key[c] = k
    return sig, key


@numba.njit(cache=True)
def _join4h(s0, k0, s1, k1, s2, k2, s3, k3, TD, kT, cap):
    """as _join4 with an open-addressing hash table on the (0, 1) pair keys (the caller puts the smaller pair side
    there); duplicate keys are chained"""
    n0 = k0.shape[0]; n1 = k1.shape[0]; n2 = k2.shape[0]; n3 = k3.shape[0]; L = TD.shape[0]
    nl = n0 * n1
    bits = 1
    while (1 << bits) < 2 * nl: bits += 1
    capt = 1 << bits; mask = capt - 1
    hk = np.zeros(capt, np.int64); hh = -np.ones(capt, np.int64); nxt = -np.ones(nl, np.int64)
    MUL = np.int64(-7046029254386353131)                    # 0x9E3779B97F4A7C15 as int64
    for a in range(n0):
        for b in range(n1):
            e = a * n1 + b; key = k0[a] + k1[b]
            h = ((key * MUL) >> (64 - bits)) & mask
            while True:
                if hh[h] < 0:
                    hk[h] = key; hh[h] = e; break
                if hk[h] == key:
                    nxt[e] = hh[h]; hh[h] = e; break
                h = (h + 1) & mask
    out = np.zeros((cap, 4), np.int64); nout = 0; ncand = 0
    for c in range(n2):
        for d in range(n3):
            need = kT - k2[c] - k3[d]
            h = ((need * MUL) >> (64 - bits)) & mask
            while hh[h] >= 0 and hk[h] != need: h = (h + 1) & mask
            e = hh[h] if hh[h] >= 0 else -1
            while e >= 0:
                a = e // n1; b = e % n1
                ncand += 1
                ok = True
                for l in range(L):
                    if s0[a, l] + s1[b, l] + s2[c, l] + s3[d, l] != TD[l]:
                        ok = False; break
                if ok and nout < cap:
                    out[nout, 0] = a; out[nout, 1] = b; out[nout, 2] = c; out[nout, 3] = d; nout += 1
                e = nxt[e]
    return out[:nout], ncand


@numba.njit(cache=True)
def _join4(s0, k0, s1, k1, s2, k2, s3, k3, TD, kT, cap):
    """all (a, b, c, d) with s0[a] + s1[b] + s2[c] + s3[d] == TD (exact), found through the int64 keys"""
    n0 = k0.shape[0]; n1 = k1.shape[0]; n2 = k2.shape[0]; n3 = k3.shape[0]; L = TD.shape[0]
    left = np.empty(n0 * n1, np.int64)
    for a in range(n0):
        for b in range(n1): left[a * n1 + b] = k0[a] + k1[b]
    order = np.argsort(left)
    ls = left[order]
    out = np.zeros((cap, 4), np.int64); nout = 0; ncand = 0
    for c in range(n2):
        for d in range(n3):
            need = kT - k2[c] - k3[d]
            lo = np.searchsorted(ls, need)
            while lo < ls.shape[0] and ls[lo] == need:
                ab = order[lo]; a = ab // n1; b = ab % n1
                ncand += 1
                ok = True
                for l in range(L):
                    if s0[a, l] + s1[b, l] + s2[c, l] + s3[d, l] != TD[l]:
                        ok = False; break
                if ok and nout < cap:
                    out[nout, 0] = a; out[nout, 1] = b; out[nout, 2] = c; out[nout, 3] = d; nout += 1
                lo += 1
    return out[:nout], ncand


class DCheck:
    def __init__(self, S, seed=9, max_cfg=1 << 22, max_left=1 << 27):
        """S: an ab_stage2.Solver (or subclass) after setup"""
        self.S = S
        no = len(S.S2.pairs); self.no = no
        self.npart = S.npart; self.L = len(S.lr)
        self.Ct = np.ascontiguousarray(S.Ct, dtype=np.int64)
        self.TD = np.ascontiguousarray(S.TD, dtype=np.int64)
        rng = np.random.default_rng(seed)
        self.w = rng.integers(1, 2 ** 62, size=self.L, dtype=np.int64)
        self.kT = np.int64(0)
        with np.errstate(over='ignore'):
            self.kT = np.int64((self.TD * self.w).sum())
        self.max_cfg, self.max_left = max_cfg, max_left
        self.dval, self.dab, self.dst, self.dcnt = [], [], [], []
        for k in range(4):
            vals = []; abl = []
            nS = max(len(S.Sv[k][o]) for o in range(no))
            st = np.zeros((no, nS), np.int64); ct = np.zeros((no, nS), np.int64)
            for o in range(no):
                P0 = S.blocks[k][o]['P0']; P1 = S.blocks[k][o]['P1']
                for i, sv in enumerate(S.Sv[k][o]):
                    ab = S.blocks[k][o]['byS'][tuple(int(x) for x in sv)]
                    st[o, i] = len(vals); ct[o, i] = len(ab)
                    for (a, b) in ab:
                        vals.append((P0[a].astype(np.int64) - P1[b]).tolist()); abl.append((a, b))
            self.dval.append(np.array(vals, np.int64).reshape(-1, self.npart)); self.dab.append(np.array(abl, np.int64))
            self.dst.append(st); self.dcnt.append(ct)

    def check(self, cfgs, cap=1 << 16):
        """exact D-check of one S-solution; returns (solutions, list sizes, number of key candidates)"""
        sigs = []; keys = []; ranges = []
        for k in range(4):
            idx = np.asarray(cfgs[k], np.int64)
            st = self.dst[k][np.arange(self.no), idx].copy(); cnt = self.dcnt[k][np.arange(self.no), idx].copy()
            n = int(np.prod(cnt.astype(np.float64)))
            if n > self.max_cfg: raise RuntimeError("D-list too large")
            s, kk = _block_sigs(self.dval[k], st, cnt, self.Ct, self.w)
            sigs.append(s); keys.append(kk); ranges.append((st, cnt))
        sizes = [len(x) for x in keys]
        # pairing with the least work (build + probe); the smaller pair side is hashed
        best = None
        for (p, q, r, s_) in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
            A = sizes[p] * sizes[q]; B = sizes[r] * sizes[s_]
            if A > B: p, q, r, s_, A, B = r, s_, p, q, B, A
            if best is None or A + B < best[0]: best = (A + B, (p, q, r, s_), A)
        _, perm, nbuild = best
        if nbuild > self.max_left: raise RuntimeError("D-join too large")
        p, q, r, s_ = perm
        hits, ncand = _join4h(sigs[p], keys[p], sigs[q], keys[q], sigs[r], keys[r], sigs[s_], keys[s_], self.TD, self.kT, cap)
        sols = []
        for row in hits:
            byblock = [0, 0, 0, 0]
            for slot, k in enumerate(perm): byblock[k] = int(row[slot])
            sol = []
            for k, ci in zip(range(4), byblock):
                st, cnt = ranges[k]
                rem = int(ci); dig = [0] * self.no
                for o in range(self.no - 1, -1, -1):
                    dig[o] = rem % int(cnt[o]); rem //= int(cnt[o])
                sol.append([tuple(int(t) for t in self.dab[k][st[o] + dig[o]]) for o in range(self.no)])
            sols.append(sol)
        return sols, sizes, int(ncand)


def _canon(sols):
    return sorted(json.dumps(s) for s in sols)


def cfgs_of_x(S, X):
    """S-value indices (per block, per orbit pair) and the (a, b) pattern choices of a lift X of the solver's y"""
    sl = S.S2.sl; cfs = []; abs_ = []
    rep = {}                                                # orbit of G -> one element
    for g in range(sl.G.n): rep.setdefault(int(sl.orb[g]), g)
    for k in range(4):
        cf = []; ab = []
        for o, (j0, j1) in enumerate(S.S2.pairs):
            p0 = np.array([X[k, rep[orb]] for orb, _ in sl.parts[j0]], np.int8)
            p1 = np.array([X[k, rep[orb]] for orb, _ in sl.parts[j1]], np.int8)
            h0, perm = S.S2.part_half[o]
            P0 = S.blocks[k][o]['P0']; P1 = S.blocks[k][o]['P1']
            a = int(np.nonzero((P0 == p0).all(1))[0][0]); b = int(np.nonzero((P1 == p1[perm]).all(1))[0][0])
            sv = P0[a].astype(np.int64) + P1[b]
            i = int(np.nonzero((S.Sv[k][o] == sv).all(1))[0][0])
            cf.append(i); ab.append((a, b))
        cfs.append(np.array(cf, np.int64)); abs_.append(ab)
    return cfs, abs_


def d_target_of(S, abs_):
    """D-level signature (at the solver's lag classes) of the lift with pattern choices abs_"""
    Hh = S.S2.Hh; tot = np.zeros(Hh.n, np.int64)
    for k in range(4):
        xd = np.zeros(Hh.n, np.int64)
        for o, (a, b) in enumerate(abs_[k]):
            d = S.blocks[k][o]['P0'][a].astype(np.int64) - S.blocks[k][o]['P1'][b]
            h0, perm = S.S2.part_half[o]
            for p, els in enumerate(h0): xd[els] = d[p]
        tot += Hh.paf(xd)
    return tot[S.lr]


def positive_tests(S, DC, cores, nplant=20, seed=5):
    """known cores must be returned from their own S-solutions; planted random lifts (D-target = their own D
    signature) must be returned; the numpy reference must agree in every case"""
    rng = np.random.default_rng(seed); res = []
    TD0, kT0 = DC.TD.copy(), DC.kT
    for name, X in cores:
        cf, ab = cfgs_of_x(S, X)
        new, _, _ = DC.check(cf); ref, _ = S.d_check_fast(cf)
        res.append(dict(case=name, found=[list(map(list, b)) for b in ab] in [[list(map(list, b)) for b in s] for s in new],
                        agree=_canon(ref) == _canon(new), D_solutions=len(new)))
    TDs = S.TD.copy()
    for t in range(nplant):
        ab = [[(int(rng.integers(len(S.blocks[k][o]['P0']))), int(rng.integers(len(S.blocks[k][o]['P1']))))
               for o in range(DC.no)] for k in range(4)]
        cf = []
        for k in range(4):
            row = []
            for o, (a, b) in enumerate(ab[k]):
                sv = S.blocks[k][o]['P0'][a].astype(np.int64) + S.blocks[k][o]['P1'][b]
                row.append(int(np.nonzero((S.Sv[k][o] == sv).all(1))[0][0]))
            cf.append(np.array(row, np.int64))
        tdp = d_target_of(S, ab)
        DC.TD = np.ascontiguousarray(tdp, dtype=np.int64)
        with np.errstate(over='ignore'): DC.kT = np.int64((DC.TD * DC.w).sum())
        S.TD = tdp
        try:
            new, _, _ = DC.check(cf); ref, _ = S.d_check_fast(cf)
        finally:
            DC.TD, DC.kT, S.TD = TD0, kT0, TDs
        want = [list(map(list, b)) for b in ab]
        res.append(dict(case=f'plant{t}', found=want in [[list(map(list, b)) for b in s] for s in new],
                        agree=_canon(ref) == _canon(new), D_solutions=len(new)))
    return res


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == 'test347':
    import ab_stage2 as S2m, ab_sweep as W, ab_assemble as AS
    from ab_s1sd import cfg_from_json
    from ab_verify_published import recoordinate
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    a = sys.argv[2:]
    nmax = int(a[a.index('--n') + 1]) if '--n' in a else 200
    cfg = cfg_from_json('@' + os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')); sl = W.build(cfg)
    T = AS.load_template347()
    XF = recoordinate(T['F'], (2, 3, 57), (2, 3, 3, 19), [(0, 2), (1, 3), (2, 3), (2, 19)])
    Y = sl.compress(XF)
    t0 = time.time()
    Sv = S2m.FastSolver2(cfg, Y, Laxis=2, log=lambda s: None)
    Sv.prepare_rows(R=[800, 800], K=[4, 128], seed=1)
    cands, nc = Sv.choose_middles(top=12, seed=101)
    print(json.dumps(dict(setup_secs=round(time.time() - t0, 1), middles=len(cands))), flush=True)
    DC = DCheck(Sv)
    cores = [('published 347 core', XF)]
    p_a1 = os.path.join(ROOT, 'logs/ab/prod347/core_a1_438c58a201c4084b_1.npy')
    if os.path.exists(p_a1):
        from ab_verify_published import recoordinate as rc
        Xa1 = rc(np.load(p_a1), (2, 3, 57), (2, 3, 3, 19), [(0, 2), (1, 3), (2, 3), (2, 19)])
        assert sl.is_invariant(Xa1) and (sl.G.paf(Xa1).sum(0) == sl.T).all()
        cores.append(('core a1 (theory/ab_prod347.py, 2026-10-14)', Xa1))
    pos = positive_tests(Sv, DC, cores)
    for r in pos: print(json.dumps(r), flush=True)
    assert all(r['found'] and r['agree'] for r in pos)
    tested = 0; agree = 0; tref = tnew = 0.0; nsol_ref = nsol_new = 0
    for ci, (c, sc) in enumerate(cands):
        sols, info = Sv.s_join(c)
        for cf in sols:
            t1 = time.time(); ref, sz_ref = Sv.d_check_fast(cf); t2 = time.time()
            new, sz_new, nc_ = DC.check(cf); t3 = time.time()
            tref += t2 - t1; tnew += t3 - t2
            assert sz_ref == sz_new, (sz_ref, sz_new)
            ok = _canon(ref) == _canon(new); agree += ok; tested += 1
            nsol_ref += len(ref); nsol_new += len(new)
            if tested >= nmax: break
        print(json.dumps(dict(middle=ci, S_solutions=len(sols), tested=tested, agree=agree, D_solutions_ref=nsol_ref,
                              D_solutions_new=nsol_new, ms_ref=round(1000 * tref / tested, 2),
                              ms_new=round(1000 * tnew / tested, 3))), flush=True)
        if tested >= nmax: break
    assert agree == tested
