r"""Exhaustive joint stage 1 of a doubly-staged slice [M] (647 x183 and x353; small controls 239, 443).

Setting (theory/nd_doubly.py): x_k on Z_m constant on the orbits of M = <mu>; qA, qB odd coprime Hall factors of m, so
Z_m = Z_2 x ... and both compressions keep the Z_2 factor.  y_B = compression along Z_qB (stage 1B, on Z_FB),
y_A = compression along Z_qA (stage 1A, on Z_FA).  A core needs a pair (y_A, y_B) that comes from one labelling and
solves both stage-1 systems.  Pipeline, every step exact:
  1. the symmetry group of the slice acting on stage-1B data [P]: multipliers u in Z_m^* (they act on Z_FB through
     Z_FB^*, modulo the image of M), per-block shifts by the M-fixed elements of Z_m, negation of row-0 blocks, and
     permutations of row-0 blocks.  Each acts on the joint problem (both compressions and the labelling), so the joint
     test of y_B depends only on the orbit of y_B;
  2. the alpha-quadruples of stage 1B (theory/nd_s1_sd.py), one per orbit of the group;
  3. every stage-1B solution over each canonical alpha-quadruple (the S/D pipeline of nd_s1_sd, streamed), reduced to
     canonical form and deduplicated.  Restart-safe and parallel: one JSONL line per alpha-quadruple;
  4. per y_B orbit the joint test: per block and Z_2-layer (M fixes Z_2, so the labels, y_B and y_A split by layer) every
     achievable layer of y_A by a vectorised DP over the M-orbits of that layer (cached on the layer's y_B values),
     then the stage-1A equations by the S/D split on Z_FA = Z_2 x Z_n:  sigma = (L0 + L1)/2, delta = (L0 - L1)/2,
     sum_k PAF_n(sigma_k) = TS and sum_k PAF_n(delta_k) = TD.  The S-level is joined first (2+2 over distinct sigma),
     then for every S-solution the D-level over the delta values compatible with it.  Restart-safe JSONL;
  5. every joint solution is handed to the cross level and completion (Doubly.cross, bordered_gs_cyclic).
Controls: --plant N (random labellings at the slice's structure; the joint test must return the planted y_A).
usage: python3 theory/nd_joint_exhaust.py m mu qA qB --dir DIR [--workers 4] [--phase enum|joint|all] [--plant N]"""
import sys, os, json, math, time, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import Problem, lag_reps, paf_lags, pair_join
from nd_doubly import Doubly, bordered_gs_cyclic

# ---------------------------------------------------------------- symmetry group of the stage-1B data
class Sym:
    def __init__(self, m, Mg, F):
        self.F = F
        MF = sorted({a % F for a in Mg})
        units = [u for u in range(1, F) if math.gcd(u, F) == 1]
        seen = set(); self.ureps = []
        for u in units:
            if u in seen: continue
            self.ureps.append(u)
            for a in MF: seen.add(u * a % F)
        self.uidx = {u: (np.arange(F) * u) % F for u in self.ureps}
        fixed = [t for t in range(m) if all((t * a - t) % m == 0 for a in Mg)]
        self.shifts = sorted({t % F for t in fixed})

    def canon_y(self, Y):
        """canonical form of a stage-1B quadruple (rows (2, 0, 0, 0)) under the group"""
        Y = np.asarray(Y, np.int64); best = None
        for u in self.ureps:
            Z = Y[:, self.uidx[u]]
            b0 = min(tuple(np.roll(Z[0], -t).tolist()) for t in self.shifts)
            rest = sorted(min(tuple((sg * np.roll(z, -t)).tolist()) for t in self.shifts for sg in (1, -1)) for z in Z[1:])
            key = (b0,) + tuple(rest)
            if best is None or key < best: best = key
        return best

    def alpha_maps(self, P, A):
        """index permutations of the alpha list A (all M-invariant +-1 sequences on Z_F, index = orbit mask)"""
        reps = P.reps
        def idx(a): return int(sum((1 << j) for j, r in enumerate(reps) if a[r] == 1))
        assert all(idx(A[i]) == i for i in range(0, len(A), max(1, len(A) // 97)))
        U = {u: np.array([idx(A[i][self.uidx[u]]) for i in range(len(A))]) for u in self.ureps}
        S = {t: np.array([idx(np.roll(A[i], -t)) for i in range(len(A))]) for t in self.shifts}
        N = np.array([idx(-A[i]) for i in range(len(A))])
        return U, S, N

    def canon_quad(self, q, U, S, N):
        best = None
        for u in self.ureps:
            a = [int(U[u][i]) for i in q]
            b0 = min(int(S[t][a[0]]) for t in self.shifts)
            rest = sorted(min(min(int(S[t][x]), int(N[S[t][x]])) for t in self.shifts) for x in a[1:])
            key = (b0,) + tuple(rest)
            if best is None or key < best: best = key
        return best

# ---------------------------------------------------------------- stage 1B, streamed
def iter_stage1(P, quad, A, rows=(2, 0, 0, 0), cache=None):
    """every stage-1B solution with the given alpha-quadruple (same steps as Problem.solve)"""
    cache = {} if cache is None else cache
    blocks = []
    for k, ai in enumerate(quad):
        key = (ai, rows[k])
        if key not in cache: cache[key] = P.block(A[ai], rows[k], int(P.TS[0]), int(P.TD[0]))
        blocks.append(cache[key])
    for (i, j, k, l) in pair_join(*[b['Ssig'] for b in blocks], P.TS):
        us = [blocks[0]['U'][i], blocks[1]['U'][j], blocks[2]['U'][k], blocks[3]['U'][l]]
        Ds = [P.dlist(b, u) for b, u in zip(blocks, us)]
        for (a, b_, c, d) in pair_join(*[x[1] for x in Ds], P.TD):
            ws = [Ds[0][0][a], Ds[1][0][b_], Ds[2][0][c], Ds[3][0][d]]
            yield np.array([P.y_of(B, u, w) for B, u, w in zip(blocks, us, ws)], np.int64)

_G = {}
def _enum_worker(args):
    m, mu, qA, qB, i, quad = args
    key = (m, mu, qA, qB)
    if _G.get('key') != key:
        P = Problem(m, mu, qB); A = P.alphas(); sym = Sym(m, P.Mg, P.F)
        _G.update(key=key, P=P, A=A, sym=sym, cache={})
    P, A, sym = _G['P'], _G['A'], _G['sym']
    t0 = time.time(); n = 0; ys = set()
    for Y in iter_stage1(P, quad, A, cache=_G['cache']):
        n += 1
        assert P.verify(Y)
        ys.add(sym.canon_y(Y))
    if len(_G['cache']) > 4000: _G['cache'].clear()
    return dict(i=i, q=[int(x) for x in quad], nsol=n, ys=[list(map(list, y)) for y in sorted(ys)], secs=round(time.time() - t0, 2))

def enumerate_yB(m, mu, qA, qB, path, workers=4, log=print):
    P = Problem(m, mu, qB); A, Q, info = P.alpha_quads()
    sym = Sym(m, P.Mg, P.F); U, S, N = sym.alpha_maps(P, A)
    canon = sorted({sym.canon_quad(q, U, S, N) for q in Q})
    log(json.dumps(dict(stage="alpha", alpha_quads=len(Q), canonical=len(canon), unit_reps=sym.ureps, shifts=sym.shifts, **info)))
    done = set()
    if os.path.exists(path):
        for line in open(path):
            try: done.add(json.loads(line)['i'])
            except Exception: pass
    todo = [(m, mu, qA, qB, i, q) for i, q in enumerate(canon) if i not in done]
    log(json.dumps(dict(stage="enum", todo=len(todo), done=len(done))))
    t0 = time.time()
    with open(path, "a") as fo:
        if workers > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(workers) as pool:
                for r in pool.imap_unordered(_enum_worker, todo, chunksize=4):
                    fo.write(json.dumps(r) + "\n"); fo.flush()
        else:
            for t in todo:
                fo.write(json.dumps(_enum_worker(t)) + "\n"); fo.flush()
    allys = {}; nsol = 0; nq = 0
    for line in open(path):
        r = json.loads(line); nq += 1; nsol += r['nsol']
        for y in r['ys']: allys.setdefault(json.dumps(y), r['i'])
    log(json.dumps(dict(stage="enum_done", canonical_quads=len(canon), quads_logged=nq, raw_solutions=nsol,
                        yB_orbits=len(allys), secs=round(time.time() - t0, 1))))
    return [json.loads(k) for k in allys]

# ---------------------------------------------------------------- joint test
class JointTest:
    def __init__(self, m, mu, qA, qB):
        D = Doubly(m, mu, qA, qB); self.D = D
        assert m % 2 == 0 and qA % 2 == 1 and qB % 2 == 1
        FA, FB = D.FA, D.FB
        assert FA % 2 == 0 and (FA // 2) % 2 == 1
        reps = np.array([int(np.nonzero(D.orb == o)[0][0]) for o in range(D.nor)])
        self.layer_orbs = [np.nonzero(reps % 2 == e)[0] for e in (0, 1)]
        self.Bcols = [[i for i, b in enumerate(D.kB) if b % 2 == e] for e in (0, 1)]
        self.Acols = [[i for i, b in enumerate(D.kA) if b % 2 == e] for e in (0, 1)]
        CB = D.cB[:, D.kB]; CA = D.cA[:, D.kA]
        for e in (0, 1):            # the layer structure [P]: orbits of layer e touch only layer-e positions
            other = self.layer_orbs[1 - e]
            assert not CB[np.ix_(other, self.Bcols[e])].any() and not CA[np.ix_(other, self.Acols[e])].any()
        self.CB, self.CA = CB, CA
        # S/D structure on Z_FA = Z_2 x Z_n
        n = FA // 2; self.n = n
        self.e0 = np.array([next(g for g in range(FA) if g % 2 == 0 and g % n == b) for b in range(n)])
        self.e1 = np.array([next(g for g in range(FA) if g % 2 == 1 and g % n == b) for b in range(n)])
        Mn = sorted({a % n for a in D.MA}); self.Mn = Mn
        On = -np.ones(n, np.int64); nreps = []
        for b in range(n):
            if On[b] >= 0: continue
            for a in Mn: On[b * a % n] = len(nreps)
            nreps.append(b)
        self.On, self.nreps = On, nreps; self.nw = np.bincount(On)
        self.nlags = lag_reps(n, Mn)
        # Z_n position b -> index into the layer's A-columns (layer 0 positions e0[b], layer 1 positions e1[b])
        kA = list(D.kA)
        def colidx(e, g):
            r = D.repA[g]; assert kA[r] % 2 == e; return self.Acols[e].index(r)
        self.pos0 = np.array([colidx(0, int(self.e0[b])) for b in range(n)])
        self.pos1 = np.array([colidx(1, int(self.e1[b])) for b in range(n)])
        TA0 = 4 * m - 4 * (qA - 1)
        self.TS = np.array([(TA0 - 4 * qA) // 4] + [-2 * qA] * (len(self.nlags) - 1), np.int64)
        self.TD = np.array([(TA0 + 4 * qA) // 4] + [0] * (len(self.nlags) - 1), np.int64)
        self.cache = {}

    def layer_options(self, e, yb_layer):
        """all (A-partial vector, count) for layer e with the layer's y_B values (at Bcols[e]) prescribed"""
        key = (e, tuple(yb_layer))
        if key in self.cache: return self.cache[key]
        orbs = self.layer_orbs[e]; CB = self.CB[np.ix_(orbs, self.Bcols[e])]; CA = self.CA[np.ix_(orbs, self.Acols[e])]
        order = sorted(range(len(orbs)), key=lambda o: (np.nonzero(CB[o])[0].max() if CB[o].any() else -1))
        CB, CA = CB[order], CA[order]
        nb, na = CB.shape[1], CA.shape[1]
        maxrem = np.zeros((len(orbs) + 1, nb), np.int64)
        for p in range(len(orbs) - 1, -1, -1): maxrem[p] = maxrem[p + 1] + np.abs(CB[p])
        boundA = np.abs(CA).sum(0)
        rad = np.concatenate([2 * maxrem[0] + 1, 2 * boundA + 1]).astype(object)
        off = np.concatenate([maxrem[0], boundA])
        mult = [1]
        for r in rad[:-1]: mult.append(mult[-1] * int(r))
        enc_ok = int(np.prod(rad)) < 2 ** 62
        mult = np.array(mult, np.int64) if enc_ok else None
        st = np.concatenate([np.array(yb_layer, np.int64), np.zeros(na, np.int64)])[None, :]
        if (np.abs(st[0, :nb]) > maxrem[0]).any(): self.cache[key] = (np.zeros((0, na), np.int64), np.zeros(0)); return self.cache[key]
        cnt = np.ones(1)
        for p in range(len(orbs)):
            d = np.concatenate([-CB[p], CA[p]])
            new = np.concatenate([st + d, st - d]); nc = np.concatenate([cnt, cnt])
            ok = (np.abs(new[:, :nb]) <= maxrem[p + 1]).all(1)
            new, nc = new[ok], nc[ok]
            if len(new) == 0: break
            if enc_ok:
                k = (new + off) @ mult
                u, inv = np.unique(k, return_inverse=True)
                first = np.zeros(len(u), np.int64); first[inv] = np.arange(len(new))
                st = new[first]; cnt = np.bincount(inv, weights=nc)
            else:
                st, inv = np.unique(new, axis=0, return_inverse=True); cnt = np.bincount(inv.ravel(), weights=nc)
        if len(st) == 0 or len(new) == 0: res = (np.zeros((0, na), np.int64), np.zeros(0))
        else:
            z = (st[:, :nb] == 0).all(1); res = (st[z, nb:], cnt[z])
        self.cache[key] = res
        return res

    def block_options(self, yb_k):
        yb_k = np.asarray(yb_k)[self.D.kB]
        L = []
        for e in (0, 1):
            vals, cnt = self.layer_options(e, yb_k[self.Bcols[e]])
            pos = self.pos0 if e == 0 else self.pos1
            L.append((vals[:, pos][:, self.nreps] if len(vals) else np.zeros((0, len(self.nreps)), np.int64), cnt))
        return L      # per layer: (values at the Z_n orbit reps, label counts)

    def paf_n(self, V):
        """PAF on Z_n at nlags of vectors given at the orbit reps"""
        X = V[:, self.On]
        return np.stack([(X * np.roll(X, -h, axis=1)).sum(1) for h in self.nlags], 1)

    def test(self, yb, TS=None, TD=None, maxsol=10 ** 6, maxpairs=6 * 10 ** 7, log=None):
        """exact joint test of one y_B: every y_A (as per-block layer pairs) consistent with y_B and solving stage 1A"""
        from nd_ejoin import ejoin4
        TS = self.TS if TS is None else np.asarray(TS, np.int64); TD = self.TD if TD is None else np.asarray(TD, np.int64)
        t0 = time.time()
        Ls = [self.block_options(yb[k]) for k in range(4)]
        sizes = [(len(a[0]), len(b[0])) for a, b in Ls]
        rec = dict(layer_options=sizes)
        if min(min(s_) for s_ in sizes) == 0: rec.update(S_solutions=0, joint=0, secs=round(time.time() - t0, 2)); return rec, []
        # exact energy filters [P]: the 8 layer energies sum to 2 (TS0 + TD0); per block ||sigma||^2 and ||delta||^2
        # sum to TS0 and TD0 over the blocks; every quantity is at least its per-block minimum
        En = [[(V * V * self.nw).sum(1) for V, c in L] for L in Ls]
        slack = 2 * int(TS[0] + TD[0]) - sum(int(e.min()) for E in En for e in E)
        rec['energy_slack'] = slack
        if slack < 0: rec.update(S_solutions=0, joint=0, secs=round(time.time() - t0, 2)); return rec, []
        P = []
        for k, ((V0, c0), (V1, c1)) in enumerate(Ls):
            E0, E1 = En[k]; bound = int(E0.min()) + int(E1.min()) + slack
            o1 = np.argsort(E1, kind='stable'); E1s = E1[o1]
            cntj = np.searchsorted(E1s, bound - E0, 'right')
            if int(cntj.sum()) > maxpairs:
                rec.update(pairs_block=k, pairs_needed=int(cntj.sum()), joint=None, note='pair list too large'); return rec, []
            I0 = np.repeat(np.arange(len(V0)), cntj)
            I1 = o1[np.concatenate([np.arange(c) for c in cntj])] if cntj.sum() else np.zeros(0, np.int64)
            sg = (V0[I0] + V1[I1]); dl = (V0[I0] - V1[I1])
            assert not (sg % 2).any(); sg //= 2; dl //= 2
            P.append(dict(I0=I0, I1=I1, sg=sg, dl=dl, se=(sg * sg * self.nw).sum(1), de=(dl * dl * self.nw).sum(1)))
        for it in range(6):
            ms = [int(p_['se'].min()) if len(p_['se']) else 10 ** 9 for p_ in P]
            md = [int(p_['de'].min()) if len(p_['de']) else 10 ** 9 for p_ in P]
            changed = False
            for k, p_ in enumerate(P):
                ok = (p_['se'] <= TS[0] - (sum(ms) - ms[k])) & (p_['de'] <= TD[0] - (sum(md) - md[k]))
                if not ok.all():
                    changed = True
                    for key in p_: p_[key] = p_[key][ok]
            if any(len(p_['se']) == 0 for p_ in P):
                rec.update(pairs=[len(p_['se']) for p_ in P], S_solutions=0, joint=0, secs=round(time.time() - t0, 2)); return rec, []
            if not changed: break
        rec['pairs'] = [int(len(p_['se'])) for p_ in P]
        # distinct sigma per block, grouped pairs
        U = []; G = []
        for p_ in P:
            u, inv = np.unique(p_['sg'], axis=0, return_inverse=True); inv = inv.ravel()
            o = np.argsort(inv, kind='stable'); bd = np.searchsorted(inv[o], np.arange(len(u) + 1))
            U.append(u); G.append((o, bd))
        rec['sigma_distinct'] = [len(u) for u in U]
        Ssig = [self.paf_n(u) for u in U]
        try:
            Sh, sst = ejoin4(Ssig, TS)
        except RuntimeError as ex:
            rec.update(S_solutions=None, joint=None, note=str(ex)[:80], secs=round(time.time() - t0, 2)); return rec, []
        rec['S_solutions'] = len(Sh); rec['S_join_log2'] = sst['work_log2']
        sols = []
        if Sh:
            # full join (S and D signatures together) over the pairs whose sigma occurs in some S-solution
            Fl = []; keep = []
            for k in range(4):
                used = np.zeros(len(U[k]), bool); used[[t[k] for t in Sh]] = True
                o, bd = G[k]
                ii = np.concatenate([o[bd[j]:bd[j + 1]] for j in np.nonzero(used)[0]])
                keep.append(ii)
                Fl.append(np.concatenate([self.paf_n(P[k]['sg'][ii]), self.paf_n(P[k]['dl'][ii])], 1))
            rec['full_lists'] = [len(x) for x in keep]
            Fh, _ = ejoin4(Fl, np.concatenate([TS, TD]))
            for t in Fh:
                ya = []; mult = 1
                for k in range(4):
                    j = keep[k][t[k]]; i0, i1 = P[k]['I0'][j], P[k]['I1'][j]
                    (V0, c0), (V1, c1) = Ls[k]
                    ya.append((V0[i0].tolist(), V1[i1].tolist())); mult *= int(c0[i0]) * int(c1[i1])
                sols.append(dict(layers=ya, mult_log2=round(math.log2(mult), 2)))
                if len(sols) >= maxsol: break
        rec['joint'] = len(sols); rec['secs'] = round(time.time() - t0, 2)
        return rec, sols

    def ya_full(self, layers):
        """y_A on Z_FA (all positions) from the per-block layer values at the Z_n orbit reps"""
        FA = self.D.FA; out = []
        for L0, L1 in layers:
            y = np.zeros(FA, np.int64)
            y[self.e0] = np.array(L0)[self.On]; y[self.e1] = np.array(L1)[self.On]
            out.append(y.tolist())
        return out

    def planted_real(self, yb, seed, per_layer_slack=None):
        """plant on a real y_B: per block and layer a random achievable layer of y_A near the layer's minimum energy
        (total slack as in the real problem); targets = the plant's own S and D sums; the test must return it"""
        rng = np.random.default_rng(seed)
        Ls = [self.block_options(yb[k]) for k in range(4)]
        En = [[(V * V * self.nw).sum(1) for V, c in L] for L in Ls]
        slack = 2 * int(self.TS[0] + self.TD[0]) - sum(int(e.min()) for E in En for e in E)
        ps = max(0, slack // 8) if per_layer_slack is None else per_layer_slack
        pick = []
        for k in range(4):
            ch = []
            for e in (0, 1):
                ok = np.nonzero(En[k][e] <= En[k][e].min() + ps)[0]; ch.append(Ls[k][e][0][rng.choice(ok)])
            pick.append(ch)
        sg = np.array([(a + b) // 2 for a, b in pick]); dl = np.array([(a - b) // 2 for a, b in pick])
        TS = self.paf_n(sg).sum(0); TD = self.paf_n(dl).sum(0)
        rec, sols = self.test(yb, TS=TS, TD=TD)
        hit = any(all(list(s_['layers'][k][0]) == pick[k][0].tolist() and list(s_['layers'][k][1]) == pick[k][1].tolist()
                      for k in range(4)) for s_ in sols)
        return dict(seed=seed, **rec, planted_found=hit)

    def planted(self, seed):
        D = self.D; rng = np.random.default_rng(seed)
        lab = rng.choice([-1, 1], (4, D.nor)); X = lab[:, D.orb]
        yb = [[int(X[k][b::D.FB].sum()) for b in range(D.FB)] for k in range(4)]
        ya = np.array([[int(X[k][b::D.FA].sum()) for b in range(D.FA)] for k in range(4)])
        sg = (ya[:, self.e0] + ya[:, self.e1]) // 2; dl = (ya[:, self.e0] - ya[:, self.e1]) // 2
        TS = paf_lags(sg, self.nlags).sum(0); TD = paf_lags(dl, self.nlags).sum(0)
        rec, sols = self.test(yb, TS=TS, TD=TD)
        hit = any(self.ya_full(s_['layers']) == ya.tolist() for s_ in sols)
        return dict(seed=seed, **rec, planted_found=hit)

def joint_phase(m, mu, qA, qB, yBs, path, log=print, cross=True):
    J = JointTest(m, mu, qA, qB)
    done = set()
    if os.path.exists(path):
        for line in open(path):
            try: done.add(json.loads(line)['k'])
            except Exception: pass
    hits = []
    with open(path, "a") as fo:
        for k, yb in enumerate(yBs):
            if k in done: continue
            rec, sols = J.test(yb)
            out = dict(k=k, **rec)
            if sols: out['first'] = sols[0]
            fo.write(json.dumps(out) + "\n"); fo.flush()
            if sols:
                hits.append((k, yb, sols)); log(json.dumps(dict(stage="JOINT_HIT", k=k, n=len(sols))))
                if cross:
                    for js, s_ in enumerate(sols[:20]):
                        YA = J.ya_full(s_['layers'])
                        rc, cores = J.D.cross(YA, yb)
                        log(json.dumps(dict(stage="cross", k=k, joint=js, **rc)))
                        if cores:
                            H = bordered_gs_cyclic(cores[0])
                            json.dump(dict(m=m, mu=mu, core=cores[0].tolist(), hadamard=H is not None),
                                      open(path + f".core{k}_{js}.json", "w"))
                            log(json.dumps(dict(stage="CORE", k=k, joint=js, hadamard=H is not None)))
    return hits

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--dir": "", "--workers": "4", "--phase": "all", "--plant": "0", "--cross": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    m, mu, qA, qB = pos
    log = lambda s: print(s, flush=True)
    if int(opt["--plant"]):
        J = JointTest(m, mu, qA, qB)
        for sd in range(int(opt["--plant"])): log(json.dumps(J.planted(sd)))
        sys.exit()
    d = opt["--dir"]; os.makedirs(d, exist_ok=True)
    tag = f"{m}_{mu}_{qA}_{qB}"
    ys = enumerate_yB(m, mu, qA, qB, os.path.join(d, f"enum_{tag}.jsonl"), int(opt["--workers"]), log)
    if opt["--phase"] in ("all", "joint"):
        hits = joint_phase(m, mu, qA, qB, ys, os.path.join(d, f"joint_{tag}.jsonl"), log, cross=bool(int(opt["--cross"])))
        log(json.dumps(dict(stage="joint_done", yB_orbits=len(ys), hits=len(hits))))
