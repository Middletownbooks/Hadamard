r"""Cross level of a doubly-staged slice, given a joint stage-1 solution (y_A, y_B) [M for the lists, E for estimates].

Setting as in theory/nd_joint_exhaust.py (m = 2 m', M fixes Z_2, qA and qB odd).  Given y_B and y_A of every block,
the labels of the M-orbits of each Z_2-layer e are constrained only by that layer's y_B and y_A values.  Per block and
layer this script lists every labelling (forward DP over the layer's orbits, then backward expansion), so the remaining
problem is exact: choose per block one labelling per layer (x on Z_m) with
    sum_k PAF_{m'}(sigma_k) = TS = (m - 1, -2, ..., -2),   sum_k PAF_{m'}(delta_k) = TD = (m + 1, 0, ..., 0),
sigma = (x|layer0 + x|layer1)/2, delta = (x|layer0 - x|layer1)/2 on Z_{m'} (the S/D split of Z_m = Z_2 x Z_{m'}),
at one lag per <M, -1>-orbit of Z_{m'}.
Estimates [E]: H = sum_k log2(|L_0,k| |L_1,k|) bits of freedom; C = lattice Gaussian (q523_gauss) of the summed
signatures under uniform per-block choices; log2E = H - C per joint solution.  Also the S and D levels separately and
the number of distinct sigma (resp. delta) per block, the input of an exact S-then-D join.
usage: python3 theory/nd_cross.py m mu qA qB JOINTFILE [--k K] [--n 20] [--N 4000]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_joint_exhaust import JointTest
from nd_s1_sd import lag_reps
from q523_gauss import gauss_bits

class Cross:
    def __init__(self, m, mu, qA, qB):
        self.J = J = JointTest(m, mu, qA, qB); D = J.D; self.D = D
        self.m = m; mp = m // 2; self.mp = mp
        self.E = [np.array([next(g for g in (b, b + mp) if g % 2 == e) for b in range(mp)]) for e in (0, 1)]
        self.Mp = sorted({a % mp for a in D.Mg})
        self.lags = lag_reps(mp, self.Mp)
        self.TS = np.array([m - 1] + [-2] * (len(self.lags) - 1), np.int64)
        self.TD = np.array([m + 1] + [0] * (len(self.lags) - 1), np.int64)
        # per layer: orbit ids, and for each orbit its y_B / y_A column contributions
        self.lo = J.layer_orbs

    def layer_labellings(self, e, yb_full, ya_full, cap=2 * 10 ** 6):
        """all labellings (rows = orbits of layer e, in self.lo[e] order) with the layer's y_B and y_A values"""
        D = self.D; J = self.J; orbs = self.lo[e]
        CB = D.cB[np.ix_(orbs, D.kB)][:, J.Bcols[e]]; CA = D.cA[np.ix_(orbs, D.kA)][:, J.Acols[e]]
        tb = np.asarray(yb_full)[D.kB][J.Bcols[e]]; ta = np.asarray(ya_full)[D.kA][J.Acols[e]]
        C = np.concatenate([CB, CA], 1); T = np.concatenate([tb, ta])
        n = len(orbs); nb = CB.shape[1]
        # orbits row by row (by their y_B column), so that each y_B coordinate closes before the next opens
        order = sorted(range(n), key=lambda o: (np.nonzero(CB[o])[0].min() if CB[o].any() else nb, np.nonzero(CA[o])[0].min() if CA[o].any() else -1))
        C = C[order]
        maxrem = np.zeros((n + 1, C.shape[1]), np.int64)
        for p in range(n - 1, -1, -1): maxrem[p] = maxrem[p + 1] + np.abs(C[p])
        touched = np.zeros((n + 1, C.shape[1]), bool)
        for p in range(n): touched[p + 1] = touched[p] | (C[p] != 0)
        def encoder(p):
            # untouched columns still hold their target (constant): radix 1; finished ones are 0 (maxrem = 0): radix 1
            off = np.where(touched[p], maxrem[p], -T); rad = np.where(touched[p], 2 * maxrem[p] + 1, 1)
            assert float(np.prod(rad.astype(float))) < 2 ** 62, "state encoding overflow"
            mult = np.cumprod(np.concatenate([[1], rad[:-1]])).astype(np.int64)
            return lambda S: (S + off) @ mult
        encs = [encoder(p) for p in range(n + 1)]
        if (np.abs(T) > maxrem[0]).any(): return np.zeros((0, n), np.int8)
        levels = [np.array([encs[0](T[None, :])[0]])]; S = T[None, :].astype(np.int64)
        for p in range(n):
            new = np.concatenate([S - C[p], S + C[p]])       # label +1 removes C[p] from the residual, -1 adds it
            ok = (np.abs(new) <= maxrem[p + 1]).all(1); new = new[ok]
            if len(new) == 0: return np.zeros((0, n), np.int8)
            k = encs[p + 1](new); k, first = np.unique(k, return_index=True); S = new[first]; levels.append(k)
        if not (S == 0).all(1).any(): return np.zeros((0, n), np.int8)
        # backward expansion from the zero residual
        Sf = np.zeros((1, C.shape[1]), np.int64); Lab = np.zeros((1, 0), np.int8)
        for p in range(n - 1, -1, -1):
            cand = []; labs = []
            for lab in (1, -1):
                prev = Sf + lab * C[p]                        # residual before orbit p
                okb = (np.abs(prev) <= maxrem[p]).all(1)
                ok = okb.copy(); ok[okb] = np.isin(encs[p](prev[okb]), levels[p])
                if ok.any():
                    cand.append(prev[ok]); labs.append(np.concatenate([np.full((ok.sum(), 1), lab, np.int8), Lab[ok]], 1))
            Sf = np.concatenate(cand); Lab = np.concatenate(labs)
            if len(Lab) > cap: raise RuntimeError(f"layer list larger than {cap}")
        assert (Sf == T[None, :]).all()
        out = np.zeros((len(Lab), n), np.int8); out[:, order] = Lab
        return out

    def layer_function(self, e, L):
        """x restricted to layer e as functions on Z_{m'} (rows), from labellings L (orbits of layer e)"""
        D = self.D; orbs = self.lo[e]
        pos = np.full(D.nor, -1, np.int64); pos[orbs] = np.arange(len(orbs))
        idx = pos[D.orb[self.E[e]]]; assert (idx >= 0).all()
        return L[:, idx].astype(np.int64)

    def paf(self, X):
        f = np.fft.rfft(X.astype(np.float64), axis=1)
        return np.rint(np.fft.irfft(f * f.conj(), n=self.mp, axis=1)[:, self.lags]).astype(np.int64)

    def lists(self, yb, ya):
        out = []
        for k in range(4):
            L0 = self.layer_labellings(0, yb[k], ya[k]); L1 = self.layer_labellings(1, yb[k], ya[k])
            out.append((L0, L1))
        return out

    def estimate(self, yb, ya, N=4000, seed=1):
        rng = np.random.default_rng(seed)
        Ls = self.lists(yb, ya)
        cnt = [(len(a), len(b)) for a, b in Ls]
        rec = dict(layer_counts=cnt)
        if min(min(c) for c in cnt) == 0: rec.update(H=None, note='no labelling'); return rec, Ls
        H = sum(math.log2(a) + math.log2(b) for a, b in cnt); rec['H'] = round(H, 2)
        VS = 0; VD = 0; nsd = []
        for k, (L0, L1) in enumerate(Ls):
            X0 = self.layer_function(0, L0); X1 = self.layer_function(1, L1)
            i0 = rng.integers(len(X0), size=N); i1 = rng.integers(len(X1), size=N)
            sg = (X0[i0] + X1[i1]) // 2; dl = (X0[i0] - X1[i1]) // 2
            VS = VS + self.paf(sg); VD = VD + self.paf(dl)
        gS = gauss_bits(VS, self.TS); gD = gauss_bits(VD, self.TD)
        gF = gauss_bits(np.concatenate([VS, VD], 1), np.concatenate([self.TS, self.TD]))
        rec.update(C=round(gF['bits'], 2) if gF.get('in_lattice') else None, rank=gF.get('rank'),
                   CS=round(gS['bits'], 2) if gS.get('in_lattice') else None, CD=round(gD['bits'], 2) if gD.get('in_lattice') else None,
                   mahal=gF.get('mahal'), log2E=round(H - gF['bits'], 2) if gF.get('in_lattice') else None)
        return rec, Ls

def load_joint(path, k):
    for line in open(path):
        r = json.loads(line)
        if r['k'] == k: return r
    return None

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--k": "0", "--n": "10", "--N": "4000"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    m, mu, qA, qB = map(int, pos[:4]); jf = pos[4]
    X = Cross(m, mu, qA, qB)
    from nd_joint_exhaust import enumerate_yB
    tag = f"{m}_{mu}_{qA}_{qB}"; d = os.path.dirname(jf)
    ys = enumerate_yB(m, mu, qA, qB, os.path.join(d, f"enum_{tag}.jsonl"), 1, log=lambda s: None)
    k = int(opt["--k"]); yb = ys[k]
    rec, sols = X.J.test(yb)
    print(json.dumps(dict(k=k, joint=rec.get('joint'))), flush=True)
    rng = np.random.default_rng(5)
    pick = rng.choice(len(sols), size=min(int(opt["--n"]), len(sols)), replace=False) if sols else []
    for j in pick:
        ya = X.J.ya_full(sols[j]['layers'])
        t0 = time.time(); est, _ = X.estimate(yb, ya, int(opt["--N"]))
        print(json.dumps(dict(k=k, joint_index=int(j), **est, secs=round(time.time() - t0, 1))), flush=True)

# ---------------------------------------------------------------- exact cross level
def orbit_maps(X):
    """the M-orbits of Z_{m'} (reps, sizes, index of every point) and, per layer, the layer-orbit -> Z_{m'}-orbit map"""
    D = X.D; mp = X.mp
    oi = -np.ones(mp, np.int64); reps = []
    for b in range(mp):
        if oi[b] >= 0: continue
        for a in X.Mp: oi[b * a % mp] = len(reps)
        reps.append(b)
    maps = []
    for e in (0, 1):
        orbs = X.lo[e]; pos = np.full(D.nor, -1, np.int64); pos[orbs] = np.arange(len(orbs))
        lay = pos[D.orb[X.E[e][reps]]]                     # layer-orbit index of each Z_{m'} orbit rep
        assert (lay >= 0).all() and len(set(lay.tolist())) == len(reps) == len(orbs)
        maps.append(lay)
    return np.array(reps), oi, maps

def expand_paf(X, V, oi, chunk=20000):
    """V: values at the Z_{m'} orbit reps (rows) -> PAF on Z_{m'} at X.lags"""
    out = np.empty((len(V), len(X.lags)), np.int64)
    for st in range(0, len(V), chunk):
        out[st:st + chunk] = X.paf(V[st:st + chunk][:, oi].astype(np.int64))
    return out

def block_pairs(X, L0, L1, maps):
    """for one block: distinct sigma (at orbit reps, int8) and the grouped (i0, i1) pairs"""
    A0 = L0[:, maps[0]].astype(np.int8); A1 = L1[:, maps[1]].astype(np.int8)
    I0 = np.repeat(np.arange(len(A0)), len(A1)); I1 = np.tile(np.arange(len(A1)), len(A0))
    sg = ((A0[I0].astype(np.int16) + A1[I1]) // 2).astype(np.int8)
    u, inv = np.unique(sg, axis=0, return_inverse=True); inv = inv.ravel()
    o = np.argsort(inv, kind='stable'); bd = np.searchsorted(inv[o], np.arange(len(u) + 1))
    return dict(A0=A0, A1=A1, I0=I0[o], I1=I1[o], bd=bd, U=u)

def solve_cross(X, yb, ya, maxcores=10, log=None, plant_T=None, Ls=None, s_only=False):
    """exact: every core (x_0..x_3 on Z_m) with the given compressions, by an S-join over the distinct sigma per block
    (theory/nd_ejoin.py, energy = ||sigma||^2) followed by a full (S and D) join over the pairs whose sigma occurs in
    some S-solution.  Every core is checked by the full cyclic PAF.  plant_T: (TS, TD) to replace the slice target."""
    from nd_ejoin import ejoin4
    t0 = time.time()
    TS, TD = (X.TS, X.TD) if plant_T is None else plant_T
    if Ls is None: Ls = X.lists(yb, ya)
    rec = dict(layer_counts=[(len(a), len(b)) for a, b in Ls])
    if min(min(c) for c in rec['layer_counts']) == 0: rec.update(cores=0, note='no labelling'); return rec, []
    reps, oi, maps = orbit_maps(X)
    P = [block_pairs(X, L0, L1, maps) for L0, L1 in Ls]
    rec['sigma_distinct'] = [len(p_['U']) for p_ in P]
    Ssig = [expand_paf(X, p_['U'], oi) for p_ in P]
    Sh, st = ejoin4(Ssig, TS, log=log)
    rec.update(S_solutions=len(Sh), S_work_log2=st['work_log2'], S_secs=st['secs'])
    cores = []
    if Sh and not s_only:
        Fl = []; keep = []
        for k, p_ in enumerate(P):
            used = np.zeros(len(p_['U']), bool); used[[t[k] for t in Sh]] = True
            ii = np.concatenate([np.arange(p_['bd'][j], p_['bd'][j + 1]) for j in np.nonzero(used)[0]])
            keep.append(ii)
            a0 = p_['A0'][p_['I0'][ii]].astype(np.int64); a1 = p_['A1'][p_['I1'][ii]].astype(np.int64)
            Fl.append(np.concatenate([expand_paf(X, (a0 + a1) // 2, oi), expand_paf(X, (a0 - a1) // 2, oi)], 1))
        rec['full_lists'] = [len(x) for x in keep]
        Fh, _ = ejoin4(Fl, np.concatenate([TS, TD]))
        for t in Fh:
            xs = []
            for k, p_ in enumerate(P):
                j = keep[k][t[k]]; a0 = p_['A0'][p_['I0'][j]].astype(np.int64); a1 = p_['A1'][p_['I1'][j]].astype(np.int64)
                x = np.zeros(X.m, np.int64); x[X.E[0]] = a0[oi]; x[X.E[1]] = a1[oi]; xs.append(x)
            xs = np.array(xs)
            Pf = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x.astype(float))) ** 2).real).astype(np.int64) for x in xs)
            ok = plant_T is not None or bool(Pf[0] == 4 * X.m and (Pf[1:] == -4).all())
            if ok: cores.append(xs)
            if len(cores) >= maxcores: break
    rec.update(cores=len(cores), secs=round(time.time() - t0, 1))
    return rec, cores

def planted_cross(X, yb, ya, seed, Ls=None):
    """plant: one random labelling per block and layer from the exact lists; targets = its own S and D sums"""
    rng = np.random.default_rng(seed)
    if Ls is None: Ls = X.lists(yb, ya)
    reps, oi, maps = orbit_maps(X)
    pick = [(int(rng.integers(len(a))), int(rng.integers(len(b)))) for a, b in Ls]
    TS = 0; TD = 0; xs = []
    for (L0, L1), (i, j) in zip(Ls, pick):
        a0 = L0[i][maps[0]].astype(np.int64); a1 = L1[j][maps[1]].astype(np.int64)
        TS = TS + expand_paf(X, ((a0 + a1) // 2)[None], oi)[0]; TD = TD + expand_paf(X, ((a0 - a1) // 2)[None], oi)[0]
        x = np.zeros(X.m, np.int64); x[X.E[0]] = a0[oi]; x[X.E[1]] = a1[oi]; xs.append(x)
    rec, cores = solve_cross(X, yb, ya, maxcores=10 ** 6, plant_T=(TS, TD), Ls=Ls)
    found = any((c == np.array(xs)).all() for c in cores)
    return dict(seed=seed, **rec, planted_found=found)
