r"""Typical-y stage-2 level budget of a GP-type slice [E]: H2, C2, the subgroup levels and the S/D deficit.

For slices whose stage 1 has no enumerated solution yet: y is drawn with the stage-1 lag-0 energy and the row
equation exact (letters i.i.d. per G'-orbit from the energy-tilted uniform law, quadruples accepted when
sum_k E_k = T1(0) and sum_k r_k^2 = R4), and every lattice Gaussian is read at the sample mean
(ab_gauss.gauss_bits_at_mean, as in ab_gauss.stage2_typical), so the y need not solve stage 1.
Levels: characters of G' trivial on a subgroup L (every nontrivial character of Z_q).  S/D split on the Z_2 axis 0 of
G' (ab_levels.sd_freedom): h_D = D-freedom once S is fixed; D_def = C(D|S) - h_D = S-solutions needed per core
(report/paley_523_answer.md 3.4; 347: 16.2, 523: 26.1).
usage: python3 theory/ab_typical_levels.py CFGJSON [--seed 3] [--N 5000] [--memgb 3.2]"""
import sys, os, json, math, resource
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json
from ab_levels import subgroup_mask_axes, sd_freedom


def typical_y_sampler(sl, seed=1):
    rng = np.random.default_rng(seed)
    T0 = int(sl.T1[0]); R4 = int(sl.T1.sum())
    w = sl.bsz.astype(np.int64); L = [np.array(sorted(x), np.int64) for x in sl.letters]; nb = len(sl.breps)
    def law(lam):
        P = []
        for j in range(nb):
            e = w[j] * L[j] ** 2; p = np.exp(-lam * (e - e.min())); P.append(p / p.sum())
        return P
    def mean_E(lam):
        return sum(float((p * w[j] * L[j] ** 2).sum()) for j, p in enumerate(law(lam)))
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if mean_E(mid) > T0 / 4: lo = mid
        else: hi = mid
    P = law((lo + hi) / 2)
    def draw(n=200000, tries=1000):
        for _ in range(tries):
            Yo = np.stack([rng.choice(L[j], size=4 * n, p=P[j]) for j in range(nb)], 1)
            E = (Yo ** 2 * w).sum(1).reshape(4, n); r = (Yo * w).sum(1).reshape(4, n)
            idx = np.nonzero((E.sum(0) == T0) & ((r ** 2).sum(0) == R4))[0]
            if len(idx): return np.array([Yo[k * n + idx[0]][sl.borb] for k in range(4)])
        return None
    return draw


def level_bits(sl, Y, levels, N=5000, seed=1, chunk=1000):
    rng = np.random.default_rng(seed); G = sl.G
    pats = []; H2 = 0.0
    for k in range(4):
        pk = []
        for j, bb in enumerate(sl.breps):
            P, _ = AGs.fibre_patterns(sl, j, int(Y[k, bb]))
            pk.append(P); H2 += math.log2(len(P))
        pats.append(pk)
    reps, _ = AGs.lag_reps(G, sl.perms)
    V = np.zeros((N, G.n), np.int64); done = 0
    while done < N:
        c = min(chunk, N - done); acc = np.zeros((c, G.n), np.int64)
        for k in range(4):
            lab = np.zeros((c, sl.r), np.int8)
            for j in range(len(sl.breps)):
                P = pats[k][j]; pk_ = P[rng.integers(len(P), size=c)]
                for i, (o, _) in enumerate(sl.parts[j]): lab[:, o] = pk_[:, i]
            acc += G.paf(lab[:, sl.orb])
        V[done:done + c] = acc; done += c
    gb = AGs.gauss_bits_at_mean(V[:, reps])
    out = dict(H2=round(H2, 2), C2=gb.get('bits_at_mean'), rank=gb.get('rank'), levels={})
    for name, Lmask_p in levels.items():
        Lmask = Lmask_p[sl.bidx] & (sl.zco == 0)
        cos, nC = AGs.subgroup_compress_index(G, Lmask)
        VL = np.zeros((N, nC), np.int64); np.add.at(VL.T, cos, V.T)
        gl = AGs.gauss_bits_at_mean(VL)
        out['levels'][name] = dict(bits=gl.get('bits_at_mean'), rank=gl.get('rank'))
    return out


def budget(cfg, seed=3, N=5000):
    sl = W.build(cfg); Y = typical_y_sampler(sl, seed)()
    Gp = sl.Gp; nax = len(Gp.mods)
    levels = {'S(axis0)': subgroup_mask_axes(Gp, [0])}
    for ax in range(1, nax): levels[f'axis{ax}'] = subgroup_mask_axes(Gp, [ax])
    r = level_bits(sl, Y, levels, N=N, seed=seed)
    sd = sd_freedom(sl, Y); r.update(h_D=sd['h_D'], H_S_distinct=sd['H_S_distinct'])
    CS = r['levels']['S(axis0)']['bits']
    r['C_D_given_S'] = round(r['C2'] - CS, 2); r['D_def'] = round(r['C2'] - CS - sd['h_D'], 2)
    r['stage2_log2E'] = round(r['H2'] - r['C2'], 2)
    r.update(v=cfg['v'], c=cfg['c'], q=cfg['q'], s=cfg['s'], H1=list(cfg['H1']), U=[list(u) for u in cfg['U']],
             orbits=int(sl.r), L1=len(sl.breps), seed=seed, N=N)
    return r


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    memgb = float(opt('--memgb', '3.2'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    print(json.dumps(budget(cfg_from_json(a[0]), int(opt('--seed', '3')), int(opt('--N', '5000')))))
