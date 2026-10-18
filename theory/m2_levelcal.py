r"""Stage-2 level budget of a REAL stage-1 solution at the true target and at the sample mean [E, M for H2 and h_D].

The typical-y budget (theory/ab_typical_levels.py, used for every b >= 3 cost in report/bge3_stage1_answer.md 4 and 6)
reads each lattice Gaussian at the sample mean (ab_gauss.gauss_bits_at_mean), so it charges every equation as if its
target sat at the centre of the lift distribution.  For 599 x4 that hides an empty slice (theory/b599x4_parseval.py) and
for the Z3^2 variants about 10 bits (theory/b599x4_z3a_rigidity.py).  Here both readings are taken on the same lift
samples of a real y: the 523 production y (logs/q523/y523_best.json, mapped to Z2 x Z9 by CRT) and the y of the
published 347 core (its compression), the two slices on which the cost model was calibrated.
Levels: S = characters trivial on the Z2 axis 0 of G' (the S/D split), and S + L'' for the U-stable subgroups L'' of
G'' of prime order (candidate middle levels).  D_def = C2 - C_S - h_D (report/paley_523_answer.md 3.4).
usage: python3 theory/m2_levelcal.py 523|347 [--N 6000] [--seed 1] [--out FILE]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json
from ab_levels import sd_freedom
from q523_gauss import gauss_bits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(which):
    if which == '523':
        cfg = dict(v=523, c=1, b=1, h=None, m=522, Q=[], H1=[2, 9], q=29, s=4, f=7, mods=(2, 9, 29), nQ=0, U=[])
        sl = W.build(cfg)
        d = json.load(open(os.path.join(ROOT, 'logs/q523/y523_best.json')))
        Yc = np.array(d['y'], np.int64)                         # on cyclic Z18
        Y = np.zeros((4, sl.Gp.n), np.int64)
        for g in range(18): Y[:, sl.Gp.idx(np.array([g % 2, g % 9]))] = Yc[:, g]
    else:
        import ab_assemble as AS
        from ab_verify_published import recoordinate
        cfg = cfg_from_json('@' + os.path.join(ROOT, 'logs/ab/cfg_347_s6.json')); sl = W.build(cfg)
        T = AS.load_template347()
        XF = recoordinate(T['F'], (2, 3, 57), (2, 3, 3, 19), [(0, 2), (1, 3), (2, 3), (2, 19)])
        assert sl.is_invariant(XF) and (sl.G.paf(XF).sum(0) == sl.T).all()
        Y = sl.compress(XF)
    assert (sl.Gp.paf(Y).sum(0) == sl.T1).all(), "not a stage-1 solution"
    return cfg, sl, Y


def subgroup_mask(Gp, gens):
    S = {0}; fr = [0]; gi = [int(Gp.idx(np.array(g))) for g in gens]
    while fr:
        x = fr.pop()
        for g in gi:
            y = int(Gp.add(x, g))
            if y not in S: S.add(y); fr.append(y)
    m = np.zeros(Gp.n, bool); m[list(S)] = True; return m


def candidate_levels(sl):
    """S (axis 0) and S + each U-stable subgroup of prime order of G'' = axes 1.. of G'"""
    Gp = sl.Gp; nax = len(Gp.mods)
    e0 = tuple([1] + [0] * (nax - 1))
    levels = {'S': subgroup_mask(Gp, [e0])}
    Up = [Gp.act(u) for u in sl.Sp]
    seen = set()
    for g in range(Gp.n):
        if Gp.co[g, 0] != 0 or g == 0: continue
        o = Gp.order(g)
        if any(o % p == 0 and o != p for p in range(2, o)): continue      # prime order only
        m = subgroup_mask(Gp, [tuple(Gp.co[g])])
        if not all((m[P] == m).all() for P in Up): continue
        key = tuple(np.nonzero(m)[0])
        if key in seen: continue
        seen.add(key)
        levels[f'S+<{tuple(int(x) for x in Gp.co[g][1:])}>'] = subgroup_mask(Gp, [e0, tuple(Gp.co[g])])
    return levels


def run(which, N=6000, seed=1):
    cfg, sl, Y = load(which)
    G = sl.G; rng = np.random.default_rng(seed)
    pats = []; H2 = 0.0
    for k in range(4):
        pk = []
        for j, bb in enumerate(sl.breps):
            P, _ = AGs.fibre_patterns(sl, j, int(Y[k, bb])); pk.append(P); H2 += math.log2(len(P))
        pats.append(pk)
    V = np.zeros((N, G.n), np.int64); done = 0
    while done < N:
        c = min(1000, N - done); acc = np.zeros((c, G.n), np.int64)
        for k in range(4):
            lab = np.zeros((c, sl.r), np.int8)
            for j in range(len(sl.breps)):
                P = pats[k][j]; pk_ = P[rng.integers(len(P), size=c)]
                for i, (o, _) in enumerate(sl.parts[j]): lab[:, o] = pk_[:, i]
            acc += G.paf(lab[:, sl.orb])
        V[done:done + c] = acc; done += c
    reps, _ = AGs.lag_reps(G, sl.perms)
    def both(VL, TL):
        gt = gauss_bits(VL, TL); gm = AGs.gauss_bits_at_mean(VL)
        return dict(true=gt.get('bits'), at_mean=gm.get('bits_at_mean'), rank=gt.get('rank'), mahal=gt.get('mahal'),
                    in_lattice=gt.get('in_lattice'))
    out = dict(slice=which, y_letters=sorted({int(a) for a in np.abs(Y).ravel()}), H2=round(H2, 2), N=N, seed=seed)
    out['C2'] = both(V[:, reps], sl.T[reps])
    out['levels'] = {}
    for name, m in candidate_levels(sl).items():
        Lmask = m[sl.bidx] & (sl.zco == 0)
        cos, nC = AGs.subgroup_compress_index(G, Lmask)
        VL = np.zeros((N, nC), np.int64); TL = np.zeros(nC, np.int64)
        np.add.at(VL.T, cos, V.T); np.add.at(TL, cos, sl.T)
        out['levels'][name] = both(VL, TL)
    sd = sd_freedom(sl, Y); out['h_D'] = sd['h_D']; out['H_S'] = round(H2 - sd['h_D'], 2)
    for key in ('true', 'at_mean'):
        C2 = out['C2'][key]; CS = out['levels']['S'][key]
        out[f'D_def_{key}'] = round(C2 - CS - sd['h_D'], 2) if C2 is not None and CS is not None else None
        out[f'stage2_per_y_{key}'] = round(H2 - C2, 2) if C2 is not None else None
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    r = run(a[0], int(opt('--N', '6000')), int(opt('--seed', '1')))
    print(json.dumps(r, indent=1))
    if '--out' in a: json.dump(r, open(a[a.index('--out') + 1], 'w'), indent=1)
