r"""c = 6 slices of 599 with uniform cyclotomy over F_9 (the F_{p^2} family of report/next_design_answer.md section 4,
here with the b = 5 targets) [P, M, E].

G = Q x Z_11 x F_9 (Q = Z_2 x Z_3 on axes 0, 1; Z_11 on axis 2; F_9 = F_3[i]/(i^2 + 1) as Z_3^2 on axes 3, 4),
H = Z_11 x F_9 (index 6).  Slice: M = the index-s subgroup of F_9^* (s | 4: uniform cyclotomy, integer Gauss periods)
acting on F_9 by multiplication, times a group U of per-axis multipliers on Z_11 (and optionally on Q).
Staging along F_9: stage 1 on G' = Q x Z_11 with letters alpha + (8/s) * (sum of s signs);
T1 = 4m - 4b*8 = 2216 at 0, -4b*9 = -180 on H' \ 0 (H' = Z_11), +36 off H'; the Q-level is the lead's.
Stage 2 [P] (uniform cyclotomy, s | p + 1 = 4): for every class C_i of M,  sum_k PAF_{G'}(w_{k,i}) = 4v delta,
w_{k,i} = alpha_k + B sigma_k + 3 c_{k,i}  (periods: s = 2: (1, -2), so B = -2; s = 4: (2, -1, -1, -1), B = -1).
PermSlice builds the slice from automorphism permutations so that theory/b599_cond.py and the stage-2 tools of
theory/ab_typical_levels.py (level_bits, which need G, perms, orb, parts, letters, T) work unchanged.
usage: python3 theory/b599_f9.py --s 2 --u11 sq [--stage1] [--stage2 SEEDS]"""
import sys, os, json, math, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import AG, target


def f9_mul_matrix(c, d):
    """multiplication by c + d i on F_9 = F_3[i]/(i^2 + 1), acting on (a, b) = a + b i"""
    return np.array([[c, -d], [d, c]]) % 3


def f9_elements():
    return [(c, d) for c in range(3) for d in range(3) if (c, d) != (0, 0)]


def f9_order(c, d):
    M = f9_mul_matrix(c, d); P = np.eye(2, dtype=int); o = 0
    while True:
        P = (P @ M) % 3; o += 1
        if (P == np.eye(2, dtype=int)).all(): return o


class PermSlice:
    """slice given by automorphisms of G (functions on coordinate arrays); staging along the axes in 'stage_axes'"""
    def __init__(self, mods, v, b, autos, stage_axes, nQ):
        G = AG(mods); self.G = G; self.v, self.b = v, b
        Hmask = (G.co[:, :nQ] == 0).all(1)
        self.Hmask = Hmask; self.T = target(G, v, b, Hmask)
        # group closure of the permutations
        gens = [G.idx(f(G.co)) for f in autos]
        idp = tuple(range(G.n)); S = {idp}; fr = [np.arange(G.n)]
        while fr:
            P = fr.pop()
            for g in gens:
                Q_ = P[g]; t = tuple(Q_.tolist())
                if t not in S: S.add(t); fr.append(Q_)
        self.perms = [np.array(p) for p in S]
        for P in self.perms: assert (Hmask[P] == Hmask).all()
        orb = -np.ones(G.n, np.int64); reps = []; sizes = []
        for g in range(G.n):
            if orb[g] >= 0: continue
            img = {int(P[g]) for P in self.perms}
            for x in img: orb[x] = len(reps)
            reps.append(g); sizes.append(len(img))
        self.orb, self.reps, self.sizes = orb, reps, np.array(sizes); self.r = len(reps)
        keep = [a for a in range(len(mods)) if a not in stage_axes]
        self.Gp = AG(tuple(mods[a] for a in keep)); Gp = self.Gp
        self.bco = G.co[:, keep]; self.bidx = Gp.idx(self.bco)
        self.zco = Gp.idx(np.zeros_like(self.bco)) * 0 + AG(tuple(mods[a] for a in stage_axes)).idx(G.co[:, stage_axes])
        # induced permutations on G' (the automorphisms preserve the fibre structure)
        Pp = set()
        for P in self.perms:
            img = np.zeros(Gp.n, np.int64)
            img[self.bidx] = self.bidx[P]
            Pp.add(tuple(img.tolist()))
        self.permsp = [np.array(p) for p in Pp]
        borb = -np.ones(Gp.n, np.int64); breps = []
        for x in range(Gp.n):
            if borb[x] >= 0: continue
            img = {int(P[x]) for P in self.permsp}
            for y in img: borb[y] = len(breps)
            breps.append(x)
        self.borb, self.breps = borb, breps; self.bsz = np.bincount(borb)
        self.parts = []
        for bb in breps:
            Pd = {}
            for g in np.nonzero(self.bidx == bb)[0]: Pd[int(orb[g])] = Pd.get(int(orb[g]), 0) + 1
            self.parts.append(sorted(Pd.items()))
        self.letters = []
        for pp in self.parts:
            L = {0}
            for _, sz in pp: L = {x + e * sz for x in L for e in (1, -1)}
            self.letters.append(sorted(L))
        self.T1 = np.zeros(Gp.n, np.int64); np.add.at(self.T1, self.bidx, self.T)
        self.Hpmask = np.zeros(Gp.n, bool); self.Hpmask[self.bidx[Hmask]] = True
        # Sp: for b599_cond (lag classes of G') use the induced permutations directly
        self.Sp = None
        self.qaxis = None; self.q = None

    def compress(self, X):
        X = np.asarray(X); Y = np.zeros(X.shape[:-1] + (self.Gp.n,), np.int64)
        for k in np.ndindex(X.shape[:-1]): np.add.at(Y[k], self.bidx, X[k])
        return Y


def build(s=2, u11='sq', uq=None):
    """s in {2, 4}: index of M in F_9^*; u11 in {'sq', 'neg', '1'}: U on Z_11 (squares / -1 / trivial)"""
    mods = (2, 3, 11, 3, 3)
    # primitive element of F_9: 1 + i (order 8)
    assert f9_order(1, 1) == 8
    gam = f9_mul_matrix(1, 1)
    Ms = np.linalg.matrix_power(gam, s) % 3                  # generator of the index-s subgroup
    autos = []
    def mulF9(A):
        def f(co):
            out = co.copy(); out[:, 3:5] = (co[:, 3:5] @ A.T) % 3; return out
        return f
    autos.append(mulF9(Ms))
    def mul11(u):
        def f(co):
            out = co.copy(); out[:, 2] = (co[:, 2] * u) % 11; return out
        return f
    if u11 == 'sq': autos.append(mul11(4))                  # 4 generates the squares mod 11 (order 5)
    elif u11 == 'neg': autos.append(mul11(10))
    elif u11 == 'all': autos.append(mul11(2))
    if uq == 'neg3':
        def f(co):
            out = co.copy(); out[:, 1] = (-co[:, 1]) % 3; return out
        autos.append(f)
    sl = PermSlice(mods, 599, 5, autos, (3, 4), 2)
    cfg = dict(v=599, c=6, b=5, h=99, m=594, Q=(2, 3), H1=('Z11',), q=9, s=s, f=8 // s, mods=mods, nQ=2,
               U=[('u11', u11)] + ([('uq', uq)] if uq else []))
    return cfg, sl


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    s = int(opt('--s', '2')); u11 = opt('--u11', 'sq'); uq = opt('--uq', None)
    cfg, sl = build(s, u11, uq)
    print(json.dumps(dict(s=s, u11=u11, uq=uq, orbits=int(sl.r), Gp_orbits=len(sl.breps), letters=sorted({tuple(L) for L in sl.letters}),
                          T1=sorted(set(sl.T1.tolist())), T1_0=int(sl.T1[0]), group=len(sl.perms))), flush=True)


def stage1(s, u11, N=2000, outp=''):
    import b599_cond as BC
    cfg, sl = build(s, u11)
    qc = AG((2, 3)).idx(sl.Gp.co[:, :2])
    C = BC.Cond(cfg, sl=sl, Qmods=(2, 3), qc=qc)
    fo = open(outp, 'a') if outp else None
    classes = C.qlevel_quads(); tot = 0.0; n = 0; best = None
    for ci, (key, us_) in enumerate(sorted(classes.items())):
        for zm in C.rows_for(us_):
            Z = np.array([[zm[k] + 2 * us_[k][q] for q in range(C.c)] for k in range(4)], np.int64)
            r = C.budget(Z, N=N); n += 1
            if r.get('log2E') is not None and np.isfinite(r['log2E']):
                tot += 2 ** r['log2E']
                if best is None or r['log2E'] > best[0]: best = (r['log2E'], ci, list(zm))
            rec = dict(s=s, u11=u11, **{'class': ci}, shapes=[list(u) for u in us_], zmin=list(zm), rows=[int(x) for x in Z.sum(1)], **r)
            if fo: fo.write(json.dumps(rec) + "\n"); fo.flush()
    return dict(s=s, u11=u11, zquads=n, sum_log2E=round(math.log2(tot), 2) if tot else None, best=best)


def stage2(s, u11, seed=3, N=3000):
    from ab_typical_levels import typical_y_sampler, level_bits
    from ab_levels import sd_freedom, subgroup_mask_axes
    cfg, sl = build(s, u11)
    Y = typical_y_sampler(sl, seed)()
    lv = {'S(axis0)': subgroup_mask_axes(sl.Gp, [0])}
    r = level_bits(sl, Y, lv, N=N, seed=seed)
    sd = sd_freedom(sl, Y)
    CS = r['levels']['S(axis0)']['bits']; HS = r['H2'] - sd['h_D']; Dd = r['C2'] - CS - sd['h_D']
    return dict(s=s, u11=u11, seed=seed, H2=r['H2'], C2=r['C2'], rank=r['rank'], per_y=round(r['H2'] - r['C2'], 2),
                C_S=CS, h_D=sd['h_D'], D_def=round(Dd, 2), cost=round(Dd + CS - HS / 2, 2))
