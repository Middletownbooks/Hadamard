r"""Gaussian estimates for slices of theory/ab_lib.Slice (any abelian G, b = 0 / 1 / >= 3) [E, counts exact].

whole(sl)          plain whole-slice: uniform S-invariant x, lattice Gaussian (theory/q523_gauss.py) of
                   sum_k PAF at the <S, -1>-lag-class representatives;  log2E = 4r - bits.
stage1(sl)         exact count N_y of S-invariant y-quadruples on G' with the stage-1 lag-0 energy T1(0) and
                   sum_k r_k^2 = 4v (DP over G'-orbits, exact integers), times the lattice Gaussian of the remaining
                   stage-1 lag classes under uniform sampling from that set;  also the x-weighted count
                   (configurations x whose compression meets lag 0 and rows).
stage2(sl, Y)      given a stage-1 quadruple Y: H2 = log2 (number of x with compression Y) (exact), and the lattice
                   Gaussian of the full core equations under uniform x given Y (stage-1 combinations drop out as
                   rank deficiency);  stage2_log2E = H2 - bits.  levels: the same restricted to the characters
                   trivial on a subgroup L <= G' (compression along L), e.g. the S level of an S/D split.
These are the estimators of theory/nd_compress_est.py / nd_stage2_est.py generalised; the 523 and 647 calibrations
(report/paley_523_answer.md 3.4-3.5, report/paley_647_answer.md 2) apply to the same estimator family."""
import sys, os, math, json
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits
from ab_lib import AG


def lag_reps(G, perms):
    negp = G.idx(-G.co)
    seen = np.zeros(G.n, bool); reps = []; cls = -np.ones(G.n, np.int64)
    for g in range(G.n):
        if seen[g]: continue
        ci = len(reps); reps.append(g)
        for P in perms:
            for x in (P[g], negp[P[g]]):
                seen[x] = True; cls[x] = ci
    return reps, cls


def whole(sl, N=20000, seed=1, chunk=2000):
    rng = np.random.default_rng(seed)
    reps, _ = lag_reps(sl.G, sl.perms)
    V = np.zeros((N, len(reps)), np.int64); done = 0
    while done < N:
        c = min(chunk, N - done)
        acc = np.zeros((c, sl.G.n), np.int64)
        for k in range(4):
            lab = rng.choice(np.array([-1, 1], np.int8), size=(c, sl.r))
            acc += sl.G.paf(lab[:, sl.orb])
        V[done:done + c] = acc[:, reps]; done += c
    gb = gauss_bits(V, sl.T[reps])
    return dict(r=int(sl.r), lag_classes=len(reps), bits=gb.get('bits'), rank=gb.get('rank'),
                in_lattice=gb.get('in_lattice'), log2E=round(4 * sl.r - gb['bits'], 2) if gb.get('in_lattice') else None)


def fibre_patterns(sl, j, y):
    """all sign vectors over the parts of fibre-orbit j with signed size-sum y; returns (patterns, part sizes)"""
    sizes = [s_ for _, s_ in sl.parts[j]]
    n = len(sizes)
    if n > 24: raise ValueError("too many parts")
    # DP-free enumeration with meet in the middle for moderate n
    out = []
    if n <= 20:
        idx = np.arange(1 << n, dtype=np.int64)
        sg = (2 * ((idx[:, None] >> np.arange(n)) & 1) - 1).astype(np.int8)
        tot = sg.astype(np.int64) @ np.array(sizes, np.int64)
        return sg[tot == y], sizes
    raise ValueError("fibre too large")


def stage1_counts(sl):
    """per block DP over G'-orbits: (E, r) -> (number of distinct y, x-weighted count)"""
    w = sl.bsz; L = sl.letters
    # x-multiplicity of letter y at fibre-orbit j: number of part-sign patterns with that sum
    mult = []
    for j in range(len(sl.breps)):
        sizes = [s_ for _, s_ in sl.parts[j]]
        cnt = defaultdict(int); cnt[0] = 1
        for sz in sizes:
            new = defaultdict(int)
            for t, c in cnt.items(): new[t + sz] += c; new[t - sz] += c
            cnt = new
        mult.append(dict(cnt))
    T0 = int(sl.T1[0])
    tab = {(0, 0): (1, 1)}
    for j in range(len(sl.breps)):
        new = defaultdict(lambda: [0, 0])
        for (E, r), (n1, nx) in tab.items():
            for y in L[j]:
                E2 = E + int(w[j]) * y * y
                if E2 > T0: continue
                key = (E2, r + int(w[j]) * y)
                new[key][0] += n1; new[key][1] += nx * mult[j][y] ** 1
        tab = {k: tuple(vv) for k, vv in new.items()}
    return tab, mult


def stage1(sl, N=20000, seed=1):
    rng = np.random.default_rng(seed)
    T0 = int(sl.T1[0]); v = sl.v
    R4 = int(sl.T1.sum())                      # sum_k r_k^2 (= 4v for b = 0 and b >= 3, 4 for b = 1)
    tab, mult = stage1_counts(sl)
    # combine 4 blocks: sum E = T0, sum r^2 = 4v  (ordered blocks)
    byE = defaultdict(list)
    for (E, r), (n1, nx) in tab.items(): byE[E].append((r, n1, nx))
    Es = sorted(byE)
    # pair tables: (E, r^2-sum) -> counts
    def pair_table():
        P = defaultdict(lambda: [0, 0])
        for E1 in Es:
            for E2 in Es:
                if E1 + E2 > T0: break
                for (r1, a1, x1) in byE[E1]:
                    for (r2, a2, x2) in byE[E2]:
                        P[(E1 + E2, r1 * r1 + r2 * r2)][0] += a1 * a2
                        P[(E1 + E2, r1 * r1 + r2 * r2)][1] += x1 * x2
        return P
    P = pair_table()
    Ny = 0; Nx = 0; splits = []
    for (E, R), (a, x) in P.items():
        k2 = (T0 - E, R4 - R)
        if k2 in P:
            Ny += a * P[k2][0]; Nx += x * P[k2][1]
            splits.append((E, R, a * P[k2][0]))
    out = dict(T0=T0, log2Ny=round(math.log2(Ny), 2) if Ny else None, log2Nx=round(math.log2(Nx), 2) if Nx else None)
    if not Ny: out['stage1_log2E'] = -math.inf; return out
    # uniform sampler over the set: choose (block energies, rows) by weight, then per block a uniform y given (E, r)
    # backward DP tables for per-block sampling given (E, r)
    w = sl.bsz; L = sl.letters; nb = len(sl.breps)
    back = [None] * (nb + 1); back[nb] = {(0, 0): 1}
    for j in range(nb - 1, -1, -1):
        nd = defaultdict(int)
        for (E, r), c in back[j + 1].items():
            for y in L[j]:
                E2 = E + int(w[j]) * y * y
                if E2 <= T0: nd[(E2, r + int(w[j]) * y)] += c
        back[j] = dict(nd)
    def sample_block(E, r):
        ys = []
        for j in range(nb):
            opts = []; wts = []
            for y in L[j]:
                e2 = E - int(w[j]) * y * y; r2 = r - int(w[j]) * y
                c = back[j + 1].get((e2, r2), 0)
                if c: opts.append(y); wts.append(c)
            p = np.array(wts, float); p /= p.sum()
            y = opts[rng.choice(len(opts), p=p)]; ys.append(y)
            E -= int(w[j]) * y * y; r -= int(w[j]) * y
        return np.array(ys)[sl.borb]
    # block-level (E, r) options with counts
    items = [(E, r, n1) for (E, r), (n1, nx) in tab.items()]
    # sample quadruples: first (E1,r1,E2,r2 | E3,r3,E4,r4) by weight -- do by rejection over a weighted list of pair
    # splits: enumerate all ordered 4-tuples of (E, r) with sums matching (could be large; cap)
    quads = []
    byER = defaultdict(int)
    for (E, r, n1) in items: byER[(E, r)] += n1
    keys = list(byER)
    pair = defaultdict(list)
    for i, k1 in enumerate(keys):
        for k2 in keys:
            if k1[0] + k2[0] > T0: continue
            pair[(k1[0] + k2[0], k1[1] ** 2 + k2[1] ** 2)].append((k1, k2, byER[k1] * byER[k2]))
    for key, lst in pair.items():
        k2 = (T0 - key[0], R4 - key[1])
        if k2 not in pair: continue
        for (a1, a2, wa) in lst:
            for (b1, b2, wb) in pair[k2]:
                quads.append((a1, a2, b1, b2, wa * wb))
        if len(quads) > 2 * 10 ** 5: break
    Wq = np.array([float(q_[4]) for q_ in quads]); Wq /= Wq.sum()
    Gp = sl.Gp
    Up = [Gp.act(u) for u in sl.Sp] if len(Gp.mods) else [np.zeros(1, np.int64)]
    reps, _ = lag_reps(Gp, Up) if len(Gp.mods) else ([0], None)
    V = np.zeros((N, len(reps)), np.int64)
    pick = rng.choice(len(quads), size=N, p=Wq)
    for i, qi in enumerate(pick):
        Y = np.array([sample_block(E, r) for (E, r) in quads[qi][:4]])
        V[i] = Gp.paf(Y).sum(0)[reps]
    gb = gauss_bits(V, sl.T1[reps])
    out.update(lag_classes=len(reps), bits=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'),
               mahal=gb.get('mahal'), stage1_log2E=round(math.log2(Ny) - gb['bits'], 2) if gb.get('in_lattice') else -math.inf)
    return out


def subgroup_compress_index(G, Lmask):
    """for a subgroup L (mask) of G: the coset index of every element and the quotient size"""
    cos = -np.ones(G.n, np.int64); n = 0
    Lel = np.nonzero(Lmask)[0]
    for g in range(G.n):
        if cos[g] >= 0: continue
        for l in Lel: cos[G.add(g, l)] = n
        n += 1
    return cos, n


def stage2(sl, Y, N=20000, seed=1, levels=None, chunk=1000):
    """levels: dict name -> mask of a subgroup L of G' (characters trivial on L)"""
    rng = np.random.default_rng(seed); Y = np.asarray(Y)
    G = sl.G
    pats = []; H2 = 0.0
    for k in range(4):
        pk = []
        for j, bb in enumerate(sl.breps):
            P, sizes = fibre_patterns(sl, j, int(Y[k, bb]))
            if len(P) == 0: return dict(H2=None, note="a fibre value is not realisable")
            pk.append(P); H2 += math.log2(len(P))
        pats.append(pk)
    reps, cls = lag_reps(G, sl.perms)
    V = np.zeros((N, G.n), np.int64); done = 0
    while done < N:
        c = min(chunk, N - done)
        acc = np.zeros((c, G.n), np.int64)
        for k in range(4):
            lab = np.zeros((c, sl.r), np.int8)
            for j in range(len(sl.breps)):
                P = pats[k][j]; pk_ = P[rng.integers(len(P), size=c)]
                for i, (o, s_) in enumerate(sl.parts[j]): lab[:, o] = pk_[:, i]
            acc += G.paf(lab[:, sl.orb])
        V[done:done + c] = acc; done += c
    out = dict(H2=round(H2, 2))
    gb = gauss_bits(V[:, reps], sl.T[reps])
    out.update(bits=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'),
               stage2_log2E=round(H2 - gb['bits'], 2) if gb.get('in_lattice') else -math.inf)
    if levels:
        out['levels'] = {}
        for name, Lmask_p in levels.items():
            # lift the G' subgroup to G (Z_q coordinate 0)
            Lmask = Lmask_p[sl.bidx] & (sl.zco == 0)
            cos, nC = subgroup_compress_index(G, Lmask)
            # compressed PAF: sum over l in L of PAF(g + l) -> function of the coset of g; one value per coset rep
            VL = np.zeros((N, nC), np.int64); TL = np.zeros(nC, np.int64)
            np.add.at(VL.T, cos, V.T); np.add.at(TL, cos, sl.T)
            # representatives: one per coset (lattice handles repeats/symmetry)
            gl = gauss_bits(VL, TL)
            out['levels'][name] = dict(bits=gl.get('bits'), rank=gl.get('rank'), in_lattice=gl.get('in_lattice'))
    return out


def stage1_sampler(sl, seed=1):
    """returns (draw(), info): draw() gives a uniform y-quadruple (4 x |G'|) among U-invariant quadruples with the
    stage-1 lag-0 energy and sum_k r_k^2 = R4 (exact weights)"""
    rng = np.random.default_rng(seed)
    T0 = int(sl.T1[0]); R4 = int(sl.T1.sum())
    w = sl.bsz; L = sl.letters; nb = len(sl.breps)
    back = [None] * (nb + 1); back[nb] = {(0, 0): 1}
    for j in range(nb - 1, -1, -1):
        nd = defaultdict(int)
        for (E, r), c in back[j + 1].items():
            for y in L[j]:
                E2 = E + int(w[j]) * y * y
                if E2 <= T0: nd[(E2, r + int(w[j]) * y)] += c
        back[j] = dict(nd)
    byER = back[0]
    keys = list(byER)
    pair = defaultdict(list)
    for k1 in keys:
        for k2 in keys:
            if k1[0] + k2[0] > T0: continue
            pair[(k1[0] + k2[0], k1[1] ** 2 + k2[1] ** 2)].append((k1, k2, byER[k1] * byER[k2]))
    quads = []
    for key, lst in pair.items():
        k2 = (T0 - key[0], R4 - key[1])
        if k2 not in pair: continue
        for (a1, a2, wa) in lst:
            for (b1, b2, wb) in pair[k2]:
                quads.append((a1, a2, b1, b2, wa * wb))
        if len(quads) > 2 * 10 ** 5: break
    if not quads: return None, dict(log2Ny=None)
    tot = sum(q_[4] for q_ in quads)
    Wq = np.array([float(q_[4]) for q_ in quads]); Wq /= Wq.sum()
    def sample_block(E, r):
        ys = []
        for j in range(nb):
            opts = []; wts = []
            for y in L[j]:
                e2 = E - int(w[j]) * y * y; r2 = r - int(w[j]) * y
                c = back[j + 1].get((e2, r2), 0)
                if c: opts.append(y); wts.append(c)
            p = np.array(wts, float); p /= p.sum()
            y = opts[rng.choice(len(opts), p=p)]; ys.append(y)
            E -= int(w[j]) * y * y; r -= int(w[j]) * y
        return np.array(ys)[sl.borb]
    def draw():
        qi = rng.choice(len(quads), p=Wq)
        return np.array([sample_block(E, r) for (E, r) in quads[qi][:4]])
    return draw, dict(log2Ny=round(math.log2(tot), 2))


def gauss_bits_at_mean(V):
    """-log2 of the lattice-Gaussian density at the sample mean (log-det part only): the stage-2 constraint of a
    'typical' target; uses the same lattice coordinates as gauss_bits"""
    V = np.asarray(V, np.int64)
    T = np.rint(V.mean(0)).astype(np.int64)
    # project T into V[0] + lattice: use gauss_bits on T' = nearest sample (in-lattice), and correct the Mahalanobis
    d = ((V - V.mean(0)) ** 2).sum(1)
    Tn = V[int(np.argmin(d))]
    gb = gauss_bits(V, Tn)
    if not gb.get('in_lattice') or gb.get('bits') is None or not np.isfinite(gb['bits']): return gb
    m = gb.get('mahal') or 0.0
    gb['bits_at_mean'] = round(gb['bits'] - 0.5 * m * m / math.log(2), 2)
    return gb


def stage2_typical(sl, Y, N=6000, seed=1, chunk=1000):
    """stage-2 Gaussian for a y that need not solve stage 1: H2 exact, constraint = log-det bits at the mean"""
    rng = np.random.default_rng(seed); Y = np.asarray(Y); G = sl.G
    pats = []; H2 = 0.0
    for k in range(4):
        pk = []
        for j, bb in enumerate(sl.breps):
            P, sizes = fibre_patterns(sl, j, int(Y[k, bb]))
            if len(P) == 0: return dict(H2=None)
            pk.append(P); H2 += math.log2(len(P))
        pats.append(pk)
    reps, cls = lag_reps(G, sl.perms)
    V = np.zeros((N, len(reps)), np.int64); done = 0
    while done < N:
        c = min(chunk, N - done)
        acc = np.zeros((c, G.n), np.int64)
        for k in range(4):
            lab = np.zeros((c, sl.r), np.int8)
            for j in range(len(sl.breps)):
                P = pats[k][j]; pk_ = P[rng.integers(len(P), size=c)]
                for i, (o, s_) in enumerate(sl.parts[j]): lab[:, o] = pk_[:, i]
            acc += G.paf(lab[:, sl.orb])
        V[done:done + c] = acc[:, reps]; done += c
    gb = gauss_bits_at_mean(V)
    return dict(H2=round(H2, 2), C2=gb.get('bits_at_mean'), rank=gb.get('rank'),
                stage2_log2E=round(H2 - gb['bits_at_mean'], 2) if gb.get('bits_at_mean') is not None else None)


def _pattern_sampler(sizes, y):
    """uniform sampler of sign vectors e (len(sizes)) with sum e_i sizes_i = y, by DP counts (any number of parts);
    returns (count, draw(rng, c) -> c x n int8)"""
    n = len(sizes)
    suf = [defaultdict(int) for _ in range(n + 1)]; suf[n][0] = 1          # suf[i][t]: patterns of parts i.. with sum t
    for i in range(n - 1, -1, -1):
        for t, m in suf[i + 1].items():
            suf[i][t + sizes[i]] += m; suf[i][t - sizes[i]] += m
    total = suf[0].get(y, 0)
    def draw(rng, c):
        out = np.zeros((c, n), np.int8)
        rem = np.full(c, y, np.int64)
        for i in range(n):
            ap = np.array([suf[i + 1].get(int(r - sizes[i]), 0) for r in rem], np.float64)
            am = np.array([suf[i + 1].get(int(r + sizes[i]), 0) for r in rem], np.float64)
            plus = rng.random(c) * (ap + am) < ap
            out[:, i] = np.where(plus, 1, -1); rem = rem - out[:, i] * sizes[i]
        return out
    return total, draw


def stage2_sampled(sl, Y, N=12000, seed=1, chunk=1000):
    """as stage2 (no levels) for fibres with any number of parts: uniform fibre sampling by DP counts"""
    rng = np.random.default_rng(seed); Y = np.asarray(Y); G = sl.G
    samp = []; H2 = 0.0
    for k in range(4):
        sk = []
        for j, bb in enumerate(sl.breps):
            sizes = [s_ for _, s_ in sl.parts[j]]
            cnt, draw = _pattern_sampler(sizes, int(Y[k, bb]))
            if cnt == 0: return dict(H2=None, note="a fibre value is not realisable")
            sk.append(draw); H2 += math.log2(cnt)
        samp.append(sk)
    reps, cls = lag_reps(G, sl.perms)
    V = np.zeros((N, len(reps)), np.int64); done = 0
    while done < N:
        c = min(chunk, N - done)
        acc = np.zeros((c, G.n), np.int64)
        for k in range(4):
            lab = np.zeros((c, sl.r), np.int8)
            for j in range(len(sl.breps)):
                pk_ = samp[k][j](rng, c)
                for i, (o, s_) in enumerate(sl.parts[j]): lab[:, o] = pk_[:, i]
            acc += G.paf(lab[:, sl.orb])
        V[done:done + c] = acc[:, reps]; done += c
    gb = gauss_bits(V, sl.T[reps])
    return dict(H2=round(H2, 2), bits=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'),
                mahal=gb.get('mahal'), stage2_log2E=round(H2 - gb['bits'], 2) if gb.get('in_lattice') else -math.inf)
