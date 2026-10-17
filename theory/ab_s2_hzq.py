r"""Exact H' x Z_q-level screens of stage 2 for a staged c = 4 slice (G = Q x H' x Z_q, Q = Z_2^2) [P, M].

For a character psi of Q put xbar_k^psi(h, t) = sum_q psi(q) x_k(q, h, t), a (U x K)-invariant function on
H' x Z_q.  Fourier-transforming sum_k PAF_G(x_k) = T over Q at psi gives, for every core,                     [P]
        sum_k PAF_{H' x Z_q}(xbar_k^psi) = T^psi,    T^psi(h, t) = sum_q psi(q) T(q, h, t).
(Compressing further along H' gives the Z_q-level condition of theory/ab_s2_zq_screen.py, so this one is stronger.)
For a stage-1 solution y: at an H'-orbit o and a K-orbit r of Z_q, xbar_k^psi takes the value
v_k[o, r] = sum_q psi(q) e_{(q,o)}(r), where e_{(q,o)} is the sign pattern chosen in the fibre of the letter y_k at
the G'-orbit (q, o) (its value on the part of the Z_q-fibre over r).  Per block: the distinct vectors v[o, .] per o
(exact multiplicities over the 4 Q-cosets' patterns), then the states (v[o, .])_o with block lag-0 energy within
budget (energies add over o), their PAF signatures on H' x Z_q at the lag classes of U x K and -1 (the last class is
implied by the fixed psi-sum of y), and the exact number of lifts by a 2+2 convolution (products of exact counts
summed in float64: zero/non-zero is exact).
usage (library): screen_y(S, y, psi) with S = setup(cfg); python3 theory/ab_s2_hzq.py CFG.json YBIN [--n N]"""
import sys, os, json, time, math, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s1_chars as AC, ab_s2_s4 as S4
from ab_lib import AG

CHARS = AC.CHARS


def setup(cfg):
    sl = W.build(cfg); st = AC.structure(cfg, sl); no = st['no']; L1 = len(sl.breps)
    Z = S4.zq_setup(sl); kr = Z['kreps']; nk = len(kr); q = Z['q']
    ksz = np.bincount(Z['kor'], minlength=nk)                                    # |K-orbit| of each class r
    # per (G'-orbit j, letter): pattern values at the K-orbit classes, per point of the G'-orbit
    vals = {}
    for j in range(L1):
        os_ = [o for o, _ in sl.parts[j]]; wj = int(sl.bsz[j])
        for y in sl.letters[j]:
            C = AGs.fibre_patterns(sl, j, int(y))[0].astype(np.int64) @ Z['Mo'][os_][:, kr]
            assert (C % wj == 0).all(); vals[(j, int(y))] = C // wj
    # H'' = H' x Z_q with U x K and -1: orbits (o, r), lag classes, structure matrices
    Hm = tuple(sl.G.mods[2:]); Hg = AG(Hm); n2 = Hg.n
    co = Hg.co; hpart = co[:, :-1]; tpart = co[:, -1]
    # H'-orbit of each h: via the stage-1 structure (orbit of (q=0, h) in G')
    Gp = sl.Gp; base = {}
    for x in range(Gp.n):
        if Gp.co[x, 0] == 0 and Gp.co[x, 1] == 0: base[tuple(Gp.co[x, 2:].tolist())] = x
    jo = {}                                                                     # G'-orbit index -> H'-orbit index
    for (qq, o), j in st['j_of'].items(): jo[j] = o
    horb = np.array([jo[int(sl.borb[base[tuple(h)]])] for h in hpart.tolist()]) if hpart.shape[1] else np.zeros(n2, np.int64)
    orb2 = horb * nk + Z['kor'][tpart]                                          # orbit index of (h, t)
    nst = no * nk
    gens = [u for u in sl.S]
    perms = [Hg.act(tuple(u[2:])) for u in gens]
    reps, _ = AGs.lag_reps(Hg, perms)
    A = np.zeros((len(reps), nst, nst))
    for li, c in enumerate(reps):
        xc = Hg.add(np.arange(n2), np.full(n2, c)); np.add.at(A[li], (orb2, orb2[xc]), 1.0)
    w2 = np.bincount(orb2, minlength=nst)                                       # orbit sizes on H''
    S = dict(sl=sl, st=st, no=no, L1=L1, nk=nk, ksz=ksz, vals=vals, A=A, reps=reps, w2=w2, nst=nst, Hg=Hg)
    set_target(S, sl.T)
    return S


def set_target(S, T):
    """T^psi at the lag-class representatives of H'' for every character psi (a test may plant its own T)"""
    sl, Hg = S['sl'], S['Hg']; Tg = np.asarray(T, np.int64); qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]
    gidx = Hg.idx(sl.G.co[:, 2:]); S['Tpsi'] = {}
    for p in range(4):
        tt = np.zeros(Hg.n, np.int64); np.add.at(tt, gidx, np.array(CHARS[p])[qg] * Tg); S['Tpsi'][p] = tt[S['reps']]


def block_dist(S, yk, p, budget):
    """distinct states (v[o, r]) of one block with lag-0 energy <= budget: (signatures, counts, min energy)"""
    st, no, nk = S['st'], S['no'], S['nk']; w2 = S['w2'].reshape(no, nk)
    per_o = []
    for o in range(no):
        opts = [S['vals'][(st['j_of'][(qq, o)], int(yk[st['j_of'][(qq, o)]]))] * CHARS[p][qq] for qq in range(4)]
        acc = {(0,) * nk: 1}
        for V in opts:
            nacc = defaultdict(int)
            for v0, c0 in acc.items():
                for v in map(tuple, V.tolist()): nacc[tuple(a + b for a, b in zip(v0, v))] += c0
            acc = nacc
        vecs = np.array(list(acc.keys()), np.int64); cnt = np.array(list(acc.values()), np.float64)
        en = (w2[o][None, :] * vecs ** 2).sum(1)
        o_ = np.argsort(en); per_o.append((vecs[o_], cnt[o_], en[o_]))
    emin = sum(float(x[2][0]) for x in per_o)
    if emin > budget: return None, None, emin
    # combine orbits with the energy budget (the minimal energies of the remaining orbits reserved)
    rest = [sum(float(x[2][0]) for x in per_o[i + 1:]) for i in range(no)]
    cur_v, cur_c, cur_e = np.zeros((1, 0), np.int64), np.ones(1), np.zeros(1)
    for o in range(no):
        vecs, cnt, en = per_o[o]
        tot = cur_e[:, None] + en[None, :]
        ia, ib = np.nonzero(tot <= budget - rest[o])
        cur_v = np.concatenate([cur_v[ia], vecs[ib]], 1); cur_c = cur_c[ia] * cnt[ib]; cur_e = tot[ia, ib]
        if len(cur_c) == 0: return None, None, emin
    X = cur_v.astype(np.float64)
    P = np.rint(np.einsum('ni,lij,nj->nl', X, S['A'], X)).astype(np.int64)
    u, inv = np.unique(P[:, :-1], axis=0, return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), cur_c)
    return u, c, emin


def screen_y(S, y, p):
    """exact number of lifts of y meeting the H' x Z_q-level condition for character p (float64 sum of exact counts)"""
    T = S['Tpsi'][p]; T0 = int(T[0])
    mins = []
    for k in range(4):
        _, _, e = block_dist(S, y[k], p, -1); mins.append(e)
    if sum(mins) > T0: return 0.0, 'lag-0 energy'
    D = []
    for k in range(4):
        u, c, _ = block_dist(S, y[k], p, T0 - (sum(mins) - mins[k]))
        if u is None: return 0.0, 'lag-0 energy'
        D.append((u, c))
    Tt = T[:-1]; d = len(Tt); dt = np.dtype([('f%d' % i_, np.int64) for i_ in range(d)])
    sv = lambda R: np.ascontiguousarray(R, np.int64).view(dt).ravel()

    def conv(a, b, chunk=1 << 22):
        us, cs = [], []
        step = max(1, chunk // max(1, len(b[0])))
        for s0 in range(0, len(a[0]), step):
            R = (a[0][s0:s0 + step][:, None, :] + b[0][None, :, :]).reshape(-1, d)
            v = (a[1][s0:s0 + step][:, None] * b[1][None, :]).ravel()
            u, inv = np.unique(sv(R), return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), v)
            us.append(u); cs.append(c)
        u, inv = np.unique(np.concatenate(us), return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), np.concatenate(cs))
        return u, c
    A_ = conv(D[0], D[1]); B_ = conv(D[2], D[3])
    Bq = sv(Tt[None, :] - B_[0].view(np.int64).reshape(-1, d))
    pos = np.minimum(np.searchsorted(A_[0], Bq), len(A_[0]) - 1); hit = A_[0][pos] == Bq
    return float((A_[1][pos[hit]] * B_[1][hit]).sum()), 'exact convolution'


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read()); S = setup(cfg); L1 = S['L1']
    n = int(a[a.index('--n') + 1]) if '--n' in a else 20
    Y = np.fromfile(a[1], np.int8).reshape(-1, 4, L1)
    idx = np.random.default_rng(1).choice(len(Y), size=min(n, len(Y)), replace=False)
    t0 = time.time()
    for i in idx:
        y = Y[i].astype(np.int64); out = []
        for p in range(4):
            c, by = screen_y(S, y, p); out.append((None if c is None else (round(math.log2(c), 2) if c > 0 else 0), by))
        print(json.dumps(dict(i=int(i), psi=out, secs=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    main()
