r"""Library for the v = 599 lead slice (b = 5, c = 6, q = 11, s = 2, H' = Z_9, U = -1 on Z_9) [P, M].

Group and slice (report/b599_answer.md section 1):
  G = Z_2 x Z_3 x Z_9 x Z_11 (mods (2, 3, 9, 11)), H = Z_9 x Z_11 (index 6), Q = G/H = Z_2 x Z_3.
  S = <squares on Z_11> x <-1 on Z_9>: 90 orbits per block; stage 1 lives on G' = Z_2 x Z_3 x Z_9.
Block vectors (stage 1): v in {+-1, +-9, +-11}^30, index o = 5 q + j with q = 3 e + t in Q (e in Z_2, t in Z_3) and
  j in {0, 1, 2, 3, 4} the representative of {j, -j} in Z_9;  weight w_o = 1 (j = 0) or 2.
  Expansion to G' (ravel index g = 27 e + 9 t + j9):  y[g] = v[5 q + min(j9, 9 - j9)].
Stage-1 target T1 on G': 2176 at 0, -220 on H' \ 0 = {(0, 0, j), j != 0}, +44 off H'.
Lag classes of G' under negation and j -> -j: (e, t, j) with e in {0, 1}, t in {0, 1}, j in {0..4}: 20 classes;
  class index L = 10 e + 5 t + j.  sig(v) = PAF_{G'}(y) at the 20 class representatives (exact integers).
S/D split on the Z_2 axis [P] (report/alpoge_design_answer.md 1.2):
  sigma = (y(0, .) + y(1, .))/2, delta = (y(0, .) - y(1, .))/2 on G'' = Z_3 x Z_9 (15 cells per block, cell = 5 t + j),
  sum_k PAF(sigma_k) = T_S = (555 at 0, -44 on (0, j != 0), 22 on (t != 0, .)),
  sum_k PAF(delta_k) = T_D = (533 at 0, -66 on (0, j != 0), 0 on (t != 0, .)).
"""
import os, sys, math, itertools, json
import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LET = np.array([-11, -9, -1, 1, 9, 11], np.int64)
W30 = np.array([1, 2, 2, 2, 2] * 6, np.int64)
NG = 54

# ---------------------------------------------------------------- group bookkeeping
def _g(e, t, j9):
    return 27 * e + 9 * (t % 3) + (j9 % 9)

EXP = np.zeros(NG, np.int64)            # g -> orbit index o
for e in range(2):
    for t in range(3):
        for j9 in range(9):
            EXP[_g(e, t, j9)] = 5 * (3 * e + t) + min(j9, 9 - j9)

LAGREP = np.array([_g(e, t, j) for e in range(2) for t in range(2) for j in range(5)], np.int64)
SHIFT = np.zeros((20, NG), np.int64)     # SHIFT[L, g] = g + d_L
for L in range(20):
    d = LAGREP[L]; de, dt, dj = d // 27, (d // 9) % 3, d % 9
    for g in range(NG):
        e, t, j = g // 27, (g // 9) % 3, g % 9
        SHIFT[L, g] = _g((e + de) % 2, t + dt, j + dj)
# class of every lag g (for full PAF -> class vector)
LAGCLS = np.zeros(NG, np.int64)
for g in range(NG):
    e, t, j = g // 27, (g // 9) % 3, g % 9
    LAGCLS[g] = 10 * e + 5 * (0 if t == 0 else 1) + min(j, 9 - j)
CLSIZE = np.bincount(LAGCLS, minlength=20)

T1 = np.full(NG, 44, np.int64)
for j9 in range(1, 9): T1[_g(0, 0, j9)] = -220
T1[0] = 2176
T1SIG = T1[LAGREP].copy()               # target signature (20 classes)


def expand(v):
    """block vector(s) (..., 30) -> y on G' (..., 54)"""
    return np.asarray(v)[..., EXP]


def paf_full(y):
    """exact PAF over G' of y (..., 54) by direct sums"""
    y = np.asarray(y, np.int64)
    out = np.zeros(y.shape, np.int64)
    for g in range(NG):
        d_e, d_t, d_j = g // 27, (g // 9) % 3, g % 9
        idx = np.array([_g((h // 27 + d_e) % 2, (h // 9) % 3 + d_t, h % 9 + d_j) for h in range(NG)])
        out[..., g] = (y * y[..., idx]).sum(-1)
    return out


@njit(cache=True)
def sig_of(V, SH):
    """V: (n, 54) int -> (n, 20) PAF at the class representatives"""
    n = V.shape[0]; out = np.zeros((n, 20), np.int64)
    for i in range(n):
        for L in range(20):
            acc = 0
            for g in range(54): acc += V[i, g] * V[i, SH[L, g]]
            out[i, L] = acc
    return out


def sig(v):
    v = np.atleast_2d(np.asarray(v, np.int64))
    return sig_of(expand(v), SHIFT)


def qcomp(v):
    """Q-compression z(q) = sum over the coset q + H' (orbit weights)"""
    v = np.atleast_2d(np.asarray(v, np.int64))
    return (v * W30).reshape(v.shape[0], 6, 5).sum(2)


def energy(v):
    v = np.atleast_2d(np.asarray(v, np.int64))
    return (v * v * W30).sum(1)


def check_stage1(Y):
    """exact check of a stage-1 quadruple (4 x 30 block vectors): letters and sum_k PAF = T1 (full PAF)"""
    Y = np.asarray(Y, np.int64)
    if Y.shape != (4, 30) or not np.isin(Y, LET).all(): return False
    return bool((paf_full(expand(Y)).sum(0) == T1).all())


# ---------------------------------------------------------------- S/D split
# cell c = 5 t + j of G'' = Z_3 x Z_9 (/ j -> -j); orbit o0 = 5 t + j (e = 0), o1 = 15 + 5 t + j (e = 1)
W15 = np.array([1, 2, 2, 2, 2] * 3, np.int64)


def sd_of(v):
    v = np.atleast_2d(np.asarray(v, np.int64))
    a, b = v[:, :15], v[:, 15:]
    return (a + b) // 2, (a - b) // 2


def _gpp(t, j9):
    return 9 * (t % 3) + (j9 % 9)

EXP15 = np.zeros(27, np.int64)
for t in range(3):
    for j9 in range(9): EXP15[_gpp(t, j9)] = 5 * t + min(j9, 9 - j9)
LAGREP15 = np.array([_gpp(t, j) for t in range(2) for j in range(5)], np.int64)     # 10 classes, L = 5 t + j
SHIFT15 = np.zeros((10, 27), np.int64)
for L in range(10):
    d = LAGREP15[L]; dt, dj = d // 9, d % 9
    for g in range(27):
        t, j = g // 9, g % 9
        SHIFT15[L, g] = _gpp(t + dt, j + dj)
TS = np.array([555] + [-44] * 4 + [22] * 5, np.int64)
TD = np.array([533] + [-66] * 4 + [0] * 5, np.int64)


@njit(cache=True)
def sig15_of(V, SH):
    n = V.shape[0]; out = np.zeros((n, 10), np.int64)
    for i in range(n):
        for L in range(10):
            acc = 0
            for g in range(27): acc += V[i, g] * V[i, SH[L, g]]
            out[i, L] = acc
    return out


def sig15(c):
    c = np.atleast_2d(np.asarray(c, np.int64))
    return sig15_of(c[:, EXP15], SHIFT15)


# ---------------------------------------------------------------- borders (347 template class)
def template_class():
    """the 347 template's per-block compression shapes: z_k = c_k 1 + 2 1_{S_k} (template order), and the kernel"""
    S = [[], [0, 3], [1, 2], [0, 4, 5]]
    return S


def bordered_z_quads(smax=40):
    """all z-quadruples z_k = comp_k + s_k 1 (template block order, s_k even) with sum_k row_k^2 = 2396, translated
    by every t_k in Q, that pass the exact border fit (theory/ab_borderfit.fits); returns list of (z (4x6), fit)"""
    import ab_border_transfer2 as BT, ab_borderfit as BF
    from ab_lib import AG
    T = BT.load(1388); comp = T['comp']; Qg = AG((2, 3))
    a = [int(comp[k].sum()) for k in range(4)]
    sols = []
    rng = range(-smax, smax + 1, 2)
    for s in itertools.product(rng, repeat=4):
        r = [a[k] + 6 * s[k] for k in range(4)]
        if sum(x * x for x in r) == 2396: sols.append(s)
    out = []
    for s in sols:
        base = [comp[k] + s[k] for k in range(4)]
        for tr in itertools.product(range(6), repeat=4):
            z = np.array([base[k][Qg.idx(Qg.co + Qg.co[tr[k]])] for k in range(4)])
            f = BF.fits(z, 1388)
            if f: out.append((z, s, tr, f[0]))
    return out


if __name__ == "__main__":
    # self-test against the generic Slice of theory/ab_sweep.py
    import ab_sweep as Wp, ab_s1sd as SD
    cfg = dict(v=599, c=6, b=5, h=99, m=594, Q=[2, 3], H1=[9], q=11, s=2, f=5, mods=[2, 3, 9, 11], nQ=2, U=[[2, 2, 8]])
    sl = Wp.build(SD.cfg_from_json(json.dumps(cfg)))
    assert (sl.T1 == T1).all()
    # orbit partition of G' agrees
    for o in range(30):
        els = np.nonzero(EXP == o)[0]
        assert len(set(sl.borb[els].tolist())) == 1 and len(els) == W30[o]
    rng = np.random.default_rng(1)
    V = LET[rng.integers(6, size=(50, 30))]
    P = sl.Gp.paf(expand(V))
    assert (P[:, LAGREP] == sig(V)).all() and (paf_full(expand(V)) == P).all()
    assert (P == P[:, np.arange(54)][:, :]).all()
    # class structure: PAF is constant on each class
    for L in range(20): assert (P[:, LAGCLS == L] == P[:, [LAGREP[L]]]).all()
    # S/D identities
    sg, dl = sd_of(V)
    PS = sig15(sg); PD = sig15(dl)
    for L in range(10):
        d = LAGREP15[L]
        assert (4 * PS[:, L] == P[:, d] + P[:, 27 + d]).all() and (4 * PD[:, L] == P[:, d] - P[:, 27 + d]).all()
    TSc = np.array([(T1[d] + T1[27 + d]) // 4 for d in LAGREP15]); TDc = np.array([(T1[d] - T1[27 + d]) // 4 for d in LAGREP15])
    assert (TSc == TS).all() and (TDc == TD).all()
    print(json.dumps(dict(selftest="ok", T1SIG=T1SIG.tolist(), CLSIZE=CLSIZE.tolist(), TS=TS.tolist(), TD=TD.tolist())))
