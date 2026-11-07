r"""H(4892) slice, step 2: exact stage 2 (the lift from Z47 x Z2 to G = Z47 x Z2 x Z13) for a stage-1 solution [M].

G = Z47 x Z2 x Z13, element (a, b, z) at index 13 (2a + b) + z.  S = U x M (U = squares mod 47 on a, M = {1,3,9}
on z).  S-orbits: (O, Z) = O x {0} and (O, C_j) = O x C_j, C_j = 2^j M, for the six U-orbits O of Z47 x Z2 (in the
order of theory/h4892_stage1.ORB) -- orbit id 5 o + (0 for Z, 1 + j for C_j).
Given y (four blocks, 6 orbit values each), alpha(O) = the +-1 with alpha = y mod 3 and sigma(O) = (y - alpha)/3;
the unknowns are the class signs c_j(O) in {+-1}, sum_j c_j(O) = sigma(O).  A block list is the product over O of
the C(4, (sigma + 4)/2) sign patterns.  PAF_G of a block is a quadratic form in its 30 orbit values: PAF(h) =
sum_{p,q} u_p u_q N_pq(h), N_pq(h) = #{g in orbit p : g + h in orbit q}, counted by brute force on G; it is constant
on S-orbits of lags and PAF(-h) = PAF(h), so it is recorded on one lag per {orbit, -orbit} class (15 nonzero
classes).  Target: -4 on every nonzero lag (b = 1 family, sum_k PAF_k = 4*1222 at 0).

Symmetry used [P]: the multiplier z -> 2z of Z13 fixes 0, maps C_j to C_{j+1}, preserves S-orbits and the target,
and fixes y; so if (x_0..x_3) is a stage-2 solution, so is its image, and some power of it brings any one chosen
block's pattern to the least of its <2>-orbit.  The search restricts ONE block (the 'pivot') to orbit minima.

Join: blocks split {a, b} | {c, d}; each list is bucketed by its S-level vector s = (PAF of the compression along
the subgroup L) -- linear in the PAF vector; a left pair (a, b) is formed only in bucket pairs whose s-sum has a
right partner with s-sum T_S - s.  Inside a matched bucket pair, a 64-bit linear hash of the 15-vector (exact mod
2^64) is sorted on the right, the left streamed against it, and every hash hit re-checked on the exact vectors.
usage (library); driver: theory/h4892_stage2_run.py"""
import sys, os, itertools, math
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h4892_stage1 as S1

P47, Q, n1 = 47, 13, 94
NG = n1 * Q
M = [1, 3, 9]
C = [sorted((pow(2, j, Q) * m) % Q for m in M) for j in range(4)]
assert sorted(sum(C, [])) == list(range(1, Q))
zcls = [0] * Q
for j, c in enumerate(C):
    for z in c: zcls[z] = 1 + j
gid = lambda g1, z: Q * g1 + z
orbit = np.zeros(NG, np.int64)
for o, O in enumerate(S1.ORB):
    for g1 in O:
        for z in range(Q): orbit[gid(g1, z)] = 5 * o + zcls[z]
NO = 30
SIZE = np.bincount(orbit, minlength=NO)
def addG(g, h):
    g1, z = divmod(g, Q); h1, t = divmod(h, Q)
    return gid(S1.add(g1, h1), (z + t) % Q)
def neg(g):
    g1, z = divmod(g, Q); a, b = divmod(g1, 2)
    return gid(2 * ((-a) % P47) + b, (-z) % Q)
# lag classes: orbit reps, merged with the orbit of -rep
rep = [int(np.nonzero(orbit == p)[0][0]) for p in range(NO)]
cls = {}; LAGREP = []
for p in range(NO):
    q = int(orbit[neg(rep[p])])
    if p in cls: continue
    cls[p] = cls[q] = len(LAGREP); LAGREP.append(rep[p])
LAGREP_NZ = [h for h in LAGREP if h != 0]                     # 15 nonzero classes
NL = len(LAGREP_NZ)
Ncnt = np.zeros((NL, NO, NO), np.int64)
for l, h in enumerate(LAGREP_NZ):
    for g in range(NG): Ncnt[l, orbit[g], orbit[addG(g, h)]] += 1
# class weights (number of lags in each class), for the linear S-level maps
LAGW = np.zeros(NL, np.int64)
for p in range(NO):
    h = rep[p]
    if h == 0: continue
    LAGW[LAGREP_NZ.index(LAGREP[cls[p]])] += SIZE[p]
assert LAGW.sum() == NG - 1


def decode(v):
    al = 1 if (v - 1) % 3 == 0 else -1
    return al, (v - al) // 3


PAT = {s: [p for p in itertools.product((1, -1), repeat=4) if sum(p) == s] for s in (-4, -2, 0, 2, 4)}


def block_list(v6):
    """all orbit-value vectors (n x 30) of the block with U-orbit values v6"""
    per = []
    for o in range(6):
        al, sg = decode(int(v6[o])); per.append([(al,) + p for p in PAT[sg]])
    L = [sum(c, ()) for c in itertools.product(*per)]
    return np.array(L, np.int8)


def paf_vec(Uv):
    U = Uv.astype(np.int64)
    return np.stack([((U @ Ncnt[l]) * U).sum(1) for l in range(NL)], 1)


def shift_pattern(Uv):
    """image under z -> 2z: the value on C_j moves to C_{j+1}"""
    X = Uv.reshape(len(Uv), 6, 5)
    Y = X.copy(); Y[:, :, 2:] = X[:, :, 1:4]; Y[:, :, 1] = X[:, :, 4]
    return Y.reshape(len(Uv), 30)


def orbit_min_mask(Uv):
    """True for patterns that are the lexicographic least of their <z -> 2z> orbit"""
    def key(A):
        A = A.astype(np.int64) + 1                            # 0 / 2
        return A @ (3 ** np.arange(30, dtype=np.int64)[::-1])
    k0 = key(Uv); best = k0.copy(); X = Uv
    for _ in range(3):
        X = shift_pattern(X); best = np.minimum(best, key(X))
    return k0 == best


def expand(Uv_row):
    return np.array([int(Uv_row[orbit[g]]) for g in range(NG)], np.int64)
