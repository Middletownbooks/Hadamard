import sys, json, math, resource
resource.setrlimit(resource.RLIMIT_AS, (3 * 2 ** 30,) * 2)
sys.path.insert(0, '/home/user/Hadamard/theory')
import numpy as np
import ab_s1sd as SD
from s1z_lib import ZOrb
from q523_gauss import echelon, coords
cfg = SD.cfg_from_json('{"v":919,"c":8,"b":7,"h":114,"m":912,"Q":[2,2,2],"H1":[3,2],"q":19,"s":2,"f":9,"mods":[2,2,2,3,2,19],"nQ":3,"U":[[3,2,2]]}')
Z = ZOrb(cfg); C = Z.compressions(); cls, idx = Z.classes(C); qd = Z.class_quads(cls)
cache = {}
for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
qd = [q for q in qd if q['n0'] > 0]
print('multisets', len(qd), 'rest lags', len(Z.rest))
rng = np.random.default_rng(1)
for t, q in enumerate(qd[:4]):
    zq = [cls[i]['rep'] for i in q['m']]
    Y, n0 = Z.sample_quads(zq, 600, rng)
    Pf = Z.Gp.paf(Y.reshape(-1, Z.L1)[:, Z.borb]).reshape(600, 4, -1).sum(1)
    V = Pf[:, Z.reps]; T = Z.T1[Z.reps]
    D = (V - V[0])
    # per-coordinate gcd and residues
    g = [math.gcd(*[int(x) for x in D[:, c]]) for c in range(D.shape[1])]
    bad = [(c, int(Z.reps[c]), g[c], int((T[c] - V[0, c]) % g[c]) if g[c] else None) for c in range(D.shape[1]) if g[c] and (T[c] - V[0, c]) % g[c]]
    print(t, q['m'], 'coordinate gcds', g[:12], '... obstructed coords (index, lag, gcd, residue):', bad[:6])
    basis, piv = echelon(D[1:200].tolist())
    print('   lattice rank', len(basis), 'target in lattice:', coords(basis, piv, (T - V[0]).tolist()) is not None)
print('--- span test')
q = qd[0]; zq = [cls[i]['rep'] for i in q['m']]
Y, n0 = Z.sample_quads(zq, 800, rng)
Pf = Z.Gp.paf(Y.reshape(-1, Z.L1)[:, Z.borb]).reshape(800, 4, -1).sum(1)
V = Pf[:, Z.reps].astype(float); T = Z.T1[Z.reps].astype(float)
D = V - V[0]
U_, S_, Vt = np.linalg.svd(D, full_matrices=False)
r = int((S_ > 1e-6 * S_[0]).sum())
P = Vt[:r]
resid = (T - V[0]) - P.T @ (P @ (T - V[0]))
print('reps', len(Z.reps), 'real rank', r, 'target residual outside the span', float(np.abs(resid).max()))
# null-space direction(s)
N = Vt[r:]
for nv in N:
    nv = nv / np.abs(nv).max()
    print('invariant direction (lags, coeff):', [(int(Z.reps[i]), round(float(c), 3)) for i, c in enumerate(nv) if abs(c) > 1e-6][:12])
    print('   value on samples', float(nv @ V[0]), float(np.std(V @ nv)), 'on target', float(nv @ T))
print('--- congruence')
from fractions import Fraction
import sympy
Vi = Pf[:, Z.reps].astype(np.int64); Ti = Z.T1[Z.reps].astype(np.int64)
basis, piv = echelon((Vi[1:400] - Vi[0]).tolist())
B = sympy.Matrix(basis)            # rows = basis vectors (r x n)
tgt = sympy.Matrix([int(x) for x in (Ti - Vi[0])])
# solve x B = tgt over Q (least squares exact since tgt in span)
sol = (B * B.T).inv() * B * tgt
fr = [sympy.nsimplify(x) for x in sol]
print('target coords in the lattice basis (non-integral ones):', [(i, str(x)) for i, x in enumerate(fr) if x.q != 1][:10])
# find a dual vector detecting it: rows of the inverse Gram-related matrix
den = sympy.ilcm(*[x.q for x in fr])
print('denominator', den)
print('--- dual functional mod 3')
from sympy import GF
Bm = sympy.Matrix(basis) % 3
# null space of B^T over GF(3): lambda with B lambda = 0 mod 3
import itertools
n = Bm.cols
M3 = np.array(Bm.tolist(), dtype=np.int64) % 3
# gaussian elimination mod 3 to get null space of M3 (r x n)
A = M3.copy(); r, c = A.shape; pivcols = []; row = 0
for col in range(c):
    pr = next((i for i in range(row, r) if A[i, col] % 3), None)
    if pr is None: continue
    A[[row, pr]] = A[[pr, row]]
    inv = 1 if A[row, col] % 3 == 1 else 2
    A[row] = (A[row] * inv) % 3
    for i in range(r):
        if i != row and A[i, col] % 3: A[i] = (A[i] - A[i, col] * A[row]) % 3
    pivcols.append(col); row += 1
    if row == r: break
free = [j for j in range(c) if j not in pivcols]
tg = np.array([int(x) for x in (Ti - Vi[0])]) % 3
for fcol in free:
    lam = np.zeros(c, np.int64); lam[fcol] = 1
    for i, pc in enumerate(pivcols): lam[pc] = (-A[i, fcol]) % 3
    val = int((lam @ tg) % 3)
    if val:
        print('lambda (lag rep -> coeff mod 3):', [(int(Z.reps[j]), int(lam[j])) for j in range(c) if lam[j]])
        print('   lambda.(T - V0) mod 3 =', val, '; on samples:', set(int(x) for x in ((Vi - Vi[0]) @ lam) % 3))
        lagvec = {int(Z.reps[j]): int(lam[j]) for j in range(c) if lam[j]}
        print('   lag coordinates:', {g: Z.Gp.co[g].tolist() for g in lagvec})
        print('   T1 at those lags:', {g: int(Z.T1[g]) for g in lagvec}, 'sample value lam.V0 mod 3', int((lam @ Vi[0]) % 3), 'lam.T mod 3', int((lam @ Ti) % 3))
        break
