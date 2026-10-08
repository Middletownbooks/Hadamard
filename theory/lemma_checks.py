"""Numerical checks of the algebraic lemmas used in report/theory_arm_2000_3000.md.
Lemma A: Williamson array W(A,B,C,D) = [[A,B,C,D],[-B,A,-D,C],[-C,D,A,-B],[-D,-C,B,A]] satisfies
  W W^T = I4 (x) S + (cross terms) where the cross terms vanish iff
  D_AB = D_DC, D_AC = D_BD, D_AD = D_CB  (D_XY = X Y^T - Y X^T)   [rows]
Lemma B: Cooper-Wallis pattern X_r = sum_k T_k (x) N_rk with circulant T(3) = (I, S, S^2, 0) and commuting
  plug-ins: sum_r X_r X_r^T = 3I (x) S  +  (S^{-1} - S) (x) (G12 + G23 + G31),  G_kl = sum_r N_rk N_rl^T.
Lemma C: no amicable T-matrices of order 3; constant-coset quaternion rows of odd length are never orthogonal.
"""
import itertools
import numpy as np
rng = np.random.default_rng(1)

def W(A, B, C, D):
    return np.block([[A, B, C, D], [-B, A, -D, C], [-C, D, A, -B], [-D, -C, B, A]])

def Dl(X, Y): return X @ Y.T - Y @ X.T

n = 5
for trial in range(200):
    A, B, C, D = [rng.integers(-3, 4, size=(n, n)).astype(float) for _ in range(4)]
    M = W(A, B, C, D) @ W(A, B, C, D).T
    S = A @ A.T + B @ B.T + C @ C.T + D @ D.T
    # predicted off-diagonal blocks
    pred_zero = np.allclose(Dl(A, B), Dl(D, C)) and np.allclose(Dl(A, C), Dl(B, D)) and np.allclose(Dl(A, D), Dl(C, B))
    actual_diag = all(np.allclose(M[i*n:(i+1)*n, i*n:(i+1)*n], S) for i in range(4))
    off = M.copy()
    for i in range(4): off[i*n:(i+1)*n, i*n:(i+1)*n] = 0
    assert actual_diag
    # build matrices satisfying the Delta conditions: random A, B, C and solve D? just check the forward implication
    # on random data the conditions fail and off-diagonal must be nonzero
    assert (np.abs(off).max() > 0) == (not pred_zero) or pred_zero
print('Lemma A diag blocks = sum XX^T: ok; now check iff on constructed examples')
# construct examples with the Delta conditions: take A,B,C,D symmetric commuting (Williamson-like): conditions hold
def circ(v):
    n = len(v); return np.array([[v[(j - i) % n] for j in range(n)] for i in range(n)], dtype=float)
for trial in range(50):
    # commuting circulants with D_AB = D_DC etc: choose A,B arbitrary circulant, D = B^T?..  test identity directly:
    a, b, c, d = [rng.integers(-2, 3, size=7) for _ in range(4)]
    A, B, C, D = map(circ, (a, b, c, d))
    M = W(A, B, C, D) @ W(A, B, C, D).T
    blk = lambda i, j: M[i*7:(i+1)*7, j*7:(j+1)*7]
    # (1,2) block predicted: -(A B^T - B A^T) ... compare against combination of Deltas
    cands = {
        'DAB+DDC': Dl(A, B) + Dl(D, C), 'DAB-DDC': Dl(A, B) - Dl(D, C), 'DBA+DDC': Dl(B, A) + Dl(D, C), 'DBA-DDC': Dl(B, A) - Dl(D, C)}
    m12 = [k for k, v in cands.items() if np.allclose(blk(0, 1), v)]
    cands13 = {'DCA+DBD': Dl(C, A) + Dl(B, D), 'DAC+DBD': Dl(A, C) + Dl(B, D), 'DCA-DBD': Dl(C, A) - Dl(B, D), 'DAC-DBD': Dl(A, C) - Dl(B, D)}
    m13 = [k for k, v in cands13.items() if np.allclose(blk(0, 2), v)]
    cands14 = {'DDA+DCB': Dl(D, A) + Dl(C, B), 'DAD+DCB': Dl(A, D) + Dl(C, B), 'DDA-DCB': Dl(D, A) - Dl(C, B), 'DAD-DCB': Dl(A, D) - Dl(C, B)}
    m14 = [k for k, v in cands14.items() if np.allclose(blk(0, 3), v)]
    if trial == 0: print('row-1 off-diagonal blocks of W W^T equal:', m12, m13, m14)
    assert m12 and m13 and m14
# Lemma B
w = 7
S3 = circ([0, 1, 0]); I3 = np.eye(3)
T = [I3, S3, S3 @ S3, np.zeros((3, 3))]
for trial in range(50):
    a, b, c, d = [rng.integers(-2, 3, size=w) for _ in range(4)]
    A, B, C, D = map(circ, (a, b, c, d))
    N = [[A, B, C, D], [-B, A, -D, C], [-C, D, A, -B], [-D, -C, B, A]]
    X = [sum(np.kron(T[k], N[r][k]) for k in range(4)) for r in range(4)]
    lhs = sum(x @ x.T for x in X)
    Ssum = A @ A.T + B @ B.T + C @ C.T + D @ D.T
    G = lambda k, l: sum(N[r][k] @ N[r][l].T for r in range(4))
    rhs = np.kron(3 * I3, Ssum) + np.kron(np.linalg.inv(S3) - S3, G(0, 1) + G(1, 2) + G(2, 0))
    assert np.allclose(lhs, rhs), 'Lemma B mismatch'
print('Lemma B: sum X_r X_r^T = 3I(x)S + (S^-1 - S)(x)(G12+G23+G31): ok')
# Lemma C(i): all T-matrices of order 3 are signed shifts; none pairwise amicable
cnt = 0; amic = 0
vals = [-1, 0, 1]
for Ts in itertools.product(itertools.product(vals, repeat=3), repeat=4):
    Ts = np.array(Ts)
    if not (np.abs(Ts).sum(0) == 1).all(): continue          # disjoint supports covering Z3
    P = sum(np.array([sum(t[i] * t[(i + s) % 3] for i in range(3)) for s in range(3)]) for t in Ts)
    if not (P[0] == 3 and P[1] == 0 and P[2] == 0): continue
    cnt += 1
    Ms = [circ(t) for t in Ts]
    if all(np.allclose(Ms[i] @ Ms[j].T, Ms[j] @ Ms[i].T) for i in range(4) for j in range(i + 1, 4)): amic += 1
print(f'Lemma C(i): T-matrices of order 3: {cnt}, pairwise amicable: {amic}')
# Lemma C(ii): rows of odd length n with entries a+bi+cj+dk (a,b,c,d=+-1) of constant product abcd never orthogonal
def qmul(p, q):
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    return (a1*a2 - b1*b2 - c1*c2 - d1*d2, a1*b2 + b1*a2 + c1*d2 - d1*c2, a1*c2 - b1*d2 + c1*a2 + d1*b2, a1*d2 + b1*c2 - c1*b2 + d1*a2)
def conj(q): return (q[0], -q[1], -q[2], -q[3])
units = [q for q in itertools.product([1, -1], repeat=4) if q[0]*q[1]*q[2]*q[3] == 1]
ok = True
for n in (1, 3, 5):
    for row1 in itertools.product(units, repeat=n):
        if n == 5 and rng.random() > 0.02: continue
        for row2 in itertools.product(units, repeat=n):
            if n == 5 and rng.random() > 0.02: continue
            s = np.zeros(4)
            for x, y in zip(row1, row2): s += qmul(x, conj(y))
            if n > 1 and not s.any(): ok = False
print('Lemma C(ii) (no orthogonal constant-coset rows of odd length, n=3 exhaustive, n=5 sampled):', ok)
