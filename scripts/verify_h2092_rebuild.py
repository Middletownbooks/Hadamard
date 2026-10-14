# Independent check of H(2092) written outside the 523 session: (1) H H^T = N I for data/H2092_packed.npz; (2) the core
# data/core_v523_b1_q523_p523f_w0.txt rebuilt into a bordered Goethals-Seidel matrix by a brute-force search over the
# border S, T in H(4) (corner solved linearly), with no code shared with theory/ or scripts/verify_h2092.py.
import numpy as np, hashlib, itertools as it
R = '/home/user/Hadamard/'
# (1) packed matrix
z = np.load(R + 'data/H2092_packed.npz'); N = int(z['N'])
H = np.unpackbits(z['bits'], axis=1)[:, :N].astype(np.int64) * 2 - 1
G = H @ H.T
print('packed: N', N, 'HH^T == N I:', bool((G == N * np.eye(N, dtype=np.int64)).all()))
print('packed sha256', hashlib.sha256(z['bits'].tobytes()).hexdigest())
# (2) core -> my own bordered GS
L = open(R + 'data/core_v523_b1_q523_p523f_w0.txt').read().split('\n')
X = [np.array([1 if c == '+' else -1 for c in l.split('=')[1].strip()], np.int64) for l in L[1:] if '=' in l]
m = len(X[0]); assert m == 522 and len(X) == 4
P = sum(np.array([int(np.dot(x, np.roll(x, -h))) for h in range(m)]) for x in X)
print('core: PAF sum at 0', P[0], 'off-zero values', sorted(set(P[1:].tolist())), 'rows', [int(x.sum()) for x in X])
# multiplier invariance x181
print('x181 invariant:', all((x[(181 * np.arange(m)) % m] == x).all() for x in X))
def circ(x): return np.array([np.roll(x, i) for i in range(len(x))])     # row i = x shifted right by i
A, B, C, D = [circ(x) for x in X]
Rm = np.fliplr(np.eye(m, dtype=np.int64))
GS = np.block([[A, B @ Rm, C @ Rm, D @ Rm],
               [-B @ Rm, A, D.T @ Rm, -C.T @ Rm],
               [-C @ Rm, -D.T @ Rm, A, B.T @ Rm],
               [-D @ Rm, C.T @ Rm, -B.T @ Rm, A]])
gg = GS @ GS.T
print('GS core: diag', set(np.diag(gg).tolist()), 'GS GS^T == 4m I + ... checking border completion')
# border: find S,T in H(4) and corner c making the (4m+4) matrix Hadamard. Rows of GS have block row sums:
one = np.ones(m, np.int64)
H4s = []
base = np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])
for p in it.permutations(range(4)):
    for sg in it.product((1,-1), repeat=4):
        H4s.append(np.diag(sg) @ base[list(p)])
found = None
for S in H4s[:]:
    # border columns: block i row-block gets column vector S[i, j] * 1 for border column j
    Bc = np.kron(S, one[:, None])            # (4m x 4)
    if not ((GS @ GS.T + Bc @ Bc.T) == (4 * m + 4) * np.eye(4 * m, dtype=np.int64)).all(): continue
    for T in H4s:
        Br = np.kron(T, one[None, :])       # (4 x 4m)
        # rows of border must be orthogonal to GS rows: Br @ GS^T + c @ Bc^T = 0  -> solve c
        M1 = Br @ GS.T                       # 4 x 4m
        # c Bc^T = -M1  ; Bc^T Bc = m * S^T S = 4m I  -> c = -M1 Bc / (4m)
        c = -(M1 @ Bc) / (4 * m)
        if not np.allclose(c, np.round(c)) or not set(np.round(c).ravel().tolist()) <= {1, -1}: continue
        c = np.round(c).astype(np.int64)
        Hf = np.block([[GS, Bc], [Br, c]])
        if ((Hf @ Hf.T) == (4 * m + 4) * np.eye(4 * m + 4, dtype=np.int64)).all():
            found = (S, T, c); break
    if found: break
print('rebuilt bordered GS from core: Hadamard', found is not None, 'order', 4 * m + 4)
