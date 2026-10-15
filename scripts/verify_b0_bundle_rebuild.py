# Independent check of a b = 0 Goethals-Seidel bundle: own parser, own PAF, own GS array, BLAS Gram (exact: |entries| < 2^53).
import sys, gzip, numpy as np
seqf, matf = sys.argv[1], sys.argv[2]
X = np.array([[1 if c == '+' else -1 for c in l.split('=', 1)[1].strip()] for l in open(seqf) if l[:2] in ('A=', 'B=', 'C=', 'D=')], np.int64)
v = X.shape[1]
F = np.fft.fft(X, axis=1); P = np.rint(np.fft.ifft((np.abs(F) ** 2).sum(0)).real).astype(np.int64)
print('v', v, 'sum PAF at 0', P[0], 'off-zero values', sorted(set(P[1:].tolist())), 'row sums', X.sum(1).tolist())
circ = lambda x: np.array([np.roll(x, i) for i in range(v)], np.float64)
A, B, C, D = map(circ, X); R = np.fliplr(np.eye(v))
H = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R], [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
N = 4 * v
print('own GS build: entries', sorted(set(np.unique(H).tolist())), 'H H^T == N I:', np.array_equal(H @ H.T, N * np.eye(N)))
M = np.array([[1.0 if c == '+' else -1.0 for c in r.strip()] for r in gzip.open(matf, 'rt') if r.strip()])
print('published normalized matrix', M.shape, 'first row/col +1:', bool((M[0] == 1).all() and (M[:, 0] == 1).all()), 'Hadamard:', np.array_equal(M @ M.T, M.shape[0] * np.eye(M.shape[0])))
