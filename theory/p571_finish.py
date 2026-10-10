"""Bordered Goethals-Seidel completion and exact verification for a cyclic b = 1 core (written independently of
theory/paley_b1pipe.py; same construction as scripts/verify_h2292.py).

Input: four +-1 sequences x_k on Z_m with sum_k PAF_k = 4m at 0 and -4 elsewhere and row sums (+-2, 0, 0, 0).
Core = GS array of circulants A, B, C, D with the back-identity R; border rows S (x) 1_m^T and columns T^T (x) 1_m with
S, T in H(4); corner forced as -S Sigma^T T / 4 (Sigma = block row sums of the core array).  Checks H H^T = N I exactly
(float64 BLAS on integers of size <= N is exact) and every entry +-1.
usage: python3 theory/p571_finish.py core_file.txt   (lines A=..., B=..., C=..., D=... over cyclic Z_m)"""
import sys, os, json, itertools, hashlib
import numpy as np

def h4_list():
    H = [np.array(t).reshape(4, 4) for t in itertools.product((1, -1), repeat=16)]
    return [h for h in H if np.array_equal(h @ h.T, 4 * np.eye(4, dtype=np.int64))]

def bordered_gs(X):
    m = X.shape[1]; N = 4 * m + 4
    P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in X)
    if not (P[0] == 4 * m and (P[1:] == -4).all()): return None, dict(core_ok=False)
    circ = lambda x: np.array([np.roll(x, j) for j in range(m)])            # circ[r, c] = x[c - r]
    R = np.fliplr(np.eye(m, dtype=np.int64))
    A, B, C, D = [circ(x).astype(np.int64) for x in X]
    core = np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])
    Sig = np.array([[core[I * m:(I + 1) * m, J * m:(J + 1) * m].sum(1)[0] for J in range(4)] for I in range(4)])
    for S in h4_list():
        for T in h4_list():
            num = -(S @ Sig.T @ T)
            if (num % 4).any() or not (np.abs(num // 4) == 1).all(): continue
            H = np.zeros((N, N), np.int8)
            H[:4, :4] = num // 4
            H[:4, 4:] = np.repeat(S, m, axis=1); H[4:, :4] = np.repeat(T.T, m, axis=0); H[4:, 4:] = core
            Hf = H.astype(np.float64)
            if np.array_equal(Hf @ Hf.T, N * np.eye(N)) and set(np.unique(H).tolist()) == {-1, 1}:
                return H, dict(core_ok=True, hadamard=True, N=N, corner=H[:4, :4].tolist(),
                               sha256_packed=hashlib.sha256(np.packbits(H > 0, axis=1).tobytes()).hexdigest())
    return None, dict(core_ok=True, hadamard=False)

def read_core(path):
    L = open(path).read().split('\n')
    return np.array([[1 if c == '+' else -1 for c in l.split('=', 1)[1].strip()] for l in L if l[:2] in ('A=', 'B=', 'C=', 'D=')])

def write_core(path, X, header):
    with open(path, 'w') as fh:
        fh.write(header + '\n')
        for nm, row in zip('ABCD', X): fh.write(nm + '=' + ''.join('+' if t > 0 else '-' for t in row) + '\n')

if __name__ == "__main__":
    X = read_core(sys.argv[1]); H, rec = bordered_gs(X)
    print(json.dumps(dict(file=sys.argv[1], m=int(X.shape[1]), rowsums=X.sum(1).tolist(), **rec)))
