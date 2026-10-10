r"""Independent completion and verification of a b = 1 core (report/paley_523_answer.md) [M].

Input: a core file in the repo format (header line 'v=.. N=.. b=1 m=.. cyclic Z_m ...', then A=, B=, C=, D= lines of
+/- over cyclic Z_m).  Steps (no code shared with theory/paley_b1pipe.py or the solvers):
 1. exact core check: sum_k PAF_k = 4m at 0 and -4 elsewhere (FFT, rounded, integer compare), row sums (+-2, 0, 0, 0);
 2. bordered Goethals-Seidel array over cyclic Z_m: circulants A..D, back-identity R, border rows S (x) 1^T and columns
    T^T (x) 1 (S, T in H(4)), corner forced to -S Sigma^T T / 4; search the H(4) pairs until H H^T = N I;
 3. H H^T = N I checked exactly (float64 BLAS; all entries are integers of absolute value <= N < 2^53);
 4. writes the matrix as packed bits (.npz) and prints its SHA-256.
usage: python3 theory/q523_finish.py core.txt [--save H.npz]"""
import sys, os, json, hashlib, itertools
import numpy as np

def read_core(path):
    L = [l for l in open(path).read().split('\n') if l.strip()]
    hd = dict(t.split('=', 1) for t in L[0].split() if '=' in t and t.split('=')[0] in ('v', 'N', 'm'))
    v, N, m = int(hd['v']), int(hd['N']), int(hd['m'])
    X = np.array([[1 if c == '+' else -1 for c in l.split('=', 1)[1].strip()] for l in L[1:5]], np.int64)
    assert X.shape == (4, m) and N == 4 * v and m == v - 1
    return v, N, m, X

def core_check(X):
    m = X.shape[1]
    P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in X)
    return bool(P[0] == 4 * m and (P[1:] == -4).all()), sorted(np.abs(X.sum(1)).tolist())

def bordered(X):
    m = X.shape[1]; N = 4 * m + 4
    circ = lambda x: np.array([np.roll(x, j) for j in range(m)])      # circ[r, c] = x[c - r]
    R = np.fliplr(np.eye(m, dtype=np.int64))
    A, B, C, D = [circ(x).astype(np.int64) for x in X]
    core = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    Sig = np.array([[core[I * m:(I + 1) * m, J * m:(J + 1) * m].sum(1)[0] for J in range(4)] for I in range(4)])
    H4 = [np.array(t).reshape(4, 4) for t in itertools.product((1, -1), repeat=16)]
    H4 = [h for h in H4 if np.array_equal(h @ h.T, 4 * np.eye(4, dtype=np.int64))]
    for S in H4:
        for T in H4:
            num = -(S @ Sig.T @ T)
            if (num % 4).any() or not (np.abs(num // 4) == 1).all(): continue
            H = np.zeros((N, N), np.int8)
            H[:4, :4] = num // 4; H[:4, 4:] = np.repeat(S, m, axis=1); H[4:, :4] = np.repeat(T.T, m, axis=0); H[4:, 4:] = core
            Hf = H.astype(np.float64)
            if np.array_equal(Hf @ Hf.T, N * np.eye(N)): return H, S, T
    return None, None, None

def main(path, save=None):
    v, N, m, X = read_core(path)
    ok, rows = core_check(X)
    out = dict(file=os.path.basename(path), v=v, N=N, core_ok=ok, abs_rowsums=rows)
    if ok and rows == [0, 0, 0, 2]:
        H, S, T = bordered(X)
        out["hadamard"] = H is not None
        if H is not None:
            assert set(np.unique(H).tolist()) == {-1, 1}
            bits = np.packbits(H > 0, axis=1)
            out["sha256_packed_rows"] = hashlib.sha256(bits.tobytes()).hexdigest()
            out["corner"] = H[:4, :4].tolist()
            if save: np.savez_compressed(save, bits=bits, N=N)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; save = None
    if '--save' in a: i = a.index('--save'); save = a[i + 1]; a = a[:i] + a[i + 2:]
    for p in a: print(json.dumps(main(p, save)))
