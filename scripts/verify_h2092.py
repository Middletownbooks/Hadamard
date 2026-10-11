"""Independent verification of a q523 b = 1 core and construction of H(4v) (pattern of scripts/verify_h2292.py; shares
no code with theory/).  For the 523 run: python3 scripts/verify_h2092.py data/core_v523_b1_q523_<tag>_w<w>.txt

1. Read the core: four +-1 sequences A, B, C, D on cyclic Z_m (the q523 writer already maps Z18 x Zq to Z_m by CRT).
2. Check sum_k PAF_k = 4m at 0 and -4 elsewhere, and row sums (+-2, 0, 0, 0) in some order (b = 1 core).
3. Bordered Goethals-Seidel: core = GS array of circulants A, B, C, D with the back-identity R; border rows
   S (x) 1_m^T and columns T^T (x) 1_m with S, T in H(4); corner forced as -S Sigma^T T / 4 (Sigma = block row sums
   of the core array).  Search all H(4) pairs S, T until the corner is +-1 and H H^T = N I.
4. Check H H^T = N I exactly (float64 BLAS; every entry is an integer of size <= N, so float64 is exact) and save H
   as packed bits with its SHA-256.
usage: python3 scripts/verify_h2092.py CORE.txt [OUT.npz]   (default OUT: data/H<N>_packed.npz)"""
import itertools, hashlib, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main(core_path, out=None):
    L = open(core_path).read().split('\n')
    X = [np.array([1 if c == '+' else -1 for c in l.split('=', 1)[1].strip()]) for l in L if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    assert len(X) == 4 and len({len(x) for x in X}) == 1
    m = len(X[0]); v = m + 1; N = 4 * v
    Xc = np.array(X)
    P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in Xc)
    core_ok = bool(P[0] == 4 * m and (P[1:] == -4).all())
    rows = Xc.sum(1).tolist()
    print(json.dumps(dict(core=os.path.basename(core_path), m=m, N=N, paf_ok=core_ok,
                          paf_offzero_values=sorted(set(P[1:].tolist())), rowsums=rows)), flush=True)
    assert core_ok and sorted(np.abs(rows).tolist()) == [0, 0, 0, 2]
    # GS array over cyclic Z_m
    circ = lambda x: np.array([np.roll(x, j) for j in range(m)])            # circ[r, c] = x[c - r]
    R = np.fliplr(np.eye(m, dtype=np.int64))
    A, B, C, D = [circ(x).astype(np.int64) for x in Xc]
    core = np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])
    Sig = np.array([[core[I * m:(I + 1) * m, J * m:(J + 1) * m].sum(1)[0] for J in range(4)] for I in range(4)])
    H4 = [np.array(t).reshape(4, 4) for t in itertools.product((1, -1), repeat=16)]
    H4 = [h for h in H4 if np.array_equal(h @ h.T, 4 * np.eye(4, dtype=np.int64))]
    for S in H4:
        for T in H4:
            num = -(S @ Sig.T @ T)
            if (num % 4).any() or not (np.abs(num // 4) == 1).all(): continue
            H = np.zeros((N, N), np.int8)
            H[:4, :4] = num // 4
            H[:4, 4:] = np.repeat(S, m, axis=1); H[4:, :4] = np.repeat(T.T, m, axis=0); H[4:, 4:] = core
            Hf = H.astype(np.float64)
            if np.array_equal(Hf @ Hf.T, N * np.eye(N)):
                assert set(np.unique(H).tolist()) == {-1, 1}
                bits = np.packbits(H > 0, axis=1)
                out = out or os.path.join(ROOT, 'data', f'H{N}_packed.npz')
                np.savez_compressed(out, bits=bits, N=N)
                print(json.dumps(dict(hadamard=True, N=N, saved=os.path.relpath(out, ROOT),
                                      sha256_of_packed_rows=hashlib.sha256(bits.tobytes()).hexdigest(),
                                      corner=H[:4, :4].tolist())), flush=True)
                return True
    print(json.dumps(dict(hadamard=False)))
    return False

if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None) else 1)
