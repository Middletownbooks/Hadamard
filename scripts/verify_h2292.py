"""Independent verification of the v = 573 core and construction of H(2292) (written separately from
theory/paley_b1pipe.py; shares no code with it).

1. Read data/core_v573_b1_paley_FOUND.txt: four +-1 sequences on G = Z52 x Z11, entry (g, z) at index g*11 + z.
2. Map to cyclic Z572 by CRT (i <-> (i mod 52, i mod 11)) and check sum_k PAF_k = 4*572 at 0 and -4 elsewhere,
   row sums (+-2, 0, 0, 0)  (b = 1 core).
3. Bordered Goethals-Seidel: core = GS array of circulants A, B, C, D with the back-identity R; border rows
   S (x) 1_m^T and columns T^T (x) 1_m with S, T in H(4); corner forced as -S Sigma^T T / 4 (Sigma = block row sums
   of the core array).  Search all H(4) pairs S, T until the corner is +-1 and H H^T = N I.
4. Check H H^T = 2292 I exactly (float64 BLAS; every entry is an integer of size <= 2292, so float64 is exact) and
   save H as packed bits to data/H2292_packed.npz with its SHA-256.
usage: python3 scripts/verify_h2292.py"""
import itertools, hashlib, json, os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    L = open(os.path.join(ROOT, 'data', 'core_v573_b1_paley_FOUND.txt')).read().split('\n')
    X = [np.array([1 if c == '+' else -1 for c in l.split('=', 1)[1].strip()]) for l in L if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    n, p = 52, 11; m = n * p; v = m + 1; N = 4 * v
    assert len(X) == 4 and all(len(x) == m for x in X)
    # CRT to cyclic Z572
    i = np.arange(m); src = (i % n) * p + (i % p)
    Xc = np.array([x[src] for x in X])
    P = sum(np.rint(np.fft.ifft(np.abs(np.fft.fft(x)) ** 2).real).astype(np.int64) for x in Xc)
    core_ok = bool(P[0] == 4 * m and (P[1:] == -4).all())
    rows = Xc.sum(1).tolist()
    print(json.dumps(dict(paf_ok=core_ok, paf_offzero_values=sorted(set(P[1:].tolist())), rowsums=rows)), flush=True)
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
                np.savez_compressed(os.path.join(ROOT, 'data', 'H2292_packed.npz'), bits=bits, N=N)
                print(json.dumps(dict(hadamard=True, N=N, sha256_of_packed_rows=hashlib.sha256(bits.tobytes()).hexdigest(),
                                      corner=H[:4, :4].tolist())), flush=True)
                return
    print(json.dumps(dict(hadamard=False)))

if __name__ == "__main__":
    main()
