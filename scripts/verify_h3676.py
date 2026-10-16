"""Independent verification of the b = 0 Goethals-Seidel core of length 919 and of H(3676) (pattern of
scripts/verify_h3732.py; shares no code with the construction scripts).
    python3 scripts/verify_h3676.py data/core_v919_b0_twist51_0.txt [OUT.npz]

1. Read the core: four +-1 sequences A, B, C, D on cyclic Z_919.
2. Check sum_k PAF_k = 3676 at 0 and 0 elsewhere, exactly (integer dot products, no FFT), and the row sums.
3. Check the slice: with u = 8 (order 51 mod 919), A(8x) = A(x), B(8^3 x) = B(x), C(x) = B(8^-1 x), D(x) = B(8^-2 x).
4. Goethals-Seidel array of the circulants A, B, C, D with the back-identity R:
       [[ A,    BR,    CR,    DR  ],
        [-BR,   A,     D^T R, -C^T R],
        [-CR,  -D^T R, A,     B^T R],
        [-DR,   C^T R, -B^T R, A   ]]
5. Check H H^T = 3676 I exactly (int64, row blocks); save H as packed bits with its SHA-256.
usage: python3 scripts/verify_h3676.py SEQ.txt [OUT.npz]   (default OUT: data/H3676_packed.npz)"""
import hashlib, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(path):
    rows = [l.split('=', 1)[1].strip() for l in open(path) if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    assert len(rows) == 4 and len({len(r) for r in rows}) == 1 and set(''.join(rows)) <= {'+', '-'}
    return np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)


def paf(x):
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(len(x))], dtype=np.int64)


def circ(x):
    v = len(x)
    return np.array([np.roll(x, j) for j in range(v)], dtype=np.int64)       # circ[r, c] = x[(c - r) mod v]


def gs_array(X):
    v = X.shape[1]
    R = np.fliplr(np.eye(v, dtype=np.int64))
    A, B, C, D = (circ(x) for x in X)
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def is_hadamard(H):
    N = H.shape[0]; H = H.astype(np.int64)
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E):
            return False
    return True


def main():
    path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'data', 'H3676_packed.npz')
    X = load(path); v = X.shape[1]
    assert v == 919
    P = sum(paf(x) for x in X)
    rec = {'v': v, 'paf0': int(P[0]), 'paf_off_zero': sorted(set(P[1:].tolist())), 'row_sums': X.sum(1).tolist()}
    ok_paf = bool(P[0] == 4 * v and (P[1:] == 0).all())
    u = 8; ui = pow(u, -1, v); xs = np.arange(v)
    A, B, C, D = X
    rec['order_of_u'] = next(k for k in range(1, v) if pow(u, k, v) == 1)
    rec['A_invariant_under_u'] = bool(np.array_equal(A[(u * xs) % v], A))
    rec['B_invariant_under_u3'] = bool(np.array_equal(B[(pow(u, 3, v) * xs) % v], B))
    rec['C_is_B_of_u^-1_x'] = bool(np.array_equal(C, B[(ui * xs) % v]))
    rec['D_is_B_of_u^-2_x'] = bool(np.array_equal(D, B[(ui * ui * xs) % v]))
    H = gs_array(X)
    rec['order'] = int(H.shape[0]); rec['entries'] = sorted(set(np.unique(H).tolist()))
    rec['hadamard'] = is_hadamard(H)
    bits = np.packbits(H > 0, axis=1)
    rec['sha256_packed_rows'] = hashlib.sha256(bits.tobytes()).hexdigest()
    np.savez_compressed(out, bits=bits, N=H.shape[0])
    rec['saved'] = os.path.relpath(out, ROOT)
    rec['ok'] = bool(ok_paf and rec['hadamard'] and rec['A_invariant_under_u'] and rec['B_invariant_under_u3']
                     and rec['C_is_B_of_u^-1_x'] and rec['D_is_B_of_u^-2_x'] and rec['order_of_u'] == 51)
    print(json.dumps(rec, indent=1))
    sys.exit(0 if rec['ok'] else 1)


if __name__ == '__main__':
    main()
