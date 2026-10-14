"""Independent verification of a b = 0 Goethals-Seidel core and construction of H(4v) (pattern of
scripts/verify_h2092.py; shares no code with theory/).  For the v = 933 run:
    python3 scripts/verify_h3732.py logs/ab/b0z3/seq_<name>_Z933.txt [OUT.npz]
Works for any length v (e.g. the same-shape controls v = 93, 213).

1. Read the core: four +-1 sequences A, B, C, D on cyclic Z_v (the search writes them already mapped from
   Z_3 x Z_q to Z_v by CRT).
2. Check sum_k PAF_k = 4v at 0 and 0 elsewhere, exactly (integer dot products, no FFT).
3. Goethals-Seidel array of the circulants A, B, C, D with the back-identity R:
       [[ A,    BR,    CR,    DR  ],
        [-BR,   A,     D^T R, -C^T R],
        [-CR,  -D^T R, A,     B^T R],
        [-DR,   C^T R, -B^T R, A   ]]
4. Check H H^T = 4v I exactly (int64, row blocks), save H as packed bits with its SHA-256.
usage: python3 scripts/verify_h3732.py SEQ.txt [OUT.npz]   (default OUT: data/H<N>_packed.npz)"""
import hashlib, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(path):
    rows = [l.split('=', 1)[1].strip() for l in open(path) if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    assert len(rows) == 4 and len({len(r) for r in rows}) == 1 and set(''.join(rows)) <= {'+', '-'}
    return np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)


def paf(x):
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(len(x))], dtype=np.int64)


def exact_hadamard(H):
    N = H.shape[0]; H = H.astype(np.int64)
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E): return False
    return True


def main(path, out=None):
    X = load(path); v = X.shape[1]; N = 4 * v
    P = sum(paf(x) for x in X)
    ok = bool(P[0] == 4 * v and (P[1:] == 0).all())
    print(json.dumps(dict(seq=os.path.basename(path), v=v, N=N, paf_ok=ok, paf_offzero_values=sorted(set(P[1:].tolist())),
                          rowsums=X.sum(1).tolist(), sum_of_squared_rowsums=int((X.sum(1) ** 2).sum()))), flush=True)
    assert ok and int((X.sum(1) ** 2).sum()) == 4 * v
    circ = lambda x: np.array([np.roll(x, j) for j in range(v)], dtype=np.int64)     # circ[r, c] = x[(c - r) mod v]
    R = np.fliplr(np.eye(v, dtype=np.int64))
    A, B, C, D = [circ(x) for x in X]
    H = np.block([[A, B @ R, C @ R, D @ R],
                  [-B @ R, A, D.T @ R, -C.T @ R],
                  [-C @ R, -D.T @ R, A, B.T @ R],
                  [-D @ R, C.T @ R, -B.T @ R, A]])
    assert set(np.unique(H).tolist()) == {-1, 1}
    had = exact_hadamard(H)
    rec = dict(hadamard=had, N=N)
    if had:
        bits = np.packbits(H > 0, axis=1)
        out = out or os.path.join(ROOT, 'data', f'H{N}_packed.npz')
        np.savez_compressed(out, bits=bits, N=N)
        rec.update(saved=os.path.relpath(out, ROOT), sha256_of_packed_rows=hashlib.sha256(bits.tobytes()).hexdigest())
    print(json.dumps(rec), flush=True)
    return had


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None) else 1)
