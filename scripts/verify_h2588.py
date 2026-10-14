"""Independent verification of a b = 1 core on cyclic Z_m and construction of H(4(m+1)); shares no code with theory/
(numpy only; written for the 647 x353 run, works for any such core, e.g. the H(2092) core as a test).

1. Read the core: four +-1 sequences A, B, C, D on cyclic Z_m (lines "A=+-+...", ...).
2. Check sum_k PAF_k = 4m at lag 0 and -4 at every other lag (direct O(m^2) integer PAF, no FFT), and row sums
   (2, 0, 0, 0) up to order and sign.
3. Bordered Goethals-Seidel: core = the GS array of the circulants A, B, C, D with the back-identity R; border rows
   S (x) 1_m^T and columns T^T (x) 1_m with S, T in H(4); the corner is forced as -S Sigma^T T / 4 (Sigma = block row
   sums of the core array).  First the border of pub/H2292/border.json is tried (it fits every core whose block A has
   row sum +2), then all H(4) pairs.
4. Check H H^T = N I exactly in int64 (row blocks) and save H as packed bits with the SHA-256 of the packed rows.
usage: python3 scripts/verify_h2588.py CORE.txt [OUT.npz]   (default OUT: data/H<N>_packed.npz)"""
import itertools, hashlib, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def read_core(path):
    rows = [l.split('=', 1)[1].strip() for l in open(path) if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    X = np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)
    assert X.shape[0] == 4 and len(set(len(r) for r in rows)) == 1 and set(''.join(rows)) <= {'+', '-'}
    return X

def paf_direct(x):
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(len(x))], dtype=np.int64)

def gs_core(X):
    m = X.shape[1]
    circ = lambda x: np.array([np.roll(x, j) for j in range(m)])            # circ[r, c] = x[(c - r) mod m]
    R = np.fliplr(np.eye(m, dtype=np.int64))
    A, B, C, D = [circ(x).astype(np.int64) for x in X]
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])

def hadamard_exact(H):
    N = H.shape[0]; H64 = H.astype(np.int64)
    for i in range(0, N, 128):
        G = H64[i:i + 128] @ H64.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E): return False
    return True

def main(core_path, out=None):
    X = read_core(core_path); m = X.shape[1]; N = 4 * (m + 1)
    P = sum(paf_direct(x) for x in X)
    rows = X.sum(1).tolist()
    ok = bool(P[0] == 4 * m and (P[1:] == -4).all()) and sorted(abs(r) for r in rows) == [0, 0, 0, 2]
    print(json.dumps(dict(core=os.path.relpath(core_path, ROOT) if core_path.startswith(ROOT) else core_path, m=m, N=N,
                          paf_lag0=int(P[0]), paf_offzero_values=sorted(set(P[1:].tolist())), rowsums=rows, core_ok=ok)), flush=True)
    if not ok: return False
    core = gs_core(X)
    Sig = np.array([[core[I * m:(I + 1) * m, J * m:(J + 1) * m].sum(1)[0] for J in range(4)] for I in range(4)])
    H4 = [np.array(t).reshape(4, 4) for t in itertools.product((1, -1), repeat=16)]
    H4 = [h for h in H4 if np.array_equal(h @ h.T, 4 * np.eye(4, dtype=np.int64))]
    cands = []
    bj = os.path.join(ROOT, 'pub', 'H2292', 'border.json')
    if os.path.exists(bj):
        b = json.load(open(bj)); cands.append((np.array(b['S']), np.array(b['T']), 'pub/H2292/border.json'))
    cands += [(S, T, 'search') for S in H4 for T in H4]
    for S, T, src in cands:
        num = -(S @ Sig.T @ T)
        if (num % 4).any() or not (np.abs(num // 4) == 1).all(): continue
        H = np.zeros((N, N), np.int8)
        H[:4, :4] = num // 4; H[:4, 4:] = np.repeat(S, m, axis=1); H[4:, :4] = np.repeat(T.T, m, axis=0); H[4:, 4:] = core
        if hadamard_exact(H):
            assert set(np.unique(H).tolist()) == {-1, 1}
            bits = np.packbits(H > 0, axis=1)
            out = out or os.path.join(ROOT, 'data', f'H{N}_packed.npz')
            np.savez_compressed(out, bits=bits, N=N)
            print(json.dumps(dict(hadamard=True, N=N, border=src, saved=os.path.relpath(out, ROOT),
                                  sha256_of_packed_rows=hashlib.sha256(bits.tobytes()).hexdigest(),
                                  S=S.tolist(), T=T.tolist(), corner=H[:4, :4].tolist())), flush=True)
            return True
    print(json.dumps(dict(hadamard=False)))
    return False

if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None) else 1)
