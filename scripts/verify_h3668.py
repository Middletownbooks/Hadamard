"""Independent verification of a b = 0 Goethals-Seidel core of length v and of H(4v) (pattern of
scripts/verify_h3732.py; shares no code with theory/ or with the construction).  For v = 917:
    python3 scripts/verify_h3668.py logs/b917/seq_H3668_Z917.txt [OUT.npz] [--slice q s] [--umult u]

1. Read the core: four +-1 sequences A, B, C, D on cyclic Z_v.
2. Check sum_k PAF_k = 4v at 0 and 0 elsewhere, exactly (integer dot products, no FFT); sum of squared row sums = 4v.
3. Goethals-Seidel array of the circulants A, B, C, D with the back-identity R:
       [[ A,    BR,    CR,    DR  ],
        [-BR,   A,     D^T R, -C^T R],
        [-CR,  -D^T R, A,     B^T R],
        [-DR,   C^T R, -B^T R, A   ]]
4. Check H H^T = 4v I exactly (int64, row blocks); save H as packed bits with its SHA-256.
5. Optional --slice q s [--umult u]: with F = v / q, check that every sequence is invariant under the multipliers
   lambda of Z_v with lambda = m (mod q), lambda = 1 (mod F) for m in the index-s subgroup M of F_q^*, and, with
   --umult u, also under lambda = 1 (mod q), lambda = u (mod F); i.e. that the core is constant on the orbits of
   M (on Z_q) x <u> (on Z_F).
usage: python3 scripts/verify_h3668.py SEQ.txt [OUT.npz] [--slice q s] [--umult u]   (default OUT: data/H<N>_packed.npz)"""
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


def crt(a, m, b, n):
    """the x mod m*n with x = a (mod m), x = b (mod n), gcd(m, n) = 1"""
    return next(x for x in range(m * n) if x % m == a % m and x % n == b % n)


def multipliers(v, q, s, u):
    F = v // q; assert F * q == v and (q - 1) % s == 0
    prime_divs = [d for d in range(2, q) if (q - 1) % d == 0 and all(d % e for e in range(2, d))]
    g = next(x for x in range(2, q) if all(pow(x, (q - 1) // r, q) != 1 for r in prime_divs))
    M = sorted({pow(g, s * i, q) for i in range((q - 1) // s)})
    lam = [crt(m, q, 1, F) for m in M]
    if u is not None:
        lam.append(crt(1, q, u, F))
    return g, M, lam


def slice_check(X, q, s, u):
    v = X.shape[1]
    g, M, lam = multipliers(v, q, s, u)
    n = np.arange(v)
    ok = all(bool((x[(l * n) % v] == x).all()) for x in X for l in lam)
    return dict(primitive_root_mod_q=g, M=M, multipliers_mod_v=lam, invariant=ok)


def main(path, out=None, sl=None, u=None):
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
    if sl: rec['slice_q_s_u'] = [sl[0], sl[1], u]; rec['slice'] = slice_check(X, sl[0], sl[1], u)
    print(json.dumps(rec), flush=True)
    return had and (not sl or rec['slice']['invariant'])


if __name__ == "__main__":
    a = sys.argv[1:]; sl = None; u = None
    if '--slice' in a:
        i = a.index('--slice'); sl = (int(a[i + 1]), int(a[i + 2])); a = a[:i] + a[i + 3:]
    if '--umult' in a:
        i = a.index('--umult'); u = int(a[i + 1]); a = a[:i] + a[i + 2:]
    sys.exit(0 if main(a[0], a[1] if len(a) > 1 else None, sl, u) else 1)
