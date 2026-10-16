"""Standalone verification of H(3532).  Needs only Python 3 and numpy.

  python3 verify.py            # (1) the four sequences: sum of periodic autocorrelations = 0 off zero
                               # (2) build H(3532) by the Goethals-Seidel array of the four circulants and check
                               #     H H^T = 3532 I exactly (integer arithmetic)
                               # (3) check the published normalized matrix H3532_normalized.txt.gz the same way, and
                               #     that it is the normalization of the matrix built in (2)
"""
import gzip, numpy as np

def load_seqs(path='sequences_Z883.txt'):
    rows = [l.split('=', 1)[1].strip() for l in open(path) if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    return np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)

def paf(x):
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(len(x))])

def is_hadamard(H):
    """exact: H H^T computed in int64 in row blocks"""
    N = H.shape[0]; H = H.astype(np.int64)
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E): return False
    return True

def build(X):
    v = X.shape[1]
    circ = lambda x: np.array([np.roll(x, j) for j in range(v)])          # circ[r, c] = x[(c - r) mod v]
    R = np.fliplr(np.eye(v, dtype=np.int64))                             # back-identity
    A, B, C, D = [circ(x) for x in X]
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])

if __name__ == '__main__':
    X = load_seqs(); v = X.shape[1]
    P = sum(paf(x) for x in X)
    print('length', v, ' PAF sum at 0:', P[0], ' off-zero values:', sorted(set(P[1:].tolist())), ' row sums:', X.sum(1).tolist())
    assert P[0] == 4 * v and (P[1:] == 0).all()
    H = build(X); N = H.shape[0]
    print('built H of order', N, ' entries', sorted(set(np.unique(H).tolist())), ' Hadamard:', is_hadamard(H))
    rows = [l.strip() for l in gzip.open('H3532_normalized.txt.gz', 'rt') if l.strip()]
    Hn = np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)
    Hb = H * H[:, [0]]; Hb = Hb * Hb[[0], :]
    print('published normalized matrix: order', Hn.shape, ' first row/col all +1:', bool((Hn[0] == 1).all() and (Hn[:, 0] == 1).all()),
          ' Hadamard:', is_hadamard(Hn), ' equals normalization of built H:', bool(np.array_equal(Hn, Hb)))
