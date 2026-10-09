"""Standalone verification of H(2292).  Needs only Python 3 and numpy.

  python3 verify.py            # (1) the four sequences: sum of periodic autocorrelations = -4 off zero
                               # (2) build H(2292) by the bordered Goethals-Seidel construction with the recorded
                               #     border (border.json) and check H H^T = 2292 I exactly (integer arithmetic)
                               # (3) check the published normalized matrix H2292_normalized.txt.gz the same way, and
                               #     that it is the normalization of the matrix built in (2)
"""
import gzip, json, numpy as np

def load_seqs(path='sequences_Z572.txt'):
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

def build(X, S, T, Cn):
    m = X.shape[1]; N = 4 * m + 4
    circ = lambda x: np.array([np.roll(x, j) for j in range(m)])          # circ[r, c] = x[(c - r) mod m]
    R = np.fliplr(np.eye(m, dtype=np.int64))                             # back-identity
    A, B, C, D = [circ(x) for x in X]
    core = np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])
    H = np.zeros((N, N), np.int64)
    H[:4, :4] = Cn; H[:4, 4:] = np.repeat(S, m, axis=1); H[4:, :4] = np.repeat(np.array(T).T, m, axis=0); H[4:, 4:] = core
    return H

if __name__ == '__main__':
    X = load_seqs(); m = X.shape[1]
    P = sum(paf(x) for x in X)
    print('length', m, ' PAF sum at 0:', P[0], ' off-zero values:', sorted(set(P[1:].tolist())), ' row sums:', X.sum(1).tolist())
    assert P[0] == 4 * m and (P[1:] == -4).all()
    b = json.load(open('border.json')); S, T, Cn = (np.array(b[k]) for k in ('S', 'T', 'corner'))
    H = build(X, S, T, Cn); N = H.shape[0]
    print('built H of order', N, ' entries', sorted(set(np.unique(H).tolist())), ' Hadamard:', is_hadamard(H))
    rows = [l.strip() for l in gzip.open('H2292_normalized.txt.gz', 'rt') if l.strip()]
    Hn = np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)
    Hb = H * H[:, [0]]; Hb = Hb * Hb[[0], :]
    print('published normalized matrix: order', Hn.shape, ' first row/col all +1:', bool((Hn[0] == 1).all() and (Hn[:, 0] == 1).all()),
          ' Hadamard:', is_hadamard(Hn), ' equals normalization of built H:', bool(np.array_equal(Hn, Hb)))
