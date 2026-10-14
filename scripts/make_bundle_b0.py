"""Write a pub/H<N>/ bundle for a b = 0 Goethals-Seidel Hadamard matrix (format of pub/H2092/; no border for b = 0).
Shares no code with theory/.
    python3 scripts/make_bundle_b0.py SEQ.txt "slice description"
writes pub/H<N>/sequences_Z<v>.txt, pub/H<N>/H<N>_normalized.txt.gz, pub/H<N>/verify.py and runs the verifier."""
import gzip, os, sys, subprocess
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

VERIFY = '''"""Standalone verification of H(%(N)d).  Needs only Python 3 and numpy.

  python3 verify.py            # (1) the four sequences: sum of periodic autocorrelations = 0 off zero
                               # (2) build H(%(N)d) by the Goethals-Seidel array of the four circulants and check
                               #     H H^T = %(N)d I exactly (integer arithmetic)
                               # (3) check the published normalized matrix H%(N)d_normalized.txt.gz the same way, and
                               #     that it is the normalization of the matrix built in (2)
"""
import gzip, numpy as np

def load_seqs(path='sequences_Z%(v)d.txt'):
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
    rows = [l.strip() for l in gzip.open('H%(N)d_normalized.txt.gz', 'rt') if l.strip()]
    Hn = np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)
    Hb = H * H[:, [0]]; Hb = Hb * Hb[[0], :]
    print('published normalized matrix: order', Hn.shape, ' first row/col all +1:', bool((Hn[0] == 1).all() and (Hn[:, 0] == 1).all()),
          ' Hadamard:', is_hadamard(Hn), ' equals normalization of built H:', bool(np.array_equal(Hn, Hb)))
'''


def main(seq, desc):
    rows = [l.split('=', 1)[1].strip() for l in open(seq) if l[:2] in ('A=', 'B=', 'C=', 'D=')]
    X = np.array([[1 if c == '+' else -1 for c in r] for r in rows], dtype=np.int64)
    v = X.shape[1]; N = 4 * v
    d = os.path.join(ROOT, 'pub', f'H{N}'); os.makedirs(d, exist_ok=True)
    rs = X.sum(1).tolist()
    with open(os.path.join(d, f'sequences_Z{v}.txt'), 'w') as fh:
        fh.write(f"# Four +-1 circulant sequences of length {v} (cyclic Z_{v}, natural order), b=0 Goethals-Seidel core of H({N})\n")
        fh.write(f"# sum_k PAF_k(s) = 0 for s != 0; row sums (A,B,C,D) = ({','.join(map(str, rs))})\n")
        for line in desc.split('\n'): fh.write(f"# {line}\n")
        for lab, r in zip('ABCD', rows): fh.write(f"{lab}={r}\n")
    open(os.path.join(d, 'verify.py'), 'w').write(VERIFY % dict(N=N, v=v))
    # normalized matrix
    circ = lambda x: np.array([np.roll(x, j) for j in range(v)], dtype=np.int8)
    A, B, C, D = [circ(x) for x in X]
    R = lambda M: M[:, ::-1]
    H = np.block([[A, R(B), R(C), R(D)], [-R(B), A, R(D.T), -R(C.T)], [-R(C), -R(D.T), A, R(B.T)], [-R(D), R(C.T), -R(B.T), A]])
    Hn = H * H[:, [0]]; Hn = Hn * Hn[[0], :]
    with gzip.open(os.path.join(d, f'H{N}_normalized.txt.gz'), 'wt') as fh:
        for r in Hn: fh.write(''.join('+' if t > 0 else '-' for t in r) + "\n")
    print(subprocess.run([sys.executable, 'verify.py'], cwd=d, capture_output=True, text=True).stdout)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "")
