"""Assemble full Hadamard matrices from slice solutions over an abelian group G and verify H H^T = N I exactly.

lp : two G-developed cores A, B with PAF_A + PAF_B = -2, row sums 1  ->  order 2|G|+2
     H = [[-1,-1, e', e'], [-1, 1, e', -e'], [e, e, A, B], [e, -e, B', -A']]
b1 : four sequences with sum PAF = -4  ->  bordered blocks W_k = [[alpha, beta e'], [beta e, X_k]] in the
     Goethals-Seidel array, order 4(|G|+1)  (row sums (2,0,0,0) after sign normalisation, or all +-1)
b0 : four sequences with sum PAF = 0   ->  Goethals-Seidel array, order 4|G|
Group-developed matrix: X[x, y] = f(y - x); R[x, y] = [x + y = 0].
"""
import numpy as np


def developed(G, f):
    D = G.idx(G.coords[None, :, :] - G.coords[:, None, :])  # D[x,y] = y - x
    return np.asarray(f)[D]


def reversal(G):
    n = G.n
    R = np.zeros((n, n), dtype=np.int64)
    R[np.arange(n), G.neg_perm()] = 1
    return R


def gs_array(A, B, C, D, R):
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def build_lp(G, a, b):
    a = np.asarray(a) * (1 if np.sum(a) == 1 else -1)
    b = np.asarray(b) * (1 if np.sum(b) == 1 else -1)
    assert a.sum() == 1 and b.sum() == 1
    A = developed(G, a)
    B = developed(G, b)
    n = G.n
    H = np.zeros((2 * n + 2, 2 * n + 2), dtype=np.int64)
    H[0, :2] = [-1, -1]
    H[0, 2:] = 1
    H[1, :2] = [-1, 1]
    H[1, 2:2 + n] = 1
    H[1, 2 + n:] = -1
    H[2:2 + n, 0] = 1
    H[2:2 + n, 1] = 1
    H[2:2 + n, 2:2 + n] = A
    H[2:2 + n, 2 + n:] = B
    H[2 + n:, 0] = 1
    H[2 + n:, 1] = -1
    H[2 + n:, 2:2 + n] = B.T
    H[2 + n:, 2 + n:] = -A.T
    return H


def build_b1(G, seqs):
    n = G.n
    seqs = [np.asarray(s) for s in seqs]
    sums = [int(s.sum()) for s in seqs]
    R0 = reversal(G)
    R = np.zeros((n + 1, n + 1), dtype=np.int64)
    R[0, 0] = 1
    R[1:, 1:] = R0
    Ws = []
    if n % 2 == 0:
        # normalise: the block with |sum| = 2 gets sum +2
        k2 = [k for k in range(4) if abs(sums[k]) == 2]
        assert len(k2) == 1 and all(sums[k] == 0 for k in range(4) if k != k2[0])
        for k in range(4):
            x = seqs[k] * (1 if sums[k] >= 0 else -1)
            X = developed(G, x)
            alpha, beta = (1, -1) if k == k2[0] else (1, 1)
            W = np.zeros((n + 1, n + 1), dtype=np.int64)
            W[0, 0] = alpha
            W[0, 1:] = beta
            W[1:, 0] = beta
            W[1:, 1:] = X
            Ws.append(W)
    else:
        # all sums +-1: normalise to +1; alpha_k, beta_k chosen by search below
        xs = [s * (1 if s.sum() > 0 else -1) for s in seqs]
        import itertools
        for signs in itertools.product([1, -1], repeat=8):
            Ws = []
            for k in range(4):
                X = developed(G, xs[k])
                W = np.zeros((n + 1, n + 1), dtype=np.int64)
                W[0, 0] = signs[2 * k]
                W[0, 1:] = signs[2 * k + 1]
                W[1:, 0] = signs[2 * k + 1]
                W[1:, 1:] = X
                Ws.append(W)
            H = gs_array(*Ws, R)
            if np.array_equal(H @ H.T, 4 * (n + 1) * np.eye(4 * (n + 1), dtype=np.int64)):
                return H
        return None
    return gs_array(*Ws, R)


def build_b0(G, seqs):
    R = reversal(G)
    Xs = [developed(G, s) for s in seqs]
    return gs_array(*Xs, R)


def verify(H):
    N = H.shape[0]
    return bool(np.all(np.abs(H) == 1) and np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64)))


def build(G, kind, sol):
    if kind == 'lp':
        return build_lp(G, *sol)
    if kind == 'b1':
        return build_b1(G, sol)
    return build_b0(G, sol)
