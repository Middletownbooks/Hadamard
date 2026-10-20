r"""Border problem and assembly for bordered Goethals-Seidel arrays with a character set T [P for the reduction,
M for the search] (Kumar, search/NOTES.md section 18; Apache-2.0; this code is new).

H = [[X, R], [C, M]], M the GS array over G (Kumar's layout):

    [ A    B2    C2    D2  ]      A[g, h] = a(h - g),  Y2[g, h] = y(g + h),  Y2'[g, h] = y(-g - h)
    [-B2   A     D2'  -C2' ]
    [-C2  -D2'   A     B2' ]
    [-D2   C2'  -B2'   A   ]

R[x, (i, g)] = Rh[x, i, g + K], C[(i, g), x] = Ch[i, g + K, x].  With the core condition, H H^T = 4n I iff
X X^T = 4t I, Rh Rh^T = 4q I, Ch Ch^T = I_4 (x) Gamma_Q and X Ch^T + Rh Mc^T = 0, Mc the same GS array over Q of the
coset sums y_i [P: R M^T = Rh Mc^T because R is constant on cosets].  For rows of Rh spanning V_T (the 4t-dimensional
space of block vectors with Fourier support in T), Ch = -Mc Rh^T X / (4t) and Ch Ch^T = I_4 (x) Gamma_Q follow
(Mc preserves V_T and has eigenvalue 4t there); so a border exists for the frame Rh iff
S_R = {eps in {+-1}^{4t} : Mc Rh^T eps in {+-4t}^{4q}} contains 4t mutually orthogonal vectors (the columns of X).
"""
import itertools, math, time, json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def gs_array(seqs, add, neg):
    """GS array (4m x 4m) of four sequences on a group with addition table add (m x m) and negation neg"""
    m = len(seqs[0])
    one = add[neg][:, :]                     # placeholder (unused)
    idx_diff = add[neg[:, None], np.arange(m)[None, :]]      # idx_diff[g, h] = h - g = (-g) + h
    idx_sum = add                                             # g + h
    a, b, c, d = [np.asarray(s, np.int64) for s in seqs]
    A = a[idx_diff]
    B2, C2, D2 = b[idx_sum], c[idx_sum], d[idx_sum]
    B2r, C2r, D2r = b[neg][idx_sum], c[neg][idx_sum], d[neg][idx_sum]   # y(-g-h)
    return np.block([[A, B2, C2, D2],
                     [-B2, A, D2r, -C2r],
                     [-C2, -D2r, A, B2r],
                     [-D2, C2r, -B2r, A]])


class QData:
    """Q = G/K with cosets numbered by smallest element (bt_engine.BStruct.coset), its addition table, Gamma"""

    def __init__(self, st):
        G = st.G
        self.q = st.q
        reps = st.coset_rep
        self.add = st.coset[G.add[np.ix_(reps, reps)]]
        self.neg = st.coset[G.neg()[reps]]
        self.gamma = st.gamma
        self.t = st.t
        self.Gmat = self.gamma[self.add[self.neg[:, None], np.arange(self.q)[None, :]]]   # Gamma(h - g)


def compressions(st, xs):
    return np.array([np.bincount(st.coset, weights=np.asarray(x, float), minlength=st.q).astype(np.int64) for x in xs])


def balanced_rows(q):
    return np.array([v for v in itertools.product((1, -1), repeat=q) if sum(v) == 0], np.int64)


def frames_k4(rng=None, nrand=0):
    """frames for q = 4, t = 3 (V_T = block-balanced vectors): Sylvester frames H4 (x) M with the columns of each
    block permuted (24^4 / symmetries, yielded in a fixed order), then random completions"""
    H4 = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])
    base = np.kron(H4, H4[1:])                                   # 12 x 16
    perms = list(itertools.permutations(range(4)))
    for p1, p2, p3 in itertools.product(range(24), repeat=3):
        cols = list(range(4))
        for J, pi in ((1, p1), (2, p2), (3, p3)):
            cols += [4 * J + perms[pi][k] for k in range(4)]
        yield base[:, cols]


def admissible(Mc, Rh, t):
    """all eps (first entry +1) in {+-1}^{4t} with Mc Rh^T eps in {+-4t}^{4q}"""
    n = Rh.shape[0]
    E = np.array(list(itertools.product((1, -1), repeat=n - 1)), np.int64)
    E = np.concatenate([np.ones((len(E), 1), np.int64), E], 1)
    W = (Mc @ Rh.T)                                              # 4q x 4t
    V = E @ W.T                                                  # (2^(n-1)) x 4q
    ok = np.all(np.abs(V) == 4 * t, axis=1)
    return E[ok]


def clique(S, size, tmax=10.0):
    if len(S) < size:
        return None
    Gm = (S @ S.T == 0)
    t0 = time.time(); best = [None]

    def rec(cands, chosen):
        if len(chosen) == size:
            best[0] = list(chosen); return True
        if time.time() - t0 > tmax or len(chosen) + len(cands) < size:
            return False
        for k, u in enumerate(cands):
            nc = [x for x in cands[k + 1:] if Gm[u, x]]
            chosen.append(u)
            if rec(nc, chosen):
                return True
            chosen.pop()
            if len(chosen) + len(cands) - k - 1 < size:
                return False
        return False
    rec(list(range(len(S))), [])
    return None if best[0] is None else S[best[0]]


def solve_border(st, xs, frames=None, max_frames=20000, tmax=5.0):
    """(X, Rh, Ch) for the core xs, or None.  Only q = 4, t = 3 (k = 4 family) frames are built in by default."""
    Q = QData(st)
    z = compressions(st, xs)
    Mc = gs_array(z, Q.add, Q.neg)
    t = st.t
    if frames is None:
        assert st.q == 4 and t == 3, "default frames are for the k = 4 family"
        frames = frames_k4()
    for k, Rh in enumerate(frames):
        if k >= max_frames:
            break
        S = admissible(Mc, Rh, t)
        if len(S) < 4 * t:
            continue
        Xc = clique(S, 4 * t, tmax=tmax)
        if Xc is None:
            continue
        X = Xc.T.copy()                                            # columns of X are the chosen eps
        num = -(Mc @ Rh.T @ X)
        if (num % (4 * t)).any():
            continue
        Ch = num // (4 * t)
        if not (np.abs(Ch) == 1).all():
            continue
        if not (Ch @ Ch.T == np.kron(np.eye(4, dtype=np.int64), Q.Gmat)).all():
            continue
        return X, Rh, Ch, k
    return None


def assemble(st, xs, X, Rh, Ch):
    """full H (4n x 4n) in Kumar's layout"""
    G = st.G
    M = gs_array([np.asarray(x, np.int64) for x in xs], G.add, G.neg())
    q = st.q; m = G.n
    R = np.zeros((X.shape[0], 4 * m), np.int64)
    C = np.zeros((4 * m, X.shape[0]), np.int64)
    for i in range(4):
        R[:, i * m:(i + 1) * m] = Rh[:, i * q + st.coset]
        C[i * m:(i + 1) * m, :] = Ch[i * q + st.coset, :]
    return np.block([[X, R], [C, M]])


def is_hadamard(H):
    n = H.shape[0]
    Hf = H.astype(np.float64)
    P = Hf @ Hf.T
    return bool(np.array_equal(P, n * np.eye(n)))
