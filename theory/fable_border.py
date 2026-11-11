#!/usr/bin/env python3
"""Borders for EKLO Theorem-7 cores with a tiny subgroup: an exact experiment.

G = Q x S (Q cyclic of order k or Z_2 x Z_2 for k = 4, S = Z_s), four +-1 sequences with
sum PAF = -4(k-1) on S \\ 0 and +4 off S; n = 4ks + 4(k-1).
1. cores by simulated annealing (they are plentiful at this size);
2. M = Goethals-Seidel array of the compressions (order 4k);
3. border: +-1 matrices E, Z (b x 4k), b = 4(k-1), and a Hadamard matrix K of order b with
       K Z = E M^T,   E block-balanced,   E E^T = 4k I
   (report/fable_fresh_answer.md section 6.3: then F = -Z^T and H = [[K, E'],[F', core]] is Hadamard).
   For a FIXED K this is a linear system in the 2 b 4k signs plus the orthogonality of E; CP-SAT decides it.
   K runs over representatives of the equivalence classes of Hadamard matrices of order b (given on input)
   and their transposes; E and Z absorb the equivalence.
4. a found border is assembled and H H^T = n I is checked on the full matrix.  [M]
Usage: fable_border.py k s ncores time_limit_s [seed]
"""
import sys, time, itertools, json
import numpy as np
from ortools.sat.python import cp_model

def paley1(q):
    chi = np.zeros(q, dtype=int)
    for x in range(1, q): chi[(x * x) % q] = 1
    chi = np.where(chi == 1, 1, -1); chi[0] = 0
    Qm = np.array([[chi[(j - i) % q] for j in range(q)] for i in range(q)])
    H = np.ones((q + 1, q + 1), dtype=int)
    H[1:, 1:] = Qm + np.eye(q, dtype=int) * 1
    H[1:, 0] = -1
    # H = [[1, e],[-e^T, Q + I]]
    assert np.array_equal(H @ H.T, (q + 1) * np.eye(q + 1, dtype=int))
    return H

def sylv(m):
    H = np.array([[1]])
    for _ in range(m): H = np.block([[H, H], [H, -H]])
    return H

def paley2_9():
    # F_9 = F_3[i], i^2 = -1 ; symmetric conference matrix of order 10, then [[C+I, C-I],[C-I, -C-I]]
    els = [(a, b) for a in range(3) for b in range(3)]
    def mul(x, y): return ((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)
    sq = set(mul(x, x) for x in els if x != (0, 0))
    def chi(x): return 0 if x == (0, 0) else (1 if x in sq else -1)
    Qm = np.array([[chi(((y[0] - x[0]) % 3, (y[1] - x[1]) % 3)) for y in els] for x in els])
    C = np.zeros((10, 10), dtype=int); C[0, 1:] = 1; C[1:, 0] = 1; C[1:, 1:] = Qm
    I = np.eye(10, dtype=int)
    H = np.block([[C + I, C - I], [C - I, -C - I]])
    assert np.array_equal(H @ H.T, 20 * np.eye(20, dtype=int))
    return H

def gs20():
    # Goethals-Seidel array of four circulants of order 5: the third equivalence class of H(20)
    # (classes told apart by the number of Hall sets inside 6-subsets of rows: see profile())
    seqs = [np.array(x) for x in itertools.product([1, -1], repeat=5)]
    def circ(r): return np.array([[r[(j - i) % 5] for j in range(5)] for i in range(5)])
    R = np.array([[1 if (i + j) % 5 == 0 else 0 for j in range(5)] for i in range(5)])
    A, B, C, D = [circ(seqs[i]) for i in (1, 1, 14, 20)]
    H = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R], [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
    assert np.array_equal(H @ H.T, 20 * np.eye(20, dtype=int))
    return H

def profile(H):
    """Equivalence invariant: distribution over 6-subsets of rows of the number of Hall sets (4 rows with
    |sum of products| = 12 at order 20; = n - 8 in general) they contain."""
    n = H.shape[0]
    if n > 20: return ()
    hall = set(q for q in itertools.combinations(range(n), 4) if abs(int(np.prod(H[list(q)], axis=0).sum())) == n - 8)
    cnt = {}
    for s in itertools.combinations(range(n), 6):
        c = sum(1 for q in itertools.combinations(s, 4) if q in hall)
        cnt[c] = cnt.get(c, 0) + 1
    return tuple(sorted(cnt.items()))

class Grp:
    def __init__(self, k, s, klein=False):
        self.k, self.s = k, s
        self.qmods = (2, 2) if klein else (k,)
        self.Qel = list(itertools.product(*[range(m) for m in self.qmods]))
        self.el = [(q, a) for q in self.Qel for a in range(s)]
        self.idx = {e: i for i, e in enumerate(self.el)}
        self.v = len(self.el)
        add = lambda x, y: (tuple((x[0][i] + y[0][i]) % self.qmods[i] for i in range(len(self.qmods))), (x[1] + y[1]) % s)
        neg = lambda x: (tuple((-x[0][i]) % self.qmods[i] for i in range(len(self.qmods))), (-x[1]) % s)
        self.addt = np.array([[self.idx[add(x, y)] for y in self.el] for x in self.el])
        self.negt = np.array([self.idx[neg(x)] for x in self.el])
        self.inS = np.array([all(c == 0 for c in e[0]) for e in self.el])
        self.qof = np.array([self.Qel.index(e[0]) for e in self.el])
    def paf(self, x):
        # PAF(g) = sum_h x(h) x(h+g)
        return (x[None, :] * x[self.addt]).sum(1)     # addt[g][h] = g + h
    def target(self):
        k, s, v = self.k, self.s, self.v
        T = np.where(self.inS, -4 * (k - 1), 4); T[self.idx[(tuple(0 for _ in self.qmods), 0)]] = 4 * v
        return T
    def dev(self, x):   # group-developed matrix X[g,h] = x(h-g)
        return x[self.addt[self.negt][:, :].T[:, :]] if False else np.array([[x[self.addt[self.negt[g]][h]] for h in range(self.v)] for g in range(self.v)])

def find_core(G, rng, iters=400000):
    T = G.target(); v = G.v
    X = rng.choice([-1, 1], size=(4, v))
    def energy(X):
        P = sum(G.paf(x) for x in X)
        return int(((P - T)[1:] ** 2).sum())
    e = energy(X); temp = 30.0
    for it in range(iters):
        b = rng.integers(4); i = rng.integers(v)
        X[b, i] *= -1
        e2 = energy(X)
        if e2 <= e or rng.random() < np.exp((e - e2) / temp):
            e = e2
            if e == 0: return X
        else:
            X[b, i] *= -1
        temp = max(1.0, temp * 0.99997)
    return None

def gs_array(G, X):
    A, B, C, D = [G.dev(x) for x in X]
    v = G.v
    R = np.zeros((v, v), dtype=int)
    for g in range(v): R[g, G.negt[g]] = 1
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])

def compress(G, Cm):
    k, v = G.k, G.v
    M = np.zeros((4 * k, 4 * k), dtype=int)
    rep = [int(np.nonzero((G.qof == q) & np.array([e[1] == 0 for e in G.el]))[0][0]) for q in range(k)]
    for i in range(4):
        for q in range(k):
            row = Cm[i * v + rep[q]]
            for j in range(4):
                for q2 in range(k):
                    M[i * k + q, j * k + q2] = row[j * v:(j + 1) * v][G.qof == q2].sum()
    return M

def solve_border(M, K, k, tl, workers=3):
    b = 4 * (k - 1); N = 4 * k
    m = cp_model.CpModel()
    E = [[m.NewBoolVar('') for c in range(N)] for r in range(b)]
    Z = [[m.NewBoolVar('') for c in range(N)] for r in range(b)]
    pm = lambda x: 2 * x - 1
    for r in range(b):
        for blk in range(4):
            m.Add(sum(E[r][blk * k + q] for q in range(k)) == k // 2)
        for c in range(N):
            m.Add(sum(int(K[r, j]) * pm(Z[j][c]) for j in range(b)) == sum(int(M[c, c2]) * pm(E[r][c2]) for c2 in range(N)))
    for r in range(b):
        for r2 in range(r + 1, b):
            xs = []
            for c in range(N):
                x = m.NewBoolVar('')
                m.AddBoolXOr([E[r][c], E[r2][c], x.Not()])     # x = E[r][c] XOR E[r2][c]
                xs.append(x)
            m.Add(sum(xs) == N // 2)
    m.Add(E[0][0] == 1)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    r = s.Solve(m); st = s.StatusName(r)
    if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        Ev = np.array([[pm(s.Value(E[r][c])) for c in range(N)] for r in range(b)])
        Zv = np.array([[pm(s.Value(Z[r][c])) for c in range(N)] for r in range(b)])
        return st, Ev, Zv
    return st, None, None

def assemble(G, Cm, K, E, Z):
    k, v = G.k, G.v; b = 4 * (k - 1)
    lift = np.array([j * k + G.qof[h] for j in range(4) for h in range(v)])
    E1 = E[:, lift]                 # b x 4v
    F1 = (-Z.T)[lift, :]            # 4v x b
    return np.block([[K, E1], [F1, Cm]])

def main():
    k = int(sys.argv[1]); s = int(sys.argv[2]); ncores = int(sys.argv[3]); tl = float(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    klein = (k == 4)
    G = Grp(k, s, klein); rng = np.random.default_rng(seed)
    b = 4 * (k - 1); n = 4 * G.v + b
    if b == 12: Ks = {'H12': paley1(11)}
    elif b == 20: Ks = {'PaleyI': paley1(19), 'PaleyII': paley2_9(), 'GS5': gs20()}
    elif b == 4: Ks = {'H4': sylv(2)}
    else: raise SystemExit('no K library for b = %d' % b)
    print('K invariants:', {nm: profile(K) for nm, K in Ks.items()}, flush=True)
    for ci in range(ncores):
        t0 = time.time(); X = None
        while X is None: X = find_core(G, rng)
        Cm = gs_array(G, X)
        M = compress(G, Cm)
        comp = [M[0, j * k:(j + 1) * k].tolist() for j in range(4)]
        okM = np.array_equal(M @ M.T, np.kron(np.eye(4, dtype=int), b * np.eye(k, dtype=int) + 4 * s * np.ones((k, k), dtype=int)))
        print(f'core {ci}: found in {time.time()-t0:.0f}s; first-row compressions {comp}; M M^T as predicted: {okM}', flush=True)
        for nm, K0 in Ks.items():
            for tr, K in (('', K0), ('^T', K0.T)):
                t1 = time.time()
                st, E, Z = solve_border(M, K, k, tl)
                msg = f'   K = {nm}{tr}: {st} ({time.time()-t1:.0f}s)'
                if E is not None:
                    H = assemble(G, Cm, K, E, Z)
                    good = np.array_equal(H @ H.T, n * np.eye(n, dtype=int))
                    msg += f'; assembled H({n}): H H^T = n I is {good}'
                print(msg, flush=True)
                if E is not None: break
            else: continue
            break

if __name__ == '__main__':
    main()
