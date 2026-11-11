#!/usr/bin/env python3
"""Partial character borders: an exact experiment.

G = Q x S with Q = Z_2^a (every character of Q is a +-1 vector) and S = Z_s.  T = a set of characters of Q.
Core: four +-1 sequences on G with  sum_k |x_k^(chi)|^2 = n  at every character chi of G outside T and = 4|T| on T,
n = 4|G| + 4|T|; equivalently  sum PAF(g) = -4 sum_{rho in T} rho(g)  for g != 0.
(T = {1}: EKLO Theorem 5; T = nontrivial characters: EKLO Theorem 7; T = all: Kumar's all-character family.)
Border: K Hadamard of order b = 4|T|, E and Z (b x 4k) +-1 with rows in V_T^4 (Fourier support in T on each block),
E E^T = 4k I, K Z = E M^T (M the Goethals-Seidel array of the compressions); then F = -Z^T and
H = [[K, E'],[F', core]].  For fixed K: CP-SAT.  A found border is assembled and H H^T = n I checked.  [M]
Usage: fable_border2.py a s "T as list of bitmasks, e.g. 1,2,4" ncores time_limit [seed]
"""
import sys, time, itertools
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, 'theory')
import fable_border as FB

class Grp2(FB.Grp):
    def __init__(self, a, s):
        self.k, self.s = 2 ** a, s
        self.qmods = (2,) * a
        self.Qel = list(itertools.product(*[range(2) for _ in range(a)]))
        self.el = [(q, x) for q in self.Qel for x in range(s)]
        self.idx = {e: i for i, e in enumerate(self.el)}
        self.v = len(self.el)
        add = lambda x, y: (tuple((x[0][i] + y[0][i]) % 2 for i in range(a)), (x[1] + y[1]) % s)
        neg = lambda x: (x[0], (-x[1]) % s)
        self.addt = np.array([[self.idx[add(x, y)] for y in self.el] for x in self.el])
        self.negt = np.array([self.idx[neg(x)] for x in self.el])
        self.qof = np.array([self.Qel.index(e[0]) for e in self.el])
        self.a = a
    def chars(self):
        # character table of Q: rows indexed by bitmask m, value (-1)^{<m,q>}
        a = self.a
        return np.array([[(-1) ** sum(((m >> i) & 1) * q[a - 1 - i] for i in range(a)) for q in self.Qel] for m in range(2 ** a)])

def find_core(G, T, rng, iters=600000):
    X0 = G.chars()
    v = G.v
    Tg = -4 * sum(X0[m][G.qof] for m in T)
    Tg = Tg.copy(); Tg[G.idx[(tuple(0 for _ in range(G.a)), 0)]] = 4 * v
    X = rng.choice([-1, 1], size=(4, v))
    def energy(X):
        P = sum(G.paf(x) for x in X)
        return int(((P - Tg)[1:] ** 2).sum())
    e = energy(X); temp = 40.0
    for it in range(iters):
        b = rng.integers(4); i = rng.integers(v)
        X[b, i] *= -1
        e2 = energy(X)
        if e2 <= e or rng.random() < np.exp((e - e2) / temp):
            e = e2
            if e == 0: return X
        else:
            X[b, i] *= -1
        temp = max(1.0, temp * 0.99998)
    return None

def solve_border(M, K, G, T, tl, workers=3):
    k = G.k; b = 4 * len(T); N = 4 * k
    ch = G.chars()
    notT = [m for m in range(k) if m not in T]
    m_ = cp_model.CpModel()
    E = [[m_.NewBoolVar('') for c in range(N)] for r in range(b)]
    Z = [[m_.NewBoolVar('') for c in range(N)] for r in range(b)]
    pm = lambda x: 2 * x - 1
    for r in range(b):
        for blk in range(4):
            for mm in notT:
                m_.Add(sum(int(ch[mm][q]) * pm(E[r][blk * k + q]) for q in range(k)) == 0)
        for c in range(N):
            m_.Add(sum(int(K[r, j]) * pm(Z[j][c]) for j in range(b)) == sum(int(M[c, c2]) * pm(E[r][c2]) for c2 in range(N)))
    for r in range(b):
        for r2 in range(r + 1, b):
            xs = []
            for c in range(N):
                x = m_.NewBoolVar('')
                m_.AddBoolXOr([E[r][c], E[r2][c], x.Not()])
                xs.append(x)
            m_.Add(sum(xs) == N // 2)
    m_.Add(E[0][0] == 1)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    r = s.Solve(m_); st = s.StatusName(r)
    if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        Ev = np.array([[pm(s.Value(E[r][c])) for c in range(N)] for r in range(b)])
        Zv = np.array([[pm(s.Value(Z[r][c])) for c in range(N)] for r in range(b)])
        return st, Ev, Zv
    return st, None, None

def main():
    a = int(sys.argv[1]); s = int(sys.argv[2]); T = [int(x) for x in sys.argv[3].split(',')]
    ncores = int(sys.argv[4]); tl = float(sys.argv[5]); seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    G = Grp2(a, s); rng = np.random.default_rng(seed)
    k = G.k; b = 4 * len(T); n = 4 * G.v + b
    if b == 4: Ks = {'H4': FB.sylv(2)}
    elif b == 8: Ks = {'H8': FB.sylv(3)}
    elif b == 12: Ks = {'H12': FB.paley1(11)}
    elif b == 16: Ks = {'Sylv16': FB.sylv(4)}
    elif b == 20: Ks = {'PaleyI': FB.paley1(19), 'PaleyII': FB.paley2_9(), 'GS5': FB.gs20()}
    else: raise SystemExit('no K library')
    print(f'Q = Z_2^{a}, S = Z_{s}, T = {T}, |G| = {G.v}, b = {b}, n = {n}', flush=True)
    for ci in range(ncores):
        t0 = time.time(); X = None; tries = 0
        while X is None and tries < 20:
            X = find_core(G, T, rng); tries += 1
        if X is None:
            print(f'core {ci}: none found in {tries} annealing runs ({time.time()-t0:.0f}s)', flush=True); continue
        Cm = FB.gs_array(G, X); M = FB.compress(G, Cm)
        comp = [M[0, j * k:(j + 1) * k].tolist() for j in range(4)]
        print(f'core {ci}: found in {time.time()-t0:.0f}s ({tries} runs); compressions {comp}', flush=True)
        for nm, K0 in Ks.items():
            done = False
            for tr, K in (('', K0), ('^T', K0.T)):
                t1 = time.time()
                st, E, Z = solve_border(M, K, G, T, tl)
                msg = f'   K = {nm}{tr}: {st} ({time.time()-t1:.0f}s)'
                if E is not None:
                    H = FB.assemble(G, Cm, K, E, Z)
                    good = np.array_equal(H @ H.T, n * np.eye(n, dtype=int))
                    msg += f'; assembled H({n}): H H^T = n I is {good}'; done = True
                print(msg, flush=True)
                if done or np.array_equal(K0, K0.T): break
            if done: break

if __name__ == '__main__':
    main()
