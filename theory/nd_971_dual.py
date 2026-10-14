r"""Exact decision of the dual compression of 971 x469 (along Z_10 onto Z_97) [P for the reduction, M].

Slice 971 x469: M = <469> on Z_970 = Z_2 x Z_5 x Z_97 (469 = -1 mod 5, order 12 mod 97 with -1 in <81>).  The
compression along Z_10, X_k(z) = sum_{g in Z_10} x_k(g, z), is invariant under <81> (8 classes C_j of size 12 and z = 0)
and satisfies  sum_k PAF_{Z_97}(X_k) = (4m - 4*9, -40, ..., -40),  rows (2, 0, 0, 0).
Reduction [P]:
  * X_k(0) is a signed sum of the parts {1, 1, 2, 2, 2, 2} (even, |.| <= 10) and X_k on C_j a sum of 10 signs (even);
  * rows: X_k(0) + 12 sum_j X_k(C_j) = R_k with |X_k(0)| <= 10 forces X_k(0) = R_k (block 0: 2; others: 0; the value
    -10 for block 0 would need sum_j X = 1, odd) and sum_j X_k(C_j) = 0;
  * with X = 2Y on the classes: lag 0 gives sum_{k,j} Y_{k,j}^2 = (3844 - 4) / 48 = 80, and every other lag h gives
    sum_k [ X_k(0) (X_k(h) + X_k(-h)) + sum_{z != 0, -h} X_k(z) X_k(z + h) ] = -40, i.e. integer equations in Y.
Search: per block every Y in [-5, 5]^8 with sum 0 and energy e = sum Y^2 <= 80 (exact list), signature = PAF of the
expanded sequence at one lag per <81, -1>-class (8 lags); exhaustive 2+2 join by energy split (e0 + e1 = E,
e2 + e3 = 80 - E), up to the symmetries (row-0 blocks: negation; blocks 1..3 sorted).  No solution => the dual
compression, hence the slice 971 x469, is empty.
usage: python3 theory/nd_971_dual.py [--plant N]"""
import sys, os, json, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import lag_reps

q = 97; g81 = 81
H = sorted({pow(81, k, q) for k in range(12)})                     # <81>, order 12, contains -1

def classes():
    """class index of every nonzero z (cosets of <81> in Z_97^*), -1 at 0"""
    cl = -np.ones(q, np.int64); j = 0
    for z in range(1, q):
        if cl[z] >= 0: continue
        for h in H: cl[z * h % q] = j
        j += 1
    assert j == 8
    return cl

CL = classes()
LAGS = lag_reps(q, H)[1:]                                           # 8 nonzero lag classes

def expand(x0, Y):
    """X on Z_97 from X(0) = x0 and class values 2*Y (rows of Y: vectors of 8)"""
    Y = np.atleast_2d(Y); X = np.zeros((len(Y), q), np.int64)
    X[:, 0] = x0; X[:, 1:] = 2 * Y[:, CL[1:]]
    return X

def block_list(x0, emax=80):
    vals = range(-5, 6); out = []
    for t in itertools.product(vals, repeat=7):
        y8 = -sum(t)
        if -5 <= y8 <= 5:
            e = sum(v * v for v in t) + y8 * y8
            if e <= emax: out.append(t + (y8,))
    Y = np.array(out, np.int64)
    X = expand(x0, Y)
    S = np.stack([(X * np.roll(X, -h, axis=1)).sum(1) for h in LAGS], 1)
    e = (Y * Y).sum(1)
    return Y, S, e

def canon_neg(Y):
    """True for the lexicographically least of {Y, -Y} (row-0 blocks are defined up to negation)"""
    d = (-Y) - Y; first = np.argmax(d != 0, axis=1); has = (d != 0).any(1)
    return ~has | (d[np.arange(len(Y)), first] > 0)

def solve(T=None, x0s=(2, 0, 0, 0), E=80, canon=True, log=print):
    T = np.full(len(LAGS), -40, np.int64) if T is None else np.asarray(T, np.int64)
    t0 = time.time()
    lists = {}
    for x0 in sorted(set(x0s)): lists[x0] = block_list(x0, E)
    Y0, S0, e0 = lists[x0s[0]]
    Y1, S1, e1 = lists[x0s[1]]
    if canon:
        keep = canon_neg(Y1); Y1, S1, e1 = Y1[keep], S1[keep], e1[keep]
    log(json.dumps(dict(stage="lists", block0=len(Y0), row0_canonical=len(Y1), secs=round(time.time() - t0, 1))))
    w = np.random.default_rng(5).integers(1, 2 ** 61, len(LAGS), dtype=np.int64)
    h0 = S0 @ w; h1 = S1 @ w; hT = T @ w
    sols = []; npairs = 0
    for EL in range(E + 1):
        ER = E - EL
        # left: block 0 (energy a) + block 1 (energy EL - a)
        Lk = []; Li = []
        for a in range(EL + 1):
            I = np.nonzero(e0 == a)[0]; J = np.nonzero(e1 == EL - a)[0]
            if len(I) and len(J):
                with np.errstate(over='ignore'):
                    Lk.append((h0[I][:, None] + h1[J][None, :]).ravel())
                Li.append(np.stack(np.meshgrid(I, J, indexing='ij'), -1).reshape(-1, 2))
        if not Lk: continue
        Lk = np.concatenate(Lk); Li = np.concatenate(Li); o = np.argsort(Lk); ls = Lk[o]
        # right: blocks 2, 3 (row-0, canonical, block2 index <= block3 index), energies summing to ER
        for a in range(ER + 1):
            I = np.nonzero(e1 == a)[0]; J = np.nonzero(e1 == ER - a)[0]
            if not (len(I) and len(J)): continue
            for i2 in I:
                Jj = J[J >= i2] if canon else J
                if not len(Jj): continue
                with np.errstate(over='ignore'):
                    r = hT - h1[i2] - h1[Jj]
                npairs += len(Jj)
                lo = np.searchsorted(ls, r, 'left'); hi = np.searchsorted(ls, r, 'right')
                for jj in np.nonzero(hi > lo)[0]:
                    for t in range(lo[jj], hi[jj]):
                        i0, i1 = Li[o[t]]
                        if ((S0[i0] + S1[i1] + S1[i2] + S1[Jj[jj]]) == T).all():
                            sols.append([Y0[i0].tolist(), Y1[i1].tolist(), Y1[i2].tolist(), Y1[Jj[jj]].tolist()])
    out = dict(stage="done", solutions=len(sols), right_pairs=int(npairs), secs=round(time.time() - t0, 1))
    log(json.dumps(out))
    return out, sols

def planted(seed):
    """random Y with the slice's structure and energy 80 split over the blocks; target = its own lag sums"""
    rng = np.random.default_rng(seed)
    Y0 = block_list(2, 80)[0]; Y1 = block_list(0, 80)[0]
    while True:
        Ys = [Y0[rng.integers(len(Y0))]] + [Y1[rng.integers(len(Y1))] for _ in range(3)]
        if sum(int((y * y).sum()) for y in Ys) == 80: break
    X = [expand(x0, y)[0] for x0, y in zip((2, 0, 0, 0), Ys)]
    T = sum(np.array([(x * np.roll(x, -h)).sum() for h in LAGS]) for x in X)
    out, sols = solve(T=T, log=lambda s: None)
    def variants(y, neg): return {tuple(y), tuple(-y)} if neg else {tuple(y)}
    V = [variants(Ys[0], False)] + [variants(Ys[k], True) for k in (1, 2, 3)]
    hit = any(tuple(s[0]) in V[0] and any(all(tuple(s[1 + i]) in V[1 + p[i]] for i in range(3))
                                          for p in itertools.permutations(range(3))) for s in sols)
    return dict(seed=seed, solutions=len(sols), planted_found=hit, secs=out['secs'])

if __name__ == "__main__":
    if "--plant" in sys.argv:
        n = int(sys.argv[sys.argv.index("--plant") + 1])
        for s in range(n): print(json.dumps(planted(s)), flush=True)
        sys.exit()
    out, sols = solve(log=lambda s: print(s, flush=True))
    if sols: print(json.dumps(dict(first=sols[0])))

def cpsat(time_limit=600, workers=2, T=None, seed=0):
    """CP-SAT on the reduced model: Y_{k,j} in [-5, 5], sum_j Y_{k,j} = 0, sum Y^2 = 80, and for each lag class h:
    sum_k [X_k(0) Y_{k,cl(h)} + sum_{j,j'} N_h(j, j') Y_{k,j} Y_{k,j'}] = T_h / 4  (T_h = -40), N_h the cyclotomic
    numbers #{z in C_j : z + h in C_j'}; X_k(0) = (2, 0, 0, 0)"""
    from ortools.sat.python import cp_model
    T = np.full(len(LAGS), -40, np.int64) if T is None else np.asarray(T, np.int64)
    x0s = (2, 0, 0, 0)
    M = cp_model.CpModel()
    Y = [[M.NewIntVar(-5, 5, f"y{k}_{j}") for j in range(8)] for k in range(4)]
    for k in range(4): M.Add(sum(Y[k]) == 0)
    sq = {}
    def prod(k, j1, j2):
        key = (k, min(j1, j2), max(j1, j2))
        if key not in sq:
            p = M.NewIntVar(-25, 25, ""); M.AddMultiplicationEquality(p, [Y[k][key[1]], Y[k][key[2]]]); sq[key] = p
        return sq[key]
    M.Add(sum(prod(k, j, j) for k in range(4) for j in range(8)) == 80)
    for hi, h in enumerate(LAGS):
        N = np.zeros((8, 8), np.int64)
        for z in range(1, q):
            w = (z + h) % q
            if w != 0: N[CL[z], CL[w]] += 1
        terms = []
        for k in range(4):
            if x0s[k]: terms.append(x0s[k] * Y[k][CL[h]])
            for j1 in range(8):
                for j2 in range(8):
                    if N[j1, j2]: terms.append(int(N[j1, j2]) * prod(k, j1, j2))
        assert T[hi] % 4 == 0
        M.Add(sum(terms) == int(T[hi]) // 4)
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = time_limit; S.parameters.num_workers = workers
    S.parameters.random_seed = seed
    t0 = time.time(); st = S.StatusName(S.Solve(M))
    out = dict(model="971 x469 dual compression, reduced", status=st, secs=round(time.time() - t0, 1))
    if st in ("FEASIBLE", "OPTIMAL"):
        Yv = np.array([[S.Value(v) for v in row] for row in Y])
        X = [expand(x0, y)[0] for x0, y in zip(x0s, Yv)]
        P = sum(np.array([(x * np.roll(x, -hh)).sum() for hh in range(q)]) for x in X)
        out.update(Y=Yv.tolist(), verified=bool(P[0] == 3844 and (P[1:] == -40).all()) if T is None or (T == -40).all() else None)
    return out
