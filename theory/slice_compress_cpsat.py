r"""Exact compressed stage 1 of ANY multiplier slice (b = 1) by CP-SAT [M when FEASIBLE/INFEASIBLE].

Slice: G = Z_m, x_k constant on the orbits of M = <mu>.  Compress along the subgroup of order K (gcd(K, m/K) = 1):
y_k(b) = sum_{g = b mod n} x_k(g), n = m/K.  Exact facts [P]:
  * y_k is M-invariant on Z_n (y(mu b) = y(b));
  * the fibre over b is a disjoint union of the sets O cap fibre over the M-orbits O that meet it, each a single
    M_b-orbit (M_b = stabiliser of b) and they belong to distinct M-orbits, so y_k(b) ranges over the letters
    {sum_parts eps_part |part|} independently per M-orbit of Z_n;
  * sum_k PAF_{Z_n}(y_k)(h) = 4m - 4(K - 1) (h = 0), -4K (h != 0);  rows sum_b y_k(b) = (+-2, 0, 0, 0).
The slice is empty if this compressed problem is INFEASIBLE.  Pure generalised-Paley slices are the case M_b = M.
usage: python3 theory/slice_compress_cpsat.py m mu K [--time 600] [--workers 4] [--seed 0]"""
import sys, json, time, math
import numpy as np
from ortools.sat.python import cp_model

def structure(m, mu, K):
    n = m // K; assert m % K == 0 and math.gcd(K, n) == 1
    orb = -np.ones(m, np.int64); sizes = []
    for g in range(m):
        if orb[g] >= 0: continue
        x = g; i = len(sizes); c = 0
        while orb[x] < 0: orb[x] = i; c += 1; x = x * mu % m
        sizes.append(c)
    # orbits of M on Z_n and letters per orbit
    borb = -np.ones(n, np.int64); reps = []
    for b in range(n):
        if borb[b] >= 0: continue
        x = b; j = len(reps)
        while borb[x] < 0: borb[x] = j; x = x * mu % n
        reps.append(b)
    letters = []
    for b in reps:
        fib = [g for g in range(b, m, n)]
        parts = {}
        for g in fib: parts[orb[g]] = parts.get(orb[g], 0) + 1
        L = {0}
        for p in parts.values(): L = {a + e * p for a in L for e in (1, -1)}
        letters.append((sorted(L), sorted(parts.values())))
    return n, borb, reps, letters

def solve(m, mu, K, time_limit=600, workers=4, seed=0):
    n, borb, reps, letters = structure(m, mu, K)
    T = [4 * m - 4 * (K - 1)] + [-4 * K] * (n // 2)
    M = cp_model.CpModel()
    Y = [[M.NewIntVarFromDomain(cp_model.Domain.FromValues(letters[j][0]), f"y{k}_{j}") for j in range(len(reps))]
         for k in range(4)]
    y = [[Y[k][borb[b]] for b in range(n)] for k in range(4)]
    prods = {}
    def prod(k, a, b):
        ja, jb = int(borb[a]), int(borb[b]); key = (k,) + tuple(sorted((ja, jb)))
        if key not in prods:
            La, Lb = letters[ja][0], letters[jb][0]
            p = M.NewIntVarFromDomain(cp_model.Domain.FromValues(sorted({u * w for u in La for w in Lb})), "")
            M.AddMultiplicationEquality(p, [Y[k][ja], Y[k][jb]]); prods[key] = p
        return prods[key]
    for h in range(n // 2 + 1):
        M.Add(sum(prod(k, b, (b + h) % n) for k in range(4) for b in range(n)) == T[h])
    for k, r in enumerate((2, 0, 0, 0)): M.Add(sum(y[k]) == r)
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = time_limit; S.parameters.num_workers = workers
    S.parameters.random_seed = seed
    t0 = time.time(); st = S.StatusName(S.Solve(M))
    out = dict(m=m, v=m + 1, mu=mu, K=K, n=n, orbits_Zn=len(reps),
               letter_sets=sorted({tuple(l[0]) for l in letters}, key=len)[:3], status=st, secs=round(time.time() - t0, 1))
    if st in ("FEASIBLE", "OPTIMAL"):
        yy = [[S.Value(v) for v in row] for row in y]
        P = [sum(sum(yy[k][b] * yy[k][(b + h) % n] for b in range(n)) for k in range(4)) for h in range(n)]
        out.update(y=yy, verified=bool(P[0] == T[0] and all(p == -4 * K for p in P[1:])))
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--time": "600", "--workers": "4", "--seed": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    m, mu, K = pos
    print(json.dumps(solve(m, mu, K, float(opt["--time"]), int(opt["--workers"]), int(opt["--seed"]))), flush=True)
