"""EXACT stage 2 by CP-SAT for a generalised-Paley slice (b = 1): G = Z_F x Z_q, Z_F fixed pointwise, Z_q split into
s cyclotomic classes C_0..C_{s-1} (cosets of the index-s subgroup; C_j = g0^j H).
Given stage-1 coset sums y (F-vectors per block), find per coset (k, g): alpha = +-1 at z = 0 and c_j = +-1 on C_j with
    alpha + f * sum_j c_j = y[k][g]      (f = (q - 1)/s)
such that the full core has sum_k PAF = -4 at every nonzero lag of G.  Exact reduction: for s in class t and lag h,
    R(h, t) = sum_k sum_g [ alpha_g c_{g+h, t} + c_{g, t'} alpha_{g+h} + sum_{i,j} N_t(i, j) c_{g, i} c_{g+h, j} ] = -4,
t' = class of -s, N_t(i, j) = #{z in C_i : z + s in C_j} (s in C_t; cyclotomic numbers).  The lags (h, 0) then follow from
stage 1.  +-1 products are encoded by equality literals.  The result is expanded to G, checked by FFT on Z_F x Z_q, and
written as a cyclic Z_{Fq} core (CRT) for completion.
usage: cyclo_stage2_cpsat.py F q s y.json [--time 1800] [--workers 4] [--out core.txt]"""
import sys, json, time, itertools as it
import numpy as np
from ortools.sat.python import cp_model

def classes(q, s):
    g0 = next(g for g in range(2, q) if all(pow(g, (q - 1) // r, q) != 1 for r in range(2, q) if (q - 1) % r == 0 and all(r % d for d in range(2, r))))
    cls = np.full(q, -1, np.int64)
    x = 1
    for e in range(q - 1):
        cls[x] = e % s; x = x * g0 % q
    return cls, g0

def cyclo_numbers(q, s, cls):
    """N[t][i][j] = #{z in C_i: z + s_t in C_j}, s_t = a representative of C_t"""
    reps = [next(z for z in range(1, q) if cls[z] == t) for t in range(s)]
    N = np.zeros((s, s, s), np.int64)
    for t, st in enumerate(reps):
        for z in range(1, q):
            w = (z + st) % q
            if w: N[t, cls[z], cls[w]] += 1
    return N, reps

def solve(F, q, s, Y, time_limit=1800, workers=4, seed=0):
    f = (q - 1) // s; cls, g0 = classes(q, s); N, reps = cyclo_numbers(q, s, cls)
    tneg = [int(cls[(q - reps[t]) % q]) for t in range(s)]          # class of -s for s in C_t
    M = cp_model.CpModel()
    # sign bits: alpha b[k][g][0], c_j b[k][g][1 + j];  value = 2b - 1
    b = [[[M.NewBoolVar(f"b{k}_{g}_{j}") for j in range(s + 1)] for g in range(F)] for k in range(4)]
    val = lambda B: 2 * B - 1
    for k in range(4):
        for g in range(F):
            M.Add(val(b[k][g][0]) + f * sum(val(b[k][g][1 + j]) for j in range(s)) == int(Y[k][g]))
    eqcache = {}
    def prod(B1, B2):
        """literal e with e <=> (B1 == B2); product of the +-1 values = 2e - 1"""
        key = (B1.Index(), B2.Index()) if B1.Index() <= B2.Index() else (B2.Index(), B1.Index())
        if key in eqcache: return eqcache[key]
        if B1.Index() == B2.Index():
            e = M.NewConstant(1)
        else:
            e = M.NewBoolVar("")
            M.AddBoolOr([B1, B2, e]); M.AddBoolOr([B1.Not(), B2.Not(), e])
            M.AddBoolOr([B1, B2.Not(), e.Not()]); M.AddBoolOr([B1.Not(), B2, e.Not()])
        eqcache[key] = e; return e
    n_eq = 0
    for h in range(F):
        for t in range(s):
            # use symmetry R(-h, t) = R(h, t'): keep (h, t) with h < F - h, or h in {0, F/2} with t <= t'
            hn = (-h) % F
            if hn < h: continue
            if hn == h and tneg[t] < t: continue
            terms = []
            for k in range(4):
                for g in range(F):
                    g2 = (g + h) % F
                    terms.append(2 * prod(b[k][g][0], b[k][g2][1 + t]) - 1)
                    terms.append(2 * prod(b[k][g][1 + tneg[t]], b[k][g2][0]) - 1)
                    for i in range(s):
                        for j in range(s):
                            if N[t, i, j]:
                                terms.append(int(N[t, i, j]) * (2 * prod(b[k][g][1 + i], b[k][g2][1 + j]) - 1))
            M.Add(sum(terms) == -4); n_eq += 1
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = time_limit; S.parameters.num_workers = workers
    S.parameters.random_seed = seed
    t0 = time.time(); st = S.Solve(M); name = S.StatusName(st)
    out = dict(F=F, q=q, s=s, g0=int(g0), equations=n_eq, status=name, secs=round(time.time() - t0, 1))
    if name in ("FEASIBLE", "OPTIMAL"):
        X = np.zeros((4, F, q), np.int64)
        for k in range(4):
            for g in range(F):
                X[k, g, 0] = 2 * S.Value(b[k][g][0]) - 1
                for z in range(1, q): X[k, g, z] = 2 * S.Value(b[k][g][1 + cls[z]]) - 1
        P = sum(np.rint(np.fft.ifft2(np.abs(np.fft.fft2(X[k])) ** 2).real).astype(np.int64) for k in range(4))
        ok = bool(P[0, 0] == 4 * F * q and (P.ravel()[1:] == -4).all())
        assert (X.reshape(4, F, q).sum(2) == np.array(Y)).all()
        out.update(core_ok=ok, X=X.reshape(4, -1).tolist())
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--time": "1800", "--workers": "4", "--seed": "0", "--out": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3]))
    Y = Y['y'] if isinstance(Y, dict) else Y
    r = solve(F, q, s, Y, float(opt["--time"]), int(opt["--workers"]), int(opt["--seed"]))
    X = r.pop('X', None)
    print(json.dumps(r), flush=True)
    if X is not None and r.get('core_ok') and opt["--out"]:
        m = F * q; idx = np.arange(m); src = (idx % F) * q + (idx % q)            # cyclic Z_m by CRT
        Xc = np.array(X)[:, src]
        with open(opt["--out"], 'w') as fh:
            fh.write(f"v={m + 1} N={4 * (m + 1)} b=1 m={m} cyclic Z{m} (CRT of Z{F} x Z{q}); generalised Paley slice s={s}; scripts/cyclo_stage2_cpsat.py\n")
            for nm, row in zip('ABCD', Xc): fh.write(nm + '=' + ''.join('+' if t > 0 else '-' for t in row) + '\n')
