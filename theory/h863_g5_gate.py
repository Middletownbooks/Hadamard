"""[M if decided] Exact gate: does the codex Z43-compressed solution y (H = Z2^2 x Z5) admit a compatible G/Z5 level
(Z2^2 x Z43) compression?  Independent model (no repo code).  Variables: class signs c[i,q,u,j] (j = 6 classes of
<4> in Z43*), e fixed by y mod 14.  w[i,q,j] = sum_u c, w0[i,q] = sum_u e.  Equations: summed PAF on Z2^2 x Z43
(targets 3392 at 0, -60 at (0, t != 0), 20 at q != 0).  2026-11-03: UNKNOWN after 1200 s on 4 workers (codex certificate).
usage: python3 theory/h863_g5_gate.py CERT.json [secs]"""
import json, sys, itertools as it, time
from ortools.sat.python import cp_model
d = json.load(open(sys.argv[1])); secs = float(sys.argv[2]) if len(sys.argv) > 2 else 600
Y = d['y']
def idx(q1, q2, u): return (q1 * 2 + q2) * 5 + u      # assumed mixed-radix order, checked below
def hpaf(y, a, b, s):
    return sum(y[idx(q1, q2, u)] * y[idx(q1 ^ a, q2 ^ b, (u + s) % 5)] for q1 in (0, 1) for q2 in (0, 1) for u in range(5))
for a, b, s in it.product((0, 1), (0, 1), range(5)):
    v = sum(hpaf(y, a, b, s) for y in Y)
    tgt = 2936 if (a, b, s) == (0, 0, 0) else (-516 if (a, b) == (0, 0) else 172)
    assert v == tgt, (a, b, s, v, tgt)
p = 43; H = [pow(4, k, p) for k in range(7)]
cls = {}; reps = []
for t in range(1, p):
    if t not in cls:
        j = len(reps); reps.append(t)
        for h in H: cls[t * h % p] = j
assert len(reps) == 6
E = [[1 if y % 14 == 1 else -1 for y in row] for row in Y]
assert all((y - e) % 14 == 0 for row, er in zip(Y, E) for y, e in zip(row, er))
M = cp_model.CpModel()
c = {}; w = {}
for i in range(4):
    for q in range(4):
        for u in range(5):
            k = q * 5 + u; s = (Y[i][k] - E[i][k]) // 7
            vs = [M.NewIntVarFromDomain(cp_model.Domain.FromValues([-1, 1]), '') for j in range(6)]
            for j in range(6): c[i, q, u, j] = vs[j]
            M.Add(sum(vs) == s)
        for j in range(6):
            w[i, q, j] = M.NewIntVar(-5, 5, ''); M.Add(w[i, q, j] == sum(c[i, q, u, j] for u in range(5)))
def val(i, q, t):        # G/Z5 sequence value at (q, t): int or var
    return sum(E[i][q * 5 + u] for u in range(5)) if t == 0 else w[i, q, cls[t]]
prods = {}
eqs = 0
for qa in range(4):
    for t0 in [0] + reps:
        terms = []
        for i in range(4):
            for q in range(4):
                for t in range(p):
                    A = (i, q, t); B = (i, q ^ qa, (t + t0) % p)
                    terms.append((A, B))
        # build linear expression with product vars
        lin = []; const = 0
        for A, B in terms:
            va, vb = val(*A), val(*B)
            ia, ib = isinstance(va, int), isinstance(vb, int)
            if ia and ib: const += va * vb
            elif ia: lin.append(va * vb)
            elif ib: lin.append(vb * va)
            else:
                key = tuple(sorted([(A[0], A[1], cls[A[2]]), (B[0], B[1], cls[B[2]])]))
                if key not in prods:
                    pv = M.NewIntVar(-25, 25, ''); M.AddMultiplicationEquality(pv, [w[key[0]], w[key[1]]]); prods[key] = pv
                lin.append(prods[key])
        tgt = 3392 if (qa, t0) == (0, 0) else (-60 if qa == 0 else 20)
        M.Add(sum(lin) + const == tgt); eqs += 1
S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = secs; S.parameters.num_workers = 4
t = time.time(); st = S.StatusName(S.Solve(M))
out = dict(status=st, secs=round(time.time() - t, 1), eqs=eqs, prods=len(prods))
if st in ('OPTIMAL', 'FEASIBLE'):
    W = [[[sum(E[i][q * 5 + u] for u in range(5))] + [S.Value(w[i, q, j]) for j in range(6)] for q in range(4)] for i in range(4)]
    # independent recheck on explicit length-172 sequences
    seq = [[W[i][q][0 if t == 0 else cls[t] + 1] for q in range(4) for t in range(p)] for i in range(4)]
    ok = True
    for qa in range(4):
        for t0 in range(p):
            v = sum(seq[i][q * p + t] * seq[i][(q ^ qa) * p + (t + t0) % p] for i in range(4) for q in range(4) for t in range(p))
            tgt = 3392 if (qa, t0) == (0, 0) else (-60 if qa == 0 else 20)
            ok &= v == tgt
    out.update(W=W, recheck=ok)
print(json.dumps(out))
