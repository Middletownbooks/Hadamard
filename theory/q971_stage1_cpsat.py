"""Independent check of theory/q971_stage1.py: CP-SAT on y = alpha + 24 j over Z_10 (971 x81 slice, composition
{23: 3, 25: 3, 1: 34}), with no alpha reduction, no symmetry breaking and no joins: alpha_k(g) = +-1, j_k(g) in {-1, 0, 1},
sum |j| = 6, #(alpha j = +1) = 3, rows (2, 0, 0, 0), sum_k PAF(y_k)(h) = -388 (h = 1..5).  INFEASIBLE is a proof.
usage: python3 theory/q971_stage1_cpsat.py [time] [workers]"""
import sys, json, time
from ortools.sat.python import cp_model
F, q = 10, 97
tl = float(sys.argv[1]) if len(sys.argv) > 1 else 1800; w = int(sys.argv[2]) if len(sys.argv) > 2 else 4
M = cp_model.CpModel()
y = [[M.NewIntVarFromDomain(cp_model.Domain.FromValues([-25, -23, -1, 1, 23, 25]), f"y{k}{g}") for g in range(F)]
     for k in range(4)]
big = [[M.NewBoolVar("") for g in range(F)] for k in range(4)]           # |y| > 1
p25 = [[M.NewBoolVar("") for g in range(F)] for k in range(4)]           # |y| = 25
for k in range(4):
    for g in range(F):
        M.AddAllowedAssignments([y[k][g], big[k][g], p25[k][g]],
                                [(-25, 1, 1), (-23, 1, 0), (-1, 0, 0), (1, 0, 0), (23, 1, 0), (25, 1, 1)])
M.Add(sum(big[k][g] for k in range(4) for g in range(F)) == 6)
M.Add(sum(p25[k][g] for k in range(4) for g in range(F)) == 3)
for k, r in enumerate((2, 0, 0, 0)): M.Add(sum(y[k]) == r)
dom = cp_model.Domain.FromValues(sorted({a * b for a in (-25, -23, -1, 1, 23, 25) for b in (-25, -23, -1, 1, 23, 25)}))
for h in range(1, F // 2 + 1):
    terms = []
    for k in range(4):
        for g in range(F):
            p = M.NewIntVarFromDomain(dom, ""); M.AddMultiplicationEquality(p, [y[k][g], y[k][(g + h) % F]]); terms.append(p)
    M.Add(sum(terms) == -4 * q)
S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = tl; S.parameters.num_workers = w
t0 = time.time(); st = S.StatusName(S.Solve(M))
print(json.dumps(dict(v=971, slice="x81", status=st, secs=round(time.time() - t0, 1))))
