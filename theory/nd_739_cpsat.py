r"""Independent CP-SAT check of the 739 x215 stage 1 (symmetric y on Z_18, letters +-1 + 10j, |j| <= 4) [M].

Model (no use of theory/nd_sym_stage1.py): y_k(g) = alpha_k(g) + 10 t_k(g) on the 10 negation orbits of Z_18,
alpha in {-1, 1}, t in [-4, 4];
    sum_k PAF(y_k)(0) = 2792, sum_k PAF(y_k)(h) = -164 (h = 1..9), rows (2, 0, 0, 0) (block 0 carries +2).
Redundant constraints (all implied, see report section 2):
    sum_k PAF(alpha_k)(h) = 16 + 20 z_h (h = 1..9)   [alpha-congruence mod 20]
    S and D levels on Z_9: sum_k PAF_9(sigma_k) = (657, -82 x 4), sum_k PAF_9(delta_k) = (739, 0 x 4).
Symmetry breaking: none beyond block 0 having row +2 (so INFEASIBLE is a statement about the whole stage 1).
usage: python3 theory/nd_739_cpsat.py [--time 3600] [--workers 4] [--redundant 1]"""
import sys, json, time
from ortools.sat.python import cp_model

def build(redundant=True):
    F, n = 18, 9
    M = cp_model.CpModel()
    reps = list(range(10))                                  # orbit representatives 0..9 of negation on Z_18
    orb = [min(g, (-g) % F) for g in range(F)]
    A = [[M.NewBoolVar(f"a{k}_{r}") for r in reps] for k in range(4)]          # alpha = 2a - 1
    T = [[M.NewIntVar(-4, 4, f"t{k}_{r}") for r in reps] for k in range(4)]
    Y = [[M.NewIntVar(-41, 41, f"y{k}_{r}") for r in reps] for k in range(4)]
    for k in range(4):
        for r in reps: M.Add(Y[k][r] == 2 * A[k][r] - 1 + 10 * T[k][r])
    y = [[Y[k][orb[g]] for g in range(F)] for k in range(4)]
    prods = {}
    def prod(k, i, j):
        key = (k,) + tuple(sorted((i, j)))
        if key not in prods:
            p = M.NewIntVar(-1681, 1681, ""); M.AddMultiplicationEquality(p, [Y[k][key[1]], Y[k][key[2]]]); prods[key] = p
        return prods[key]
    for h in range(10):
        tgt = 2792 if h == 0 else -164
        M.Add(sum(prod(k, orb[g], orb[(g + h) % F]) for k in range(4) for g in range(F)) == tgt)
    for k, r in enumerate((2, 0, 0, 0)): M.Add(sum(y[k]) == r)
    if redundant:
        # alpha-congruence: alpha products via booleans: alpha_i alpha_j = 1 - 2 (a_i xor a_j)
        xors = {}
        def axor(k, i, j):
            key = (k,) + tuple(sorted((i, j)))
            if key[1] == key[2]: return None
            if key not in xors:
                x = M.NewBoolVar(""); M.AddBoolXOr([A[k][key[1]], A[k][key[2]], x.Not()]); xors[key] = x
            return xors[key]
        for h in range(1, 10):
            terms = []; const = 0
            for k in range(4):
                for g in range(F):
                    x = axor(k, orb[g], orb[(g + h) % F])
                    if x is None: const += 1
                    else: const += 1; terms.append(x)
            z = M.NewIntVar(-10, 10, f"z{h}")
            M.Add(const - 2 * sum(terms) == 16 + 20 * z)
        # S/D levels on Z_9: sigma(b) = (y(b') + y(b' + 9))/2, delta(b) = (-1)^b' (y(b') - y(b' + 9))/2, b' = 0..8
        S = [[M.NewIntVar(-41, 41, "") for b in range(n)] for k in range(4)]
        D = [[M.NewIntVar(-41, 41, "") for b in range(n)] for k in range(4)]
        for k in range(4):
            for b in range(n):
                M.Add(2 * S[k][b] == y[k][b] + y[k][b + n])
                M.Add(2 * D[k][b] == (1 if b % 2 == 0 else -1) * (y[k][b] - y[k][b + n]))
        def lvl(V, T0, Th):
            pp = {}
            for h in range(5):
                acc = []
                for k in range(4):
                    for b in range(n):
                        i, j = b, (b + h) % n
                        key = (k, min(i, j), max(i, j))
                        if key not in pp:
                            p = M.NewIntVar(-1681, 1681, ""); M.AddMultiplicationEquality(p, [V[k][i], V[k][j]]); pp[key] = p
                        acc.append(pp[key])
                M.Add(sum(acc) == (T0 if h == 0 else Th))
        lvl(S, 657, -82); lvl(D, 739, 0)
    return M

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--time": "3600", "--workers": "4", "--redundant": "1", "--seed": "0"}; i = 0
    while i < len(a): opt[a[i]] = a[i + 1]; i += 2
    M = build(bool(int(opt["--redundant"])))
    S = cp_model.CpSolver(); S.parameters.max_time_in_seconds = float(opt["--time"])
    S.parameters.num_workers = int(opt["--workers"]); S.parameters.random_seed = int(opt["--seed"])
    t0 = time.time(); st = S.StatusName(S.Solve(M))
    print(json.dumps(dict(v=739, slice="x215", model="y + alpha/t + S/D redundant" if int(opt["--redundant"]) else "y only",
                          status=st, secs=round(time.time() - t0, 1), workers=int(opt["--workers"]))), flush=True)
