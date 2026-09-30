"""Independent brute force for Williamson quadruples (W0). No compression, no multipliers, no product condition.
Counts normalised quadruples: row sums s_A>=s_B>=s_C>=s_D>0. Usage: wil_brute.py n"""
import sys, json, itertools, numpy as np
from wil_common import expand, paf, verify, rowsum_cases, epsilon

n = int(sys.argv[1]); hl = (n + 1) // 2
H = np.array(list(itertools.product([1, -1], repeat=hl)), dtype=np.int64)
X = expand(H, n); S = X.sum(1); P = paf(X)[:, 1:hl]            # lags 1..(n-1)/2
by = {s: np.where(S == s)[0] for s in set(S.tolist())}
total = 0; sols = []; prod_ok = True
for case in rowsum_cases(n):
    sa, sb, sc, sd = case
    A, B, C, D = (by.get(s, np.array([], int)) for s in case)
    tab = {}
    for i in C:
        for j in D:
            tab.setdefault((-(P[i] + P[j])).tobytes(), []).append((i, j))
    for i in A:
        for j in B:
            for (k, l) in tab.get((P[i] + P[j]).tobytes(), []):
                Q = X[[i, j, k, l]]; assert verify(Q)
                total += 1; sols.append(Q[:, :hl].tolist())
                pr = Q[0, 1:] * Q[1, 1:] * Q[2, 1:] * Q[3, 1:]
                prod_ok &= bool(np.all(pr == epsilon(n, case)))
print(json.dumps({"n": n, "cases": rowsum_cases(n), "N": total, "product_condition_holds": prod_ok}))
json.dump({"n": n, "N": total, "solutions": sols}, open(f"logs/wil_brute_n{n}.json", "w"))
