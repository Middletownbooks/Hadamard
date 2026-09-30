"""Stage-C exactness: summed over all margin pairs, lift counts must equal 2^((n+1)/2); also each lift has the requested margins."""
import sys, itertools, json, numpy as np
import wil_search as ws
n, p, q = map(int, sys.argv[1:4]); hp, hq = (p + 1) // 2, (q + 1) // 2
def halves(m, k):
    H = np.array(list(itertools.product(range(-k, k + 1, 2), repeat=(m + 1) // 2)))
    return H, H[:, 0] + 2 * H[:, 1:].sum(1)
Rs, SR = halves(p, q); Cs, SC = halves(q, p)
total = 0; margins_ok = True; i = np.arange(n)
for R, sr in zip(Rs, SR):
    for C in Cs[SC == sr]:
        X = ws.lifts(n, p, q, R, C); total += len(X)
        if len(X):
            rr = np.stack([X[:, i % p == a].sum(1) for a in range(hp)], 1)
            cc = np.stack([X[:, i % q == b].sum(1) for b in range(hq)], 1)
            margins_ok &= bool(np.all(rr == R) and np.all(cc == C))
            margins_ok &= len(np.unique(X, axis=0)) == len(X)
print(json.dumps({"n": n, "lifts_total": total, "expected": 2 ** ((n + 1) // 2), "margins_ok": margins_ok, "pass": total == 2 ** ((n + 1) // 2) and margins_ok}))
