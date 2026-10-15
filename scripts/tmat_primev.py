"""T(v) (w = 1) screen at prime open v: every odd k | v-1 with >= 4 orbits; exact row-sum condition, then Gaussian.
Usage: tmat_primev.py v1,v2,... [N]          odd-k K-invariant slices
       tmat_primev.py --asym t1,t2,... [N]   no symmetry (k = 1): Gaussian only (any row sums are admissible; the
                                             row-sum check of tmat_estimate assumes k >= 3)"""
import os, sys, json, math, itertools, resource
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tmat_estimate import estimate

def rowsum_ok(t, k):
    r = math.isqrt(t)
    return any(sum(x * x for x in q) == t and sorted(x % k for x in q).count(0) == 3 and any(x % k in (1, k - 1) for x in q)
               for q in itertools.product(range(-r, r + 1), repeat=4))

def asym(t, N):
    import numpy as np
    from tmat_estimate import sample_f
    from tmat_search import subgroup, setup
    orbits, oid, M, w = setup(t, subgroup(t, 1)); m = len(orbits)
    F = sample_f(t, orbits, N, np.random.default_rng(1)); C = F.shape[1]
    mu = F.mean(0); S = np.cov(F.T); d = F - F[0]
    steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, c])))) for c in range(C)])
    sign, logdet = np.linalg.slogdet(S); q = mu @ np.linalg.solve(S, mu)
    log2P = (np.log(steps).sum() - 0.5 * C * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q) / math.log(2)
    return {"t": t, "k": 1, "orbits": m, "classes": C, "log2E": round(3 * m + log2P, 1), "n_bits": round(-log2P, 1), "N": N}

if __name__ == "__main__":
    if sys.argv[1] == "--asym":
        N = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
        for t in [int(x) for x in sys.argv[2].split(",")]:
            print(json.dumps(asym(t, N)), flush=True)
        sys.exit()
    vs = [int(x) for x in sys.argv[1].split(",")]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
    for v in vs:
        for k in range(3, v, 2):
            if (v - 1) % k or (v - 1) // k + 1 < 4: continue
            m = (v - 1) // k + 1
            if not rowsum_ok(v, k):
                print(json.dumps({"t": v, "k": k, "orbits": m, "rowsum_feasible": False}), flush=True); continue
            rec = estimate(v, k, N)
            rec["n_bits"] = round(3 * m - rec["log2E"], 1)
            print(json.dumps(rec), flush=True)
