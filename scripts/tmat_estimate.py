"""Gaussian estimate of the number of K-invariant T-matrices of prime order t (K <= (Z/t)* of order k).

Random labelling: each K-orbit of Z_t gets a uniform label in 0..7 (row = l & 3, sign = l & 4).
f_c = (sum of the four PAFs) at a representative lag of each nonzero orbit class c; T-matrices <=> f = 0.
log2 E = m*3 + log2 P(f = 0), P from a multivariate normal with Monte-Carlo mean/covariance times the lattice
cell volume. The lattice is estimated as the per-coordinate gcd of sampled differences, plus joint mod-2^j checks.
Reported per symmetry orbit: log2(E / (4! * 2^4 * (t-1)/k)).
"""
import sys, json, math
import numpy as np
from tmat_search import subgroup, setup

def sample_f(t, orbits, N, rng, batch=20000):
    m = len(orbits)
    reps, seen = [], set()
    for o in orbits[1:]:                                 # PAF(s) = PAF(-s): one rep per class pair {c, -c}
        if o[0] in seen: continue
        reps.append(o[0]); seen.update(o); seen.update((-x) % t for x in o)
    oid = np.zeros(t, dtype=np.int64)
    for i, o in enumerate(orbits): oid[o] = i
    out = []
    for s in range(0, N, batch):
        n = min(batch, N - s)
        L = rng.integers(0, 8, (n, m))
        Lx = L[:, oid]                                   # n x t labels
        sgn = np.where(Lx & 4, -1.0, 1.0)
        tot = np.zeros((n, t))
        for r in range(4):
            x = np.where((Lx & 3) == r, sgn, 0.0)
            X = np.fft.rfft(x, axis=1)
            tot += np.fft.irfft(X * X.conj(), n=t, axis=1)
        out.append(np.rint(tot[:, reps]).astype(np.int64))
    return np.vstack(out)

def estimate(t, k, N=200000, seed=1):
    K = subgroup(t, k)
    orbits, oid, M, w = setup(t, K)
    m = len(orbits)
    F = sample_f(t, orbits, N, np.random.default_rng(seed)); C = F.shape[1]
    mu = F.mean(0); S = np.cov(F.T)
    # lattice: per-coordinate step = gcd of deviations from first sample; parity of f_c
    d = F - F[0]
    steps = np.array([np.gcd.reduce(np.abs(d[:, c])) for c in range(C)])
    steps[steps == 0] = 1
    # check 0 is on the lattice coset
    on = all((0 - F[0, c]) % steps[c] == 0 for c in range(C))
    sign, logdet = np.linalg.slogdet(S)
    q = mu @ np.linalg.solve(S, mu)
    log2P = (np.sum(np.log(steps)) - 0.5 * C * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q) / math.log(2)
    log2E = 3 * m + log2P
    sym = math.log2(24 * 16 * (t - 1) / k)
    # row sums: each block's sum is +-1 + k*j with j = (m-1) mod 2, |j| <= m-1; need a^2+b^2+c^2+d^2 = 4t
    ach = {e + k * j for e in (1, -1) for j in range(-(m - 1), m) if (j - (m - 1)) % 2 == 0}
    from tmat_sat import rowsum_cases
    feas = [c for c in rowsum_cases(t) if all(x in ach or -x in ach for x in c)]
    if not feas:
        log2E = float("-inf")
    return {"rowsum_feasible": bool(feas), "t": t, "k": k, "orbits": m, "classes": C, "steps": sorted(set(steps.tolist())), "zero_on_lattice": on,
            "log2E": round(log2E, 1), "log2E_per_symclass": round(log2E - sym, 1), "N": N}

if __name__ == "__main__":
    t = int(sys.argv[1]); ks = [int(x) for x in sys.argv[2].split(",")]
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
    for k in ks:
        print(json.dumps(estimate(t, k, N)), flush=True)
