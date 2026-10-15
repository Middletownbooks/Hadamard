"""Exact counts for the t = 2k^2 + 1 family (row sums (1, k, k, 0)), to calibrate tmat_estimate on the shape of 883/k=21.
t = 19, k = 3: exhaustive over all 8^7 K-invariant labellings.
t = 163, k = 9: complete 4-list joins (all solutions) on random balanced supports S, as in tmat_kt.py; labelled
solutions per S = 16 * (normalised count); predicted = E / 2^m."""
import os, sys, json, math, time, resource
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tmat_kt import orbit_setup, side_list, tm_check
from tmat_estimate import estimate

def exhaustive(t, k):
    K, orbits, reps = orbit_setup(t, k); m = len(orbits)
    oid = np.zeros(t, np.int64)
    for i, o in enumerate(orbits): oid[o] = i
    n = 8 ** m; cnt = 0
    for s in range(0, n, 1 << 18):
        idx = np.arange(s, min(n, s + (1 << 18)), dtype=np.int64)
        L = np.stack([(idx >> (3 * j)) & 7 for j in range(m)], 1)[:, oid]
        tot = np.zeros((len(idx), t))
        for r in range(4):
            x = np.where((L & 3) == r, np.where(L & 4, -1.0, 1.0), 0.0)
            X = np.fft.rfft(x, axis=1); tot += np.fft.irfft(np.abs(X) ** 2, n=t, axis=1)
        tot = np.rint(tot).astype(np.int64)
        cnt += int(((tot[:, 0] == t) & (tot[:, 1:] == 0).all(1)).sum())
    return m, cnt

def pair_hash(P, r):
    H = (P.astype(np.int64) @ r)                          # random linear hash of each PAF vector
    return (H[:, None] + H[None, :]).ravel()

def census(t, k, nS, seed):
    K, orbits, reps = orbit_setup(t, k); m = len(orbits)
    rng = np.random.default_rng(seed); r = rng.integers(1, 2 ** 40, len(reps)).astype(np.int64)
    tot, rows, t0 = 0, [], time.time()
    for i in range(nS):
        S = sorted(rng.choice(m, m // 2 + int(rng.integers(0, 2)) * (m % 2), replace=False).tolist())
        Sc = [o for o in range(m) if o not in S]
        sgS, PS = side_list(t, orbits, S, reps); sgC, PC = side_list(t, orbits, Sc, reps)
        hA, hB = pair_hash(PS, r), pair_hash(PC, r)
        uA, cA = np.unique(hA, return_counts=True); uB, cB = np.unique(-hB, return_counts=True)
        common, ia, ib = np.intersect1d(uA, uB, return_indices=True)
        c = int((cA[ia] * cB[ib]).sum())
        if c:                                             # verify every hashed match exactly
            good = 0
            for h in common:
                for a in np.flatnonzero(hA == h):
                    for b in np.flatnonzero(-hB == h):
                        i1, i2 = divmod(int(a), len(PS)); i3, i4 = divmod(int(b), len(PC))
                        T = np.zeros((4, t), np.int64)
                        for j, o in enumerate(S):
                            T[0, orbits[o]] = (int(sgS[i1, j]) + int(sgS[i2, j])) // 2; T[1, orbits[o]] = (int(sgS[i1, j]) - int(sgS[i2, j])) // 2
                        for j, o in enumerate(Sc):
                            T[2, orbits[o]] = (int(sgC[i3, j]) + int(sgC[i4, j])) // 2; T[3, orbits[o]] = (int(sgC[i3, j]) - int(sgC[i4, j])) // 2
                        good += tm_check(T)
            assert good == c, (good, c)
        tot += c; rows.append(c)
    return m, tot, rows, time.time() - t0

out = {}
m, c = exhaustive(19, 3); e = estimate(19, 3, 200000)
out["19/3"] = {"orbits": m, "exact_labelled": c, "exact_log2": round(math.log2(c), 2) if c else None, "est_log2E": e["log2E"]}
print(json.dumps(out), flush=True)
nS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
e = estimate(163, 9, 200000)
m, tot, rows, secs = census(163, 9, nS, 7)
lab = 16 * tot / nS
print(json.dumps({"163/9": {"orbits": m, "supports": nS, "normalised_solutions": tot, "per_support": round(tot / nS, 4),
                  "predicted_per_support": round(2 ** e["log2E"] / 2 ** (m + 4), 4),
                  "census_log2E": round(math.log2(lab * 2 ** m), 2) if tot else None, "est_log2E": e["log2E"],
                  "supports_with_solutions": sum(1 for x in rows if x), "max_per_support": max(rows), "secs": round(secs)}}), flush=True)
