"""Exact counts for the 3-cycle-twisted T-matrix ansatz (tmat_twist.py, g = 3-cycle on three rows, fourth row fixed),
to calibrate its Gaussian on the shape of 919/d = 27.

Per support S (the nonzero U-orbits whose base label lies in the three cycled rows): list A = all base labels on S
(6 per orbit, first orbit's label fixed: the 3 row rotations x global sign of the block act freely), list B = all signs
of the fixed row on S^c (position 0 carries +, global sign of that row removed). Exact 2-list join on the class sums,
every match re-checked. Labelled solutions on S = 12 * matches.
  all supports:      tmat_twist3_census.py t d all
  importance sample: tmat_twist3_census.py t d sample NS p seed   (each orbit in S with probability p; unbiased E)"""
import os, sys, json, math, time, itertools, resource
resource.setrlimit(resource.RLIMIT_AS, (8 << 30, 8 << 30))
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tmat_twist as tw

G3 = [g for ct, g in tw.class_reps().items() if ct == ((1, 1), (3, 1))][0]
ROWS, SGN = tw.ROWS, tw.SGN

def class_sums(rowsX, reps):
    """rowsX: n x 4 x t (float32) -> n x len(reps) PAF sums."""
    F = np.fft.rfft(rowsX, axis=2)
    return np.rint(np.fft.irfft((np.abs(F) ** 2).sum(1), n=rowsX.shape[2], axis=1)[:, reps]).astype(np.int64)

def count_S(t, orbits, G, reps, S, Sc, frow, cyc_labels, fix_label, rng_hash, chunk=1 << 15):
    pos = np.array(orbits)
    # list B: fixed row on S^c plus position 0 (+)
    nB = 2 ** len(Sc); B = np.empty((nB, len(reps)), np.int64)
    for s0 in range(0, nB, chunk):
        idx = np.arange(s0, min(nB, s0 + chunk))
        X = np.zeros((len(idx), 4, t), np.float32); X[:, frow, 0] = 1
        for j, o in enumerate(Sc):
            sg = np.where((idx >> j) & 1, -1.0, 1.0).astype(np.float32)
            X[:, frow, pos[o]] = sg[:, None]
        B[s0:s0 + len(idx)] = class_sums(X, reps)
    hB = B @ rng_hash
    # list A: base labels on S (first fixed)
    nA = 6 ** (len(S) - 1) if S else 1; total = 0; matches = []
    uB, cB = np.unique(-hB, return_counts=True)
    for s0 in range(0, nA, chunk):
        idx = np.arange(s0, min(nA, s0 + chunk))
        X = np.zeros((len(idx), 4, t), np.float32)
        for j, o in enumerate(S):
            if j == 0: lab = np.full(len(idx), fix_label)
            else: lab = cyc_labels[(idx // 6 ** (j - 1)) % 6]
            for q in range(G.shape[0]):
                L = G[q][lab]; X[np.arange(len(idx))[:, None], ROWS[L][:, None], pos[o][q]] = SGN[L][:, None].astype(np.float32)
        A = class_sums(X, reps); hA = A @ rng_hash
        k = np.searchsorted(uB, hA); k[k >= len(uB)] = 0
        hit = uB[k] == hA
        total += int(cB[k[hit]].sum())
        for a in np.flatnonzero(hit):                         # exact re-check of each match
            for b in np.flatnonzero(-hB == hA[a]):
                assert ((A[a] + B[b]) == 0).all()
                matches.append((int(idx[a]), int(b)))
    return total, matches

def setup3(t, d):
    u, orbits, G, fix0, reps = tw.setup(t, d, G3)
    mg = tw.LMAP[G3]
    frow = [r for r in range(4) if mg[r] == r and mg[r | 4] == (r | 4)][0]
    cyc = np.array([l for l in range(8) if ROWS[l] != frow]); fix_label = int(cyc[0])
    return orbits, G, reps, frow, cyc, fix_label

if __name__ == "__main__":
    t, d, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    orbits, G, reps, frow, cyc, fix_label = setup3(t, d); m = len(orbits)
    est = tw.estimate(t, d, G3, 100000)
    rngh = np.random.default_rng(12345).integers(1, 2 ** 40, len(reps)).astype(np.int64)
    t0 = time.time()
    if mode == "all":
        E = 0
        for mask in range(2 ** m):
            S = [o for o in range(m) if (mask >> o) & 1]; Sc = [o for o in range(m) if not (mask >> o) & 1]
            if not S: continue                                  # all positions in one row: impossible for t > 1
            c, _ = count_S(t, orbits, G, reps, S, Sc, frow, cyc, fix_label, rngh); E += 12 * c
        print(json.dumps({"t": t, "d": d, "g": "3-cycle", "orbits": m + 1, "mode": "all", "exact_labelled": E,
                          "exact_log2": round(math.log2(E), 2) if E else None, "est_log2E": est.get("log2E"),
                          "secs": round(time.time() - t0)}), flush=True)
    else:
        NS, p, seed = int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
        rng = np.random.default_rng(seed); vals = []
        for i in range(NS):
            inS = rng.random(m) < p
            S = [o for o in range(m) if inS[o]]; Sc = [o for o in range(m) if not inS[o]]
            if not S: vals.append(0.0); continue
            c, _ = count_S(t, orbits, G, reps, S, Sc, frow, cyc, fix_label, rngh)
            q = p ** len(S) * (1 - p) ** len(Sc); vals.append(12 * c / q)
            if (i + 1) % 10 == 0 or i + 1 == NS:
                v = np.array(vals); mean = v.mean(); se = v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float("nan")
                print(json.dumps({"t": t, "d": d, "g": "3-cycle", "orbits": m + 1, "mode": "sample", "p": p, "supports": len(v),
                                  "supports_with_solutions": int((v > 0).sum()), "E_hat": mean, "se": se,
                                  "E_hat_log2": round(math.log2(mean), 2) if mean > 0 else None, "est_log2E": est.get("log2E"),
                                  "secs": round(time.time() - t0)}), flush=True)
