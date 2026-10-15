"""Calibrate tmat_twist's d = 2 involution classes against exact counts at small t = 3 (mod 4) (as 263, 191, 179).
Exact count by vectorised enumeration of all base labellings (8^(orbits) * |fix0|)."""
import os, sys, json, math, itertools, resource
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tmat_twist as tw

def exact(t, d, g, chunk=1 << 18):
    u, orbits, G, fix0, reps = tw.setup(t, d, g)
    pos = np.array(orbits); nb = len(orbits); total = 0
    for l0 in fix0:
        for start in range(0, 8 ** nb, chunk):
            idx = np.arange(start, min(start + chunk, 8 ** nb), dtype=np.int64)
            base = np.stack([(idx >> (3 * j)) & 7 for j in range(nb)], 1)
            lab = np.empty((len(idx), t), np.int64); lab[:, 0] = l0
            for j in range(G.shape[0]): lab[:, pos[:, j]] = G[j][base]
            tot = np.zeros((len(idx), len(reps)), np.int64)
            for r in range(4):
                x = np.where(tw.ROWS[lab] == r, tw.SGN[lab], 0)
                for c, s in enumerate(reps):
                    tot[:, c] += (x * np.roll(x, -s, axis=1)).sum(1)
            total += int((tot == 0).all(1).sum())
    return total

reps = tw.class_reps()
for t in [int(x) for x in sys.argv[1].split(",")]:
    for ct, g in sorted(reps.items()):
        if 2 % tw.order_of(tw.LMAP[g]) != 0: continue
        e = tw.estimate(t, 2, g, 100000)
        u, orbits, G, fix0, rp = tw.setup(t, 2, g)
        n_conf = 8 ** len(orbits) * max(1, len(fix0))
        b = exact(t, 2, g) if (fix0 and e.get("feasible") is not False and n_conf <= 2 ** 30) else None
        print(json.dumps({"t": t, "d": 2, "g": ct, "orbits": len(orbits) + 1, "exact": b,
                          "exact_log2": round(math.log2(b), 1) if b else None, "est_log2E": e.get("log2E"),
                          "feasible": e.get("feasible"), "why": e.get("why")}), flush=True)
