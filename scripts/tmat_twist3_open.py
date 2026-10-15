"""3-cycle-twisted ansatz (the only nontrivial twist for odd d) at the open v: exact row sums (ILP), then Gaussian
when the orbit count is moderate. Odd d must be divisible by 3 for a 3-cycle twist."""
import os, sys, json, math, resource
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tmat_twist as tw
g3 = [g for ct, g in tw.class_reps().items() if ct == ((1, 1), (3, 1))][0]
for v in [int(x) for x in sys.argv[1].split(",")]:
    for d in range(3, v, 6):                      # odd multiples of 3
        if (v - 1) % d or (v - 1) // d < 3: continue
        u, orbits, G, fix0, reps = tw.setup(v, d, g3)
        rec = {"t": v, "d": d, "g": "3-cycle", "orbits": len(orbits) + 1, "classes": len(reps)}
        if not tw.rowsum_feasible(v, orbits, G, fix0):
            rec.update(feasible=False, why="row sums")
        elif len(orbits) <= 150:
            rec = tw.estimate(v, d, g3, 100000); rec["g"] = "3-cycle"
        else:
            rec.update(feasible=True, why="row sums ok; not estimated (orbits > 150)")
        print(json.dumps(rec), flush=True)
