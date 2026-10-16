"""b = 0 (Goethals-Seidel core of prime length v) multiplier-slice screen at the prime open v: every odd d | v-1 with
15-130 orbits, untwisted (g = id) and, when 3 | d, the 3-cycle twist. Exact block-sum feasibility, then the Gaussian of
gs_twist.py (n_bits = 4m - log2E constraint bits).  Usage: b0_open_screen.py v1,v2,...  [N]"""
import os, sys, json, resource
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_twist as gt

REPS = gt.twisted_reps()
ID = REPS[((1, 1), (1, 1), (1, 1), (1, 1))]
C3 = REPS[((1, 1), (3, 1))]
if __name__ == "__main__":
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
    for v in [int(x) for x in sys.argv[1].split(",")]:
        for d in range(3, v, 2):
            if (v - 1) % d or not 15 <= (v - 1) // d + 1 <= 130: continue
            for g in ([ID, C3] if d % 3 == 0 else [ID]):
                rec = gt.estimate(v, d, g, N); rec["twist"] = "id" if g == ID else "3-cycle"
                print(json.dumps(rec), flush=True)
