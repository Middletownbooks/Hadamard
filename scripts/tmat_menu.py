"""Menu of findable cyclotomic T-matrices: every prime t in [lo, hi) and odd k | t-1 with 5 <= orbits <= 60; exact
row-sum condition (three sums = 0 mod k, one = +-1 mod k), then the Gaussian of tmat_estimate (n_bits = 3m - log2E).
Usage: tmat_menu.py part/W lo hi   (job part of W)"""
import os, sys, json, math, itertools, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
resource.setrlimit(resource.RLIMIT_AS, (int(2.5 * 2**30),) * 2)
from tmat_estimate import estimate
def isprime(n): return n > 1 and all(n % p for p in range(2, int(n ** .5) + 1))
def feasible(t, k):
    R = int(math.isqrt(t))
    for s0 in range(-R, R + 1):
        if (s0 - 1) % k and (s0 + 1) % k: continue
        rest = t - s0 * s0
        for a in range(0, R + 1, k):
            for b in range(a, R + 1, k):
                for c in range(b, R + 1, k):
                    if a*a + b*b + c*c == rest: return True
    return False
part, W = (int(x) for x in sys.argv[1].split('/'))
lo, hi = int(sys.argv[2]), int(sys.argv[3])
jobs = []
for t in range(lo, hi):
    if not isprime(t): continue
    for k in range(3, t, 2):
        if (t - 1) % k: continue
        m = (t - 1) // k + 1
        if m > 60 or m < 5: continue
        jobs.append((t, k, m))
for i, (t, k, m) in enumerate(jobs):
    if i % W != part: continue
    if not feasible(t, k):
        print(json.dumps(dict(t=t, k=k, orbits=m, rowsum_feasible=False)), flush=True); continue
    r = estimate(t, k, N=100000)
    r['n_bits'] = round(3 * m - r['log2E'], 1) if r['log2E'] != float('-inf') else None
    print(json.dumps(r), flush=True)
