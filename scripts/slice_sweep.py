"""Sweep Paley(p) + multiplier orbit slices (b = 1) over the open orders; Gaussian estimate per distinct slice
(theory/paley_slices.estimate, calibrated at 359/419/573).  One JSON line per slice.
Candidates: m = v - 1 = n' p with p prime, p || m, cyclic G' = Z_n'; multipliers u in (Z/n')^*, one per distinct
orbit partition, with orbits per block in [omin, omax].
usage: slice_sweep.py v p [--b 1] [--omin 8] [--omax 40] [--N 20000]   (b >= 3: H' = c Z_n', c = b + 1)"""
import sys, os, json, math
from math import gcd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'theory'))
from paley_slices import estimate

def orbit_partition(n, u):
    seen, out = set(), []
    for g in range(n):
        if g in seen: continue
        o, x = [], g
        while x not in o: o.append(x); x = (x * u) % n
        seen.update(o); out.append(tuple(sorted(o)))
    return tuple(sorted(out))

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--omin": "8", "--omax": "40", "--N": "20000", "--b": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, p = pos; b = int(opt["--b"]); n = (v - b) // p; assert n * p == v - b and n % p
    seen = {}
    for u in range(2, n):
        if gcd(u, n) != 1: continue
        part = orbit_partition(n, u)
        if part in seen: continue
        seen[part] = u
        k = len(part)
        if not (int(opt["--omin"]) <= k <= int(opt["--omax"])): continue
        d = 1; x = u
        while x != 1: x = x * u % n; d += 1
        try:
            r = estimate(v, b, p, (n,), None, f"mult:{u}", N=int(opt["--N"]))
        except Exception as e:
            r = dict(error=repr(e))
        print(json.dumps({**dict(v=v, p=p, n=n, u=u, order=d, orbits=k), **r}), flush=True)
