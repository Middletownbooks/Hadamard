"""Stage A validation driver: every cyclic subgroup U = <u> of (Z/m)* and every twist ('id', 'BC' if d even, 'BCD' if
3 | d) on small Z_m whose row-sum groups have at most MAXBITS bits; exact count (b1tw_exact.exact) against the Gaussian
(b1tw_lib.gaussian).  Appends JSON lines to OUT; skips (m, u, tw) already in OUT (restart-safe).
usage: b1tw_validate.py OUT m [m ...]"""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b1tw_lib import Slice, gaussian, cyclic_subgroups
from b1tw_exact import exact, MAXBITS

OUT = sys.argv[1]
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        z = json.loads(l); done.add((z["m"], z["u"], z["tw"]))
for m in map(int, sys.argv[2:]):
    for d, u, U in cyclic_subgroups(m):
        if d == 1: continue
        for tw in ("id", "BC", "BCD"):
            if (tw == "BC" and d % 2) or (tw == "BCD" and d % 3) or (m, u, tw) in done: continue
            sl = Slice(m, u, tw)
            if sl.r < 5 or max(len(b) for b in sl.groups.values()) > MAXBITS: continue
            t0 = time.time(); ex = exact(sl)
            if ex is None: continue
            ga = gaussian(sl, N=40000)
            rec = dict(ga, exact=ex, exact_log2=(round(math.log2(ex), 2) if ex else None),
                       group_bits={g: len(b) for g, b in sl.groups.items()}, secs=round(time.time() - t0, 1))
            print(json.dumps(rec), flush=True)
            with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
