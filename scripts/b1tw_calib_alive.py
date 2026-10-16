"""Stage A calibration on live slices: every cyclic subgroup <u> of (Z/m)* and twist ('id', 'BC', 'BCD') on Z_m for the
given m whose row-sum groups are small: the exact count (b1tw_exact.exact); for slices with solutions, the Gaussian
(b1tw_lib.gaussian).  Restart-safe JSONL (every decided slice is logged; Gaussian only where exact > 0).
usage: b1tw_calib_alive.py OUT m [m ...]"""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b1tw_lib import Slice, gaussian, cyclic_subgroups, log2_rowsum_configs
from b1tw_lag0 import divisors
from b1tw_screen import quotient_ok
from b1tw_exact import exact, MAXBITS

OUT = sys.argv[1]; done = set()
if os.path.exists(OUT):
    for l in open(OUT): z = json.loads(l); done.add((z["m"], z["u"], z["tw"]))
for m in map(int, sys.argv[2:]):
    for d, u, U in cyclic_subgroups(m):
        if d == 1: continue
        for tw in ("id", "BC", "BCD"):
            if (tw == "BC" and d % 2) or (tw == "BCD" and d % 3) or (m, u, tw) in done: continue
            sl = Slice(m, u, tw)
            if sl.r < 6 or max(len(b) for b in sl.groups.values()) > MAXBITS: continue
            if log2_rowsum_configs(sl) == float("-inf"): continue
            t0 = time.time()
            ex = exact(sl)
            if ex is None: continue
            rec = dict(m=m, u=u, d=d, tw=tw, orbits=sl.r, classes=sl.L, minus1_in_U=sl.minus1_in_U, exact=ex,
                       exact_log2=(round(math.log2(ex), 2) if ex else None))
            if ex:
                ga = gaussian(sl, N=40000)
                rec.update(**{k: ga.get(k) for k in ("log2E", "n_cond", "why", "steps", "mahal")})
            rec["secs"] = round(time.time() - t0, 1)
            with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
            if rec['exact']: print(json.dumps(rec), flush=True)
