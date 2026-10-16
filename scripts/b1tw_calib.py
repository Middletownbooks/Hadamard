"""Stage A calibration table: for each record of a validation log (b1tw_validate.py), add the exact quotient screens
(b1tw_lag0: lag-0 energy DP and all-lags quotient test for every subgroup) and compare the Gaussian log2 E with the exact
count.  usage: b1tw_calib.py VALIDATION.jsonl OUT.jsonl"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b1tw_lib import Slice
from b1tw_lag0 import quotient_screen, divisors, screen

seen = set()
if os.path.exists(sys.argv[2]):
    for l in open(sys.argv[2]): z = json.loads(l); seen.add((z["m"], z["u"], z["tw"]))
out = open(sys.argv[2], "a")
for l in open(sys.argv[1]):
    z = json.loads(l)
    if (z["m"], z["u"], z["tw"]) in seen: continue
    sl = Slice(z["m"], z["u"], z["tw"])
    q = {K: quotient_screen(sl, K) for K in divisors(sl.m) if 1 < K < sl.m}
    lag0 = screen(sl)
    rec = dict(m=z["m"], u=z["u"], d=z["d"], tw=z["tw"], orbits=z["orbits"], classes=z["classes"],
               minus1_in_U=z["minus1_in_U"], log2E=z.get("log2E"), why=z.get("why"), exact=z["exact"],
               exact_log2=z["exact_log2"], lag0_alive=all(x["attainable"] for x in lag0),
               quotient={str(k): v for k, v in q.items()}, quotient_alive=not any(v is False for v in q.values()))
    out.write(json.dumps(rec) + "\n"); out.flush(); print(json.dumps(rec), flush=True)
