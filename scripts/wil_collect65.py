"""Collect Kaggle shard outputs for n=65 and check the W4 completeness ledger.
Usage: wil_collect65.py DIR [K]   (DIR holds wil_n65_c<case>_s<i>of<K>.json files)"""
import sys, json, glob, os, numpy as np
from wil_common import expand, verify
D = sys.argv[1]; K = int(sys.argv[2]) if len(sys.argv) > 2 else 32
cost = {tuple(c["case"]): c for c in json.load(open("logs/wil_n65_cost.json"))["cases"]}
cases = list(cost)
report = {"K": K, "cases": [], "solutions": []}
for ci, case in enumerate(cases):
    got = {}
    for f in glob.glob(os.path.join(D, f"wil_n65_c{ci}_s*of{K}.json")):
        r = json.load(open(f)); c = r["cases"][0]
        assert tuple(c["case"]) == case and r["shard"][1] == K
        got[r["shard"][0]] = c
        for s in r["solutions"]:
            Q = expand(np.array(s["halves"]), 65)
            report["solutions"].append({"case": case, "verified": verify(Q), "halves": s["halves"]})
    missing = sorted(set(range(K)) - set(got))
    ledger = sum(c["ledger"] for c in got.values()); work = sum(c["work"] for c in got.values())
    report["cases"].append({"case": case, "shards_present": len(got), "missing_shards": missing,
                            "ledger": ledger, "combos": cost[case]["combos"], "work": work, "work_expected": cost[case]["work"],
                            "complete": not missing and ledger == cost[case]["combos"] and work == cost[case]["work"],
                            "sols": sum(c["sols"] for c in got.values())})
report["all_complete"] = all(c["complete"] for c in report["cases"])
report["N_found"] = len(report["solutions"]); report["all_verified"] = all(s["verified"] for s in report["solutions"])
for c in report["cases"]:
    m = c["missing_shards"]
    print(c["case"], "shards", c["shards_present"], "/", K, "missing", (m if len(m) <= 12 else m[:12] + ["..."]),
          "ledger", c["ledger"], "/", c["combos"], "work %.3e / %.3e" % (c["work"], c["work_expected"]), "complete", c["complete"], "sols", c["sols"])
print("all_complete", report["all_complete"]); print("solutions:", report["N_found"], "verified:", report["all_verified"])
json.dump(report, open("williamson/n65_collected.json", "w"), default=int)
