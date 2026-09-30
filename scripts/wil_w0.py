"""W0 build validation: pipeline vs independent brute force, as sets, plus stage-A inclusion. Usage: wil_w0.py n p q"""
import sys, json, math, numpy as np
from wil_common import expand, rowsum_cases, epsilon
import wil_search as ws

n, p, q = map(int, sys.argv[1:4]); hn = (n + 1) // 2
brute = json.load(open(f"logs/wil_brute_n{n}.json"))
B = {np.array(h, np.int64).tobytes() for h in brute["solutions"]}
nom = json.load(open(f"logs/wil_n{n}_nomult.json")); mul = json.load(open(f"logs/wil_n{n}.json"))
Pn = {np.array(s["halves"], np.int64).tobytes() for s in nom["solutions"]}
Pm = set()
for s in mul["solutions"]:
    X = expand(np.array(s["halves"], np.int64), n)
    for u in range(1, n):
        if math.gcd(u, n) == 1: Pm.add(X[:, (u * np.arange(n)) % n][:, :hn].astype(np.int64).tobytes())
# stage-A inclusion: compressions of every brute solution are among the complementary quadruples (canonical and raw)
rng = np.random.default_rng(0); incl = True
for case in rowsum_cases(n):
    for m in (p, q):
        _, Q = ws.comp_quads(n, m, case, rng)
        raw = {x.tobytes() for x in Q.astype(np.int64)}
        can = {x.tobytes() for x in ws.canon(Q, m).astype(np.int64)} if len(Q) else set()
        for h in brute["solutions"]:
            X = expand(np.array(h), n)
            if tuple(X.sum(1)) != case: continue
            C = X.reshape(4, n // m, m).sum(1)[:, :(m + 1) // 2].astype(np.int64)    # i -> i mod m
            incl &= C.tobytes() in raw and ws.canon(C[None], m)[0].astype(np.int64).tobytes() in can
res = {"n": n, "N_brute": len(B), "N_pipeline_nomult": len(Pn), "N_pipeline_mult_expanded": len(Pm),
       "sets_equal_nomult": Pn == B, "sets_equal_mult": Pm == B, "stageA_inclusion": bool(incl),
       "product_condition_on_brute": all(np.all(np.prod(expand(np.array(h), n)[:, 1:], 0) == epsilon(n, tuple(expand(np.array(h), n).sum(1)))) for h in brute["solutions"])}
res["W0_pass"] = all(v for k, v in res.items() if k.startswith(("sets", "stageA", "product"))) and len(B) == len(Pn) == len(Pm)
print(json.dumps(res))
