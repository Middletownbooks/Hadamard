"""W2-T63 (amendment A3): targeted check of the pipeline on the constructed Turyn quadruple of order 63."""
import sys, json, math, numpy as np
sys.argv = ["x", "63", "7", "9"]
import wil_search as ws
from wil_common import expand, verify, rowsum_cases, epsilon

n, p, q = 63, 7, 9; hn = (n + 1) // 2
Q = np.array(json.load(open("williamson/turyn63.json"))["Q"], dtype=np.int64)
assert verify(Q)
Q = Q * np.sign(Q.sum(1))[:, None]                     # normalise: positive row sums, sorted
Q = Q[np.argsort(-Q.sum(1), kind="stable")]
case = tuple(int(x) for x in Q.sum(1)); assert case in rowsum_cases(n)
rng = np.random.default_rng(12345); res = {"case": case}
comp = lambda X, m: X.reshape(4, n // m, m).sum(1)[:, :(m + 1) // 2]   # i -> i mod m
# (a) stage-A inclusion
_, QP = ws.comp_quads(n, p, case, rng); _, QQ = ws.comp_quads(n, q, case, rng)
cp, cq = comp(Q, p), comp(Q, q)
res["a_Zp_in_stageA"] = any(np.array_equal(cp, x) for x in QP)
res["a_Zq_in_stageA"] = any(np.array_equal(cq, x) for x in QQ)
# (b) canonical combination present, every block liftable (after PSD filter)
QPr, QQr = ws.uniq_rows(ws.canon(QP, p)), ws.uniq_rows(ws.canon(QQ, q))
kp, kq = ws.canon(cp[None], p)[0], ws.canon(cq[None], q)[0]
res["b_canonical_in_reduced_lists"] = any(np.array_equal(kp, x) for x in QPr) and any(np.array_equal(kq, x) for x in QQr)
T = [ws.lift_table(n, ws.lifts(n, p, q, kp[b], kq[b])) for b in range(4)]
res["b_lift_counts"] = [len(t[0]) for t in T]; res["b_all_liftable"] = min(res["b_lift_counts"]) > 0
# (c) stage D on that one combination
W = rng.integers(1, 2**63, hn - 1, dtype=np.uint64) | np.uint64(1); K1 = np.uint64(rng.integers(1, 2**63, dtype=np.uint64) | 1)
sols, work = ws.join4(n, T, epsilon(n, case), W, K1)
res["c_solutions_on_combo"] = len(sols); res["c_all_verified"] = all(verify(S) for S in sols)
orbit = set()
for u in range(1, n):
    if math.gcd(u, n) == 1: orbit.add(Q[:, (u * np.arange(n)) % n].tobytes())
res["c_constructed_orbit_hit"] = any(S.astype(np.int64).tobytes() in orbit for S in sols)
res["W2_T63_pass"] = bool(res["a_Zp_in_stageA"] and res["a_Zq_in_stageA"] and res["b_canonical_in_reduced_lists"]
                          and res["b_all_liftable"] and res["c_solutions_on_combo"] > 0 and res["c_all_verified"] and res["c_constructed_orbit_hit"])
print(json.dumps(res, default=lambda o: bool(o) if isinstance(o, np.bool_) else int(o)))
json.dump(res, open("williamson/W2_T63_result.json", "w"), default=lambda o: bool(o) if isinstance(o, np.bool_) else int(o))
