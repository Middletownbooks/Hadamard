"""Stage-1 annealing for a generalised-Paley slice (b = 1; Z_F pointwise fixed, index-s classes on Z_q): find coset sums
y_k on Z_F (letters a + f * sum of s signs) with sum_k PAF(y_k) = T_J.  Row sums then follow (sum_h T_J = 4).
usage: cyclo_stage1_sa.py F q s [restarts] [steps]"""
import sys, os, json, time, itertools as it
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_staged import _s1_anneal
F, q, s = map(int, sys.argv[1:4]); R = int(sys.argv[4]) if len(sys.argv) > 4 else 10; steps = int(float(sys.argv[5])) if len(sys.argv) > 5 else 2e8
f = (q - 1) // s; m = F * q
vals = np.array(sorted({a + f * sum(c) for a in (1, -1) for c in it.product((1, -1), repeat=s)}), np.int64)
T = np.full(F, -4 * q, np.int64); T[0] = 4 * m - 4 * (q - 1)
add = np.array([[(g + h) % F for h in range(F)] for g in range(F)], np.int64); sub = np.array([[(g - h) % F for h in range(F)] for g in range(F)], np.int64)
pafs = lambda y: np.array([(y * np.roll(y, -h)).sum() for h in range(F)])
rng = np.random.default_rng(1); found = []
for r in range(R):
    Y = rng.choice(np.array([1, -1]), (4, F)).astype(np.int64)
    E = (sum(pafs(Y[k]) for k in range(4)) - T).astype(np.int64)
    t0 = time.time()
    c = _s1_anneal(Y, E, add, sub, vals, steps, float(4 * q * q), 0.5, int(rng.integers(1 << 30)), np.ones(F, np.bool_))
    rec = dict(F=F, q=q, s=s, restart=r, final_cost=int(c), secs=round(time.time() - t0, 1))
    if c == 0:
        assert (sum(pafs(Y[k]) for k in range(4)) == T).all()
        rec.update(y=Y.tolist(), rows=Y.sum(1).tolist()); found.append(Y.tolist())
    print(json.dumps(rec), flush=True)
