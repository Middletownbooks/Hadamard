r"""[M] Planted test of the exact S14 stage-1 decider (theory/h863_s14_stage1_exact.py): a random tau-invariant y with
few large orbit values (low energy, like real stage-1 solutions) is planted; its own Q-compressions and 20-lag H-level
PAF vector become z and the target; the exhaustive decider must return it (among possibly other solutions).
usage: python3 theory/h863_s14_stage1_planted.py [trials] [seed]"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_s14_stage1_exact as X
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
ORBS = sorted({tuple(sorted({k, X.TAU[k]})) for k in range(20)})
big = [v for v in X.A if abs(v) in (13, 15)]
for tr in range(trials):
    ys = []
    for i in range(4):
        y = [0] * 20
        for o in ORBS:
            v = int(rng.choice(big)) if rng.random() < 0.12 else int(rng.choice([-1, 1]))
            for k in o: y[k] = v
        ys.append(y)
    z = tuple((sum(y[0:5]), sum(y[5:10]), sum(y[15:20])) for y in ys)
    tgt = X.sigs(np.array(ys)).sum(0)
    r = X.decide(z, tgt)
    found = any(sorted(map(tuple, s['y'])) == sorted(map(tuple, ys)) or
                all(s['y'][i] == ys[i] for i in range(4)) for s in r['solutions'])
    print(json.dumps(dict(trial=tr, z=z, energy=int(tgt[0]), sizes=r['sizes'], nsol=len(r['solutions']),
                          planted_found=found, all_recheck=all(s['recheck'] for s in r['solutions']))), flush=True)
