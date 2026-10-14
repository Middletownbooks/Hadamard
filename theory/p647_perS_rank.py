r'''Per-S-solution D-lattice rank and candidate count, 647 x353 (calibration j = 11514) against the 495 x83 control
(y_B #9): the rank of the lattice spanned by the D-residual differences of each S-solution's candidates [M].
Output: logs/p647/perS_rank.log.'''
import sys, json, math, collections, itertools; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np
from nd_cross import Cross, orbit_maps
from nd_joint_exhaust import enumerate_yB
import p647_common as PC, p647_cross as PX
rng = np.random.default_rng(6)
def ranks(X, bl, S, n=600, cap=20000):
    TE = 2 * (X.TS + X.TD); out = collections.Counter(); cands = []
    for s in S[rng.choice(len(S), min(n, len(S)), replace=False)]:
        rngs = [bl[k].pairs[bl[k].bd[s[k]]:bl[k].bd[s[k] + 1]] for k in range(4)]
        tot = int(np.prod([len(r) for r in rngs])); cands.append(tot)
        if tot > cap: out['big'] += 1; continue
        V = []
        for combo in itertools.product(*rngs):
            v = -TE.copy()
            for k, p in enumerate(combo):
                i, j = divmod(int(p), bl[k].n1); v = v + bl[k].P0[i] + bl[k].P1[j]
            V.append(v)
        V = np.array(V, float)
        out[int(np.linalg.matrix_rank(V[1:] - V[0])) if len(V) > 1 else 0] += 1
    return out, cands
X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
bl = PX.build(X, ys[0], X.J.ya_full(sols[11514]['layers']), maps, NH, oi)
S = np.load(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))), 'logs') + '/p647/calib/j11514_ssols.npy')
o, c = ranks(X, bl, S)
print('647 j11514: per-S D-lattice rank (full = 11):', sorted(o.items(), key=lambda kv: str(kv[0])), 'median candidates', int(np.median(c)), flush=True)
X5 = Cross(494, 83, 13, 19)
ys5 = enumerate_yB(494, 83, 13, 19, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))), 'logs') + '/p647/control/enum_494_83_13_19.jsonl', 1, log=lambda s: None)
reps5, oi5, maps5 = orbit_maps(X5); NH5 = PX.lag_tensor(X5, oi5)
rec5, sols5 = X5.J.test(ys5[9], maxsol=10 ** 7)
tot = collections.Counter(); cc = []
for j in rng.choice(len(sols5), 20, replace=False):
    b5 = PX.build(X5, ys5[9], X5.J.ya_full(sols5[int(j)]['layers']), maps5, NH5, oi5)
    r, Ss, hits, cand = PX.solve(X5, b5, P=61)
    if Ss:
        o, c = ranks(X5, b5, np.asarray(Ss).reshape(-1, 4), n=60); tot.update(o); cc += c
print('495 yb9: per-S D-lattice rank (full = 8):', sorted(tot.items(), key=lambda kv: str(kv[0])), 'median candidates', int(np.median(cc)), flush=True)
