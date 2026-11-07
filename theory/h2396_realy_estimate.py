r"""Gate A2: the repository estimator's stage-2 prediction for a real stage-1 solution y of the plug-in analogue on
Z_n = Z_{n/2} x Z_2 (CRT), n = 2k with k odd (theory/open_screen_deficit.analyse with real = True) [E].
Input: a file with the four u-sequences on Z_n (one per line, as written by theory/h2396_s1cp.py).
usage: python3 theory/h2396_realy_estimate.py n FILE"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import open_screen_deficit as D
import ab_sweep as SW
n = int(sys.argv[1]); k = n // 2; assert n % 2 == 0 and k % 2 == 1
U = np.array([list(map(int, l.split())) for l in open(sys.argv[2]) if l.strip()])
key = ["b1", 13 * n + 1, [], [k, 2], 13, 2, []]
sl = SW.build(D.cfg_from_key(key))
Y = np.zeros((4, n), np.int64)
for t in range(n): Y[:, (t % k) * 2 + t % 2] = U[:, t]
assert (sl.Gp.paf(Y).sum(0) == sl.T1).all(), "index map or stage-1 check failed"
r = D.analyse(sl, Y, real=True)
print(json.dumps(dict(key=key, **r), default=float))
