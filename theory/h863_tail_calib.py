r"""[M vs E] Tail calibration of the lattice Gaussian on the mixed level at v = 863 (as report/paley_523_answer.md 3.4):
exact fraction of uniform (y, w)-fibre configurations hitting the target on nested subsets of the lag classes that
are not fixed by (y, w), against q523_gauss.gauss_bits fitted to the same samples.  Input: the int16 class-PAF arrays
written by theory/h863_tail_sample.py (Codex certificate y, w = data/h863_g5_solutions.json[0]).
usage: python3 theory/h863_tail_calib.py OUT.jsonl V1.npy [V2.npy ...]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h863_lib as L, q523_gauss as QG
out = sys.argv[1]
V = np.concatenate([np.load(f) for f in sys.argv[2:]]).astype(np.int64)
N = len(V); T = L.FULL_TARGET
var = [c for c in range(72) if V[:, c].min() != V[:, c].max()]
fixed_ok = all((V[:, c] == T[c]).all() for c in range(72) if c not in var)
rng = np.random.default_rng(0)
res = []
fo = open(out, 'w')
def log(d):
    s = json.dumps(d); print(s, flush=True); fo.write(s + '\n')
log(dict(N=N, varying=len(var), fixed_classes_on_target=fixed_ok))
orders = [('K-lags first', sorted(var, key=lambda c: (L.LAGREPS[c][0] != 0, c)))]
for k in range(3):
    orders.append(('random %d' % k, list(rng.permutation(var))))
for name, order in orders:
    for k in range(1, 8):
        sub = [int(c) for c in order[:k]]
        hits = int((V[:, sub] == T[sub]).all(1).sum())
        g = QG.gauss_bits(V[:200000, sub], T[sub])
        log(dict(order=name, k=k, lags=[L.LAGREPS[c] for c in sub], hits=hits,
                 exact_bits=round(-math.log2(hits / N), 2) if hits else None, gauss_bits=g['bits'], rank=g['rank']))
        if hits < 20: break
