r"""Calibrate theory/nd_stage2_est.py on the 573 Paley(11) + x29 slice: predicted stage-2 completions summed over the
y's that scripts/paley_s2join.py processed exactly (logs/s2join_573_*.jsonl) against the solutions it found [E vs M].
usage: python3 theory/nd_calib573.py [--N 20000]"""
import sys, os, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_stage2_est import Slice

N = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else 20000
S = Slice(572, 289, 11)
tot = 0.0; n = 0; found = 0
for f in sorted(glob.glob('logs/s2join_573_*.jsonl')):
    for line in open(f):
        r = json.loads(line)
        if r.get('FOUND'): continue
        d = np.load(os.path.join('theory/logs', os.path.basename(r['file'])))
        A, W = d['A'][r['y']], d['W'][r['y']]
        Y = (A + 10 * W).tolist()
        e = S.estimate(Y, N=N, seed=1)
        sol = r.get('solutions', 0) or 0
        found += sol; n += 1
        if e.get('stage2_log2E') is not None and e['stage2_log2E'] > -1e9: tot += 2.0 ** e['stage2_log2E']
        print(json.dumps(dict(file=os.path.basename(r['file']), y=r['y'], solutions=sol, **e)), flush=True)
print(json.dumps(dict(ys=n, predicted_total=tot, log2_predicted_total=float(np.log2(tot)) if tot > 0 else None,
                      found_total=found)), flush=True)
