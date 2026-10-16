"""Stage-1 and stage-2 Gaussian estimates for the six undecided Paley-type c = 4 configurations [E]
(theory/ab_gauss.stage1: exact count of energy/row-sum feasible y-quadruples x lag-class Gaussian;
 theory/ab_gauss.stage2_typical on y's drawn by stage1_sampler).  usage: estimate.py NAME"""
import sys, os, json, time, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'theory'))
import numpy as np
import ab_s1sd as SD, ab_sweep as W, ab_gauss as AGs
name = sys.argv[1]
cfg = SD.cfg_from_json(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'cfg_{name}.json')).read())
sl = W.build(cfg); t0 = time.time()
s1 = AGs.stage1(sl, N=20000)
rec = dict(name=name, orbits=int(sl.r), L1=len(sl.breps), stage1=s1, secs1=round(time.time() - t0, 1))
draw, info = AGs.stage1_sampler(sl, seed=3)
s2 = []
for i in range(4):
    Y = draw(); r = AGs.stage2_typical(sl, Y, N=6000, seed=i + 1); s2.append(r)
rec['stage2_typical'] = s2
vals = [r['stage2_log2E'] for r in s2 if r.get('stage2_log2E') is not None]
rec['stage2_median'] = float(np.median(vals)) if vals else None
if s1.get('stage1_log2E') not in (None, -math.inf) and rec['stage2_median'] is not None:
    rec['log2_cores_expected'] = round(s1['stage1_log2E'] + rec['stage2_median'], 1)
rec['secs'] = round(time.time() - t0, 1)
print(json.dumps(rec), flush=True)
