"""Planted version of z2zq_crosscheck.py (positive counts): random lifts x on a small slice, target T' = sum_k PAF_G(x_k),
compare the independent exact count with theory/ab_s2_z2zq_screen.run for r = 1, 2, 3.  usage: CFG.json N SEED"""
import sys, os, json, tempfile, importlib.util
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '../../../../theory'))
import numpy as np, ab_sweep as W, ab_s1sd as SD, ab_s2_z2zq_screen as SC
spec = importlib.util.spec_from_file_location('cc', os.path.join(here, 'z2zq_crosscheck.py'))
cfgf, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
sys.argv = [sys.argv[0], cfgf, '/dev/null', '0', '0']
src = open(os.path.join(here, 'z2zq_crosscheck.py')).read().split('Y = np.fromfile')[0]
g = {'__file__': os.path.join(here, 'z2zq_crosscheck.py')}; exec(compile(src, 'z2zq_crosscheck.py', 'exec'), g)
rng = np.random.default_rng(seed); d = tempfile.mkdtemp(); bad = 0
for trial in range(N):
    sl = W.build(SD.cfg_from_json(open(cfgf).read())); L1 = len(sl.breps)
    lab = rng.choice(np.array([-1, 1], np.int8), size=(4, sl.r)); X = lab[:, sl.orb].astype(np.int64); Tp = sl.G.paf(X).sum(0)
    sl.T = Tp; g['T'] = Tp; g['sl'] = sl
    Y = np.zeros((4, L1), np.int64)
    for j, bb in enumerate(sl.breps): Y[:, j] = X[:, np.nonzero(sl.bidx == bb)[0]].sum(1)
    yb = os.path.join(d, f'y{trial}.bin'); Y.astype(np.int8).tofile(yb); out = []
    for r in (1, 2, 3):
        ex = g['count_y'](Y, r); res = SC.run(sl, yb, os.path.join(d, f'o{trial}_{r}.jsonl'), r)
        sc = list(res['done'].values())[0]['count']; out.append((str(ex), sc)); bad += float(ex) != sc
    print(json.dumps(dict(trial=trial, exact_int_vs_screen=out)), flush=True)
print(json.dumps(dict(mismatches=bad)))
