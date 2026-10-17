"""Planted test of theory/ab_s2_zq_screen.run itself: random lifts x of random letter quadruples, target replaced by
T' = sum_k PAF_G(x_k); the screen must keep y = compress(x) alive under every character.
usage: python3 planted_screen_run.py CFG.json N SEED"""
import sys, os, json, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../theory'))
import numpy as np, ab_sweep as W, ab_s1sd as SD, ab_s2_zq_screen as SC
cfg = SD.cfg_from_json(open(sys.argv[1]).read()); N = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
fails = 0; d = tempfile.mkdtemp()
for trial in range(N):
    sl = W.build(cfg); L1 = len(sl.breps)
    lab = rng.choice(np.array([-1, 1], np.int8), size=(4, sl.r)); X = lab[:, sl.orb].astype(np.int64)
    sl.T = sl.G.paf(X).sum(0)
    Y = np.zeros((4, L1), np.int64)
    for j, bb in enumerate(sl.breps): Y[:, j] = X[:, np.nonzero(sl.bidx == bb)[0]].sum(1)
    yb = os.path.join(d, f'y{trial}.bin'); Y.astype(np.int8).tofile(yb)
    r = SC.run(sl, yb, os.path.join(d, f'o{trial}.jsonl'), [0, 1, 2, 3])
    fails += r['alive'] != 1
    print(json.dumps(dict(trial=trial, alive=r['alive'], log2_counts=[v['log2_count'] for v in r['done'].values()])), flush=True)
print(json.dumps(dict(trials=N, false_kills=fails)))
