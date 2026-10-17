"""Planted test of theory/ab_s2_hzq.screen_y: random lifts x of random letter quadruples, target T' = sum_k PAF_G(x_k);
the H' x Z_q-level count of y = compress(x) must be >= 1 for every character.  usage: planted_hzq_test.py CFG N SEED"""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../theory'))
import numpy as np, ab_s1sd as SD, ab_s2_hzq as HZ
cfg = SD.cfg_from_json(open(sys.argv[1]).read()); S = HZ.setup(cfg); sl = S['sl']; L1 = S['L1']
rng = np.random.default_rng(int(sys.argv[3])); fails = 0
for trial in range(int(sys.argv[2])):
    lab = rng.choice(np.array([-1, 1], np.int8), size=(4, sl.r)); X = lab[:, sl.orb].astype(np.int64)
    HZ.set_target(S, sl.G.paf(X).sum(0))
    Y = np.zeros((4, L1), np.int64)
    for j, bb in enumerate(sl.breps): Y[:, j] = X[:, np.nonzero(sl.bidx == bb)[0]].sum(1)
    res = [HZ.screen_y(S, Y, p) for p in range(4)]
    fails += any(c == 0 for c, _ in res)
    print(json.dumps(dict(trial=trial, log2=[round(math.log2(c), 2) if c else c for c, _ in res], by=[b for _, b in res])), flush=True)
print(json.dumps(dict(trials=int(sys.argv[2]), false_kills=fails)))
