"""Planted check of theory/ab_s2_qzq_screen.run (the Q x Z_q level): random lifts x on a small slice, target
T' = sum_k PAF_G(x_k); the screen must keep y = compress(x), and its count must equal an independent exact integer
count (per-coset pattern products enumerated directly, PAF on Q x Z_q computed pointwise, dict convolution).
usage: python3 qzq_crosscheck_planted.py CFG.json N SEED"""
import sys, os, json, tempfile, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../theory'))
import numpy as np, ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s2_qzq_screen as SC
from collections import defaultdict
cfgf, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); rng = np.random.default_rng(seed); d = tempfile.mkdtemp(); bad = 0
for trial in range(N):
    sl = W.build(SD.cfg_from_json(open(cfgf).read())); L1 = len(sl.breps); q = sl.q
    lab = rng.choice(np.array([-1, 1], np.int8), size=(4, sl.r)); X = lab[:, sl.orb].astype(np.int64); sl.T = sl.G.paf(X).sum(0)
    Y = np.zeros((4, L1), np.int64)
    for j, bb in enumerate(sl.breps): Y[:, j] = X[:, np.nonzero(sl.bidx == bb)[0]].sum(1)
    # independent count: every lift of every block (all fibre patterns), compressed along H' pointwise
    qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]; cell = qg * q + sl.zco                    # point -> (q, t) cell
    tt = np.zeros(4 * q, np.int64); np.add.at(tt, cell, sl.T); Tq = tt.reshape(4, q)
    def paf_qz(xb):                                                                    # xb: (4, q) -> (4, q) PAF on Q x Z_q
        return np.array([[sum(int(xb[a, t] * xb[a ^ dl, (t + tau) % q]) for a in range(4) for t in range(q)) for tau in range(q)] for dl in range(4)])
    dists = []
    for k in range(4):
        pats = [AGs.fibre_patterns(sl, j, int(Y[k, j]))[0] for j in range(L1)]
        dist = defaultdict(int)
        for choice in itertools.product(*[range(len(P)) for P in pats]):
            labk = np.zeros(sl.r, np.int64)
            for j, ci in enumerate(choice):
                for t_, (o, _) in enumerate(sl.parts[j]): labk[o] = pats[j][ci][t_]
            xb = np.zeros(4 * q, np.int64); np.add.at(xb, cell, labk[sl.orb])
            dist[tuple(paf_qz(xb.reshape(4, q)).ravel().tolist())] += 1
        dists.append(dist)
    A = defaultdict(int)
    for s1, c1 in dists[0].items():
        for s2, c2 in dists[1].items(): A[tuple(a + b for a, b in zip(s1, s2))] += c1 * c2
    T = tuple(Tq.ravel().tolist()); ex = 0
    for s3, c3 in dists[2].items():
        for s4, c4 in dists[3].items(): ex += A.get(tuple(t - a - b for t, a, b in zip(T, s3, s4)), 0) * c3 * c4
    yb = os.path.join(d, f'y{trial}.bin'); Y.astype(np.int8).tofile(yb)
    res = SC.run(sl, yb, os.path.join(d, f'o{trial}.jsonl')); sc = list(res['done'].values())[0]['count']
    bad += float(ex) != sc or res['alive'] != 1
    print(json.dumps(dict(trial=trial, exact=str(ex), screen=sc, alive=res['alive'])), flush=True)
print(json.dumps(dict(trials=N, mismatches_or_false_kills=bad)))
