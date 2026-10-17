"""Cross-check theory/ab_s2_z2zq_screen.run against an independent dict-DP implementation (the prototype's
count_y, reimplemented here) on real stage-1 survivors: the exact Z_2 x Z_q-level counts must be equal.
usage: python3 z2zq_crosscheck.py CFG.json YBIN N SEED"""
import sys, os, json, math, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../theory'))
import numpy as np, ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s2_s4 as S4, ab_s2_z2zq_screen as SC
from collections import defaultdict
cfg = SD.cfg_from_json(open(sys.argv[1]).read()); sl = W.build(cfg); L1 = len(sl.breps)
Z = S4.zq_setup(sl); kr = Z['kreps']; nk = len(kr); q = Z['q']; kor = Z['kor']; lreps = Z['lreps']; ksz = np.bincount(kor, minlength=nk)
contrib = {(j, int(y)): AGs.fibre_patterns(sl, j, int(y))[0].astype(np.int64) @ Z['Mo'][[o for o, _ in sl.parts[j]]][:, kr]
           for j in range(L1) for y in sl.letters[j]}
qj = [int(sl.Gp.co[b, 0] * 2 + sl.Gp.co[b, 1]) for b in sl.breps]; T = np.asarray(sl.T, np.int64); qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]
def count_y(y, r):
    a_of = np.array([0 if qq in (0, r) else 1 for qq in range(4)])[qg]
    tt = np.zeros((2, q), np.int64); np.add.at(tt, (a_of, sl.zco), T); Tt = np.concatenate([tt[0][lreps], tt[1][lreps]])[:-1]
    sig = []
    for k in range(4):
        halves = []
        for a in (0, 1):
            acc = {(0,) * nk: 1}
            for j in range(L1):
                if (qj[j] in (0, r)) != (a == 0): continue
                nacc = defaultdict(int)
                for v0, c0 in acc.items():
                    for v in map(tuple, contrib[(j, int(y[k, j]))].tolist()): nacc[tuple(x + z for x, z in zip(v0, v))] += c0
                acc = nacc
            halves.append(acc)
        dist = defaultdict(int)
        for h0, c0 in halves[0].items():
            X0 = np.array(h0)[kor]
            if (np.array(h0) ** 2 * ksz).sum() > Tt[0]: continue
            for h1, c1 in halves[1].items():
                if (np.array(h0) ** 2 * ksz).sum() + (np.array(h1) ** 2 * ksz).sum() > Tt[0]: continue
                X1 = np.array(h1)[kor]
                P0 = [int((X0 * np.roll(X0, -t)).sum() + (X1 * np.roll(X1, -t)).sum()) for t in lreps]
                P1 = [int((X0 * np.roll(X1, -t)).sum() + (X1 * np.roll(X0, -t)).sum()) for t in lreps]
                dist[tuple((P0 + P1)[:-1])] += c0 * c1
        sig.append(dist)
    A = defaultdict(int)
    for s1, c1 in sig[0].items():
        for s2, c2 in sig[1].items(): A[tuple(a + b for a, b in zip(s1, s2))] += c1 * c2
    tot = 0
    for s3, c3 in sig[2].items():
        for s4, c4 in sig[3].items():
            tot += A.get(tuple(t - a - b for t, a, b in zip(Tt.tolist(), s3, s4)), 0) * c3 * c4
    return tot
Y = np.fromfile(sys.argv[2], np.int8, count=200000 * 4 * L1).reshape(-1, 4, L1)
rng = np.random.default_rng(int(sys.argv[4])); d = tempfile.mkdtemp(); bad = 0
for i in rng.choice(len(Y), int(sys.argv[3]), replace=False):
    y = Y[i].astype(np.int64); yb = os.path.join(d, f'y{i}.bin'); Y[i].tofile(yb); out = []
    for r in (1, 2, 3):
        ex = count_y(y, r); res = SC.run(sl, yb, os.path.join(d, f'o{i}_{r}.jsonl'), r)
        sc = list(res['done'].values())[0]['count']; out.append((ex, sc)); bad += float(ex) != sc
    print(json.dumps(dict(i=int(i), exact_int_vs_screen=[(str(e), s) for e, s in out])), flush=True)
print(json.dumps(dict(mismatches=bad)))
