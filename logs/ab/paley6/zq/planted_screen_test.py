"""Planted test of the Z_q-level screens (theory/ab_s2_zq_screen.py logic): random lifts x of random letter
quadruples in a slice, target T' := sum_k PAF_G(x_k); the exact count for y = compress(x) must be >= 1 for every
character (no false kills).  usage: python3 planted_screen_test.py CFG.json N"""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../theory'))
import numpy as np, ab_sweep as W, ab_s1sd as SD, ab_gauss as AGs, ab_s2_s4 as S4
CHARS = [(1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)]
cfg = SD.cfg_from_json(open(sys.argv[1]).read()); sl = W.build(cfg); L1 = len(sl.breps)
Z = S4.zq_setup(sl); kr = Z['kreps']; nl = len(Z['lreps'])
qj = [int(sl.Gp.co[b, 0] * 2 + sl.Gp.co[b, 1]) for b in sl.breps]; qg = sl.G.co[:, 0] * 2 + sl.G.co[:, 1]
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1); bad = 0
for trial in range(int(sys.argv[2])):
    lab = rng.choice(np.array([-1, 1], np.int8), size=(4, sl.r)); X = lab[:, sl.orb].astype(np.int64)
    Tp = sl.G.paf(X).sum(0)
    Y = np.zeros((4, L1), np.int64)
    for j, bb in enumerate(sl.breps): Y[:, j] = X[:, np.nonzero(sl.bidx == bb)[0]].sum(1)
    out = []
    for p in range(4):
        tt = np.zeros(Z['q'], np.int64); np.add.at(tt, sl.zco, np.array(CHARS[p])[qg] * Tp); T = tt[Z['lreps']][:nl - 1]
        D = []
        for k in range(4):
            pats = [AGs.fibre_patterns(sl, j, int(Y[k, j]))[0] for j in range(L1)]
            st = np.zeros((1, len(kr)), np.int64); cnt = np.ones(1, np.int64)
            for j in range(L1):
                os_ = [o for o, _ in sl.parts[j]]
                V = CHARS[p][qj[j]] * (pats[j].astype(np.int64) @ Z['Mo'][os_][:, kr])
                nst = (st[:, None, :] + V[None, :, :]).reshape(-1, len(kr)); ncnt = np.repeat(cnt, len(V))
                u, inv = np.unique(nst, axis=0, return_inverse=True); c2 = np.zeros(len(u), np.int64); np.add.at(c2, inv.ravel(), ncnt); st, cnt = u, c2
            Pq = S4.zq_paf(st[:, Z['kor']], Z['lreps'])[:, :nl - 1]
            u, inv = np.unique(Pq, axis=0, return_inverse=True); c = np.zeros(len(u)); np.add.at(c, inv.ravel(), cnt.astype(float)); D.append((u, c))
        # brute 4-fold convolution on Python dicts (independent of the screen's encoding)
        A = {}
        for s1, c1 in zip(map(tuple, D[0][0].tolist()), D[0][1]):
            for s2, c2 in zip(map(tuple, D[1][0].tolist()), D[1][1]):
                k_ = tuple(a + b for a, b in zip(s1, s2)); A[k_] = A.get(k_, 0) + c1 * c2
        ex = 0.0; Tt = tuple(T.tolist())
        for s3, c3 in zip(map(tuple, D[2][0].tolist()), D[2][1]):
            for s4, c4 in zip(map(tuple, D[3][0].tolist()), D[3][1]):
                need = tuple(t - a - b for t, a, b in zip(Tt, s3, s4)); ex += A.get(need, 0) * c3 * c4
        out.append(round(math.log2(ex), 2) if ex > 0 else None); bad += ex == 0
    print(json.dumps(dict(trial=trial, log2_counts=out)), flush=True)
print(json.dumps(dict(trials=int(sys.argv[2]), false_kills=int(bad))))
