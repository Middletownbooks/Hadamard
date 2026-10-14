r"""Checks on a control slice (first: 495 x83): plants through the full pipeline, and the D-tail [M against E].

(a) plants: a random (L0, L1) pair per block of a random joint solution (y_B given); targets = its own S-sums and
    even-lag sums; the full S-join (every energy class) and the D-check must return it;
(b) D-tail: for a sample of joint solutions with the real target, every D-candidate of every S-solution is enumerated
    (theory/p647_calib.dtail); the exact fraction hitting the D-target on subsets of the cross lags (and on all of
    them = the cores) is compared with the lattice Gaussian fitted to the pooled residuals.
usage: python3 theory/p647_ctrl_check.py m mu qA qB --yb I [--plants 10] [--tail 40]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_cross import Cross, orbit_maps
from nd_joint_exhaust import enumerate_yB
import p647_common as PC
import p647_cross as PX
import p647_calib as CAL
from q523_gauss import gauss_bits

def main(m, mu, qA, qB, ib, nplants=10, ntail=40, seed=0, log=print):
    d = os.path.join(PC.LOG, "control"); tag = f"{m}_{mu}_{qA}_{qB}"
    X = Cross(m, mu, qA, qB)
    ys = enumerate_yB(m, mu, qA, qB, os.path.join(d, f"enum_{tag}.jsonl"), 1, log=lambda s: None)
    yb = ys[ib]; rec, sols = X.J.test(yb, maxsol=10 ** 7)
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    rng = np.random.default_rng(seed)
    ok = 0
    for t in range(nplants):
        j = int(rng.integers(len(sols)))
        blocks = PX.build(X, yb, X.J.ya_full(sols[j]['layers']), maps, NH, oi)
        pr, TS, TE = PX.plant(X, blocks, rng)
        r, Ss, hits, cand = PX.solve(X, blocks, TS=TS, TE=TE, P=61)
        found = tuple(p for p, s in pr) in set(tuple(int(v) for v in h[1:]) for h in hits)
        ok += found
        log(json.dumps(dict(plant=t, jsol=j, S=r['S_solutions'], cand=r['S_pairs'], D_hits=r['D_hits'], plant_found=found)))
    log(json.dumps(dict(plants=nplants, found=ok)))
    # D-tail over a sample of joint solutions, pooled
    cross = [t for t, h in enumerate(X.lags) if h % qA and h % qB]
    subsets = [[t] for t in cross[:3]] + [cross[:2], cross[:3], cross[:4], cross[:5], cross[:6]] + \
              [[t for t in cross if t != u] for u in cross] + [cross]          # leave-one-out subsets, then all
    cnt = np.zeros(len(subsets), np.int64); tot = 0; smp = []; cores = 0
    pick = rng.choice(len(sols), size=min(ntail, len(sols)), replace=False)
    for j in pick:
        blocks = PX.build(X, yb, X.J.ya_full(sols[int(j)]['layers']), maps, NH, oi)
        r, Ss, hits, cand = PX.solve(X, blocks, P=61)
        cores += r['D_hits']
        if not Ss: continue
        dt = CAL.dtail(X, blocks, Ss, subsets, nsmp=5000)
        tot += dt['candidates']
        for q, row in enumerate(dt['tail']): cnt[q] += row['hits']
        # residual sample for the pooled Gaussian: recompute a small sample directly
        S = np.asarray(Ss, np.int64).reshape(-1, 4)
        TEv = 2 * (X.TS + X.TD)
        wS = np.prod([blocks[k].mult[S[:, k]].astype(np.float64) for k in range(4)], axis=0); wS /= wS.sum()
        for _ in range(min(3000, 40 * len(S))):
            s = S[rng.choice(len(S), p=wS)]; v = -TEv.copy()     # uniform over D-candidates
            for k, B in enumerate(blocks):
                p = int(B.pairs[rng.integers(B.bd[s[k]], B.bd[s[k] + 1])]); i, jj = divmod(p, B.n1)
                v = v + B.P0[i] + B.P1[jj]
            smp.append(v)
    V = np.array(smp)
    out = []
    for q, sb in enumerate(subsets):
        g = gauss_bits(V[:, sb], np.zeros(len(sb), np.int64))
        out.append(dict(lags=[int(X.lags[t]) for t in sb], hits=int(cnt[q]), exact_bits=round(math.log2(tot / cnt[q]), 2) if cnt[q] else None,
                        gauss_bits=g.get('bits'), rank=g.get('rank'),
                        expected_hits=round(tot * 2.0 ** -g['bits'], 2) if g.get('bits') is not None and math.isfinite(g['bits']) else None))
    log(json.dumps(dict(stage="dtail", joint_solutions=len(pick), candidates=tot, candidates_log2=round(math.log2(tot), 2) if tot else None,
                        cores=cores, tail=out)))

if __name__ == "__main__":
    a = sys.argv[1:]
    pos = [int(x) for x in a[:4]]
    ib = int(a[a.index("--yb") + 1])
    npl = int(a[a.index("--plants") + 1]) if "--plants" in a else 10
    nt = int(a[a.index("--tail") + 1]) if "--tail" in a else 40
    main(*pos, ib, npl, nt, log=lambda s: print(s, flush=True))
