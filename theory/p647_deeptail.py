r"""Deep D-tail of the 647 x353 cross level on the complete calibration S-solution sets [M against E].

For each calibration joint solution j (logs/p647/calib/j{j}_ssols.npy, the complete S-solution set) every D-candidate
of every S-solution is enumerated (compiled) and the exact number whose D-residual vanishes on the first r cross lags
(nested subsets, r = 1 .. 12, ranks 1 .. 11) is counted; the lattice Gaussian fitted to a sample of the same
candidates' residuals gives the expected count.  The 495 x83 control showed agreement to rank 6 and a collapse at
rank 7-8; this measures where (if anywhere) the 647 tail departs from the Gaussian.
usage: python3 theory/p647_deeptail.py [J ...]"""
import sys, os, json, math, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as PX
from nd_cross import orbit_maps
from q523_gauss import gauss_bits

@njit(cache=True)
def nested_counts(sols, pr0, bd0, n10, pr1, bd1, n11, pr2, bd2, n12, pr3, bd3, n13, P00, P10, P01, P11, P02, P12, P03, P13,
                  TE, order, cnt, smp, every, nsmp_max):
    """cnt[r] += 1 when the residual vanishes on order[0..r] (nested); every 'every'-th candidate's residual -> smp"""
    L = TE.shape[0]; R = order.shape[0]; ns = 0; tot = 0
    r0 = np.empty(L, np.int64); r1 = np.empty(L, np.int64); r2 = np.empty(L, np.int64); v = np.empty(L, np.int64)
    for s in range(sols.shape[0]):
        for x0 in range(bd0[sols[s, 0]], bd0[sols[s, 0] + 1]):
            p = pr0[x0]
            for t in range(L): r0[t] = P00[p // n10, t] + P10[p % n10, t] - TE[t]
            for x1 in range(bd1[sols[s, 1]], bd1[sols[s, 1] + 1]):
                p = pr1[x1]
                for t in range(L): r1[t] = r0[t] + P01[p // n11, t] + P11[p % n11, t]
                for x2 in range(bd2[sols[s, 2]], bd2[sols[s, 2] + 1]):
                    p = pr2[x2]
                    for t in range(L): r2[t] = r1[t] + P02[p // n12, t] + P12[p % n12, t]
                    for x3 in range(bd3[sols[s, 3]], bd3[sols[s, 3] + 1]):
                        p = pr3[x3]
                        for q in range(R):
                            t = order[q]
                            if r2[t] + P03[p // n13, t] + P13[p % n13, t] != 0: break
                            cnt[q] += 1
                        if tot % every == 0 and ns < nsmp_max:
                            for t in range(L): smp[ns, t] = (r2[t] + P03[p // n13, t] + P13[p % n13, t]) // 2
                            ns += 1
                        tot += 1
    return tot, ns

def main(js, log=print, nsmp=300000):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    TE = 2 * (X.TS + X.TD)
    cross = np.array([t for t, h in enumerate(X.lags) if h % 17 and h % 19], np.int64)
    CNT = np.zeros(len(cross), np.int64); TOT = 0; SMP = []
    for j in js:
        path = os.path.join(PC.LOG, "calib", f"j{j}_ssols.npy")
        if not os.path.exists(path): log(json.dumps(dict(j=j, note="no S-solution file"))); continue
        S = np.ascontiguousarray(np.load(path).astype(np.int64))
        B = PX.build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
        cand = int(np.prod([B[k].mult[S[:, k]].astype(np.float64) for k in range(4)], axis=0).sum())
        every = max(1, cand // (nsmp // max(1, len(js))))
        cnt = np.zeros(len(cross), np.int64); smp = np.zeros((nsmp // max(1, len(js)) + 2, len(TE)), np.int64)
        t0 = time.time()
        tot, ns = nested_counts(S, B[0].pairs, B[0].bd, B[0].n1, B[1].pairs, B[1].bd, B[1].n1, B[2].pairs, B[2].bd, B[2].n1,
                                B[3].pairs, B[3].bd, B[3].n1, B[0].P0, B[0].P1, B[1].P0, B[1].P1, B[2].P0, B[2].P1, B[3].P0, B[3].P1,
                                TE, cross, cnt, smp, every, len(smp))
        CNT += cnt; TOT += tot; SMP.append(smp[:ns])
        log(json.dumps(dict(j=j, candidates_log2=round(math.log2(tot), 3), nested_hits=cnt.tolist(), secs=round(time.time() - t0, 1))))
    V = np.concatenate(SMP)
    rows = []
    for q in range(len(cross)):
        sb = cross[:q + 1]
        g = gauss_bits(V[:, sb], np.zeros(len(sb), np.int64))
        exp_ = TOT * 2.0 ** -g['bits'] if g.get('bits') is not None and math.isfinite(g['bits']) else None
        rows.append(dict(r=q + 1, lag=int(X.lags[cross[q]]), rank=g.get('rank'), hits=int(CNT[q]),
                         exact_bits=round(math.log2(TOT / CNT[q]), 2) if CNT[q] else None, gauss_bits=g.get('bits'),
                         expected=round(exp_, 2) if exp_ is not None else None))
        log(json.dumps(rows[-1]))
    log(json.dumps(dict(stage="summary", joint_solutions=list(js), candidates_log2=round(math.log2(TOT), 3), rows=rows)))

if __name__ == "__main__":
    js = [int(x) for x in sys.argv[1:]] or [4208, 9902, 11514, 16846]
    main(js, log=lambda s: print(s, flush=True))
