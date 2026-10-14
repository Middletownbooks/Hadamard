r"""Energy-bucket selection for the 647 x353 S-join, evaluated exactly on complete S-solution sets [M].

Within a left-energy class EL the left pairs split into buckets by (e_a, e_b = EL - e_a) and the right pairs by
(e_c, e_d = T0 - EL - e_c).  Keeping a set L of left buckets and R of right buckets costs exactly
sum_L n_a(e_a) n_b(e_b) + sum_R n_c(e_c) n_d(e_d) emitted pairs and yields exactly the D-candidates of the S-solutions
in L x R.  Over all (EL, left bucket, right bucket) cells of a joint solution this script computes the efficiency
frontier by greedy removal of the bucket whose removal raises yield / work the most.
usage: python3 theory/p647_buckets.py J [J ...]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as P
from nd_cross import orbit_maps

def frontier(blocks, S, order, E0=645):
    a, b, c, d = order
    e = [B.Ssig[:, 0] for B in blocks]
    cnt = [np.bincount(x, minlength=E0 + 1)[:E0 + 1].astype(np.float64) for x in e]
    y = np.prod([blocks[k].mult[S[:, k]].astype(np.float64) for k in range(4)], axis=0)
    ea = e[a][S[:, a]]; eb = e[b][S[:, b]]; ec = e[c][S[:, c]]
    EL = ea + eb
    # left bucket id = (EL, ea), right bucket id = (EL, ec); work of every bucket with any pair
    Lw = {}; Rw = {}
    for EL_ in range(E0 + 1):
        for x in np.nonzero(cnt[a])[0]:
            yb = EL_ - x
            if 0 <= yb <= E0 and cnt[b][yb] > 0: Lw[(EL_, int(x))] = cnt[a][x] * cnt[b][yb]
        ER = E0 - EL_
        for x in np.nonzero(cnt[c])[0]:
            yd = ER - x
            if 0 <= yd <= E0 and cnt[d][yd] > 0: Rw[(EL_, int(x))] = cnt[c][x] * cnt[d][yd]
    # drop buckets whose class has no partner side
    ELs_L = {k[0] for k in Lw}; ELs_R = {k[0] for k in Rw}
    Lw = {k: v for k, v in Lw.items() if k[0] in ELs_R}; Rw = {k: v for k, v in Rw.items() if k[0] in ELs_L}
    Y = {}
    for i in range(len(S)):
        key = (int(EL[i]), int(ea[i]), int(ec[i])); Y[key] = Y.get(key, 0.0) + y[i]
    Lk = sorted(Lw); Rk = sorted(Rw); li = {k: n for n, k in enumerate(Lk)}; ri = {k: n for n, k in enumerate(Rk)}
    M = np.zeros((len(Lk), len(Rk)))
    for (EL_, xa, xc), v in Y.items(): M[li[(EL_, xa)], ri[(EL_, xc)]] += v
    wl = np.array([Lw[k] for k in Lk]); wr = np.array([Rw[k] for k in Rk])
    keepL = np.ones(len(Lk), bool); keepR = np.ones(len(Rk), bool)
    Ytot = M.sum(); Wtot = wl.sum() + wr.sum()
    pts = []
    while True:
        Ycur = M[np.ix_(keepL, keepR)].sum(); Wcur = wl[keepL].sum() + wr[keepR].sum()
        if Ycur <= 0: break
        pts.append((Wcur / Wtot, Ycur / Ytot, (Ycur / Wcur) / (Ytot / Wtot), int(keepL.sum()), int(keepR.sum())))
        best = None
        rowY = M[:, keepR].sum(1); colY = M[keepL, :].sum(0)
        for i in np.nonzero(keepL)[0]:
            r = (Ycur - rowY[i]) / max(Wcur - wl[i], 1)
            if Ycur - rowY[i] > 0 and (best is None or r > best[0]): best = (r, 'L', i)
        for i in np.nonzero(keepR)[0]:
            r = (Ycur - colY[i]) / max(Wcur - wr[i], 1)
            if Ycur - colY[i] > 0 and (best is None or r > best[0]): best = (r, 'R', i)
        if best is None: break
        if best[1] == 'L': keepL[best[2]] = False
        else: keepR[best[2]] = False
    return pts, Lk, Rk, M, wl, wr

def main(js, log=print):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
    orders = {json.loads(l)['j']: json.loads(l)['join_order'] for l in open(os.path.join(PC.LOG, "calib", "calib.jsonl"))}
    for j in js:
        blocks = P.build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
        S = np.load(os.path.join(PC.LOG, "calib", f"j{j}_ssols.npy"))
        pts, Lk, Rk, M, wl, wr = frontier(blocks, S, orders.get(j, [0, 1, 2, 3]))
        sel = []
        for f in [0.9, 0.75, 0.5, 0.35, 0.25, 0.15, 0.1, 0.05]:
            p = min(pts, key=lambda t: abs(t[1] - f))
            sel.append(dict(yield_frac=round(p[1], 3), work_frac=round(p[0], 4), gain=round(p[2], 2), left_buckets=p[3], right_buckets=p[4]))
        log(json.dumps(dict(j=j, buckets=[len(Lk), len(Rk)], frontier=sel)))

if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]], log=lambda s: print(s, flush=True))
