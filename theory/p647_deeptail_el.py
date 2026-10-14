r"""Deep D-tail of the 647 x353 calibration runs split by energy class [M against E]: the two central left-energy classes
(the production choice, theory/p647_prod.central_els) against the others.  For each group: exact hits on nested
cross-lag subsets (theory/p647_deeptail.nested_counts) and the Gaussian fitted to that group's own residuals.
usage: python3 theory/p647_deeptail_el.py [J ...]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p647_common as PC
import p647_cross as PX
from p647_deeptail import nested_counts
from p647_prod import central_els
from nd_cross import orbit_maps
from q523_gauss import gauss_bits

def main(js, log=print, nsmp=150000):
    X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
    reps, oi, maps = orbit_maps(X); NH = PX.lag_tensor(X, oi)
    TE = 2 * (X.TS + X.TD)
    cross = np.array([t for t, h in enumerate(X.lags) if h % 17 and h % 19], np.int64)
    orders = {json.loads(l)['j']: json.loads(l)['join_order'] for l in open(os.path.join(PC.LOG, "calib", "calib.jsonl"))}
    G = {g: dict(cnt=np.zeros(len(cross), np.int64), tot=0, smp=[]) for g in ("central", "other")}
    for j in js:
        S = np.load(os.path.join(PC.LOG, "calib", f"j{j}_ssols.npy")).astype(np.int64)
        B = PX.build(X, ys[0], X.J.ya_full(sols[j]['layers']), maps, NH, oi)
        order = orders.get(j, [0, 1, 2, 3]); a, b = order[0], order[1]
        els, centre = central_els(B, [Bk.Ssig for Bk in B], order, 2)
        EL = B[a].Ssig[S[:, a], 0] + B[b].Ssig[S[:, b], 0]
        for g, mask in (("central", np.isin(EL, els)), ("other", ~np.isin(EL, els))):
            Sg = np.ascontiguousarray(S[mask])
            if not len(Sg): continue
            cand = int(np.prod([B[k].mult[Sg[:, k]].astype(np.float64) for k in range(4)], axis=0).sum())
            every = max(1, cand // (nsmp // len(js))); cnt = np.zeros(len(cross), np.int64)
            smp = np.zeros((nsmp // len(js) + 2, len(TE)), np.int64)
            tot, ns = nested_counts(Sg, B[0].pairs, B[0].bd, B[0].n1, B[1].pairs, B[1].bd, B[1].n1, B[2].pairs, B[2].bd, B[2].n1,
                                    B[3].pairs, B[3].bd, B[3].n1, B[0].P0, B[0].P1, B[1].P0, B[1].P1, B[2].P0, B[2].P1, B[3].P0, B[3].P1,
                                    TE, cross, cnt, smp, every, len(smp))
            G[g]['cnt'] += cnt; G[g]['tot'] += tot; G[g]['smp'].append(smp[:ns])
        log(json.dumps(dict(j=j, central_els=els, centre=centre)))
    for g, d in G.items():
        V = np.concatenate(d['smp']); rows = []
        for q in range(len(cross)):
            gb = gauss_bits(V[:, cross[:q + 1]], np.zeros(q + 1, np.int64))
            rows.append(dict(rank=gb.get('rank'), hits=int(d['cnt'][q]), gauss_bits=gb.get('bits'),
                             expected=round(d['tot'] * 2.0 ** -gb['bits'], 2) if gb.get('bits') is not None and math.isfinite(gb['bits']) else None))
        log(json.dumps(dict(group=g, candidates_log2=round(math.log2(d['tot']), 3), rows=rows[4:])))

if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or [4208, 9902, 11514, 16846], log=lambda s: print(s, flush=True))
