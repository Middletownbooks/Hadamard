r"""Calibration of the last level (order-18 characters given S and D3) [E from measured residuals].
Runs the S-step of theory/q523_wagner.py on a y, then for every S-solution and every D3-solution enumerates all
D-candidates (configurations inside the four D3 groups), records the full D-level signature sums, and fits the
lattice-exact Gaussian to them: log2 P(full D-level match | D3 match).  Together with the measured D3-solutions per
S-solution this gives the D-success probability per S-solution (the D-deficit) without the unconditional model.
usage: python3 theory/q523_dfine_calib.py F q s y.json [--beta 1.5] [--rows 64] [--K 24] [--logN 28] [--seed 0]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q523_wagner as W
from q523_dcheck import DChecker, hjoin4
from q523_gauss import gauss_bits
from q523_sdsolve import SDInst

def candidates_sums(dc, pidx4, cap=4000):
    B = [dc.block(k, pidx4[k]) for k in range(4)]
    h3, n3 = hjoin4(B[0]['k3'], B[1]['k3'], B[2]['k3'], B[3]['k3'], dc.k3T, 16, 1 << 20)
    out = []; ncand = 0; nd3 = 0
    for q in h3:
        if not (B[0]['sig3'][q[0]] + B[1]['sig3'][q[1]] + B[2]['sig3'][q[2]] + B[3]['sig3'][q[3]] == dc.T3).all(): continue
        nd3 += 1
        idx = [B[k]['order'][B[k]['start'][q[k]]:B[k]['start'][q[k] + 1]] for k in range(4)]
        n = int(np.prod([len(x) for x in idx])); ncand += n
        if len(out) < cap:
            rng = np.random.default_rng(len(out))
            for _ in range(min(n, 8)):
                cfg = [x[rng.integers(len(x))] for x in idx]
                out.append(sum(B[k]['sigf'][cfg[k]] for k in range(4)))
    return out, nd3, ncand, sum(math.log2(b['n']) for b in B)

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--beta": "1.5", "--rows": "64", "--K": "24", "--logN": "28", "--seed": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y)
    sols = []
    W.DChecker = None
    # capture S-solutions by monkeypatching the D-check stage: run the S-step only
    beta, nrow, K, logN, seed = float(opt["--beta"]), int(opt["--rows"]), int(opt["--K"]), float(opt["--logN"]), int(opt["--seed"])
    class Grab:
        def __init__(self, *a, **k): pass
        def check(self, pidx4): sols.append(pidx4); return [], dict(d3_solutions=0)
    W.DChecker = Grab
    W.run(inst, beta, nrow, K, logN, 1, seed, 20, None, lambda d: print(json.dumps(d), flush=True))
    dc = DChecker(inst, seed); allv = []; nd3s = []; ncs = []; hds = []
    for p4 in sols:
        v, nd3, nc, hd = candidates_sums(dc, p4); allv += v; nd3s.append(nd3); ncs.append(nc); hds.append(hd)
    V = np.array(allv); g = gauss_bits(V, dc.TF)
    tot = [math.log2(max(1, c)) for c in ncs]
    print(json.dumps(dict(stage="D-fine calibration", s_solutions=len(sols), mean_hD=round(float(np.mean(hds)), 2),
                          d3_per_S_log2=round(math.log2(max(1e-9, np.mean(nd3s))), 2), candidates_per_S_log2=round(math.log2(max(1e-9, np.mean(ncs))), 2),
                          fine_bits=g['bits'], mahal=g.get('mahal'), samples=len(V),
                          D_success_per_S_log2=round(math.log2(max(1e-9, np.mean(ncs))) - g['bits'], 2))))
