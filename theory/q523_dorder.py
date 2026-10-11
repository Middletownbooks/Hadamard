r"""Is the D-check worth ordering by D-freedom? (report/paley_523_answer.md, section 3.4) [M on real S-solutions].

Runs the S-step of theory/q523_wagner.py, then for every S-solution records
    hd    = log2 #D-configurations given S (sum over blocks of log2 of the option products),
    nd3   = D3-solutions (Z3-level matches of D), ncand = D-candidates (configurations inside matched D3 groups),
    secs  = time of the compiled D-check (theory/q523_dcheck_nb.py).
The expected number of cores from an S-solution is ncand * 2^-(order-18 bits) (section 3.4), so the yield of the top
fraction of S-solutions ranked by hd, against the D-check time they take, says whether skipping low-hd S-solutions
pays.  Prints the per-S table (JSON lines) and the capture curve.
usage: python3 theory/q523_dorder.py F q s y.json [--logN 28] [--seed 0] [--beta 1] [--rows 64] [--K 24]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q523_wagner as W
from q523_dcheck import DChecker
from q523_dcheck_nb import DCheckNB
from q523_dfine_calib import candidates_sums
from q523_sdsolve import SDInst

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--beta": "1", "--rows": "64", "--K": "24", "--logN": "28", "--seed": "0"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, q, s = map(int, pos[:3]); Y = json.load(open(pos[3])); Y = Y['y'] if isinstance(Y, dict) else Y
    inst = SDInst(F, q, s, Y); sols = []
    class Grab:
        def __init__(self, *a, **k): pass
        def check(self, pidx4): sols.append(pidx4); return [], dict(d3_solutions=0)
    W.DChecker = Grab
    W.run(inst, float(opt["--beta"]), int(opt["--rows"]), int(opt["--K"]), float(opt["--logN"]), 1, int(opt["--seed"]), 20,
          None, lambda d: print(json.dumps(d), flush=True))
    dc = DChecker(inst, 1); nb = DCheckNB(inst, 1); rows = []
    nb.check(sols[0])                                                   # compile
    for p4 in sols:
        t0 = time.time(); _, info = nb.check(p4); dt = time.time() - t0
        _, nd3, nc, hd = candidates_sums(dc, p4)
        rows.append(dict(hd=round(hd, 2), nd3=nd3, ncand=nc, secs=round(dt, 5)))
        print(json.dumps(rows[-1]), flush=True)
    R = sorted(rows, key=lambda r: -r['hd']); tc = sum(r['ncand'] for r in R); ts = sum(r['secs'] for r in R)
    hd = np.array([r['hd'] for r in rows]); lc = np.log2(np.array([max(1, r['ncand']) for r in rows]))
    curve = []
    for frac in (0.1, 0.25, 0.5, 0.75, 1.0):
        n = max(1, int(round(frac * len(R))))
        curve.append(dict(top=frac, cand_share=round(sum(r['ncand'] for r in R[:n]) / max(1, tc), 3),
                          time_share=round(sum(r['secs'] for r in R[:n]) / ts, 3)))
    print(json.dumps(dict(stage="summary", s_solutions=len(rows), hd_mean=round(float(hd.mean()), 2),
                          hd_sd=round(float(hd.std()), 2), corr_hd_log2cand=round(float(np.corrcoef(hd, lc)[0, 1]), 3),
                          mean_secs=round(ts / len(rows), 5), curve=curve)), flush=True)
