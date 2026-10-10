r"""S-level calibration AT 571: run the exact S-level join (theory/p571_hsolve.stage2) on random S-class sublists of a
real stage-1 solution y and compare the number of S-solutions with the Gaussian prediction
    E[#S-sol] = 2^{sum_k log2(sub_k) - C3}     (C3 from theory/p571_s2est.level_bits_y; marginal L3 Gaussian)
usage: python3 theory/p571_s3calib.py y.npz index sub [--seeds 3]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p571_hsolve import Inst, stage2, Stage2Block
from p571_s2est import level_bits_y

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--seeds": "3"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    Y = np.load(pos[0])['Y']; idx = int(pos[1]); sub = int(pos[2])
    inst = Inst(570, [103]); y = Y[idx]
    lb = level_bits_y(inst, y, N=12000, seed=idx); C3 = lb['C3'][0]
    ns = [Stage2Block(inst, y[k]).n_classes for k in range(4)]
    used = [min(sub, n) for n in ns]
    pred = sum(math.log2(u) for u in used) - C3
    for s in range(int(opt["--seeds"])):
        t = time.time(); rec, cores = stage2(inst, y, sub=sub, seed=100 + s, tbits=22)
        print(json.dumps(dict(y=idx, seed=s, sclass_log2=[round(math.log2(n), 1) for n in ns], used_log2=[round(math.log2(u), 1) for u in used],
                              C3=C3, predicted_S_log2=round(pred, 1), S_solutions=rec['s_solutions'],
                              actual_S_log2=(round(math.log2(rec['s_solutions']), 1) if rec['s_solutions'] else None),
                              d_solutions=rec['d_solutions'], cores=rec['cores'], secs=round(time.time() - t, 1))), flush=True)
