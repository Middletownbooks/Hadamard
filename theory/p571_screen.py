r"""Screen stage-1 solutions y of the 571 x27 slice by their predicted stage-2 cost per core [E].

For each y: the exact S-class list sizes 2^{h3_k} (theory/p571_hsolve.Stage2Block) give the S-level join cost
max(h3_a + h3_b, h3_c + h3_d) over the best pairing (log2 pairs per side); the stage-2 yield per y is the Gaussian of
theory/p571_ternary.stage2 (Z190 lags, exact index-3 lattice; calibrated at v = 79: it OVERestimates exact per-y core
counts by 2^1.5-2^3.5, theory/logs/p571/calib_v79.log).  Expected work per core:
    work = join - min(yield, 0) + (yield > 0 ? -yield / 2 : 0)      (log2; subsampling when yield > 0)
usage: python3 theory/p571_screen.py y.npz [y2.npz ...] [--N 8000] [--top 20]"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p571_ternary as T
from p571_hsolve import Inst, Stage2Block

def main(files, N=8000, top=20, m=570, gens=(103,)):
    inst = Inst(m, list(gens)); TS = T.TernSlice(m, list(gens)); q = inst.q; S = inst.S
    def to_tern(y):
        out = np.zeros((4, TS.R), np.int64)
        for j, (oa, ob, fr) in enumerate(TS.cls):
            g = oa[0]; out[:, j] = y[:, g % 2, S.oid[g % q]]
        return out
    Ys = np.concatenate([np.load(f)['Y'] for f in files])
    Ys = np.unique(Ys.reshape(len(Ys), -1), axis=0).reshape(-1, 4, 2, inst.R)
    rows = []
    for i, y in enumerate(Ys):
        h3 = sorted(Stage2Block(inst, y[k]).log2n for k in range(4))
        join = max(h3[3] + h3[0], h3[2] + h3[1])
        r2 = T.stage2(TS, to_tern(y), N=N, seed=i)
        yl = r2['stage2_log2E_per_y']
        work = join - min(yl, 0) - (yl / 2 if yl > 0 else 0)
        rows.append(dict(i=i, h3=round(sum(h3), 1), join=round(join, 1), H2=r2['H2'], C2=r2['C2'], yield_log2=yl, work_log2=round(work, 1)))
        if (i + 1) % 100 == 0: print(json.dumps(dict(screened=i + 1)), flush=True)
    rows.sort(key=lambda r: r['work_log2'])
    w = np.array([r['work_log2'] for r in rows]); yl = np.array([r['yield_log2'] for r in rows])
    print(json.dumps(dict(summary=True, ys=len(rows), work_min=float(w.min()), work_median=float(np.median(w)),
                          yield_max=float(yl.max()), yield_median=float(np.median(yl)))), flush=True)
    for r in rows[:top]: print(json.dumps(r), flush=True)
    np.savez_compressed(os.path.join(os.path.dirname(files[0]), 'screen_best.npz'), Y=Ys[[r['i'] for r in rows[:top]]],
                        work=np.array([r['work_log2'] for r in rows[:top]]))

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "8000", "--top": "20"}; files = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: files.append(a[i]); i += 1
    main(files, N=int(opt["--N"]), top=int(opt["--top"]))
