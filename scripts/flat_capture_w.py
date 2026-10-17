"""Capture-recapture flat-set size (multiplier orbits) from two independent generations, for 1-D or (N, W) masks.
Same estimate as scripts/flat_capture.py (Lincoln-Petersen N = n1 n2 / m orbits, about N e blocks; a lower estimate,
since tabu walks reach some orbits more easily).  usage: python3 scripts/flat_capture_w.py V D DIR1 DIR2"""
import os, sys, glob, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flatw as FW


def main():
    v, d, d1, d2 = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]; e = (v - 1) // d
    for f in sorted(glob.glob(os.path.join(d1, 'orbits_*.npy'))):
        s = f.split('_')[-1][:-4]; g = os.path.join(d2, f'orbits_{s}.npy')
        if not os.path.exists(g): continue
        a, b = np.load(f), np.load(g); a = a.reshape(len(a), -1); b = b.reshape(len(b), -1)
        if not len(a) or not len(b): continue
        m = int(FW.isin_w(a, b).sum()); N = len(a) * len(b) / m if m else float('inf')
        print(json.dumps(dict(rowsum=int(s), n1=len(a), n2=len(b), both=m, orbits_est=round(N) if m else None,
                              log2_blocks_est=round(float(np.log2(N * e)), 2) if m else None)))


if __name__ == '__main__':
    main()
