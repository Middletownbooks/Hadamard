"""Capture-recapture size of a flat set, at the level of multiplier orbits, from two independent probes.

n1, n2 = orbit representatives found by two independently seeded generations (scripts/flat_probe.py --seed),
m = orbits found by both; Lincoln-Petersen N = n1 n2 / m orbits, i.e. about N e blocks per row sum.  Tabu walks
do not sample the flat set uniformly (some orbits are easier to reach), which biases N low; treat it as a lower
estimate.  usage: python3 scripts/flat_capture.py V D DIR1 DIR2"""
import os, sys, glob, json
import numpy as np


def main():
    v, d, d1, d2 = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]; e = (v - 1) // d
    for f in sorted(glob.glob(os.path.join(d1, 'orbits_*.npy'))):
        s = f.split('_')[-1][:-4]; g = os.path.join(d2, f'orbits_{s}.npy')
        if not os.path.exists(g): continue
        a, b = np.load(f), np.load(g); m = int(np.isin(a, b, assume_unique=True).sum())
        N = len(a) * len(b) / m if m else float('inf')
        print(json.dumps(dict(rowsum=int(s), n1=len(a), n2=len(b), both=m, orbits_est=round(N), log2_blocks_est=round(float(np.log2(N * e)), 2) if m else None)))


if __name__ == '__main__':
    main()
