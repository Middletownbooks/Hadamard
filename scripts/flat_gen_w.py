"""Generate rho-flat blocks for a b = 0 multiplier slice with any number of orbits (W-word masks, scripts/flatw.py)
and reduce them to multiplier-orbit representatives (canonical rotations), per admissible |row sum|.

Same outputs as scripts/flat_probe.py, but masks are (N, W) uint64 arrays (W = ceil(r / 64)):
OUTDIR/orbits_<s>.npy (canonical representatives).  Restart-safe (skips if the orbit files exist).  Prints the raw
count, the number of orbits and the closure size per row sum.
usage: python3 scripts/flat_gen_w.py V D RHO MINUTES WORKERS OUTDIR [--seed S] [--rows a,b,...]"""
import os, sys, json, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import flatw as FW
KJ = FW.KJ


def gen(args):
    v, d, rho, secs, targets, seed = args
    R = KJ.Red(v, d); M = KJ.build_M(R); F0 = int(rho ** 2 * R.m * R.L); W = FW.nwords(R.r)
    rng = np.random.default_rng(seed); buf = np.zeros((1 << 21, W), np.uint64); n = 0; t0 = time.time(); walks = 0; out = []
    while time.time() - t0 < secs:
        y = rng.choice([-1, 1], size=R.r).astype(np.int64)
        n = FW.walk_w(y, M, F0, 300000, max(3, R.r // 10), int(rng.integers(1 << 31)), R.d, np.array(targets, np.int64), buf, n)
        walks += 1
        if n > len(buf) - (1 << 18): out.append(FW.unique_w(buf[:n])); n = 0
    out.append(FW.unique_w(buf[:n]))
    return FW.unique_w(np.concatenate(out)), walks


def main():
    a = sys.argv[1:]
    v, d, rho, minutes, workers, outdir = int(a[0]), int(a[1]), float(a[2]), float(a[3]), int(a[4]), a[5]
    seed0 = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    R = KJ.Red(v, d); e = (v - 1) // d; W = FW.nwords(R.r)
    targets = [int(x) for x in a[a.index('--rows') + 1].split(',')] if '--rows' in a else \
        [s for s in range(1, int((4 * v) ** .5) + 1) if s % (2 * d) in (1, 2 * d - 1)]
    os.makedirs(outdir, exist_ok=True); t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, t=round(time.time() - t0, 1))), flush=True)
    if all(os.path.exists(f'{outdir}/orbits_{s}.npy') for s in targets): log(stage='exists'); return
    with Pool(workers) as p:
        res = p.map(gen, [(v, d, rho, minutes * 60, targets, 1000 * k + int(rho * 1000) + 7919 * seed0) for k in range(workers)])
    allw = FW.unique_w(np.concatenate([x for x, _ in res])); rs = np.abs(FW.rowsum_w(allw, R.r, d))
    log(stage='generated', walks=int(sum(w for _, w in res)), raw=len(allw), W=W, r=R.r, e=e)
    for s in targets:
        x = allw[rs == s]
        if len(x) == 0: np.save(f'{outdir}/orbits_{s}.npy', np.zeros((0, W), np.uint64)); continue
        reps = FW.unique_w(FW.canonical_w(np.ascontiguousarray(x), e))
        np.save(f'{outdir}/orbits_{s}.npy', reps)
        log(stage='pool', rowsum=s, raw=len(x), orbits=len(reps), log2_closed=round(float(np.log2(len(reps) * e)), 2))


if __name__ == '__main__':
    main()
