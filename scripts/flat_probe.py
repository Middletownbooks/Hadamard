"""Probe the flat set of a b = 0 multiplier slice: generation yield, coverage, multiplier closure, constraint bits.

Blocks are U-invariant +-1 sequences on Z_v (v prime, |U| = d), packed as orbit-sign bits (bit 0 = orbit {0}, bit
1 + j = orbit g^j U, j < e = (v-1)/d; needs e + 1 <= 64).  A block is rho-flat if F = sum over lag classes of PAF^2
is <= rho^2 * v * L (the kt_join_kaggle convention).

Multiplier closure [P]: x -> x(g^s y) maps U-invariant blocks to U-invariant blocks, rotates the e orbit bits by s,
keeps the row sum, and permutes the lag classes (all of size 2d), so it keeps F: the flat set is a union of rotation
orbits.  Rotation by e/2 is the reversal y -> -y.  So (i) the fraction of rotations of pool members that are already
in the pool estimates the coverage of the flat set, and (ii) closing a pool under rotations enlarges it for free.

usage: python3 scripts/flat_probe.py V D RHO MINUTES WORKERS OUTDIR [--profile a,b,c,d] [--P 4096] [--kmax 10] [--seed S]
  writes OUTDIR/pool_<rowsum>.npy (raw, deduplicated) and OUTDIR/closed_<rowsum>.npy (rotation-closed),
  prints JSON lines: yields, coverage, closure sizes, Gaussian n' and partial-join bits for the profile."""
import os, sys, json, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kaggle'))
import kt_join_kaggle as KJ                                   # walk, Red, build_M, sparse_N, paf_from_bits (no main on import)


def rotate(w, s, e):
    """rotate the e orbit bits (bits 1..e) of packed blocks w by s: new orbit j = old orbit j + s"""
    s %= e
    if s == 0: return w.copy()
    me = np.uint64((1 << e) - 1)
    f = (w >> np.uint64(1)) & me
    g = ((f >> np.uint64(s)) | (f << np.uint64(e - s))) & me
    return (w & np.uint64(1)) | (g << np.uint64(1))


def canonical(w, e):
    """least rotation of each packed block (orbit representative under the multipliers), O(N) memory"""
    best = w.copy(); cur = w
    for s in range(1, e):
        cur = rotate(w, s, e); np.minimum(best, cur, out=best)
    return best


def closure_sample(canon_reps, e, n, rng):
    """n blocks drawn uniformly from the rotation closure of the given orbit representatives (aperiodic orbits)"""
    i = rng.integers(0, len(canon_reps), n); s = rng.integers(0, e, n)
    out = np.empty(n, np.uint64)
    for s_ in range(e):
        m = s == s_
        if m.any(): out[m] = rotate(canon_reps[i[m]], s_, e)
    return np.unique(out)


def rowsum(w, r, d):
    ones = np.zeros(len(w), np.int64)
    for q in range(1, r): ones += ((w >> np.uint64(q)) & np.uint64(1)).astype(np.int64)
    return np.where(w & np.uint64(1), 1, -1) + d * (2 * ones - (r - 1))


def gen(args):
    v, d, rho, secs, targets, seed = args
    R = KJ.Red(v, d); M = KJ.build_M(R); F0 = int(rho ** 2 * R.m * R.L)
    rng = np.random.default_rng(seed); buf = np.zeros(1 << 23, np.uint64); n = 0; t0 = time.time(); walks = 0
    out = []
    while time.time() - t0 < secs:
        y = rng.choice([-1, 1], size=R.r).astype(np.int64)
        n = KJ.walk(y, M, F0, 300000, max(3, R.r // 10), int(rng.integers(1 << 31)), R.d, np.array(targets, np.int64), buf, n)
        walks += 1
        if n > len(buf) - (1 << 20): out.append(np.unique(buf[:n])); n = 0
    out.append(np.unique(buf[:n]))
    return np.unique(np.concatenate(out)), walks


def paf_classes(w, R, Nsp):
    return KJ.paf_from_bits(np.ascontiguousarray(w), *Nsp, R.L).astype(np.int64)


def main():
    a = sys.argv[1:]
    v, d, rho, minutes, workers, outdir = int(a[0]), int(a[1]), float(a[2]), float(a[3]), int(a[4]), a[5]
    prof = [int(x) for x in a[a.index('--profile') + 1].split(',')] if '--profile' in a else None
    P = int(a[a.index('--P') + 1]) if '--P' in a else 4096
    kmax = int(a[a.index('--kmax') + 1]) if '--kmax' in a else 10
    seed0 = int(a[a.index('--seed') + 1]) if '--seed' in a else 0
    os.makedirs(outdir, exist_ok=True)
    R = KJ.Red(v, d); e = (v - 1) // d; assert e + 1 <= 64 and R.r == e + 1
    Nsp = KJ.sparse_N(R); t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, t=round(time.time() - t0, 1))), flush=True)
    targets = sorted(set(prof)) if prof else [s for s in range(1, int((4 * v) ** .5) + 1) if s % (2 * d) in (1, 2 * d - 1)]
    # 1. generation (skipped if the raw pools exist: restart-safe)
    if not all(os.path.exists(f'{outdir}/pool_{s}.npy') for s in targets):
        with Pool(workers) as p:
            res = p.map(gen, [(v, d, rho, minutes * 60, targets, 1000 * k + int(rho * 1000) + 7919 * seed0) for k in range(workers)])
        allw = np.unique(np.concatenate([x for x, _ in res])); rs = np.abs(rowsum(allw, R.r, d))
        for s in targets: np.save(f'{outdir}/pool_{s}.npy', allw[rs == s])
        log(stage='generated', walks=int(sum(w for _, w in res)), minutes=minutes, workers=workers)
    pools = {s: np.load(f'{outdir}/pool_{s}.npy') for s in targets}
    # 2. F invariance under rotation (spot check) and coverage / closure
    smp = pools[max(pools, key=lambda s: len(pools[s]))][:200]
    F = lambda w: (paf_classes(w, R, Nsp) ** 2).sum(1)
    assert all(np.array_equal(F(smp), F(rotate(smp, s_, e))) for s_ in (1, 7, e // 2)), 'F not rotation invariant'
    for s in targets:
        x = pools[s]
        if len(x) == 0: log(stage='pool', rowsum=s, size=0); continue
        sx = np.sort(x); idx = np.random.default_rng(1).choice(len(x), min(len(x), 20000), replace=False); smp = x[idx]
        cov_rev = float(np.isin(rotate(smp, e // 2, e), sx).mean())
        reps = np.unique(canonical(x, e)); np.save(f'{outdir}/orbits_{s}.npy', reps)      # one block per multiplier orbit
        closed_size = len(reps) * e                                                       # aperiodic orbits (period < e is rare)
        log(stage='pool', rowsum=s, raw=len(x), log2_raw=round(float(np.log2(len(x))), 2), coverage_reversal=round(cov_rev, 3),
            orbits=len(reps), coverage_rotations=round(len(x) / closed_size, 4), log2_closed=round(float(np.log2(closed_size)), 2))
    if not prof: return
    # 3. Gaussian n' and partial joins on the profile (rotation-closed pools when available)
    rng = np.random.default_rng(0)
    use = {}
    for s in set(prof):                                   # samples of the rotation-closed pools (uniform over the closure)
        reps = np.load(f'{outdir}/orbits_{s}.npy'); use[s] = closure_sample(reps, e, 60000, rng)
    sizes = {s: len(np.load(f'{outdir}/orbits_{s}.npy')) * e for s in set(prof)}
    paf = {s: paf_classes(use[s][rng.choice(len(use[s]), min(P, len(use[s])), replace=False)], R, Nsp) for s in set(prof)}
    big = {s: paf_classes(use[s][rng.choice(len(use[s]), min(50000, len(use[s])), replace=False)], R, Nsp) for s in set(prof)}
    res = [big[s] - big[s].mean(0) for s in big]; Cv = np.cov(np.concatenate(res).T); mu = sum(big[s].mean(0) for s in prof)
    def bits(idx):
        Cs = 4 * Cv[np.ix_(idx, idx)]; sgn, ld = np.linalg.slogdet(2 * np.pi * Cs); m_ = mu[idx]
        return float(-(len(idx) * np.log(4) - 0.5 * ld - 0.5 * m_ @ np.linalg.solve(Cs, m_)) / np.log(2))
    full = list(range(R.L - 1)); npr = bits(full)
    log(stage='nprime_gauss', rho=rho, profile=prof, n_prime=round(npr, 1), per_class=round(npr / (R.L - 1), 2),
        log2_closed_tuples=round(float(sum(np.log2(sizes[s]) for s in prof)), 1),
        log2_solutions_in_closed_pools=round(float(sum(np.log2(sizes[s]) for s in prof) - npr), 1))
    order = rng.permutation(R.L - 1)
    for k in range(1, kmax + 1):
        idx = order[:k]; W = rng.integers(1, 2 ** 62, size=k).astype(np.uint64)
        h = [(paf[s][:, idx].astype(np.uint64) * W).sum(1) for s in prof]
        AB = (h[0][:, None] + h[1][None, :]).ravel(); CD = np.sort((np.uint64(0) - (h[2][:, None] + h[3][None, :])).ravel())
        cnt = int((np.searchsorted(CD, AB, 'right') - np.searchsorted(CD, AB, 'left')).sum()); del AB, CD
        tot = float(sum(np.log2(len(paf[s])) for s in prof))
        log(stage='partial_join', k=k, matches=cnt, bits_measured=round(tot - float(np.log2(cnt)), 2) if cnt else None,
            bits_gaussian=round(bits(list(idx)), 2))


if __name__ == '__main__':
    main()
