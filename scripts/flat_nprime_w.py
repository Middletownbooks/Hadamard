"""Flat-subspace constraint bits n' for a b = 0 multiplier slice from orbit representatives (any number of orbits).

W-word version of step 3 of scripts/flat_probe.py: samples uniformly from the rotation closures of the orbit
representatives in DIR... (scripts/flat_gen_w.py or scripts/flat_probe.py output; 1-D masks are read as W = 1),
fits the class-PAF Gaussian (pooled residual covariance, per-row-sum means) and prints the Gaussian n' over the
L - 1 independent classes, then partial joins on k = 1..kmax random classes (pools of P per slot) against the
Gaussian bits of the same classes.  Also prints the closure sizes and log2 solution tuples in the closed pools.
usage: python3 scripts/flat_nprime_w.py V D a,b,c,d DIR1 [DIR2 ...] [--P 4096] [--kmax 10] [--big 50000]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import flatw as FW
KJ = FW.KJ


def load_reps(dirs, s):
    parts = []
    for dd in dirs:
        f = os.path.join(dd, f'orbits_{s}.npy')
        if os.path.exists(f):
            x = np.load(f); parts.append(x if x.ndim == 2 else x.reshape(-1, 1))
    return FW.unique_w(np.concatenate(parts)) if parts else None


def closure_sample_w(reps, e, n, rng):
    i = rng.integers(0, len(reps), n); s = rng.integers(0, e, n); out = np.empty((n, reps.shape[1]), np.uint64)
    for s_ in range(e):
        m = s == s_
        if m.any(): out[m] = FW.rotate_w(np.ascontiguousarray(reps[i[m]]), s_, e)
    return FW.unique_w(out)


def main():
    a = sys.argv[1:]
    opt = lambda nm, df: type(df)(a[a.index(nm) + 1]) if nm in a else df
    P = opt('--P', 4096); kmax = opt('--kmax', 10); NBIG = opt('--big', 50000)
    pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or not a[i - 1].startswith('--'))]
    v, d = int(pos[0]), int(pos[1]); prof = [int(x) for x in pos[2].split(',')]; dirs = pos[3:]
    R = KJ.Red(v, d); e = (v - 1) // d; Nsp = KJ.sparse_N(R); t0 = time.time()
    log = lambda **k: print(json.dumps(dict(k, t=round(time.time() - t0, 1))), flush=True)
    rng = np.random.default_rng(0)
    reps = {s: load_reps(dirs, s) for s in set(prof)}
    sizes = {s: len(reps[s]) * e for s in set(prof)}
    use = {s: closure_sample_w(reps[s], e, 3 * NBIG, rng) for s in set(prof)}
    paf = lambda x: FW.paf_w(np.ascontiguousarray(x), *Nsp, R.L).astype(np.int64)
    big = {s: paf(use[s][rng.choice(len(use[s]), min(NBIG, len(use[s])), replace=False)]) for s in set(prof)}
    small = {s: paf(use[s][rng.choice(len(use[s]), min(P, len(use[s])), replace=False)]) for s in set(prof)}
    F = {s: (big[s] ** 2).sum(1) for s in big}
    log(stage='pools', v=v, d=d, profile=prof, orbits={str(s): len(reps[s]) for s in reps},
        log2_closure={str(s): round(float(np.log2(sizes[s])), 2) for s in sizes},
        rho_max={str(s): round(float(np.sqrt(F[s].max() / (v * R.L))), 4) for s in F},
        rho_median={str(s): round(float(np.sqrt(np.median(F[s]) / (v * R.L))), 4) for s in F})
    res = [big[s] - big[s].mean(0) for s in big]; Cv = np.cov(np.concatenate(res).T); mu = sum(big[s].mean(0) for s in prof)
    def bits(idx):
        Cs = 4 * Cv[np.ix_(idx, idx)]; sgn, ld = np.linalg.slogdet(2 * np.pi * Cs); m_ = mu[idx]
        return float(-(len(idx) * np.log(4) - 0.5 * ld - 0.5 * m_ @ np.linalg.solve(Cs, m_)) / np.log(2))
    full = list(range(R.L - 1)); npr = bits(full); lt = float(sum(np.log2(sizes[s]) for s in prof))
    log(stage='nprime_gauss', n_prime=round(npr, 2), per_class=round(npr / (R.L - 1), 3),
        log2_closed_tuples=round(lt, 2), log2_solutions_in_closed_pools=round(lt - npr, 2))
    order = rng.permutation(R.L - 1); excess = []
    for k in range(1, kmax + 1):
        idx = order[:k]; Wt = rng.integers(1, 2 ** 62, size=k).astype(np.uint64)
        h = [(small[s][:, idx].astype(np.uint64) * Wt).sum(1) for s in prof]
        AB = (h[0][:, None] + h[1][None, :]).ravel(); CD = np.sort((np.uint64(0) - (h[2][:, None] + h[3][None, :])).ravel())
        cnt = int((np.searchsorted(CD, AB, 'right') - np.searchsorted(CD, AB, 'left')).sum()); del AB, CD
        tot = float(sum(np.log2(len(small[s])) for s in prof)); bg = bits(list(idx))
        bm = tot - float(np.log2(cnt)) if cnt else None
        if bm is not None and cnt >= 30: excess.append((bm - bg) / k)
        log(stage='partial_join', k=k, matches=cnt, bits_measured=round(bm, 2) if bm else None, bits_gaussian=round(bg, 2))
    if excess:
        ex = float(excess[-1])                    # the excess grows with k: take the largest k with >= 30 matches
        log(stage='nprime', n_prime_gauss=round(npr, 2), excess_per_class=round(ex, 3), n_prime=round(npr + ex * (R.L - 1), 2),
            per_class=round(npr / (R.L - 1) + ex, 3))


if __name__ == '__main__':
    main()
