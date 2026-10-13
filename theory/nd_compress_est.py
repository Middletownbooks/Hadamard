r"""Compressed stage-1 Gaussian for ANY cyclic slice (m, mu, K): exact count of quadruples with the lag-0 energy and the
rows (DP over the M-orbits of Z_F), times the lattice Gaussian of the other lag classes [E, counts exact].
Same estimator as theory/nd_fq_split.py (s1_estimate), structure from theory/nd_s1_sd.py (structure).
usage: python3 theory/nd_compress_est.py m mu K [K2 ...] [--N 20000]"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nd_s1_sd import structure, lag_reps
from nd_fq_split import s1_estimate

def compressed_cyclic(m, mu, K):
    F, Mg, Fo, reps, letters = structure(m, mu, K)
    MF = sorted({a % F for a in Mg})
    fsz = list(np.bincount(Fo).astype(int))
    lr = lag_reps(F, MF)
    T = [4 * m - 4 * (K - 1)] + [-4 * K] * (len(lr) - 1)
    return dict(F=F, m=m, q=K, gens=[[mu, 0]], Fo=Fo.tolist(), freps=reps, fsz=fsz, letters=letters, lreps=lr, T=T)

if __name__ == "__main__":
    a = sys.argv[1:]; N = 20000
    if "--N" in a: i = a.index("--N"); N = int(a[i + 1]); a = a[:i] + a[i + 2:]
    m, mu = int(a[0]), int(a[1])
    for K in map(int, a[2:]):
        C = compressed_cyclic(m, mu, K)
        r = s1_estimate(C, N)
        print(json.dumps(dict(m=m, v=m + 1, mu=mu, K=K, F=C['F'], orbits=len(C['freps']),
                              letter_sets=sorted({tuple(l) for l in C['letters']}, key=len)[:2], lag_classes=len(C['lreps']) - 1,
                              **{k: r[k] for k in r if k != 'lag_bits'}, lag_bits=r.get('lag_bits', {}).get('bits'))), flush=True)

def quad_sampler(C, seed=0):
    """returns a function drawing quadruples y (4 x F, M-invariant, letters per orbit) uniformly among those with the
    slice's exact lag-0 energy and rows (2, 0, 0, 0) (energy split drawn with its exact weight)"""
    import math
    from collections import defaultdict
    rng = np.random.default_rng(seed); F = C['F']; T0 = int(C['T'][0])
    L = C['letters']; w = C['fsz']; allL = sorted({x for l in L for x in l}); g = 0
    for y in allL: g = math.gcd(g, y * y - 1)
    unit = lambda y: (y * y - 1) // g
    E = (T0 - 4 * F) // g; nr = len(w)
    tabs = [None] * (nr + 1); tabs[nr] = {(0, 0): 1}
    for i in range(nr - 1, -1, -1):
        nd = defaultdict(int)
        for (e, r), c in tabs[i + 1].items():
            for y in L[i]:
                e2 = e + w[i] * unit(y)
                if e2 <= E: nd[(e2, r + w[i] * y)] += c
        tabs[i] = dict(nd)
    c0 = defaultdict(int); c2 = defaultdict(int)
    for (e, r), c in tabs[0].items():
        if r == 0: c0[e] += c
        if r == 2: c2[e] += c
    parts = []
    for e0 in range(E + 1):
        for e1 in range(E + 1 - e0):
            for e2 in range(E + 1 - e0 - e1):
                e3 = E - e0 - e1 - e2
                ww = c2.get(e0, 0) * c0.get(e1, 0) * c0.get(e2, 0) * c0.get(e3, 0)
                if ww: parts.append((e0, e1, e2, e3, ww))
    W = np.array([float(p_[4]) for p_ in parts]); W /= W.sum(); Fo = np.array(C['Fo'])
    def block(e, r):
        vals = []
        for i in range(nr):
            opts, wts = [], []
            for y in L[i]:
                c = tabs[i + 1].get((e - w[i] * unit(y), r - w[i] * y), 0)
                if c: opts.append(y); wts.append(float(c))
            wts = np.array(wts); y = opts[rng.choice(len(opts), p=wts / wts.sum())]
            vals.append(y); e -= w[i] * unit(y); r -= w[i] * y
        return np.array(vals)[Fo]
    def draw():
        e0, e1, e2, e3, _ = parts[rng.choice(len(parts), p=W)]
        return np.array([block(e0, 2), block(e1, 0), block(e2, 0), block(e3, 0)])
    return draw
