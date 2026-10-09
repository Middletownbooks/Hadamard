r"""Gaussian (row-sum conditioned) estimates for ORBIT SLICES of the Paley-coset ansatz, any b (0, 1, >= 3).

A slice is a partition of G' into orbits; every block's label (alpha, type, sign) is constant on each orbit.  The
natural slices: 'full' (singletons), 'even' (orbits {g, -g}: X == 0 identically, report section 3), 'mult:u' (orbits of
g -> u g for an automorphism u of a cyclic G'), 'center:t' (orbits {g, t - g}).  Same estimator as
scripts/paley_bge3.gauss_estimate (calibrated there on 359 / 347 / 479 and on exact counts), with labels constant on
orbits and row sums counted with orbit sizes.

usage: python3 theory/paley_slices.py v b p mods... [--Hgen a.b.c,...] [--slice full|even|center:i|mult:u] [--N 40000]"""
import sys, os, json, itertools as it, math
from collections import Counter
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from paley_ansatz import LABELS, LC, corr_table
import paley_ls as LSm

def orbits(P, kind):
    n = P.n
    if kind == 'full': return [[g] for g in range(n)]
    if kind == 'even' or kind.startswith('center:'):
        t = 0 if kind == 'even' else int(kind.split(':')[1])
        seen, out = set(), []
        for g in range(n):
            if g in seen: continue
            o = sorted({g, int(P.sub[t, g])}); seen.update(o); out.append(o)
        return out
    if kind.startswith('mult:'):
        u = int(kind.split(':')[1]); assert len(P.G.mods) == 1
        seen, out = set(), []
        for g in range(n):
            if g in seen: continue
            o, x = [], g
            while x not in o: o.append(x); x = (x * u) % n
            seen.update(o); out.append(sorted(o))
        return out
    raise ValueError(kind)

def reduced_paf_batch(P, lab):
    """lab (N, n) -> (N, n, 3)"""
    A = LC[lab].astype(np.int64)
    R = np.zeros((lab.shape[0], P.n, 3), np.int64)
    for h in range(P.n):
        B = A[:, P.add[:, h]]
        R[:, h, :] = np.einsum('ngi,ngj,ijc->nc', A, B, P.C)
    return R

def estimate(v, b, p, mods, Hgen=None, kind='full', N=40000, seed=1, maxprof=6):
    P = LSm.Prob(v, b, p, mods, Hgen); rng = np.random.default_rng(seed)
    orbs = orbits(P, kind); k = len(orbs); sizes = np.array([len(o) for o in orbs])
    rv = np.array([LABELS[i][0] + (p - 1) * (0 if LABELS[i][1] else LABELS[i][2]) for i in range(8)])
    cnt = Counter({0: 1})
    for s in sizes:
        E = Counter()
        for a, c in cnt.items():
            for L in range(8): E[a + s * rv[L]] += c
        cnt = E
    if b == 1: profs = [(2, 0, 0, 0)] if (2 in cnt and 0 in cnt) else []
    else:
        A_ = sorted(a for a in cnt if a >= 0)
        profs = [q for q in it.combinations_with_replacement(A_, 4) if sum(x * x for x in q) == 4 * v]
    profs = sorted(profs, key=lambda q: -sum(math.log2(cnt[x]) for x in q))[:maxprof]
    T = P.T.reshape(-1)
    def expand(lo):
        lab = np.zeros((lo.shape[0], P.n), np.int64)
        for j, o in enumerate(orbs): lab[:, o] = lo[:, [j]]
        return lab
    stats = {}
    def samp(a):
        if a not in stats:
            res = []
            for _ in range(600):
                lo = rng.integers(0, 8, (40000, k)); lo = lo[(rv[lo] * sizes).sum(1) == a]; res.append(lo)
                if sum(len(x) for x in res) >= N: break
            lo = np.concatenate(res)[:N]
            stats[a] = None if len(lo) < 1000 else reduced_paf_batch(P, expand(lo)).reshape(len(lo), -1)
        return stats[a]
    best, tot = None, 0.0
    for q in profs:
        Rs = [samp(a) for a in q]
        if any(r is None for r in Rs): continue
        mu = sum(r.mean(0) for r in Rs); Cv = sum(np.cov(r.T) for r in Rs); Rc = sum(r[:1500] for r in Rs)
        cols, seen = [], {}
        for j in range(1, len(T)):
            kk = Rc[:, j].tobytes()
            if kk in seen:
                if T[j] != T[seen[kk]]: cols = None; break
                continue
            seen[kk] = j; cols.append(j)
        if cols is None: continue
        w, V = np.linalg.eigh(Cv[np.ix_(cols, cols)]); kw = w > 1e-6 * w.max(); pr = V.T @ (T[cols] - mu[cols])
        if (~kw).any() and np.abs(pr[~kw]).max() > 1e-6 * max(1, np.abs(pr).max()): continue
        d = Rc[:, cols] - Rc[0, cols]; steps = np.array([max(1, int(np.gcd.reduce(np.abs(d[:, j])))) for j in range(len(cols))])
        if not all((T[cols][j] - Rc[0, cols][j]) % steps[j] == 0 for j in range(len(cols))): continue
        lp = (np.log(steps).sum() - 0.5 * kw.sum() * math.log(2 * math.pi) - 0.5 * np.log(w[kw]).sum()
              - 0.5 * (pr[kw] ** 2 / w[kw]).sum()) / math.log(2)
        perms = len(set(it.permutations(q))); signs = 2 ** sum(1 for x in q if x != 0)
        lE = sum(math.log2(cnt[x]) for x in q) + lp + math.log2(perms * signs)
        tot += 2 ** lE
        if best is None or lE > best[1]: best = ([int(t) for t in q], round(float(lE), 1), round(float(-lp), 1), int(kw.sum()))
    rec = dict(v=v, b=b, p=p, mods=list(mods), slice=kind, orbits=k, free_bits=round(4 * k * 3, 1))
    if best is None: rec.update(feasible=False)
    else: rec.update(log2E=round(float(math.log2(tot)), 1), profile=best[0], constraint_bits=best[2], eqs=best[3])
    return rec

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--Hgen": None, "--slice": "full", "--N": "40000", "--seed": "1"}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    v, b, p = pos[:3]; mods = tuple(pos[3:])
    Hgen = None if opt["--Hgen"] is None else [tuple(map(int, e.split('.'))) for e in opt["--Hgen"].split(',')]
    for kind in opt["--slice"].split(','):
        print(json.dumps(estimate(v, b, p, mods, Hgen, kind, N=int(opt["--N"]), seed=int(opt["--seed"]))), flush=True)
