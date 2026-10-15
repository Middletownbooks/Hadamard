r"""Calibration of the conditional stage-1 Gaussian (theory/b599_flevels.py) against exact counts [M against E].

For a family member (theory/b599_fam.py) and Q-compression quadruples z: every block vector with compression z_k is
listed exactly (theory/b599_exact.block_list); blocks 0 and 2 are subsampled independently with probability p
(seeded); every 4-tuple of the subsampled lists with sum of class signatures = T1 is found by the exact join
theory/p647_join.join4 (energy split, residue passes, exact re-check; imported unchanged) and re-verified by the full
G' PAF.  hits / p^2 is an unbiased estimate of the exact number of ordered stage-1 solutions with compressions z
(p = 1: the exact count).  The Gaussian prediction is recomputed for the same z with N samples.
usage: python3 theory/b599_calib.py --n 5 --q 7 --top 6 --p 0.0156 [--N 6000] [--out FILE] [--seed 1]"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b599_fam import Fam
import b599_exact as EX, b599_flevels as FL
import p647_join as PJ


def calib_z(F, z, p, seed=1, log=None):
    T0 = int(F.T1SIG[0]); rng = np.random.default_rng(seed)
    Ls = []; mins = []
    for k in range(4):
        V, E = EX.block_list(F, z[k], T0)
        if len(V) == 0: return dict(hits=0, sizes=[0] * 4)
        Ls.append(V); mins.append(int(E.min()))
    full = []
    for k in range(4):
        cap = T0 - (sum(mins) - mins[k]); Ls[k] = Ls[k][F.energy(Ls[k]) <= cap]; full.append(len(Ls[k]))
    for k in (0, 2):
        if p < 1: Ls[k] = Ls[k][rng.random(len(Ls[k])) < p]
    S = [F.sig(x) for x in Ls]
    t0 = time.time()
    sols, st = PJ.join4(S, F.T1SIG.astype(np.int64), log=log)
    ok = 0
    for (i0, i1, i2, i3) in sols:
        Y = np.array([Ls[0][i0], Ls[1][i1], Ls[2][i2], Ls[3][i3]])
        ok += int(F.check_stage1(Y) and (F.qcomp(Y) == z).all())
    hits = len(sols)
    est = hits / (p * p)
    return dict(hits=hits, verified=ok, p=p, count_est=est, log2_count=round(math.log2(est), 2) if hits else None,
                rel_sd=round(1 / math.sqrt(hits), 3) if hits else None, sizes=full,
                work_log2=st.get('work_log2'), join_secs=st.get('secs'), ns_per_pair=st.get('ns_per_pair'),
                secs=round(time.time() - t0, 1))


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    n, q = int(opt('--n', '5')), int(opt('--q', '7'))
    top = int(opt('--top', '6')); p = float(opt('--p', '0.0156')); N = int(opt('--N', '6000'))
    seed = int(opt('--seed', '1')); outp = opt('--out', ''); skip = int(opt('--skip', '0'))
    F = Fam(n, q)
    est = {}
    fn = f'logs/b599/flevels_n{n}_q{q}.jsonl'
    if os.path.exists(fn):
        for l in open(fn):
            d = json.loads(l); est[(d['class'], tuple(d['c']))] = d.get('log2E')
    quads = FL.class_quads(F)
    quads.sort(key=lambda t: -(est.get((t[0], tuple(t[2]))) or -1e9))
    fo = open(outp, 'a') if outp else None
    for ci, cl, c, z in quads[skip:skip + top]:
        g = FL.budget(F, z, N=N, seed=seed)
        r = calib_z(F, z, p, seed=seed)
        rec = dict(n=n, q=q, v=F.v, **{'class': ci}, c=list(c), gauss_log2E=g.get('log2E'), gauss_bits=g['full']['bits'],
                   H1=g.get('H1'), rank=g['full']['rank'], **r)
        s = json.dumps(rec); print(s, flush=True)
        if fo: fo.write(s + "\n"); fo.flush()
