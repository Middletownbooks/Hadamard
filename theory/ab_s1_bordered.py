r"""Exact stage 1 restricted to Q-compressions that a published frame borders [M].

For a c = 6 (or 4, 8) GP slice of an open order: the allowed per-block compressions are the template's compressions
moved along the corner-kernel direction w (mean shifts s with s = h' - h mod 2) and translated by Q; the exact
stage-1 solver (theory/ab_s1sd.py, zallow) enumerates every stage-1 solution whose blocks' compressions lie in that
set; each solution is then tested exactly against the frame (theory/ab_borderfit.py: +-1 forced corner for some block
order and signs).
usage: python3 theory/ab_s1_bordered.py CFGJSON TEMPLATE_N [--smax 16] [--max-sol N] [--out FILE] [--no-translate]"""
import sys, os, json, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD, ab_borderfit as BF, ab_border_transfer2 as BT
from ab_lib import AG

def allowed(cfg, N, smax=16, translate=True):
    T = BT.load(N); W = BT.kernel_dir(T)
    c = AG(T['Qmods']).n; Q = AG(T['Qmods'])
    hT = {1388: 57, 1916: 119, 1436: 44}[N]
    par = (cfg['h'] - hT) % 2
    out = set()
    for k in range(4):
        for s in range(-smax, smax + 1):
            if (s - par) % 2: continue
            z = T['comp'][k] + s * W[k]
            vs = [z]
            if translate:
                vs = [z[Q.idx(Q.co + Q.co[t])] for t in range(c)]
            for zz in vs: out.add(tuple(int(x) for x in zz))
    return [list(z) for z in out]

if __name__ == "__main__":
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(a[0]); N = int(a[1])
    smax = int(a[a.index('--smax') + 1]) if '--smax' in a else 16
    max_sol = int(a[a.index('--max-sol') + 1]) if '--max-sol' in a else 10 ** 6
    Z = allowed(cfg, N, smax, '--no-translate' not in a)
    P = SD.S1SD(cfg, zallow=Z)
    rq, Ainfo = SD.row_quads_of(cfg)
    t0 = time.time()
    res, sols = P.solve(rq, max_sol=max_sol, log=lambda s: print(s, flush=True), progress=500)
    fit = []
    for s_ in sols:
        Yb = np.array(s_['y']); c = P.zallow.shape[1]
        z = np.zeros((4, c), np.int64)
        for k in range(4):
            for q in range(c): z[k, q] = Yb[k][P.qfull == q].sum()
        f = BF.fits(z, N)
        s_['border_fit'] = [list(f[0][0]), list(f[0][1])] if f else None
        s_['zcomp'] = z.tolist()
        if f: fit.append(s_)
    out = dict(cfg=dict(v=cfg['v'], c=cfg['c'], q=cfg['q'], s=cfg['s'], U=[list(u) for u in cfg['U']]), template=N,
               allowed_vectors=len(Z), **{k: v for k, v in res.items() if k != 'stage'}, bordered_solutions=len(fit),
               secs_total=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    if '--out' in a:
        json.dump(dict(result=out, solutions=fit[:5000]), open(a[a.index('--out') + 1], 'w'))
