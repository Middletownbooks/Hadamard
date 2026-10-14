r"""Exact stage 1 of one b >= 3 configuration, optionally restricted to a template's border class [M].
usage: python3 theory/ab_s1_one.py '<cfg json>' [--tpl 1388|1436] [--maxsol N] [--out FILE]
(same code path as theory/ab_cbatch.py, for a single configuration, without the time limit)"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_s1sd as SD, ab_borderfit as BF
from ab_s1_bordered import allowed

a = sys.argv[1:]
cfg = json.loads(a[0]); tpl = int(a[a.index('--tpl') + 1]) if '--tpl' in a else None
maxsol = int(a[a.index('--maxsol') + 1]) if '--maxsol' in a else 20000
out = a[a.index('--out') + 1] if '--out' in a else None
c2 = SD.cfg_from_json(json.dumps(cfg)); sl = W.build(c2)
t0 = time.time()
Z = allowed(c2, tpl) if tpl else None
P = SD.S1SD(c2, zallow=Z, sl=sl)
rq, A = SD.row_quads_of(c2, sl=sl)
rec = dict(cfg=cfg, tpl=tpl, orbits=int(sl.r), L1=len(sl.breps), A=A.get('A'), row_quads=len(rq) if rq is not None else None)
print(json.dumps(rec), flush=True)
if A.get('A'):
    res, sols = P.solve(rq, max_sol=maxsol, log=lambda s: print(s, flush=True))
    nfit = 0; ex = []
    for s_ in sols:
        Yb = np.array(s_['y']); ok = True
        if tpl:
            cq = len(P.zallow[0]); z = np.zeros((4, cq), np.int64)
            for k in range(4):
                for qq in range(cq): z[k, qq] = Yb[k][P.qfull == qq].sum()
            ok = bool(BF.fits(z, tpl))
        if ok:
            nfit += 1
            if len(ex) < (10 ** 7 if "--all" in a else 50): ex.append(Yb.tolist())
    rec.update(result=res, bordered=nfit, secs=round(time.time() - t0, 1))
    print(json.dumps({k: v for k, v in rec.items()}), flush=True)
    if out: json.dump(dict(rec, examples=ex), open(out, 'w'))
