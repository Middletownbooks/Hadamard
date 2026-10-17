r"""Which compression multisets of a c = 6 (or c = 8) configuration carry a known border [M].

For every class multiset with a non-zero lag-0 count (theory/s1z_lib.py), every ordered member reachable by per-block
translations (by Q-projections of U-fixed elements of G') and reversals is tested, with all block orders and signs:
  * the published template frame (c = 6: 347's H(1388) template; c = 8: 359's H(1436)), by the forced corner
    C = -K Cb(z)^T Kc^T / (4c) in {+-1} (theory/ab_borderfit.fits);
  * monomial compressions against the Singer templates (theory/ab_mono_border.border_for).
Bordering depends on the actual z (not only on its class), and all members of a multiset have the same stage-1
count, so a multiset is bordered if one member is.
usage: python3 theory/s1z_border.py CFG.json OUT.jsonl"""
import sys, os, json, math, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import ZOrb
import ab_borderfit as BF
import ab_mono_border as MB

TPL = {6: 1388, 8: 1436}


def members(Z, cls, m):
    """ordered compression quadruples (one per translation/reversal choice per block) of multiset m"""
    per = []
    for i in m:
        rep = np.asarray(cls[i]['rep'], np.int64); imgs = set()
        for zz in (rep, rep[Z.qneg]):
            for t in Z.tQ: imgs.add(tuple(int(x) for x in zz[Z.qshift[t]]))
        per.append(sorted(imgs))
    for combo in itertools.product(*per): yield [list(z) for z in combo]


def main():
    a = sys.argv[1:]
    cfg = SD.cfg_from_json(open(a[0]).read() if a[0].endswith('.json') else a[0]); outp = a[1]
    Z = ZOrb(cfg); c = Z.c; Qmods = Z.Q.mods
    C = Z.compressions(); cls, idx = Z.classes(C); qd = Z.class_quads(cls)
    cache = {}
    for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
    qd = [q for q in qd if q['n0'] > 0]
    fo = open(outp, 'w'); t0 = time.time(); nb = 0
    for t, q in enumerate(qd):
        rec = dict(i=t, m=list(q['m']), ordered=q['ordered'], log2N0=round(math.log2(q['n0']), 3), template=None, singer=None)
        for zq in members(Z, cls, q['m']):
            if rec['template'] is None and c in TPL:
                f = BF.fits(np.array(zq), TPL[c])
                if f: rec['template'] = dict(z=zq, perm=list(f[0][0]), signs=list(f[0][1]))
            if rec['singer'] is None and c == 6:
                bd = MB.border_for(np.array(zq), Qmods)
                if bd is not None: rec['singer'] = dict(z=zq, template=bd['template'], a=bd['a'], q=bd['q'])
            if rec['template'] is not None and (rec['singer'] is not None or c != 6): break
        nb += rec['template'] is not None or rec['singer'] is not None
        rec['secs'] = round(time.time() - t0, 1)
        fo.write(json.dumps(rec) + "\n"); fo.flush()
        print(json.dumps({k: (v if k not in ('template', 'singer') else (v is not None)) for k, v in rec.items()}), flush=True)
    print(json.dumps(dict(stage='done', multisets=len(qd), bordered=nb, secs=round(time.time() - t0, 1))), flush=True)


if __name__ == '__main__':
    main()
