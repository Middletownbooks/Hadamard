# exact (p = 1) counts at v = 215 for compressions with small Gaussian estimates (theory/b599_calib.py)
import sys, json, os
sys.path.insert(0, 'theory')
import numpy as np
from b599_fam import Fam
import b599_flevels as FL, b599_calib as CB
part, nparts = int(sys.argv[1]), int(sys.argv[2])
F = Fam(5, 7)
est = {}
for l in open('logs/b599/flevels_n5_q7.jsonl'):
    d = json.loads(l); est[(d['class'], tuple(d['c']))] = d.get('log2E')
quads = [t for t in FL.class_quads(F) if est.get((t[0], tuple(t[2]))) is not None and -3 <= est[(t[0], tuple(t[2]))] <= 8]
quads.sort(key=lambda t: est[(t[0], tuple(t[2]))])
sel = quads[::max(1, len(quads) // 16)][part::nparts]
fo = open(f'logs/b599/calib/calib_215_low.jsonl', 'a')
for ci, cl, c, z in sel:
    g = FL.budget(F, z, N=6000, seed=1)
    r = CB.calib_z(F, z, 1.0)
    rec = dict(n=5, q=7, v=F.v, **{'class': ci}, c=list(c), gauss_log2E=g.get('log2E'), gauss_bits=g['full']['bits'], H1=g.get('H1'), rank=g['full']['rank'], **r)
    fo.write(json.dumps(rec) + "\n"); fo.flush(); print(json.dumps(rec), flush=True)
