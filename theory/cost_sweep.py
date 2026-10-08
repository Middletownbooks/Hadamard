"""Gaussian costs of all symmetric slices (Goursat subgroups of Aut(G)) with r in [RMIN, RMAX] orbits.
usage: cost_sweep.py "comps" kind RMIN RMAX [TOP]"""
import sys, ast, json, time
import numpy as np
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare
from cost import gauss_cost
comps = ast.literal_eval(sys.argv[1]); kind = sys.argv[2]; RMIN = int(sys.argv[3]); RMAX = int(sys.argv[4])
TOP = int(sys.argv[5]) if len(sys.argv) > 5 else 15
G = AbGroup(comps); AG = AutGroup(G)
subs = goursat_subgroups(AG)
seen = set(); res = []
for gens in subs:
    auts = [AG.to_aut(x) for x in gens]
    orb, lag, frq = slice_structure(G, auts)
    key = tuple(orb)
    if key in seen: continue
    seen.add(key)
    r = int(orb.max()) + 1
    if not (RMIN <= r <= RMAX): continue
    P = prepare(G, auts)
    c = gauss_cost(P, kind, ns=2000)
    if not c.get('feasible'): continue
    S = len(AG.closure(gens))
    res.append(dict(G=str(G), kind=kind, S=S, r=r, L=len(P['lrep']), **c))
res = [x for x in res if x['null_mismatch'] < 1e-6]
res.sort(key=lambda x: x['n'])
print(f'# G={G} kind={kind} slices costed={len(res)}')
for x in [y for y in res if y['log2E'] > 0][:TOP]:
    print(json.dumps(x))
print('# lowest-n with log2E<=0:')
for x in [y for y in res if y['log2E'] <= 0][:5]:
    print(json.dumps(x))
