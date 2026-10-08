"""Apply exact compression screens to the Gaussian-feasible b1 slices of a group (top-K by n).
usage: screen_top.py "comps" kind RMIN RMAX K"""
import sys, ast, json
import numpy as np
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare
from cost import gauss_cost
from compress_screen import screen
comps = ast.literal_eval(sys.argv[1]); kind = sys.argv[2]; RMIN, RMAX, K = map(int, sys.argv[3:6])
G = AbGroup(comps); AG = AutGroup(G)
seen = set(); cand = []
for gens in goursat_subgroups(AG):
    auts = [AG.to_aut(x) for x in gens]
    orb, lag, frq = slice_structure(G, auts)
    key = tuple(orb)
    if key in seen: continue
    seen.add(key)
    r = int(orb.max()) + 1
    if not (RMIN <= r <= RMAX): continue
    c = gauss_cost(prepare(G, auts), kind, ns=1500)
    if c.get('feasible') and c['null_mismatch'] < 1e-6 and c['log2E'] > 0:
        cand.append((c['n'], r, len(AG.closure(gens)), gens, c))
cand.sort(key=lambda t: t[0])
ncomp = len(comps)
for n, r, Sord, gens, c in cand[:K]:
    auts = [AG.to_aut(x) for x in gens]
    res = {}
    for Hc in [[i] for i in range(ncomp)] + [list(x) for x in __import__('itertools').combinations(range(ncomp), 2)]:
        if len(Hc) == ncomp: continue
        ok, info = screen(G, auts, kind, Hc, verbose=False)
        res[str(Hc)] = ok if ok is not None else info.get('reason')
        if ok is False: break
    verdict = 'KILLED' if any(v is False for v in res.values()) else 'survives'
    print(json.dumps(dict(G=str(G), S=Sord, r=r, n=n, log2E=c['log2E'], verdict=verdict, screens=res, gens=[list(g) for g in gens])), flush=True)
