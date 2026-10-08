"""Positive control for compress_screen: every slice with a solution (found by slicesolve) must pass every screen."""
import itertools, sys
import numpy as np
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare, solve
from compress_screen import screen
bad = 0; checked = 0
for comps, kind in [([('E',2,2),('C',7)], 'b1'), ([('E',2,3),('C',3)], 'b1'), ([('E',2,2),('C',5)], 'b1'), ([('C',3),('C',5)], 'lp'),
                    ([('C',5),('C',7)], 'lp'), ([('C',3),('C',7)], 'b0'), ([('E',2,2),('C',3),('C',5)], 'b1')]:
    G = AbGroup(comps); AG = AutGroup(G); seen = set()
    for gens in goursat_subgroups(AG):
        auts = [AG.to_aut(x) for x in gens]
        orb, _, _ = slice_structure(G, auts)
        if tuple(orb) in seen: continue
        seen.add(tuple(orb))
        r = int(orb.max()) + 1
        if r > 22: continue
        sols, info = solve(prepare(G, auts), kind, max_solutions=1, verbose=False)
        if not sols: continue
        for Hc in [[i] for i in range(len(comps))]:
            ok, inf = screen(G, auts, kind, Hc, verbose=False)
            checked += 1
            if ok is False:
                bad += 1; print('FALSE NEGATIVE', G, kind, r, Hc)
    print(G, kind, 'done', flush=True)
print(f'screens checked on solvable slices: {checked}, false negatives: {bad}')
