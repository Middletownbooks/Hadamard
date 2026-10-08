"""Positive control for the 573-type configuration: b=1 over Z2^2 x Z3 x Z5 (m = 60, v = 61, H(244)).
Find a solvable symmetric slice with nontrivial S, certify the full matrix, and run the quotient screens on it."""
import numpy as np
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare, solve
from certify import build, verify
from compress_screen import screen
G = AbGroup([('E', 2, 2), ('C', 3), ('C', 5)]); AG = AutGroup(G)
seen = set(); found = 0
for gens in sorted(goursat_subgroups(AG), key=lambda g: -len(AG.closure(g))):
    auts = [AG.to_aut(x) for x in gens]
    orb, _, _ = slice_structure(G, auts)
    if tuple(orb) in seen: continue
    seen.add(tuple(orb)); r = int(orb.max()) + 1
    S = len(AG.closure(gens))
    if S < 2 or r > 24: continue
    sols, info = solve(prepare(G, auts), 'b1', max_solutions=1, verbose=False)
    if sols:
        H = build(G, 'b1', sols[0]); ok = verify(H)
        scr = {str(Hc): screen(G, auts, 'b1', Hc, verbose=False)[0] for Hc in ([0], [1], [2], [0, 1], [0, 2], [1, 2])}
        print(f'|S|={S} r={r}: solution found, H({H.shape[0]}) certified={ok}; screens on this slice: {scr}', flush=True)
        found += 1
        if found >= 3: break
print('solvable nontrivial slices shown:', found)
