"""Positive control for the 573-shape non-cyclic b=1 pipeline: G = Z2^2 x Z3 x Z7 (v = 85), |S| = 6.
Finds the solvable slices, certifies H(340), saves the core, and runs every quotient screen on the slice."""
import json
import numpy as np
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare, solve
from certify import build, verify
from compress_screen import screen
G = AbGroup([('E', 2, 2), ('C', 3), ('C', 7)]); AG = AutGroup(G); seen = set(); out = []
for gens in goursat_subgroups(AG):
    auts = [AG.to_aut(x) for x in gens]
    orb, _, _ = slice_structure(G, auts)
    if tuple(orb) in seen: continue
    seen.add(tuple(orb)); r = int(orb.max()) + 1; S = len(AG.closure(gens))
    if S != 6 or r > 21: continue
    sols, info = solve(prepare(G, auts), 'b1', max_solutions=1, verbose=False)
    if not sols: continue
    H = build(G, 'b1', sols[0]); ok = verify(H)
    scr = {str(Hc): screen(G, auts, 'b1', Hc, verbose=False)[0] for Hc in ([0], [1], [2], [0, 1], [0, 2], [1, 2])}
    rec = dict(G=str(G), S=S, r=r, gens=[list(g) for g in gens], certified=ok, order=H.shape[0], screens=scr,
               core=[''.join('+' if x > 0 else '-' for x in s) for s in sols[0]])
    out.append(rec); print(json.dumps({k: rec[k] for k in ('G', 'S', 'r', 'certified', 'order', 'screens')}), flush=True)
    np.save(f'logs/controls/H340_b1_Z2sqZ3Z7_S{S}_r{r}.npy', (H < 0).astype(np.uint8))
json.dump(out, open('logs/controls/control_340.json', 'w'), indent=1)
