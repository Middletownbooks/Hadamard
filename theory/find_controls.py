"""Search small groups of the 573 shape (Z2^2 x Zp x Zq) and other non-cyclic shapes for b=1 slices with
nontrivial symmetry that HAVE solutions (positive controls for the non-cyclic b=1 pipeline)."""
import sys
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare, solve
from certify import build, verify
from compress_screen import screen
RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 26
for comps in ([('E',2,2),('C',3),('C',5)], [('E',2,2),('C',3),('C',7)], [('E',2,2),('C',5),('C',7)], [('E',2,2),('C',3),('C',11)],
              [('E',2,2),('C',7)], [('E',2,2),('C',11)], [('E',2,2),('C',13)], [('E',2,3),('C',5)], [('E',2,3),('C',7)], [('C',2),('E',3,2),('C',5)]):
    G = AbGroup(comps); AG = AutGroup(G); seen = set(); hits = []
    for gens in goursat_subgroups(AG):
        auts = [AG.to_aut(x) for x in gens]
        orb, _, _ = slice_structure(G, auts)
        if tuple(orb) in seen: continue
        seen.add(tuple(orb)); r = int(orb.max()) + 1
        S = len(AG.closure(gens))
        if S < 2 or r > RMAX: continue
        sols, info = solve(prepare(G, auts), 'b1', max_solutions=1, verbose=False)
        if sols:
            H = build(G, 'b1', sols[0]); ok = verify(H)
            hits.append((S, r, H.shape[0], ok))
    print(G, 'v =', G.n + 1, ' solvable slices with |S|>=2, r<=%d:' % RMAX, sorted(hits, key=lambda t: -t[0])[:6], flush=True)
