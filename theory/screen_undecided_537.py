"""Exact quotient screens for the one undecided cyclic slice: b1 over Z8 x Z67 with |S| = 22, r = 29."""
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from compress_screen import screen
G = AbGroup([('C', 8), ('C', 67)]); AG = AutGroup(G)
for gens in goursat_subgroups(AG):
    auts = [AG.to_aut(x) for x in gens]
    orb, lag, _ = slice_structure(G, auts)
    r = int(orb.max()) + 1
    if r == 29 and len(AG.closure(gens)) == 22:
        for Hc in ([1], [0]):
            print('slice r=29 |S|=22, H=component', Hc, '->', screen(G, auts, 'b1', Hc, verbose=False), flush=True)
