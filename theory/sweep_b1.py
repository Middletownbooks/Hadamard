"""b=1 (bordered Goethals-Seidel) slices over abelian groups G of order m = v-1, all symmetry subgroups
S <= Aut(G) (Goursat; GL-part up to conjugacy) with at most RMAX orbits; exhaustive per slice.
usage: sweep_b1.py "comps" RMAX [kind]   e.g.  sweep_b1.py "[('E',2,3),('C',67)]" 26"""
import sys, ast, json, time
import numpy as np
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare, solve
from certify import build, verify

comps = ast.literal_eval(sys.argv[1]); RMAX = int(sys.argv[2]); kind = sys.argv[3] if len(sys.argv) > 3 else 'b1'
LMIN = int(sys.argv[4]) if len(sys.argv) > 4 else 1
G = AbGroup(comps); AG = AutGroup(G)
t0 = time.time()
subs = goursat_subgroups(AG)
print(f'# G={G} n={G.n} |Aut|={AG.order} subgroups(GL up to conj)={len(subs)} ({time.time()-t0:.1f}s)', flush=True)
seen = set(); rows = []
for gens in subs:
    auts = [AG.to_aut(x) for x in gens]
    orb, lag, frq = slice_structure(G, auts)
    key = tuple(orb)
    if key in seen: continue
    seen.add(key)
    r = int(orb.max()) + 1
    L = len(set(lag.tolist())) - 1
    if r > RMAX or L < LMIN: continue
    S = AG.closure(gens)
    rows.append((r, L, len(S), gens))
rows.sort(key=lambda t: (t[0], t[1]))
print(f'# distinct slices with r<={RMAX}: {len(rows)}', flush=True)
tot = 0
for r, L, order, gens in rows:
    auts = [AG.to_aut(x) for x in gens]
    P = prepare(G, auts)
    try:
        sols, info = solve(P, kind, max_solutions=2, verbose=False)
    except RuntimeError as e:
        print(json.dumps(dict(G=str(G), kind=kind, S=order, r=r, L=L, error=str(e))), flush=True); continue
    cert = None
    if sols:
        H = build(G, kind, sols[0]); cert = bool(H is not None and verify(H))
        np.save(f'sol_{kind}_{G}_{order}_{r}.npy', np.array(sols[0]))
    tot += len(sols)
    print(json.dumps(dict(G=str(G), kind=kind, S=order, r=r, L=L, surv=info['survivors'], nsol=len(sols), certified=cert, secs=info['secs'], undecided=info.get('undecided'))), flush=True)
print(f'# done: total solutions found {tot}  ({time.time()-t0:.1f}s)')
