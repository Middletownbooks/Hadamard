"""Try to decide the slices that exceeded the join cap, with exact quotient screens.
Slices are re-identified by (|S|, r, L) from the sweep logs."""
import json, itertools
from abgroup import AbGroup, AutGroup, slice_structure
from subgroups import goursat_subgroups
from slicesolve import prepare
from compress_screen import screen
from cost import gauss_cost
targets = {
  "[('C',2),('C',13),('C',23)]": 'logs/b0b1_cyclic_r30.log', "[('C',2),('C',17),('C',19)]": 'logs/b0b1_cyclic_r30.log',
  "[('C',2),('E',3,2),('C',29)]": 'logs/b1_nc_523.log', "[('E',2,2),('C',11),('C',13)]": 'logs/b1_nc_573.log'}
for comps_s, logf in targets.items():
    comps = eval(comps_s); G = AbGroup(comps); AG = AutGroup(G)
    und = [json.loads(l) for l in open(logf) if l.startswith('{') and '"undecided": "' in l and json.loads(l)['G'] == str(G)]
    keys = {(d['S'], d['r'], d['L']) for d in und}
    seen = set()
    for gens in goursat_subgroups(AG):
        auts = [AG.to_aut(x) for x in gens]
        orb, lag, _ = slice_structure(G, auts)
        if tuple(orb) in seen: continue
        seen.add(tuple(orb))
        r = int(orb.max()) + 1; L = len(set(lag.tolist())) - 1; S = len(AG.closure(gens))
        if (S, r, L) not in keys: continue
        res = {}
        for k in range(1, len(comps)):
            for Hc in itertools.combinations(range(len(comps)), k):
                ok, info = screen(G, auts, 'b1', list(Hc), verbose=False)
                res[str(list(Hc))] = ok if ok is not None else 'cap'
                if ok is False: break
            if any(v is False for v in res.values()): break
        c = gauss_cost(prepare(G, auts), 'b1', ns=1500)
        verdict = 'DECIDED: no solutions' if any(v is False for v in res.values()) else 'still undecided'
        print(json.dumps(dict(G=str(G), S=S, r=r, L=L, verdict=verdict, screens=res, gauss_log2E=c.get('log2E'), gauss_n=c.get('n'))), flush=True)
