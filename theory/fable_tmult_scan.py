#!/usr/bin/env python3
"""CP-SAT scan of typed formal multiplications by Z_t (types of real dimension >= 5 by default).
Usage: fable_tmult_scan.py t time_limit_s mindim outfile   (appends one line per type; skips types already logged)."""
import sys, os, json
sys.path.insert(0, 'theory')
import fable_tmult as F, fable_tmult_types as TT, fable_tmult_sat as S
t = int(sys.argv[1]); tl = float(sys.argv[2]); mindim = int(sys.argv[3]); out = sys.argv[4]
done = set()
if os.path.exists(out):
    for l in open(out):
        try: done.add(tuple(json.loads(l)['type']))
        except Exception: pass
tys = [ty for ty in TT.types() if sum(2 if e == 'f' else 1 if e == 's' else 0 for e in ty) >= mindim]
# canonical representatives only: skip types equal up to permuting blocks
seen = set(); reps = []
def canon(ty):
    # multiset signature: free blocks with their multiplicities/reversal pattern, symmetric with multiplicities
    groups = {}
    for k, e in enumerate(ty):
        if e in ('f', 's'): groups[k] = [e, 1, 0]
    for k, e in enumerate(ty):
        if e not in ('f', 's'):
            g = groups[int(e[1])]; g[1] += 1; g[2] += (e[0] == '*')
    out = []
    for g in groups.values():
        kind, mult, star = g
        if kind == 'f': star = min(star, mult - star)
        out.append((kind, mult, star))
    return tuple(sorted(out))
for ty in tys:
    c = canon(ty)
    if c in seen: continue
    seen.add(c); reps.append(ty)
print(len(reps), 'type classes', flush=True)
for ty in sorted(reps, key=lambda ty: sum(2 if e == 'f' else 1 if e == 's' else 0 for e in ty)):
    if tuple(ty) in done: continue
    F.symbols_for_type = lambda _x, ty=ty: TT.make(ty)
    S.symbols_for_type = F.symbols_for_type
    st = S.solve(t, 'x', tl, 4)
    with open(out, 'a') as f:
        f.write(json.dumps({'t': t, 'type': ty, 'status': st, 'limit': tl}) + '\n')
