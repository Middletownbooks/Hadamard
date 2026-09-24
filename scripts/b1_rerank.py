"""Re-screen b=1 multiplier configurations (log2E > 0 in the Gaussian enumeration) with the exact
Z_t-compression obstruction (b1red.compression_feasible) for every divisor 2 <= t <= TMAX of m.
Verdict: INFEASIBLE if any t proves it; 'passes' if all decided t pass; undecided t (state blow-up) listed.
usage: b1_rerank.py logfile [logfile ...]    env TMAX (default 24), MAXST (default 4e6)"""
import sys, os, json, time, numpy as np
sys.path.insert(0,'scripts')
from b1red import RedC, compression_feasible
from b1_multiplier_enum import subgroups, orbits, lagclasses
TMAX=int(os.environ.get('TMAX','24')); MAXST=float(os.environ.get('MAXST','4e6'))
for f in sys.argv[1:]:
    for line in open(f):
        try: rec=json.loads(line)
        except Exception: continue
        v=rec['v']; m=rec['m']; want=[c for c in rec['configs'] if c.get('log2E') is not None and c['log2E']>0]
        if not want: print(json.dumps(dict(v=v,note='no config with log2E>0')),flush=True); continue
        subs=[sorted(U) for U in subgroups(m) if len(U)>1]
        for c in want:
            for U in subs:
                if len(U)!=c['d']: continue
                lab,orbs=orbits(U,m)
                if len(orbs)!=c['r'] or len(lagclasses(U,m)[0])!=c['L']: continue
                R=RedC(m,U); res={}; t0=time.time(); verdict='passes'
                for t in [t for t in range(2,TMAX+1) if m%t==0]:
                    ok,info=compression_feasible(R,t,maxstates=MAXST); res[t]=ok
                    if ok is False: verdict='INFEASIBLE'; break
                print(json.dumps(dict(v=v,d=c['d'],r=c['r'],L=c['L'],log2E=c['log2E'],n=c['n_cond'],U=U[:6],verdict=verdict,
                      tested={str(k):v_ for k,v_ in res.items()},secs=round(time.time()-t0,1))),flush=True)
