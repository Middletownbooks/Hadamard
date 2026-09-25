"""515 pipeline: stage A (compressed, t=17) -> forced-residue filter (per-class, then pairwise joint mod 16/64) -> sampled stage B.
usage: bge3_515_pipeline.py nA nbuckets_A   (collect nA compressed solutions, report how many pass the residue filter)"""
import sys, os, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageA import solve_A
from bge3_stageB_join import setup, components
from bge3_residue import residue_ok, joint_residue_ok
nA=int(sys.argv[1]); os.environ['A_BUCKETS']=sys.argv[2]
h,p,u=85,5,53
g=setup(h,p,u); comps=components(g['footB'].T); co=g['co']; m=g['m']; NEG=g['NEG']
e=int(np.lcm.reduce([2,3,h])); perm=g['idx'](co*(u%e)); cl=-np.ones(m,int); reps=[]
for x in range(1,m):
    if cl[x]<0:
        S=set(); y=x
        while y not in S: S.add(y); S.add(int(NEG[y])); y=int(perm[y])
        for z in S: cl[z]=len(reps)
        reps.append(x)
reps=np.array(reps)
sols=solve_A(h,p,u,log=lambda s: print(s,flush=True),max_solutions=nA)
good=[]
for i,(s,xc) in enumerate(sols):
    ok,bad=residue_ok(g,comps,reps,xc,seed=i); jinfo=None
    if ok:
        for M in (16,64):
            ok,jinfo=joint_residue_ok(g,comps,reps,xc,h,ns=8000,M=M,seed=i)
            if not ok: jinfo['M']=M; break
    if ok: good.append(xc)
    print(json.dumps(dict(sol=i,shift=[int(x) for x in s],residue_ok=bool(ok),bad_lags=bad[:6] if isinstance(bad,list) else str(bad),joint=jinfo)),flush=True)
print(json.dumps(dict(compressed=len(sols),joint_residue_feasible=len(good))),flush=True)
np.save('data/stageA_515_residue_ok.npy',np.array(good) if good else np.zeros((0,)))
