"""Exact pool-level feasibility of b=1 multiplier configurations.
Integer-valued frequency classes (P(w) takes few values on U-invariant blocks) are found by sampling; for those,
the exact achievable P(w) sets over the COMPLETE PSD-filtered pools (row sums 2 and 0) are computed by a
streaming scan, and we test whether P_A + P_B + P_C + P_D = 4m+4 is reachable.
usage: b1_pool_feasible.py v rmax   (all subgroups with r <= rmax)"""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from b1red import RedC
from b1_pool_enum import pool_value_sets
from b1_multiplier_enum import subgroups, orbits, count_rowsum
from math import log2
v=int(sys.argv[1]); rmax=int(sys.argv[2]); m=v-1; bound=4*m+4
cfg=[]
for U in subgroups(m):
    r=len(orbits(U,m)[1])
    if len(U)>1 and r<=rmax: cfg.append((r,sorted(U)))
for r,U in sorted(cfg):
    R=RedC(m,U); t0=time.time()
    if not (count_rowsum(list(R.size),0) and count_rowsum(list(R.size),2)):
        print(json.dumps(dict(v=v,d=len(U),r=r,L=R.L,feasible=False,reason='row sums')),flush=True); continue
    Y=np.random.default_rng(0).choice([-1,1],size=(20000,R.r)); P=np.abs(Y@R.S[1:].T)**2
    fsel=[f for f in range(P.shape[1]) if np.all(np.abs(P[:,f]-np.rint(P[:,f]))<1e-6)]
    res={}; n={}
    for a in (2,0): n[a],res[a]=pool_value_sets(R,a,bound,fsel)
    bad=[]
    for j,f in enumerate(fsel):
        V0=res[0][j]; S2=np.unique((V0[:,None]+V0[None,:]).ravel()); S3=np.unique((S2[:,None]+V0[None,:]).ravel())
        if not np.isin(bound-res[2][j],S3).any(): bad.append(int(R.freqs[f+1]))
    out=dict(v=v,d=len(U),r=r,L=R.L,U=U[:6],pool_log2={a:round(log2(n[a]),2) if n[a] else None for a in n},
             int_freqs=[int(R.freqs[f+1]) for f in fsel],feasible=(len(bad)==0 and n[0]>0 and n[2]>0),bad_freqs=bad,secs=round(time.time()-t0,1))
    print(json.dumps(out),flush=True)
