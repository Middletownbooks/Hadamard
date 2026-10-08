"""Exact solution count for small b>=3 configurations (validation of bge3_cost / bge3_cost2).
Enumerates all 2^r S-invariant blocks, PAF at lag-class reps, then counts ordered 4-tuples with
sum_k PAF_k = target (and sum a_k^2 = 4v) by pair-sum hashing.
usage: bge3_exact_small.py v c Qmods d r L [group_index]"""
import sys, os, json, itertools, numpy as np
from collections import Counter, defaultdict
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; d0=int(sys.argv[4]); r0=int(sys.argv[5]); L0=int(sys.argv[6])
gsel=int(sys.argv[7]) if len(sys.argv)>7 else 0
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(r0),str(r0)]
exec(src[:cut])
cands=[(key,val) for key,val in seen.items() if val[0]==d0 and val[1]==r0 and len(lagclasses(val[3]))==L0]
key,(d,r,lab,grp,abel)=cands[gsel]
reps=lagclasses(grp); L=len(reps)
tgt=np.array([-4*b if int(x) in Hset else 4 for x in reps],np.int64)
assert r<=22
idx_=np.arange(2**r,dtype=np.int64)
Y=(1-2*((idx_[:,None]>>np.arange(r))&1)).astype(np.int8)
P=np.empty((len(Y),L),np.int32); a=np.empty(len(Y),np.int64)
for s in range(0,len(Y),1<<14):
    Xs=Y[s:s+(1<<14)][:,lab].reshape((-1,)+tuple(mods)).astype(float)
    F=np.fft.fftn(Xs,axes=tuple(range(1,len(mods)+1)))
    Pf=np.fft.ifftn(np.abs(F)**2,axes=tuple(range(1,len(mods)+1))).real.reshape(len(Xs),-1)
    P[s:s+(1<<14)]=np.rint(Pf[:,reps]).astype(np.int32); a[s:s+(1<<14)]=Xs.reshape(len(Xs),-1).sum(1)
# linear 64-bit hash of (a^2, PAF vector); count ordered 4-tuples with hash(B1)+..+hash(B4) == hash(target)
rng_=np.random.default_rng(12345); wv=rng_.integers(1,2**62,size=L+1,dtype=np.int64)|1
Hb=(a.astype(np.int64)**2)*wv[0]+(P.astype(np.int64)*wv[1:]).sum(1)          # wraps mod 2^64
Ht=(np.array([4*v],np.int64)*wv[:1]).sum()+(tgt.astype(np.int64)*wv[1:]).sum()   # wraps mod 2^64 like Hb
uh,cnt=np.unique(Hb,return_counts=True)
# pair hashes with multiplicity
ps=[]; pw=[]
for i0 in range(0,len(uh),2048):
    S=(uh[i0:i0+2048,None]+uh[None,:]).ravel(); Wt=(cnt[i0:i0+2048,None]*cnt[None,:]).ravel()
    o=np.argsort(S,kind='stable'); S=S[o]; Wt=Wt[o]
    u2,st2=np.unique(S,return_index=True); ws=np.add.reduceat(Wt,st2)
    ps.append(u2); pw.append(ws)
S=np.concatenate(ps); Wt=np.concatenate(pw); o=np.argsort(S,kind='stable'); S=S[o]; Wt=Wt[o]
uS,stS=np.unique(S,return_index=True); wS=np.add.reduceat(Wt,stS)
need=Ht-uS
pos=np.searchsorted(uS,need); pos=np.clip(pos,0,len(uS)-1); hit=uS[pos]==need
tot=int((wS[hit].astype(object)*wS[pos[hit]].astype(object)).sum()) if hit.any() else 0
K=uh
import math
print(json.dumps(dict(v=v,b=b,c=c,Q=Q,h=h,d=d,r=r,L=L,blocks=len(Y),distinct=len(K),ordered_solutions=tot,log2E=(round(math.log2(tot),2) if tot else None))))
