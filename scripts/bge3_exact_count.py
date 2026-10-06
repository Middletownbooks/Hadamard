"""Exact solution COUNT per row-sum profile for b>=3 slices of moderate size (r <= ~26), for calibrating bge3_cost.
Blocks with row sum a (all a in the profile set) are enumerated and filtered by the exact per-block PSD bound
|X^(chi)|^2 <= T^(chi) (necessary, so counts stay exact). Ordered 4-tuples (B1..B4) with row sums exactly (a1..a4)
and sum PAF = target are counted by 64-bit linear hashing of pair sums (collisions negligible at these sizes).
Output log2E per profile, directly comparable with bge3_cost.py's per-profile log2E.
usage: bge3_exact_count.py v c Qmods d r L [group_index]"""
import sys, os, json, itertools, numpy as np
from math import log2
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; d0=int(sys.argv[4]); r0=int(sys.argv[5]); L0=int(sys.argv[6])
gsel=int(sys.argv[7]) if len(sys.argv)>7 else 0
os.environ.setdefault('GMAX','400'); os.environ.setdefault('PAIRMAX','40')
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(r0),str(r0)]
exec(src[:cut])
cands=[(key,val) for key,val in seen.items() if val[0]==d0 and val[1]==r0 and len(lagclasses(val[3]))==L0]
key,(d,r,lab,grp,abel)=cands[gsel]
reps=lagclasses(grp); L=len(reps); sz=np.bincount(lab)
tgt=np.array([-4*b if int(x) in Hset else 4 for x in reps],np.int64)
profs=profiles(sz); A=sorted({a for p in profs for a in p})
nQ=len(Q); chis=np.array(np.unravel_index(np.arange(m),mods)).T
That=np.where((chis[:,-1]==0)&(chis[:,:nQ]!=0).any(1),4*b,4*v).reshape(mods)
H=[]; Aa=[]; CH=1<<15
rng_=np.random.default_rng(7); wv=rng_.integers(1,2**62,size=L,dtype=np.int64)|1
for s in range(0,2**r,CH):
    idx_=np.arange(s,min(2**r,s+CH),dtype=np.int64)
    Y=(1-2*((idx_[:,None]>>np.arange(r))&1)).astype(np.int8)
    Xs=Y[:,lab]; a=Xs.sum(1).astype(np.int64); sel=np.isin(a,A)
    if not sel.any(): continue
    Xs=Xs[sel]; a=a[sel]
    F=np.fft.fftn(Xs.reshape((-1,)+tuple(mods)).astype(np.float64),axes=tuple(range(1,len(mods)+1)))
    Pw=np.abs(F)**2; ok=(Pw<=That[None]+1e-6).reshape(len(Xs),-1).all(1)
    if not ok.any(): continue
    Pf=np.fft.ifftn(Pw[ok],axes=tuple(range(1,len(mods)+1))).real.reshape(int(ok.sum()),-1)
    P=np.rint(Pf[:,reps]).astype(np.int64)
    H.append((P*wv).sum(1)); Aa.append(a[ok])
H=np.concatenate(H) if H else np.zeros(0,np.int64); Aa=np.concatenate(Aa) if Aa else np.zeros(0,np.int64)
Ht=(tgt*wv).sum()
def pairhist(x,y):
    S=(x[:,None]+y[None,:]).ravel(); u,cnt=np.unique(S,return_counts=True); return u,cnt
out=[]
for p in profs:
    for perm in sorted(set(itertools.permutations(p))):
        lists=[H[Aa==x] for x in perm]
        if any(len(l)==0 for l in lists): out.append((perm,0)); continue
        u12,c12=pairhist(lists[0],lists[1]); u34,c34=pairhist(lists[2],lists[3])
        need=Ht-u34; pos=np.searchsorted(u12,need); pos=np.clip(pos,0,len(u12)-1); hit=u12[pos]==need
        tot=int((c12[pos[hit]].astype(object)*c34[hit].astype(object)).sum()) if hit.any() else 0
        out.append((perm,tot))
best={}
for perm,tot in out:
    key_=tuple(sorted(perm)); best.setdefault(key_,[]).append(tot)
res=dict(v=v,c=c,Q=Q,h=h,d=d,r=r,L=L,psd_filtered_blocks=int(len(H)),
         per_profile={str(tuple(int(t) for t in k)):dict(ordered_counts=sorted(set(vv)),log2E=(round(log2(max(vv)),2) if max(vv)>0 else None)) for k,vv in best.items()})
print(json.dumps(res,default=int))
