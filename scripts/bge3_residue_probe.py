"""Configuration-level 2-adic screen for small slices: per lag class, PAF residues mod 2^j of all PSD-filtered
S-invariant blocks (by row sum); is the target reachable mod 2^j (j=3,4,5) for each row-sum profile ordering?
(exact for slices small enough to enumerate; used to explain the v=63/d4/r19 calibration zero: unreachable mod 32)
usage: bge3_residue_probe.py v c Qmods d r L"""
"""Why does a slice with predicted solutions have none? Probe: per lag class, the set of PAF residues mod 2^j of
PSD-filtered blocks (by row sum), and whether the target is reachable mod 2^j for every row-sum profile ordering.
Also: partial exact joins on the first k lag classes (count of 4-tuples matching k classes)."""
import sys, os, itertools, numpy as np, json
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; d0=int(sys.argv[4]); r0=int(sys.argv[5]); L0=int(sys.argv[6])
os.environ.setdefault('GMAX','400'); os.environ.setdefault('PAIRMAX','40')
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(r0),str(r0)]
exec(src[:cut])
cands_=[(key,val) for key,val in seen.items() if val[0]==d0 and val[1]==r0 and len(lagclasses(val[3]))==L0]
key,(d,r,lab,grp,abel)=cands_[0]
reps=lagclasses(grp); L=len(reps); sz=np.bincount(lab)
tgt=np.array([-4*b if int(x) in Hset else 4 for x in reps],np.int64)
profs=profiles(sz); A=sorted({a for p in profs for a in p})
nQ=len(Q); chis=np.array(np.unravel_index(np.arange(m),mods)).T
That=np.where((chis[:,-1]==0)&(chis[:,:nQ]!=0).any(1),4*b,4*v).reshape(mods)
Ps=[]; As=[]
for s in range(0,2**r,1<<15):
    idx_=np.arange(s,min(2**r,s+(1<<15)),dtype=np.int64)
    Y=(1-2*((idx_[:,None]>>np.arange(r))&1)).astype(np.int8)
    Xs=Y[:,lab]; a=Xs.sum(1).astype(np.int64); sel=np.isin(a,A)
    Xs=Xs[sel]; a=a[sel]
    if not len(a): continue
    F=np.fft.fftn(Xs.reshape((-1,)+tuple(mods)).astype(np.float64),axes=tuple(range(1,len(mods)+1)))
    Pw=np.abs(F)**2; ok=(Pw<=That[None]+1e-6).reshape(len(Xs),-1).all(1)
    if not ok.any(): continue
    Pf=np.fft.ifftn(Pw[ok],axes=tuple(range(1,len(mods)+1))).real.reshape(int(ok.sum()),-1)
    Ps.append(np.rint(Pf[:,reps]).astype(np.int64)); As.append(a[ok])
P=np.concatenate(Ps); Aa=np.concatenate(As)
print('profiles',profs,' filtered blocks per row sum',{int(x):int((Aa==x).sum()) for x in A})
print('class reps',reps.tolist(),'target',tgt.tolist())
for j in (3,4,5):
    md=2**j; bad=[]
    for p in profs:
        for perm in set(itertools.permutations(p)):
            res=[ {tuple(row) for row in (P[Aa==x]%md)} for x in perm]
            if any(len(rr)==0 for rr in res): continue
            s12={tuple((np.array(u)+np.array(w_))%md) for u in res[0] for w_ in res[1]}
            s34={tuple((np.array(u)+np.array(w_))%md) for u in res[2] for w_ in res[3]}
            need=tuple((-np.array(t))%md for t in s34)
            T_=tuple(tgt%md)
            need_set={tuple((np.array(tgt)-np.array(t))%md) for t in s34}
            reach=len(need_set & s12)>0
            bad.append((perm,reach))
    print('mod',md,': orderings reachable',sum(1 for _,x in bad if x),'of',len(bad))
