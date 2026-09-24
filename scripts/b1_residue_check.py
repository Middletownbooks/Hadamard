"""Congruence obstruction test: for a b=1 multiplier configuration, sample random U-invariant 4-tuples with row sums
(2,0,0,0) and record sum_k PAF_k(c) mod M per lag class; if some class never hits (-4 mod M) the config is infeasible."""
import sys, numpy as np, json
sys.path.insert(0,'scripts')
from b1_multiplier_enum import subgroups, orbits, lagclasses
rng=np.random.default_rng(1)
def check(v,d,r,M=(8,16,32),ns=3000):
    m=v-1
    for U in subgroups(m):
        lab,orbs=orbits(U,m)
        if len(U)!=d or len(orbs)!=r: continue
        sz=np.array([len(o) for o in orbs]); reps,_=lagclasses(U,m)
        def samp(a,n):
            out=[]
            while len(out)<n:
                Y=rng.choice([-1,1],size=(20000,r)); out.extend(Y[(Y@sz)==a][:n-len(out)])
            X=np.array(out)[:,lab]; F=np.fft.fft(X.astype(float),axis=1)
            return np.rint(np.fft.ifft(np.abs(F)**2,axis=1).real).astype(np.int64)[:,reps]
        A=samp(2,ns); B=samp(0,ns); C=samp(0,ns); D=samp(0,ns)
        S=A+B+C+D
        res={}
        for Mm in M:
            hit=np.mean(np.all((S-(-4))%Mm==0,axis=1)); per=[bool(((S[:,c]+4)%Mm==0).any()) for c in range(S.shape[1])]
            res[Mm]=dict(all_classes_hit_frac=round(float(hit),5),classes_ever_hit=f"{sum(per)}/{len(per)}",
                         residues_class0=sorted(set(((S[:,0])%Mm).tolist()))[:8])
        return dict(v=v,d=d,r=r,**{str(k):x for k,x in res.items()})
for args in [(179,11,18),(179,22,10),(419,15,42),(443,12,50),(537,22,30),(537,22,35),(537,22,32)]:
    print(json.dumps(check(*args)),flush=True)
