"""Exact solution counts for small-r b=1 multiplier configurations (profile (2,0,0,0)), by
convolution of PAF-class-vector histograms (MITM on distinct vectors)."""
import sys, numpy as np, json
from collections import Counter
sys.path.insert(0,'scripts')
from b1_multiplier_enum import subgroups, orbits, lagclasses
def exact(v,dwant,rwant):
    m=v-1
    for U in subgroups(m):
        lab,orbs=orbits(U,m)
        if len(U)!=dwant or len(orbs)!=rwant: continue
        sz=np.array([len(o) for o in orbs]); r=len(orbs); reps,_=lagclasses(U,m)
        Y=1-2*((np.arange(2**r)[:,None]>>np.arange(r))&1); rs=Y@sz
        def hist(a):
            Z=Y[rs==a]; X=Z[:,lab]; F=np.fft.fft(X.astype(float),axis=1)
            P=np.rint(np.fft.ifft(np.abs(F)**2,axis=1).real).astype(np.int64)[:,reps]
            return Counter(map(tuple,P)),len(Z)
        H2,N2=hist(2); H0,N0=hist(0)
        def conv(Ha,Hb):
            out=Counter()
            for a,ca in Ha.items():
                for b,cb in Hb.items(): out[tuple(x+y for x,y in zip(a,b))]+=ca*cb
            return out
        AB=conv(H2,H0); CD=conv(H0,H0)
        tgt=tuple([-4]*len(reps)); sol=0
        for vec,c in AB.items():
            need=tuple(t-x for t,x in zip(tgt,vec)); sol+=c*CD.get(need,0)
        return dict(v=v,d=dwant,r=r,N2=N2,N0=N0,distinct=(len(H2),len(H0)),exact_solutions=int(sol))
for v,d,r in [(167,41,6),(179,44,6),(179,22,10),(179,11,18)]:
    print(json.dumps(exact(v,d,r)),flush=True)
