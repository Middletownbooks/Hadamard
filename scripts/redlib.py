"""Reduced (U-orbit) setting for b=0, prime m, U = subgroup of order d of (Z/m)*."""
import numpy as np
from math import comb
def primroot(p):
    fs=[q for q in range(2,p) if (p-1)%q==0 and all(q%r for r in range(2,int(q**.5)+1))]
    for g in range(2,p):
        if all(pow(g,(p-1)//q,p)!=1 for q in fs): return g
class Red:
    def __init__(s,m,d):
        assert (m-1)%d==0
        s.m=m; s.d=d; g=primroot(m); s.g=g
        h=pow(g,(m-1)//d,m); s.U=sorted({pow(h,k,m) for k in range(d)})
        s.orb=np.zeros(m,int); s.orbits=[[0]]
        for j in range((m-1)//d):
            o=sorted({g**j*u%m for u in s.U} if False else {pow(g,j,m)*u%m for u in s.U})
            for i in o: s.orb[i]=len(s.orbits)
            s.orbits.append(o)
        s.r=len(s.orbits)
        s.size=np.array([len(o) for o in s.orbits])
        # lag classes: orbits of <U,-1> on nonzero
        seen=set(); s.lags=[]
        for j in range((m-1)//d):
            x=pow(g,j,m)
            if x in seen: continue
            cl={x*u%m for u in s.U}|{(-x*u)%m for u in s.U}; seen|=cl; s.lags.append(x)
        s.lags=np.array(s.lags); s.L=len(s.lags)
        # frequency classes are the same sets (multiplier action on Z_m)
    def expand(s,Y):  # Y: (..., r) orbit signs -> (..., m)
        return Y[...,s.orb]
    def paf(s,X):
        F=np.fft.fft(X.astype(float),axis=-1); R=np.fft.ifft(np.abs(F)**2,axis=-1).real
        return np.rint(R).astype(np.int64)
    def psd(s,X):
        return np.abs(np.fft.fft(X.astype(float),axis=-1))**2
    def rowsums(s):
        """admissible row sums: y0 + d*t"""
        out=set()
        for y0 in (1,-1):
            for k in range(s.r):  # k = #(+1) among r-1 nonzero orbits
                if k<=s.r-1: out.add(y0+s.d*(2*k-(s.r-1)))
        return sorted(out)
    def count_rowsum(s,a):
        c=0
        for y0 in (1,-1):
            t=(a-y0)
            if t%s.d: continue
            t//=s.d; k=(t+(s.r-1))
            if k%2: continue
            k//=2
            if 0<=k<=s.r-1: c+=comb(s.r-1,k)
        return c
    def sample(s,a,n,rng):
        """n uniform orbit-constant blocks with row sum a; returns (n, r) orbit signs"""
        opts=[];w=[]
        for y0 in (1,-1):
            t=a-y0
            if t%s.d: continue
            t//=s.d; k=t+(s.r-1)
            if k%2: continue
            k//=2
            if 0<=k<=s.r-1: opts.append((y0,k)); w.append(comb(s.r-1,k))
        w=np.array(w,float); w/=w.sum()
        Y=-np.ones((n,s.r),int)
        ch=rng.choice(len(opts),size=n,p=w)
        for i in range(n):
            y0,k=opts[ch[i]]; Y[i,0]=y0
            Y[i,1+rng.choice(s.r-1,size=k,replace=False)]=1
        return Y
    def profiles(s,v):
        rs=[a for a in s.rowsums() if a>=0]
        out=set()
        for a in rs:
            for b in rs:
                for c in rs:
                    dd=4*v-a*a-b*b-c*c
                    if dd<0: continue
                    e=int(round(dd**.5))
                    if e*e==dd and e in rs: out.add(tuple(sorted((a,b,c,e),reverse=True)))
        return sorted(out)
