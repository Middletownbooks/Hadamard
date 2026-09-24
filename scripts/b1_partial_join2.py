"""Deeper partial-join validation (chunked, low memory). Same statistics as b1_partial_join.py.
usage: b1_partial_join2.py v d r P seed"""
import sys, json, numpy as np
from math import log, log2, pi
sys.path.insert(0,'scripts')
from b1_multiplier_enum import subgroups, orbits, lagclasses
v,d,r,P,seed=[int(x) for x in sys.argv[1:6]]
m=v-1; rng=np.random.default_rng(seed)
for U in subgroups(m):
    lab,orbs=orbits(U,m)
    if len(U)==d and len(orbs)==r: break
sz=np.array([len(o) for o in orbs]); reps,_=lagclasses(U,m); L=len(reps)
def pool(a,n):
    out=[]
    while len(out)<n:
        Y=rng.choice([-1,1],size=(40000,r)); out.extend(Y[(Y@sz)==a][:n-len(out)])
    X=np.array(out)[:,lab]; res=[]
    for s in range(0,n,4096):
        F=np.fft.fft(X[s:s+4096].astype(float),axis=1); res.append(np.rint(np.fft.ifft(np.abs(F)**2,axis=1).real).astype(np.int16)[:,reps])
    return np.concatenate(res)
A=pool(2,P); B=pool(0,P); C=pool(0,P); D=pool(0,P)
mu=A.mean(0)+B.mean(0)+C.mean(0)+D.mean(0)+4
Cov=np.cov(A.T.astype(float))+np.cov(B.T.astype(float))+np.cov(C.T.astype(float))+np.cov(D.T.astype(float))
order=list(np.random.default_rng(5).permutation(L))
W=rng.integers(1,2**62,size=L,dtype=np.int64).astype(np.uint64)
print(json.dumps(dict(v=v,d=d,r=r,L=L,P=P,seed=seed)),flush=True)
prev=None
for k in range(1,L+1):
    idx=order[:k]; w=W[idx]
    h=[(X[:,idx].astype(np.int64).astype(np.uint64)*w).sum(1) for X in (A,B,C,D)]
    tgt=(np.full(k,-4).astype(np.uint64)*w).sum()
    CD=np.sort((tgt-(h[2][:,None]+h[3][None,:])).ravel())
    cnt=0
    for s in range(0,P,256):
        AB=(h[0][s:s+256,None]+h[1][None,:]).ravel()
        cnt+=int((np.searchsorted(CD,AB,'right')-np.searchsorted(CD,AB,'left')).sum())
    del CD
    Ck=Cov[np.ix_(idx,idx)]; mk=mu[idx]; w_,V=np.linalg.eigh(Ck); keep=w_>1e-9*w_.max(); pr=V.T@mk
    lp=keep.sum()*log(4)-0.5*np.sum(np.log(2*pi*w_[keep]))-0.5*np.sum(pr[keep]**2/w_[keep])
    pred=4*log2(P)+lp/log(2)
    print(json.dumps(dict(k=k,matches=cnt,log2_meas=round(log2(cnt),2) if cnt else None,log2_pred=round(pred,2),
                          offset=round(log2(cnt)-pred,2) if cnt else None)),flush=True)
    if cnt==0 and pred<1: break
