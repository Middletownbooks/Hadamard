"""Partial-join validation of the b=1 Gaussian constraint model for a multiplier configuration (v,d,r).
Pools of P random U-invariant blocks at row sum 2 (slot A) and 0 (B,C,D, independent pools); for k=1..L,
exact count of 4-tuples with sum_k PAF = -4 on the first k lag classes (linear 64-bit hash of AB and CD sums,
sort/searchsorted), versus Gaussian prediction from the same pools' mean/covariance restricted to those classes."""
import sys, json, numpy as np
from math import log, log2, pi
sys.path.insert(0,'scripts')
from b1_multiplier_enum import subgroups, orbits, lagclasses
v,d,r,P=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
m=v-1; rng=np.random.default_rng(5)
for U in subgroups(m):
    lab,orbs=orbits(U,m)
    if len(U)==d and len(orbs)==r: break
sz=np.array([len(o) for o in orbs]); reps,csz=lagclasses(U,m); L=len(reps)
def pool(a,n):
    out=[]
    while len(out)<n:
        Y=rng.choice([-1,1],size=(40000,r)); out.extend(Y[(Y@sz)==a][:n-len(out)])
    X=np.array(out)[:,lab]; F=np.fft.fft(X.astype(float),axis=1)
    return np.rint(np.fft.ifft(np.abs(F)**2,axis=1).real).astype(np.int64)[:,reps]
A=pool(2,P); B=pool(0,P); C=pool(0,P); D=pool(0,P)
# drop classes with zero variance (deterministic) but check them
S=A[:1000]+B[:1000]+C[:1000]+D[:1000]
var=S.var(0); det=np.where(var==0)[0]
print(json.dumps(dict(v=v,d=d,r=r,L=L,P=P,deterministic_classes=det.tolist(),det_values=S[0,det].tolist())),flush=True)
order=[c for c in rng.permutation(L) if c not in set(det)]
mu=A.mean(0)+B.mean(0)+C.mean(0)+D.mean(0)+4
Cov=np.cov(A.T)+np.cov(B.T)+np.cov(C.T)+np.cov(D.T)
W=rng.integers(1,2**62,size=L,dtype=np.int64).astype(np.uint64)
for k in range(1,len(order)+1):
    idx=order[:k]; w=W[idx]
    hA=(A[:,idx].astype(np.uint64)*w).sum(1); hB=(B[:,idx].astype(np.uint64)*w).sum(1)
    hC=(C[:,idx].astype(np.uint64)*w).sum(1); hD=(D[:,idx].astype(np.uint64)*w).sum(1)
    tgt=(np.full(k,-4).astype(np.uint64)*w).sum()
    AB=(hA[:,None]+hB[None,:]).ravel(); CD=np.sort((tgt-(hC[:,None]+hD[None,:])).ravel())
    lo=np.searchsorted(CD,AB,'left'); hi=np.searchsorted(CD,AB,'right'); cnt=int((hi-lo).sum())
    Ck=Cov[np.ix_(idx,idx)]; mk=mu[idx]
    w_,V=np.linalg.eigh(Ck); keep=w_>1e-9*w_.max(); pr=V.T@mk
    lp=keep.sum()*log(4)-0.5*np.sum(np.log(2*pi*w_[keep]))-0.5*np.sum(pr[keep]**2/w_[keep])
    pred=4*log2(P)+lp/log(2)
    print(json.dumps(dict(k=k,matches=cnt,log2_meas=round(log2(cnt),2) if cnt else None,log2_pred=round(pred,2))),flush=True)
    if cnt==0 and pred<3: break
