"""S1: constraint bits in the flat subspace.
Input: flat pools data/pool_v{v}_d{d}_rho{rho}_s*.npz (blocks with ratio<=rho, keyed by row sum).
(1) Gaussian n' = -log2 P(4 random flat blocks, profile-weighted by flat counts, form a solution):
    covariance pooled across row sums (per-row-sum mean removed); mean per class = sample mean per row sum
    (fallback: (a^2-m)/(m-1)).
(2) partial joins: with the pools of the real core's profile, count 4-tuples matching on the first k
    lag classes (random order), compare with the Gaussian marginal prediction for those k classes.
usage: s1_nprime.py v d rho profile(a,b,c,d) kmax"""
import sys, glob, json, itertools, numpy as np
sys.path.insert(0,'scripts')
from redlib import Red
v,d,rho=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3])
prof=tuple(int(x) for x in sys.argv[4].split(',')); kmax=int(sys.argv[5])
R=Red(v,d); m=v; L=R.L
pools={}
for f in glob.glob(f'data/pool_v{v}_d{d}_rho{rho}_s*.npz'):
    z=np.load(f)
    for a in z.files:
        pools.setdefault(int(a),[]).append(z[a])
pools={a:np.unique(np.concatenate(x),axis=0) for a,x in pools.items()}
MAXP=int(sys.argv[6]) if len(sys.argv)>6 else 10**9
_r=np.random.default_rng(1)
pools={a:(x if len(x)<=MAXP else x[_r.choice(len(x),MAXP,replace=False)]) for a,x in pools.items()}
paf={a:R.paf(R.expand(Y))[:,R.lags] for a,Y in pools.items()}
print('pool sizes',{a:len(x) for a,x in sorted(pools.items())})
# pooled covariance
res=[];means={}
for a,P in paf.items():
    if len(P)>=30: means[a]=P.mean(0); res.append(P-means[a])
C=np.cov(np.concatenate(res).T)
def mean_of(a): return means.get(a, np.full(L,(a*a-m)/(m-1)))
def logp_zero(profile,idx):
    mu=sum(mean_of(a) for a in profile)[idx]; Cs=4*C[np.ix_(idx,idx)]
    sgn,ld=np.linalg.slogdet(2*np.pi*Cs)
    return len(idx)*np.log(4)-0.5*ld-0.5*mu@np.linalg.solve(Cs,mu)
full=list(range(L-1))
print(json.dumps(dict(v=v,d=d,rho=rho,profile=prof,n_prime_given_profile=round(float(-logp_zero(prof,full)/np.log(2)),1),
                      per_class_bits=round(float(-logp_zero(prof,full)/np.log(2)/(L-1)),2))))
# partial joins on the real profile
Ps=[paf[a] for a in prof]
rng=np.random.default_rng(0); order=rng.permutation(L-1)
for k in range(1,kmax+1):
    idx=order[:k]
    W=rng.integers(1,2**62,size=k,dtype=np.int64).astype(np.uint64)
    h=[(P[:,idx].astype(np.uint64)*W).sum(1) for P in Ps]
    AB=(h[0][:,None]+h[1][None,:]).ravel(); CD=np.sort((np.uint64(0)-(h[2][:,None]+h[3][None,:])).ravel())
    lo=np.searchsorted(CD,AB,'left'); hi=np.searchsorted(CD,AB,'right'); cnt=int((hi-lo).sum())
    tot=np.prod([float(len(P)) for P in Ps])
    pred=-logp_zero(prof,list(idx))/np.log(2)
    meas=-np.log2(cnt/tot) if cnt else float('inf')
    print(json.dumps(dict(k=k,matches=cnt,tuples_log2=round(float(np.log2(tot)),1),bits_measured=round(float(meas),2),bits_gaussian=round(float(pred),2))),flush=True)
