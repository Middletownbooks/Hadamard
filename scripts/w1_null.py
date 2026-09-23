"""W1: unbiased (pool-)MITM sampler of b=0 d-reduced solutions and their fingerprint stats.
For each profile, pools of uniform blocks per slot (or the whole population if smaller than N:
then the enumeration is exhaustive).  Linear 64-bit hash of the lag-class PAF vector makes
the AB/CD join an outer sum.  Every solution in pool product is found -> unbiased.
usage: w1_null.py v d N seconds"""
import sys, time, json, numpy as np
sys.path.insert(0,'scripts')
from redlib import Red
from dirichlet_null import cdf
from scipy.stats import chi2
v,d,N,secs=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]); MAXHIT=int(sys.argv[5]) if len(sys.argv)>5 else 400
R=Red(v,d); m=v; rng=np.random.default_rng(int(time.time()*1e3)%2**32)
W=rng.integers(0,2**63,size=R.L,dtype=np.uint64)|np.uint64(1)
def enum_all(a):
    from itertools import combinations
    out=[]
    for y0 in (1,-1):
        t=a-y0
        if t%d: continue
        k=t//d+(R.r-1)
        if k%2: continue
        k//=2
        if not 0<=k<=R.r-1: continue
        for S in combinations(range(1,R.r),k):
            y=-np.ones(R.r,int); y[0]=y0; y[list(S)]=1; out.append(y)
    return np.array(out)
def pool(a):
    c=R.count_rowsum(a)
    if c<=N: return enum_all(a),True
    return R.sample(a,N,rng),False
def stats(X):
    P=R.psd(X); S=P[:,1:]/(4*v); sv=float(np.mean((S-0.25)**2))
    Pa=R.paf(X)[:,1:]; ratio=float(np.mean(Pa.std(1))/np.sqrt(m))
    lp=0.0; Ts=[]
    for x in R.lags:
        s=P[:,x]/P[:,x].sum(); t=float(((s-0.25)**2).sum()); Ts.append(t); lp+=float(np.log10(cdf(t)[0]))
    fisher=float(chi2.logsf(-2*lp*np.log(10),2*R.L)/np.log(10))
    return sv,ratio,lp,fisher,Ts
profiles=R.profiles(v)
import os
if os.environ.get('PROFILE'): profiles=[tuple(int(x) for x in os.environ['PROFILE'].split(','))]
OUT=os.environ.get('OUT',f'logs/w1_v{v}_d{d}.jsonl')
out=open(OUT,'a'); t0=time.time(); rnd=0
while time.time()-t0<secs:
    for prof in profiles:
        pools=[];exh=True
        for a in prof:
            Y,e=pool(a); pools.append(Y); exh&=e
        H=[]; PAFs=[]
        for Y in pools:
            pf=R.paf(R.expand(Y))[:,R.lags]; PAFs.append(pf)
            H.append((pf.astype(np.uint64)*W[None,:]).sum(1))   # wraps mod 2^64
        AB=(H[0][:,None]+H[1][None,:]).ravel()
        CD=(np.uint64(0)-(H[2][:,None]+H[3][None,:])).ravel()
        order=np.argsort(CD,kind='stable'); CDs=CD[order]
        left=np.searchsorted(CDs,AB,'left'); right=np.searchsorted(CDs,AB,'right')
        cnt=right-left; hit=np.where(cnt>0)[0]; nhit=int(cnt.sum())   # all hash matches
        cum=np.cumsum(cnt[hit])
        pick=np.sort(rng.choice(nhit,size=min(MAXHIT,nhit),replace=False)) if nhit else []
        nsol=0
        for q in pick:          # uniform over all matches
            h=int(np.searchsorted(cum,q,'right')); ab=hit[h]; off=q-(cum[h]-cnt[ab])
            i,j=divmod(int(ab),len(H[1])); k,l=divmod(int(order[left[ab]+off]),len(H[3]))
            if np.any(PAFs[0][i]+PAFs[1][j]+PAFs[2][k]+PAFs[3][l]): continue
            Y=np.array([pools[0][i],pools[1][j],pools[2][k],pools[3][l]]); X=R.expand(Y)
            tot=R.paf(X).sum(0); assert np.all(tot[1:]==0)
            sv,ratio,lp,fi,Ts=stats(X); nsol+=1
            out.write(json.dumps(dict(v=v,d=d,prof=prof,round=rnd,exh=exh,nmatch=nhit,sv=sv,ratio=ratio,lp=lp,fisher=fi,T=Ts,
                Y=[''.join('+' if t>0 else '-' for t in y) for y in Y]))+'\n')
        out.flush()
        sizes=[len(p) for p in pools]
        print(json.dumps(dict(round=rnd,prof=prof,sizes=sizes,exh=exh,nsol=nsol,nhit=nhit,pairs=len(AB),t=round(time.time()-t0,1))),flush=True)
    rnd+=1
    if all(R.count_rowsum(a)<=N for p in profiles for a in p): break
