"""S4: RRR / divide-and-concur with thin per-block spectral sets.
x: (4, m) real, orbit-constant.  P_A: orbit-constant +-1 (sign of orbit mean).
P_B: per block, rescale |X_k(w)| to sqrt(T_k(w)), phases kept.
Targets T_k updated from the binary projection: T_k = P_k + (4v - sum_j P_j)/4 (clipped >= eps), start T_k = v.
usage: s4_rrr.py v d beta maxit seeds update_every"""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from redlib import Red
v,d,beta,maxit,seeds,upd=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]); eta=float(sys.argv[7]) if len(sys.argv)>7 else 1.0
R=Red(v,d); m=v
orb=R.orb; size=np.bincount(orb)
def PA(x):
    s=np.zeros((4,R.r)); 
    for k in range(4): s[k]=np.bincount(orb,weights=x[k],minlength=R.r)
    y=np.where(s>=0,1.0,-1.0); return y[:,orb], y
def PB(x,T):
    F=np.fft.fft(x,axis=1); a=np.abs(F); a[a<1e-12]=1e-12
    G=F/a*np.sqrt(T); return np.fft.ifft(G,axis=1).real
for seed in range(seeds):
    rng=np.random.default_rng(seed)
    x=rng.choice([-1.0,1.0],size=(4,R.r))[:,orb]
    T=np.full((4,m),float(v)); t0=time.time(); solved=False
    for it in range(maxit):
        xb,y=PA(x)
        P=np.abs(np.fft.fft(xb,axis=1))**2
        tot=P.sum(0)
        err=np.abs(tot[1:]-4*v).max()
        if err<1e-6: solved=True;break
        if upd and it%upd==0:
            T=(1-eta)*T+eta*np.maximum(P+(4*v-tot)[None,:]/4,1e-3)
        pb=PB(x,T)
        x=x+beta*(PA(2*pb-x)[0]-pb)
    out=dict(v=v,d=d,alg='S4rrr',beta=beta,upd=upd,eta=eta,seed=seed,iters=it+1,solved=solved,secs=round(time.time()-t0,1))
    if solved:
        S=P[:,1:]/(4*v); pa=R.paf(xb.astype(int))[:,1:]
        out.update(sv=float(np.mean((S-0.25)**2)),ratio=float(np.mean(pa.std(1))/np.sqrt(m)))
    else:
        out.update(final_maxerr=float(err))
    print(json.dumps(out),flush=True)
