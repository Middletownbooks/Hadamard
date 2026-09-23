"""S1 measurements at a given (v,d):
 (1) Wang-Landau density of states of a single orbit-constant block over F = sum_c PAF(c)^2
     -> log2 #blocks with PAF-std ratio <= rho (any row sum), and per row sum via final histogram.
 (2) flat pools (tabu trajectories, fixed row sum profile) + partial joins on k lag classes:
     measured match rate per class vs Gaussian prediction -> n' (bits) for the flat subspace.
usage: s1_flat_measure.py v d mode[wl|join] ..."""
import sys, json, time, numpy as np, numba as nb
sys.path.insert(0,'scripts')
from engine import Red, build, init_state, flip

@nb.njit(cache=True)
def wang_landau(y,M,edges,nsweeps_per_check,lnf_final,seed,maxflips):
    np.random.seed(seed)
    r=y.shape[1]; L=M.shape[0]; nb_=len(edges)-1
    G,P=init_state(y,M)
    F=0
    for c in range(L): F+=P[0,c]*P[0,c]
    lng=np.zeros(nb_); H=np.zeros(nb_,np.int64)
    # rowsum tracking: counts of (bin, rowsum) visits in the final stage
    rs=0
    sz=np.ones(r,np.int64)
    lnf=1.0; flips=0
    def_bin=-1
    b=np.searchsorted(edges,F)-1
    while lnf>lnf_final and flips<maxflips:
        for it in range(nsweeps_per_check):
            q=np.random.randint(r); yq=y[0,q]; dF=0
            for c in range(L):
                dd=-2*yq*G[0,c,q]+2*M[c,q,q]; dF+=2*P[0,c]*dd+dd*dd
            Fn=F+dF; bn=np.searchsorted(edges,Fn)-1
            if bn<0 or bn>=nb_:
                accept=False
            else:
                accept = (b<0 or b>=nb_) or np.log(np.random.random())<lng[b]-lng[bn]
            if accept:
                flip(y,G,P,M,0,q); F=Fn; b=bn
            if 0<=b<nb_:
                lng[b]+=lnf; H[b]+=1
            flips+=1
        # flatness check
        mn=1<<62; tot=0; cnt=0
        for i in range(nb_):
            if lng[i]<=0: continue       # only bins discovered so far
            tot+=H[i]; cnt+=1
            if H[i]<mn: mn=H[i]
        if cnt>0 and mn>0.8*tot/cnt:
            lnf/=2.0; H[:]=0
    return lng,flips,lnf

def run_wl(v,d,seed,maxflips):
    R=Red(v,d); M=build(R)
    # F in units: ratio^2 * m * L ; bins from ratio 0.24 to 1.3
    lo,hi=float(sys.argv[6]) if len(sys.argv)>6 else 0.24,1.30
    edges=np.unique(np.rint(np.linspace(lo**2,hi**2,121)*R.m*R.L).astype(np.int64))
    # start from a typical state in range
    rng=np.random.default_rng(seed)
    y=rng.choice([-1,1],size=(1,R.r)).astype(np.int64)
    t=time.time(); lng,flips,lnf=wang_landau(y,M,edges,2_000_000,float(sys.argv[7]) if len(sys.argv)>7 else 1e-5,seed,maxflips)
    # normalise: total count of blocks 2^r ; high bins cover typical mass -> use sum over all bins + estimate of mass outside
    # we normalise to the upper range by matching exact typical fraction via direct sampling
    Xs=rng.choice([-1,1],size=(200000,R.r)); Ps=R.paf(R.expand(Xs))[:,R.lags]; Fs=(Ps**2).sum(1)
    hist=np.histogram(Fs,bins=edges)[0]/len(Fs)
    ok=(hist>2e-3)&(lng>0)
    ln_norm=np.log(np.sum(hist[ok]))-np.log(np.sum(np.exp(lng[ok]-lng[ok].max())))-lng[ok].max()   # ln(fraction) = lng + ln_norm
    ln_frac=lng+ln_norm
    ratio_edges=np.sqrt(edges/(R.m*R.L))
    cum=np.log2(np.cumsum(np.where(lng>0,np.exp(ln_frac),0))+1e-300)+R.r   # log2 #blocks with F <= upper edge
    return dict(v=v,d=d,seed=seed,r=R.r,flips=int(flips),lnf=float(lnf),secs=round(time.time()-t,1),
                ratio_upper=[float(x) for x in ratio_edges[1:]],log2_count_le=[float(x) for x in cum],visited=[bool(x) for x in (lng>0)])
if __name__=='__main__':
    mode=sys.argv[3]; v,d=int(sys.argv[1]),int(sys.argv[2])
    if mode=='wl':
        res=run_wl(v,d,int(sys.argv[4]),int(float(sys.argv[5])))
        print(json.dumps(res))

def run_wl_range(v,d,seed,maxflips,lo,hi,lnf_final,nbins=60):
    """restricted-range WL: returns raw ln g on [lo,hi] (ratio units); normalise offline by overlap"""
    R=Red(v,d); M=build(R)
    edges=np.unique(np.rint(np.linspace(lo**2,hi**2,nbins+1)*R.m*R.L).astype(np.int64))
    rng=np.random.default_rng(seed)
    # start inside range: tabu-descend from random until F<=upper edge
    from engine import flat_pool
    y=rng.choice([-1,1],size=(1,R.r)).astype(np.int64); out=np.zeros((1,R.r),np.int64)
    flat_pool(y,M,int(edges[-1])-1,200000,max(3,R.r//10),seed,1,out,1); y=out.copy()
    t=time.time(); lng,flips,lnf=wang_landau(y,M,edges,2_000_000,lnf_final,seed,maxflips)
    return dict(v=v,d=d,seed=seed,lo=lo,hi=hi,edges=[int(e) for e in edges],ratio_upper=[float(x) for x in np.sqrt(edges[1:]/(R.m*R.L))],
                lng=[float(x) for x in lng],flips=int(flips),lnf=float(lnf),secs=round(time.time()-t,1))
if __name__=='__main__' and sys.argv[3]=='wlr':
    print(json.dumps(run_wl_range(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[4]),int(float(sys.argv[5])),float(sys.argv[6]),float(sys.argv[7]),float(sys.argv[8]))))
