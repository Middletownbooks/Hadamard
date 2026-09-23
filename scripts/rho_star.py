"""Join-balance test: for (v,d,profile), compute across rho:
 log2 N_flat(rho, a) per row sum a  = WL total count(rho) x pool fraction of row sum a (tabu pools; approximate)
 n'(rho) = Gaussian flat-subspace bits (pool covariance at that rho), x1.02 partial-join correction
 rho_exist: sum_k log2 N(a_k) >= n'  (flat solutions expected to exist)
 rho_join : min_k log2 N(a_k) >= n'/3 (k-tree 4-list join feasible with whole flat set as lists)"""
import sys, glob, json, numpy as np
sys.path.insert(0,'scripts')
from redlib import Red
v,d=int(sys.argv[1]),int(sys.argv[2]); prof=tuple(int(x) for x in sys.argv[3].split(',')); wlfile=sys.argv[4]
R=Red(v,d); m=v; L=R.L
if wlfile!='none':
    wl=json.load(open(wlfile)); ru=np.array(wl['ratio_upper']); cl=np.array(wl['log2_count_le'])
else:
    ru=np.array([0.0,10.0]); cl=np.array([np.nan,np.nan])
pools={}
for f in glob.glob(f'data/pool_v{v}_d{d}_rho*_s*.npz'):
    z=np.load(f)
    for a in z.files: pools.setdefault(int(a),[]).append(z[a])
rng=np.random.default_rng(0)
pools={a:np.concatenate(x) for a,x in pools.items()}
pools={a:np.unique(x if len(x)<=20000 else x[rng.choice(len(x),20000,replace=False)],axis=0) for a,x in pools.items()}
counts_all={a:len(x) for a,x in pools.items()}
paf={a:R.paf(R.expand(Y))[:,R.lags] for a,Y in pools.items()}
F={a:(P**2).sum(1) for a,P in paf.items()}
ratio={a:np.sqrt(F[a]/(m*L)) for a in F}
print('pool sizes',counts_all)
for rho in np.arange(0.26,0.47,0.02):
    sel={a:ratio[a]<=rho for a in ratio}
    n_by={a:int(s.sum()) for a,s in sel.items()}
    tot=sum(n_by.values())
    if tot<200 or any(n_by.get(a,0)<30 for a in set(prof)): continue
    i=np.searchsorted(ru,rho,'right')-1; log2N=cl[i]
    res=[];means={}
    for a,P in paf.items():
        Q=P[sel[a]]
        if len(Q)>=30: means[a]=Q.mean(0); res.append(Q-means[a])
    C=np.cov(np.concatenate(res).T)[:L-1,:L-1]*4
    mu=sum(means[a] for a in prof)[:L-1]
    sg,ld=np.linalg.slogdet(2*np.pi*C); lp=(L-1)*np.log(4)-0.5*ld-0.5*mu@np.linalg.solve(C,mu)
    nprime=-lp/np.log(2)*1.02
    lN={a:log2N+np.log2(n_by[a]/tot) for a in set(prof)}
    s4=sum(lN[a] for a in prof); mn=min(lN[a] for a in prof)
    print(json.dumps(dict(v=v,rho=round(float(rho),2),log2N_total=round(float(log2N),1),log2N_per_rowsum={a:round(float(x),1) for a,x in lN.items()},
        n_prime=round(float(nprime),1),log2_expected_flat_solutions=round(float(s4-nprime),1),join_margin=round(float(mn-nprime/3),1))))
