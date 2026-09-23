import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d=int(sys.argv[1]),int(sys.argv[2]); mus=[float(x) for x in sys.argv[3].split(',')]; inner=int(sys.argv[4]); seeds=int(sys.argv[5]); maxs=int(sys.argv[6]); ten=int(sys.argv[7]); gs=int(sys.argv[8])
R=Red(v,d); M=build(R); rng=np.random.default_rng(5)
for mu in mus:
    for s in range(seeds):
        y=random_y(R,rng); t=time.time(); st,ok=run_shared(y,M,mu,maxs,inner,ten,s,gs); dt=time.time()-t
        fp=fingerprint(R,y) if ok else (None,None)
        X=R.expand(y); pa=R.paf(X)[:,1:]; rr=int((pa.sum(0)**2).sum())
        print(json.dumps(dict(v=v,d=d,alg='S2',mu=mu,inner=inner,gs=gs,tenure=ten,seed=s,steps=st,solved=bool(ok),secs=round(dt,2),sv=fp[0],ratio=float(np.mean(pa.std(1))/np.sqrt(v)),final_R2=rr)),flush=True)
