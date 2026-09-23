import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d=int(sys.argv[1]),int(sys.argv[2]); rhos=[float(x) for x in sys.argv[3].split(',')]; lam=float(sys.argv[4]); seeds=int(sys.argv[5]); maxs=int(sys.argv[6]); ten=int(sys.argv[7])
R=Red(v,d); M=build(R); rng=np.random.default_rng(11)
for rho in rhos:
    F0=int(rho**2*v*R.L)   # sum_c PAF(c)^2 over L classes ~ rho^2 * m per class
    for s in range(seeds):
        y=random_y(R,rng); t=time.time(); st,ok=run_hinge(y,M,lam,F0,maxs,ten,s); dt=time.time()-t
        fp=fingerprint(R,y) if ok else (None,None)
        print(json.dumps(dict(v=v,d=d,alg='S3hinge',rho=rho,lam=lam,tenure=ten,seed=s,steps=st,solved=bool(ok),secs=round(dt,2),sv=fp[0],ratio=fp[1])),flush=True)
