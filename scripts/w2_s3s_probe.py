import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d=int(sys.argv[1]),int(sys.argv[2]); lam0=float(sys.argv[3]); taus=[float(x) for x in sys.argv[4].split(',')]; lmin=float(sys.argv[5]); seeds=int(sys.argv[6]); maxs=int(sys.argv[7]); ten=int(sys.argv[8])
R=Red(v,d); M=build(R); rng=np.random.default_rng(7)
for tau in taus:
    for s in range(seeds):
        y=random_y(R,rng); t=time.time(); st,ok=run_joint_sched(y,M,lam0,tau,lmin,maxs,ten,s); dt=time.time()-t
        fp=fingerprint(R,y) if ok else (None,None)
        print(json.dumps(dict(v=v,d=d,alg='S3sched',lam0=lam0,tau=tau,lam_min=lmin,tenure=ten,seed=s,steps=st,solved=bool(ok),secs=round(dt,2),sv=fp[0],ratio=fp[1])),flush=True)
