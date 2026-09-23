import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d=int(sys.argv[1]),int(sys.argv[2]); lams=[float(x) for x in sys.argv[3].split(',')]; seeds=int(sys.argv[4]); maxs=int(sys.argv[5]); ten=int(sys.argv[6]) if len(sys.argv)>6 else 8
R=Red(v,d); M=build(R); rng=np.random.default_rng(1)
for lam in lams:
    for s in range(seeds):
        y=random_y(R,rng); t=time.time(); st,ok=run_joint(y,M,lam,maxs,ten,s,0.0); dt=time.time()-t
        fp=fingerprint(R,y) if ok else (None,None)
        print(json.dumps(dict(v=v,d=d,alg='S3tabu',lam=lam,tenure=ten,seed=s,steps=st,solved=bool(ok),secs=round(dt,2),sv=fp[0],ratio=fp[1])),flush=True)
