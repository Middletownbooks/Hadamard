import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d=int(sys.argv[1]),int(sys.argv[2]); lam=float(sys.argv[3]); wincs=[float(x) for x in sys.argv[4].split(',')]; seeds=int(sys.argv[5]); maxs=int(sys.argv[6]); ten=int(sys.argv[7])
R=Red(v,d); M=build(R); rng=np.random.default_rng(21)
for wi in wincs:
    for s in range(seeds):
        y=random_y(R,rng); t=time.time(); st,ok=run_breakout(y,M,lam,maxs,ten,s,wi); dt=time.time()-t
        fp=fingerprint(R,y) if ok else (None,None)
        print(json.dumps(dict(v=v,d=d,alg='S3breakout',lam=lam,winc=wi,tenure=ten,seed=s,steps=st,solved=bool(ok),secs=round(dt,2),sv=fp[0],ratio=fp[1])),flush=True)
