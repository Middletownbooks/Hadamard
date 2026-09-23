import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d=int(sys.argv[1]),int(sys.argv[2]); mu=float(sys.argv[3]); shares=[float(x) for x in sys.argv[4].split(',')]; inners=[int(x) for x in sys.argv[5].split(',')]; seeds=int(sys.argv[6]); maxs=int(sys.argv[7]); ten=int(sys.argv[8]); gs=int(sys.argv[9])
R=Red(v,d); M=build(R); rng=np.random.default_rng(31)
for sh in shares:
  for inn in inners:
    for s in range(seeds):
        y=random_y(R,rng); t=time.time(); st,ok=run_shared2(y,M,mu,sh,maxs,inn,ten,s,gs); dt=time.time()-t
        fp=fingerprint(R,y) if ok else (None,None)
        print(json.dumps(dict(v=v,d=d,alg='S2fixed',mu=mu,share=sh,inner=inn,gs=gs,tenure=ten,seed=s,steps=st,solved=bool(ok),secs=round(dt,2),sv=fp[0],ratio=fp[1])),flush=True)
