import sys,time,json,numpy as np; sys.path.insert(0,'scripts')
from engine import *
v,d,steps,seeds,mode=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])
a=int(sys.argv[6]) if len(sys.argv)>6 else None
R=Red(v,d); M=build(R); rng=np.random.default_rng(3)
for s in range(seeds):
    y=(R.sample(a,1,rng) if a is not None else random_y(R,rng)[:1]).astype(np.int64)
    t=time.time(); F=flat_tabu(y,M,steps,max(3,R.r//10),s,mode); dt=time.time()-t
    X=R.expand(y); pa=R.paf(X)[0,1:]
    print(json.dumps(dict(v=v,d=d,mode=('swap' if mode else 'flip'),rowsum=int(X.sum()),steps=steps,seed=s,paf_std=round(float(pa.std()),3),ratio=round(float(pa.std()/np.sqrt(v)),3),maxabs=int(np.abs(pa).max()),secs=round(dt,1))),flush=True)
