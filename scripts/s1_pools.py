"""build flat pools at (v,d): distinct blocks with ratio<=rho, per row sum. saves data/pool_v{v}_d{d}_rho{rho}.npz"""
import sys,time,numpy as np; sys.path.insert(0,'scripts')
from engine import Red,build,flat_pool
v,d,rho,secs,stride,seed=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6])
R=Red(v,d); M=build(R); F0=int(rho**2*R.m*R.L); rng=np.random.default_rng(seed)
store={}; t0=time.time(); it=0; out=np.zeros((20000,R.r),np.int64)
while time.time()-t0<secs:
    y=rng.choice([-1,1],size=(1,R.r)).astype(np.int64)
    n=flat_pool(y,M,F0,200000,max(3,R.r//10),seed*100000+it,stride,out,20000); it+=1
    for row in out[:n]:
        a=int(row[0]+d*row[1:].sum()); key=row.tobytes()
        if a<0: row=-row; a=-a; key=row.tobytes()
        store.setdefault(a,{})[key]=row.copy()
res={a:np.array(list(dct.values())) for a,dct in store.items()}
np.savez_compressed(f'data/pool_v{v}_d{d}_rho{rho}_s{seed}.npz',**{str(a):x for a,x in res.items()})
print({a:len(x) for a,x in sorted(res.items())}, 'walks',it,'secs',round(time.time()-t0))
