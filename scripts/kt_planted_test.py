import sys,numpy as np,json; sys.path.insert(0,'scripts')
from hadlib import all_cores
from redlib import Red
from kt_join import run
R=Red(487,9); c=[x for x in all_cores() if x.v==487][0]
row=[r for r in json.load(open('logs/census.json')) if r['v']==487][0]
X=np.array([np.roll(c.X[k],-row['shifts'][k]) for k in range(4)])
Y=[]
for k in range(4):
    x=X[k]; x=x if x.sum()>0 else -x; Y.append(np.array([x[o[0]] for o in R.orbits]))
rsr=[int(R.expand(y).sum()) for y in Y]
i17=[i for i,a in enumerate(rsr) if a==17]; order=[rsr.index(37),i17[0],i17[1],rsr.index(1)]
Ys=[Y[i] for i in order]
def pack(y): return np.uint64(sum(1<<q for q in range(R.r) if y[q]>0))
g=np.load(sys.argv[1]); bits=((g[:,None]>>np.arange(55,dtype=np.uint64))&np.uint64(1)).astype(int); yy=2*bits-1; rs=yy[:,0]+9*yy[:,1:].sum(1)
pools={a:g[rs==a] for a in (37,17,1)}
for y in Ys: a=int(R.expand(y).sum()); pools[a]=np.unique(np.append(pools[a],pack(y)))
k=int(sys.argv[2]); seed=0
rng=np.random.default_rng(seed); classes=rng.permutation(R.L); chunk=classes[:k]
P=R.paf(R.expand(np.array(Ys)))[:,R.lags]
tau=(P[0]+P[1])[chunk]; print('tau from real AB',tau, ' CD',(P[2]+P[3])[chunk])
sols,tot=run(487,9,(37,17,17,1),k,pools,npasses=2,seed=seed,tau=tau)
real=[''.join('+' if t>0 else '-' for t in y) for y in Ys]
print('found real?',any(s==real for s in sols),'solutions',len(sols))
