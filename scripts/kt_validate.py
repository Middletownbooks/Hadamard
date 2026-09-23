"""Check stage-2 selectivity on trajectory pools: for the top tau bucket, count AB/CD pairs that also
match on j further classes (exact), vs independent-pair expectation from sampled pair-sum distributions."""
import sys,glob,numpy as np,json; sys.path.insert(0,'scripts')
from redlib import Red
from kt_join import sparse_N, paf_from_bits
v,d=487,9; R=Red(v,d); Nsp=sparse_N(R); rng=np.random.default_rng(0)
g=np.unique(np.concatenate([np.load(f) for f in glob.glob(sys.argv[1])]))
ones=np.zeros(len(g),np.int64)
for q in range(1,R.r): ones+=((g>>np.uint64(q))&np.uint64(1)).astype(np.int64)
rs=np.where(g&np.uint64(1),1,-1)+d*(2*ones-(R.r-1))
P={a:paf_from_bits(g[rs==a],*Nsp,R.L,R.r).astype(np.int16) for a in (37,17,1)}
k=6; cl=rng.permutation(R.L); chunk=cl[:k]; extra=cl[k:k+6]
tau=np.full(k,2)
def side(a1,a2,t):
    X=P[a1];Y=P[a2]
    ky=[tuple(r) for r in Y[:,chunk]]; idx={}
    for j,kk in enumerate(ky): idx.setdefault(kk,[]).append(j)
    out=[]
    for i in range(len(X)):
        for j in idx.get(tuple(t-X[i,chunk]),()): out.append(X[i,extra]+Y[j,extra])
    return np.array(out)
AB=side(37,17,tau); CD=side(17,1,-tau)
print('stage1',len(AB),len(CD))
for j in range(1,7):
    from collections import Counter
    ca=Counter(map(tuple,AB[:,:j])); cc=Counter(map(tuple,(-CD[:,:j])))
    m=sum(n*cc.get(key,0) for key,n in ca.items())
    # independence expectation: shuffle AB rows' components independently? use product of marginal collision prob per class from random pairs
    ia=rng.integers(0,len(P[37]),400000); ib=rng.integers(0,len(P[17]),400000); ic=rng.integers(0,len(P[17]),400000); idd=rng.integers(0,len(P[1]),400000)
    S=(P[37][ia]+P[17][ib]+P[17][ic]+P[1][idd])[:,extra[:j]]
    p=np.mean(np.all(S==0,axis=1))
    print(json.dumps(dict(j=j,matches=m,expected_indep=round(len(AB)*len(CD)*p,1),bits_meas=round(float(np.log2(len(AB)*len(CD)/max(m,1))),2),bits_indep=round(float(-np.log2(max(p,1e-12))),2))),flush=True)
