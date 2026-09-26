"""D-optimal transfer audit, checks P1-P7 (see dopt/PREREGISTRATION.md).  Objects are verified, not function returns:
every pair is rebuilt as the 2v x 2v matrix R and checked with exact integer arithmetic."""
import itertools, math, json, numpy as np

def seqs(v): return (((np.arange(2**v)[:,None]>>np.arange(v)[None,:])&1)*2-1).astype(np.int64)
def pafs(X): v=X.shape[1]; return np.stack([(X*np.roll(X,-s,axis=1)).sum(1) for s in range(v)],1)
def circ(x): v=len(x); return np.array([np.roll(x,i) for i in range(v)])       # row i = x shifted right by i
def bareiss_det(M):
    M=[[int(x) for x in r] for r in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n): M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
def pair_ok(a,b):                                   # the defining identity, checked on the circulant matrices
    v=len(a); A=circ(a); B=circ(b); return np.array_equal(A@A.T+B@B.T,(2*v-2)*np.eye(v,dtype=np.int64)+2)
def R_of(a,b): A=circ(a); B=circ(b); return np.block([[A,B],[-B.T,A.T]])
def R_ok(a,b,ew_base_off=2):
    v=len(a); n=2*v; R=R_of(a,b); L=(2*v-2)*np.eye(v,dtype=np.int64)+2
    gram=np.array_equal(R@R.T,np.block([[L,np.zeros((v,v),np.int64)],[np.zeros((v,v),np.int64),L]]))
    ew=(2*n-2)*(n-ew_base_off)**((n-2)//2)
    return gram and abs(bareiss_det(R))==ew

res={}
# P1 + P2 (exhaustive two-list hash join on PAF vectors, then object checks)
counts={}; allP2=True; example=None
for v in (7,9,11):
    X=seqs(v); P=pafs(X)[:,1:]; key={}
    for i,p in enumerate(map(tuple,P.tolist())): key.setdefault(p,[]).append(i)
    pairs=[(i,j) for i,p in enumerate(P.tolist()) for j in key.get(tuple(2-x for x in p),[])]
    counts[v]=len(pairs)
    for i,j in pairs:
        if not (pair_ok(X[i],X[j]) and R_ok(X[i],X[j])): allP2=False
    if pairs and example is None: example=(X[pairs[0][0]].copy(),X[pairs[0][1]].copy())
res['P1']=dict(counts=counts,PASS=counts=={7:392,9:1944,11:0})
res['P2']=dict(PASS=allP2 and sum(counts.values())>0)
# P3 negative tests
a,b=example; a2=a.copy(); a2[0]*=-1
res['P3']=dict(flip_pair_ok=bool(pair_ok(a2,b)),flip_R_ok=bool(R_ok(a2,b)),wrong_EW_ok=bool(R_ok(a,b,ew_base_off=1)),
               control_true_ok=bool(pair_ok(a,b) and R_ok(a,b)))
res['P3']['PASS']=(not res['P3']['flip_pair_ok']) and (not res['P3']['flip_R_ok']) and (not res['P3']['wrong_EW_ok']) and res['P3']['control_true_ok']
# P4
bad=0
for v in (3,5,7,9,11,13):
    X=seqs(v); P=pafs(X); a=X.sum(1)
    bad+=int(((P[:,1:]-v)%4!=0).sum())+int((P[:,1:].sum(1)!=a**2-v).sum())
res['P4']=dict(counterexamples=bad,PASS=bad==0)
# P5
mins={}
for v in (7,9):
    X=seqs(v); P=pafs(X)[:,1:(v-1)//2+1]; a=X.sum(1); best=None
    for sa in np.unique(a):
        for sb in np.unique(a):
            if sa*sa+sb*sb!=4*v-2: continue
            PA=P[a==sa]; PB=P[a==sb]
            c=((PA[:,None,:]+PB[None,:,:]-2)**2).sum(2); nz=c[c>0]
            if nz.size: best=nz.min() if best is None else min(best,nz.min())
    mins[v]=int(best)
res['P5']=dict(min_nonzero=mins,PASS=all(x==32 for x in mins.values()))
# P6
def two_sq(N): return any(math.isqrt(N-x*x)**2==N-x*x for x in range(math.isqrt(N)+1))
work=[487,523,571,659,739,643]; ok6={v:two_sq(4*v-2) for v in work}
res['P6']=dict(sum_of_two_squares=ok6,PASS=[v for v in work if ok6[v]]==[643])
wide=[515,523,537,571,573,599,643,647,659,669,719,721,739]; res['info_worklist13']={v:two_sq(4*v-2) for v in wide}
# P7
tab={(4,67):3.50,(4,79):3.56,(4,103):3.76,(4,151):4.02,(2,67):2.86,(2,79):3.00,(2,103):3.19,(2,151):3.47}
bmod=lambda k,v: 0.5*math.log2(2*math.pi*k*v)-2
dev={f'{k},{v}':round(bmod(k,v)-b,3) for (k,v),b in tab.items()}
res['P7']=dict(deviations=dev,b_4_523=round(bmod(4,523),3),PASS=all(abs(x)<=0.15 for x in dev.values()) and abs(bmod(4,523)-4.87)<=0.15)
res['zero_crossing_v']={k:round(2**(4*k+4)/(2*math.pi*k)) for k in (2,4)}
print(json.dumps(res,indent=1,default=str))
