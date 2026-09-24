"""b=1 k-tree join over tau-buckets (GS cores of length m=v-1, row sums (2,0,0,0), sum_k PAF_k(s) = -4 for s != 0).
Pools are COMPLETE: every U-invariant block with the row sum and P(w) <= 4m+4 at all w (b1_pool_enum).
Slot A's PAF is shifted by +4 on every lag class, so the join target becomes 0 as in the b=0 join (kt_run.py).
usage: kt_b1.py v d r which k nbuckets seed
env: B0 (resume bucket), PLANT (core file: plant the real core, run its bucket only), HKEEP (log2 of the
     high-half fraction kept, <0 subsamples the pools; for planted tests), LOGM (bloom size)"""
import sys, os, time, json, heapq, numpy as np
from math import log, log2, pi
sys.path.insert(0,'scripts')
from b1red import find_U, RedC, canon, read_core
from b1_pool_enum import enum_pool, bits_to_Y
from kt_join import sparse_N, paf_from_bits, merge_thist, stage1_merge_t, bloom_match
import numba as nb
v,d,r,which,k,nbuck,seed=[int(x) for x in sys.argv[1:8]]
m=v-1; U,_=find_U(m,d,r,which); R=RedC(m,U); t0=time.time()
def log_(**kw): kw['secs']=round(time.time()-t0,1); print(json.dumps(kw),flush=True)
log_(stage='config',v=v,m=m,d=d,r=R.r,L=R.L,U=U[:8])
bound=4*m+4; hk=float(os.environ.get('HKEEP','0'))
pools={}
for a in (2,0):
    fn=f'data/b1pool_v{v}_d{d}_r{r}_w{which}_rs{a}.npy'
    if hk==0 and os.path.exists(fn): pools[a]=np.load(fn)
    else:
        keep=None
        if hk<0:
            n=1<<(R.r-R.r//2); keep=np.sort(np.random.default_rng(seed+a).choice(n,int(n*2**hk),replace=False))
        pools[a]=enum_pool(R,a,bound,keep)
        if hk==0: np.save(fn,pools[a])
plant=None
if os.environ.get('PLANT'):
    X=read_core(os.environ['PLANT']); rs=X.sum(1); o=np.argsort(-np.abs(rs)); X=X[o]; X[0]*=np.sign(X[0].sum())
    Y=canon(R,X); assert Y is not None
    pb=[np.uint64(int(sum(1<<q for q in range(R.r) if Y[i,q]>0))) for i in range(4)]
    for i,b_ in enumerate(pb):
        a=2 if i==0 else 0; pools[a]=np.unique(np.concatenate([pools[a],[b_]]))
    plant=pb
log_(stage='pools',log2={a:round(log2(len(x)),2) for a,x in pools.items()},sizes={a:len(x) for a,x in pools.items()})
# lag-class PAF values for every pool block
lagidx=R.lags
Nsp=sparse_N(R)
def pafs(bits): return paf_from_bits(bits,*Nsp,R.L,R.r)
rng=np.random.default_rng(seed)
cls=rng.permutation(R.L); chunk=cls[:k]; cc=cls[k]; rest=cls[k+1:]
W=rng.integers(1,2**63,size=len(rest),dtype=np.uint64)|np.uint64(1); G=np.uint64(0x9E3779B97F4A7C15)
Dd={}
for a,x in pools.items():
    P=pafs(x)
    if a==2: P+=4                                      # target -4 -> 0
    h=np.empty(len(x),np.uint64)
    for s_ in range(0,len(x),1<<20): h[s_:s_+(1<<20)]=(P[s_:s_+(1<<20)][:,rest].astype(np.int64).astype(np.uint64)*W[None,:]).sum(1)
    Dd[a]=dict(ch=np.ascontiguousarray(P[:,chunk]),c=np.ascontiguousarray(P[:,cc]),h=h,
               P=P[np.random.default_rng(1).integers(0,len(x),300000)].astype(np.int64))
    del P
# Gaussian estimate of match probability on the classes outside the chunk
def nbits(idx):
    mu=Dd[2]['P'][:,idx].mean(0)+3*Dd[0]['P'][:,idx].mean(0)
    C=np.atleast_2d(np.cov(Dd[2]['P'][:,idx].T.astype(float)))+3*np.atleast_2d(np.cov(Dd[0]['P'][:,idx].T.astype(float)))
    w,V=np.linalg.eigh(C); keep=w>1e-6*w.max(); pr=V.T@mu
    lp=keep.sum()*log(4)-0.5*np.sum(np.log(2*pi*w[keep]))-0.5*np.sum(pr[keep]**2/w[keep])
    return -lp/log(2)
n_all=nbits(np.arange(R.L)); n_chunk=nbits(chunk); n_rest=n_all-n_chunk
log_(stage='model',n_all=round(n_all,1),n_chunk=round(n_chunk,1),n_rest=round(n_rest,1),
     log2E=round(log2(len(pools[2]))+3*log2(len(pools[0]))-n_all,1),chunk=chunk.tolist(),cc=int(cc))
for a in Dd: del Dd[a]['P']
c0=m%4
@nb.njit(parallel=True,cache=True)
def keyof(ch,tau,sgn,c0):
    N=ch.shape[0]; K=ch.shape[1]; out=np.empty(N,np.uint64)
    for n in nb.prange(N):
        w=np.uint64(0)
        for i in range(K):
            val=np.int64(tau[i])+sgn*np.int64(ch[n,i])
            w|=np.uint64(((val-c0)>>2)&255)<<np.uint64(8*i)
        out[n]=w
    return out
# own keys (B, C, D all come from pool 0); partner keys: A -> tau - chA, C -> -tau - chC
z=np.zeros(k,np.int64)
sO=keyof(Dd[0]['ch'],z,1,c0); oO=np.argsort(sO,kind='stable'); sO=sO[oO]
ns=2_000_000
def pairdist(a1,a2):
    i=rng.integers(0,len(pools[a1]),ns); j=rng.integers(0,len(pools[a2]),ns)
    S=Dd[a1]['ch'][i].astype(np.int32)+Dd[a2]['ch'][j].astype(np.int32)
    return [dict(zip(*np.unique(S[:,q],return_counts=True))) for q in range(k)]
dAB=pairdist(2,0); dCD=pairdist(0,0)
opts=[]
for q in range(k):
    vals=sorted(set(dAB[q])&{-x for x in dCD[q]},key=lambda t:-(dAB[q][t]/ns)*(dCD[q].get(-t,0)/ns))
    opts.append([(int(t),log2(dAB[q][t]/ns)+log2(dCD[q][-t]/ns)) for t in vals[:24]])
score=lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
buckets=[]
if plant is not None:
    Pp=np.array([R.paf(R.expand(bits_to_Y(np.array([b_],np.uint64),R.r)))[0,lagidx] for b_ in plant]); Pp[0]+=4
    tp=(Pp[0]+Pp[1])[chunk]
    for q in range(k):
        if tp[q] not in [o[0] for o in opts[q]]: opts[q].append((int(tp[q]),-99.0))
    buckets=[tuple([o[0] for o in opts[q]].index(tp[q]) for q in range(k))]
    log_(stage='planted',tau=tp.tolist())
else:
    best=tuple(0 for _ in range(k)); heap=[(-score(best),best)]; seen={best}
    while heap and len(buckets)<nbuck:
        s,ix=heapq.heappop(heap); buckets.append(ix)
        for q in range(k):
            if ix[q]+1<len(opts[q]):
                nx=ix[:q]+(ix[q]+1,)+ix[q+1:]
                if nx not in seen: seen.add(nx); heapq.heappush(heap,(-score(nx),nx))
out=open(f'data/kt_b1_solutions_v{v}.txt','a'); nsol=0; exp_tot=0.0; B0=int(os.environ.get('B0','0'))
T=2*m+16; A=Dd[2]; O=Dd[0]; logm=int(os.environ.get('LOGM','28'))
for b,ix in enumerate(buckets):
    if b<B0: continue
    tau=np.array([opts[q][ix[q]][0] for q in range(k)],np.int64)
    nA=keyof(A['ch'],tau,-1,c0); oA=np.argsort(nA); nA=nA[oA]
    nC=keyof(O['ch'],-tau,-1,c0); oC=np.argsort(nC); nC=nC[oC]
    hAB=merge_thist(nA,oA,A['c'],sO,oO,O['c'],64,T); hCD=merge_thist(nC,oC,O['c'],sO,oO,O['c'],64,T)
    both=hAB*hCD[::-1]; cAB=int(hAB.sum()); cCD=int(hCD.sum())
    expb=float(cAB)*float(cCD)*2.0**(-n_rest)
    budget=4.0e9/24; passes=[]; cur=[]; load=0
    for t in np.arange(-T,T+1)[np.argsort(-both)]:
        if both[t+T]==0: continue
        w_=hAB[t+T]+hCD[T-t]
        if load+w_>budget and cur: passes.append(cur); cur=[]; load=0
        cur.append(int(t)); load+=w_
    if cur: passes.append(cur)
    found=0; ncand=0
    for pl in passes:
        mA=np.zeros(2*T+1,np.bool_); mC=np.zeros(2*T+1,np.bool_)
        for t in pl: mA[t+T]=True; mC[-t+T]=True
        KA,IA,JA=stage1_merge_t(nA,oA,A['c'],A['h'],sO,oO,O['c'],O['h'],1,G,64,mA,T)
        KC,IC,JC=stage1_merge_t(nC,oC,O['c'],O['h'],sO,oO,O['c'],O['h'],-1,G,64,mC,T)
        ia,ic,_=bloom_match(KA,KC,logm); ncand+=len(ia)
        for q,qc in zip(ia,ic):
            bits=np.array([pools[2][IA[q]],pools[0][JA[q]],pools[0][IC[qc]],pools[0][JC[qc]]],np.uint64)
            Y=bits_to_Y(bits,R.r); X=R.expand(Y); tot=R.paf(X).sum(0)
            if np.all(tot[1:]==-4):
                found+=1; nsol+=1
                out.write(json.dumps(dict(v=v,d=d,r=r,U=U,rowsums=X.sum(1).tolist(),
                    blocks=[''.join('+' if t_>0 else '-' for t_ in x) for x in X]))+'\n'); out.flush()
        del KA,IA,JA,KC,IC,JC
    exp_tot+=expb
    log_(stage='bucket',b=b,tau=tau.tolist(),log2p=round(score(ix),2),AB=cAB,CD=cCD,passes=len(passes),cand=int(ncand),
         solutions=found,total=nsol,exp_hits=round(expb,4),exp_cum=round(exp_tot,3))
log_(stage='done',solutions=nsol)
