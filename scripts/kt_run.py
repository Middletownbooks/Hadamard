"""Production k-tree join over tau-buckets.
usage: kt_run.py v d profile(csv) k Lmax gen_glob nbuckets seed
 - pools: union of generator outputs, split by |row sum|, subsampled to Lmax per row sum
 - per slot: chunk PAF (int8), bucket class c (int16), linear hash h of rest classes
 - buckets tau (AB chunk sum; CD uses -tau) enumerated best-first by estimated probability
Logs JSON lines; solutions to data/kt_solutions_v{v}.txt"""
import sys, glob, time, json, heapq, numpy as np
sys.path.insert(0,'scripts')
from redlib import Red
from kt_join import sparse_N, paf_from_bits, stage1, stage1_par, stage1_merge, merge_count_total, merge_thist, stage1_merge_t
v,d=int(sys.argv[1]),int(sys.argv[2]); prof=tuple(int(x) for x in sys.argv[3].split(','))
k=int(sys.argv[4]); Lmax=int(float(sys.argv[5])); files=sorted(glob.glob(sys.argv[6])); nbuck=int(sys.argv[7]); seed=int(sys.argv[8])
R=Red(v,d); Nsp=sparse_N(R); rng=np.random.default_rng(seed); t0=time.time()
def log(**kw): kw['secs']=round(time.time()-t0,1); print(json.dumps(kw),flush=True)
import os
cache=os.environ.get('KTCACHE','data/kt487/cache'); os.makedirs(cache,exist_ok=True)
def cpath(name): return os.path.join(cache,name+'.npy')
if all(os.path.exists(cpath(f'pool_{a}')) for a in set(prof)):
    pools={a:np.load(cpath(f'pool_{a}'),mmap_mode='r').view(np.ndarray) for a in set(prof)}
else:
    pools={}
    for a in set(prof):
        fs=[f for f in files if f'_rs{a}_' in f]
        x=np.unique(np.concatenate([np.load(f) for f in fs]+[np.zeros(0,np.uint64)]))
        # legacy single files (all row sums mixed)
        leg=[f for f in files if '_rs' not in f]
        for f in leg:
            g=np.load(f); ones=np.zeros(len(g),np.int64)
            for q in range(1,R.r): ones+=((g>>np.uint64(q))&np.uint64(1)).astype(np.int64)
            rs=np.where(g&np.uint64(1),1,-1)+d*(2*ones-(R.r-1)); x=np.unique(np.concatenate([x,g[rs==a]]))
        if len(x)>Lmax: x=x[rng.choice(len(x),Lmax,replace=False)]
        pools[a]=x
    
    log(stage='pools',sizes={int(a):len(x) for a,x in pools.items()},log2={int(a):round(float(np.log2(len(x))),2) for a,x in pools.items()})
    for a,x in pools.items(): np.save(cpath(f'pool_{a}'),x)
    pools={a:np.load(cpath(f'pool_{a}'),mmap_mode='r').view(np.ndarray) for a in set(prof)}
classes=rng.permutation(R.L); chunk=classes[:k]; cc=classes[k]; rest=classes[k+1:]
W=rng.integers(1,2**63,size=len(rest),dtype=np.uint64)|np.uint64(1); G=np.uint64(0x9E3779B97F4A7C15)
if all(os.path.exists(cpath(f'h_{a}_seed{seed}_k{k}')) for a in set(prof)):
    D={a:dict(ch=np.load(cpath(f'ch_{a}_seed{seed}_k{k}'),mmap_mode='r').view(np.ndarray),c=np.load(cpath(f'c_{a}_seed{seed}_k{k}')),h=np.load(cpath(f'h_{a}_seed{seed}_k{k}'),mmap_mode='r').view(np.ndarray)) for a in set(prof)}
else:
    D={}
    for a,x in pools.items():
        ch=np.empty((len(x),k),np.int8); c=np.empty(len(x),np.int16); h=np.empty(len(x),np.uint64)
        for s in range(0,len(x),1<<21):
            P=paf_from_bits(x[s:s+(1<<21)],*Nsp,R.L,R.r)
            ch[s:s+len(P)]=P[:,chunk]; c[s:s+len(P)]=P[:,cc]
            h[s:s+len(P)]=(P[:,rest].astype(np.int64).astype(np.uint64)*W[None,:]).sum(1)
        D[a]=dict(ch=ch,c=c,h=h)
    for a,X in D.items():
        for nm in ('ch','c','h'): np.save(cpath(f'{nm}_{a}_seed{seed}_k{k}'),X[nm])
    del D
    D={a:dict(ch=np.load(cpath(f'ch_{a}_seed{seed}_k{k}'),mmap_mode='r').view(np.ndarray),c=np.load(cpath(f'c_{a}_seed{seed}_k{k}')),h=np.load(cpath(f'h_{a}_seed{seed}_k{k}'),mmap_mode='r').view(np.ndarray)) for a in set(prof)}
log(stage='prep')
import numba as nb
@nb.njit(parallel=True,cache=True)
def keyof_nb(ch,tau,sgn):
    N=ch.shape[0]; K=ch.shape[1]; out=np.empty(N,np.uint64)
    for n in nb.prange(N):
        w=np.uint64(0)
        for i in range(K):
            val=np.int64(tau[i])+sgn*np.int64(ch[n,i])
            w|=np.uint64(val & 255)<<np.uint64(8*i)
        out[n]=w
    return out
KEFF=min(k,8)
def keyof_own(ch): return keyof_nb(ch[:,:KEFF],np.zeros(KEFF,np.int64),1)
def keyof_partner(ch,tau): return keyof_nb(ch[:,:KEFF],np.asarray(tau,np.int64)[:KEFF],-1)
# own keys for B, D slots
sB=keyof_own(np.asarray(D[prof[1]]['ch'])); oB=np.argsort(sB,kind='stable'); sB=sB[oB]
sD=keyof_own(np.asarray(D[prof[3]]['ch'])); oD=np.argsort(sD,kind='stable'); sD=sD[oD]
# per-class pair-sum distributions (AB and CD) from random pairs
ns=2_000_000
def pairdist(a1,a2):
    i=rng.integers(0,len(pools[a1]),ns); j=rng.integers(0,len(pools[a2]),ns)
    S=D[a1]['ch'][i].astype(np.int16)+D[a2]['ch'][j].astype(np.int16)
    return [dict(zip(*np.unique(S[:,q],return_counts=True))) for q in range(k)]
dAB=pairdist(prof[0],prof[1]); dCD=pairdist(prof[2],prof[3])
opts=[]
for q in range(k):
    vals=sorted(set(dAB[q])&{-x for x in dCD[q]},key=lambda t:-(dAB[q][t]/ns)*(dCD[q].get(-t,0)/ns))
    opts.append([(t,np.log2(dAB[q][t]/ns)+np.log2(dCD[q][-t]/ns)) for t in vals[:12]])
# best-first enumeration of tau vectors
best=tuple(0 for _ in range(k)); score=lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
heap=[(-score(best),best)]; seen={best}; buckets=[]
while heap and len(buckets)<nbuck:
    s,ix=heapq.heappop(heap); buckets.append(ix)
    for q in range(k):
        if ix[q]+1<len(opts[q]):
            nx=ix[:q]+(ix[q]+1,)+ix[q+1:]
            if nx not in seen: seen.add(nx); heapq.heappush(heap,(-score(nx),nx))
# rest-bits estimate for logging: per-class quad bits ~3.3 on (L-1-k) independent classes
out=open(f'data/kt_solutions_v{v}.txt','a'); nsol=0; exp_tot=0.0
for b,ix in enumerate(buckets):
    tau=np.array([opts[q][ix[q]][0] for q in range(k)],np.int16)
    nA=keyof_partner(np.asarray(D[prof[0]]['ch']),tau); nC=keyof_partner(np.asarray(D[prof[2]]['ch']),-tau)
    A,B,C,Dd=[D[a] for a in prof]
    oA=np.argsort(nA); nA=nA[oA]; oCc=np.argsort(nC); nC=nC[oCc]
    T=300
    hAB=merge_thist(nA,oA,A['c'],sB,oB,B['c'],64,T); hCD=merge_thist(nC,oCc,C['c'],sD,oD,Dd['c'],64,T)
    tv=np.arange(-T,T+1); both=hAB*hCD[::-1]          # AB at t, CD at -t
    cAB=int(hAB.sum()); cCD=int(hCD.sum()); exp_b=float((both.sum())*2.0**(-3.42*(R.L-1-KEFF-1)))
    log(stage='count',b=b,AB=cAB,CD=cCD,exp_hits_bucket=round(exp_b,3))
    # group t values into passes under the memory budget
    budget=5.0e9/24; passes=[]; cur=[];load=0
    for t in tv[np.argsort(-both)]:
        if both[t+T]==0: continue
        w=hAB[t+T]+hCD[T-t]
        if load+w>budget and cur: passes.append(cur); cur=[];load=0
        cur.append(int(t)); load+=w
    if cur: passes.append(cur)
    log(stage='passes',b=b,n=len(passes))
    found=0; ncand=0
    for pl in passes:
      for _once in (0,):
        mA=np.zeros(2*T+1,np.bool_); mC=np.zeros(2*T+1,np.bool_)
        for t in pl: mA[t+T]=True; mC[-t+T]=True
        KA,IA,JA=stage1_merge_t(nA,oA,A['c'],A['h'],sB,oB,B['c'],B['h'],1,G,64,mA,T)
        KC,IC,JC=stage1_merge_t(nC,oCc,C['c'],C['h'],sD,oD,Dd['c'],Dd['h'],-1,G,64,mC,T)
        oc=np.argsort(KC); KC=KC[oc]
        pos=np.searchsorted(KC,KA); pos[pos>=len(KC)]=len(KC)-1
        m=np.where(KC[pos]==KA)[0] if len(KC) and len(KA) else []; ncand+=len(m)
        for q in m:
            i,j=IA[q],JA[q]; pp=pos[q]
            while pp<len(KC) and KC[pp]==KA[q]:
                kk,l=IC[oc[pp]],JC[oc[pp]]; pp+=1
                bits=[pools[prof[0]][i],pools[prof[1]][j],pools[prof[2]][kk],pools[prof[3]][l]]
                Y=np.array([[1 if (int(bb)>>qq)&1 else -1 for qq in range(R.r)] for bb in bits])
                tot=R.paf(R.expand(Y)).sum(0)
                if np.all(tot[1:]==0):
                    found+=1; nsol+=1
                    out.write(json.dumps(dict(v=v,d=d,prof=prof,Y=[''.join('+' if t>0 else '-' for t in y) for y in Y]))+'\n'); out.flush()
        del KA,IA,JA,KC,IC,JC
    nAB=cAB; nCD=cCD; m=range(ncand)
    exp=exp_b
    exp_tot+=exp
    log(stage='bucket',b=b,tau=tau.tolist(),log2_prob=round(float(-score(ix)),1) if False else round(float(score(ix)),2),AB=int(nAB),CD=int(nCD),cand=int(len(m)),solutions=found,total=nsol,exp_hits_est=round(exp,3),exp_cum=round(exp_tot,2))
log(stage='done',solutions=nsol)
