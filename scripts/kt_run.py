"""Production k-tree join over tau-buckets.
usage: kt_run.py v d profile(csv) k Lmax gen_glob nbuckets seed
 - pools: union of generator outputs, split by |row sum|, subsampled to Lmax per row sum
 - per slot: chunk PAF (int8), bucket class c (int16), linear hash h of rest classes
 - buckets tau (AB chunk sum; CD uses -tau) enumerated best-first by estimated probability
Logs JSON lines; solutions to data/kt_solutions_v{v}.txt"""
import sys, glob, time, json, heapq, numpy as np
sys.path.insert(0,'scripts')
from redlib import Red
from kt_join import sparse_N, paf_from_bits, stage1
v,d=int(sys.argv[1]),int(sys.argv[2]); prof=tuple(int(x) for x in sys.argv[3].split(','))
k=int(sys.argv[4]); Lmax=int(float(sys.argv[5])); files=sorted(glob.glob(sys.argv[6])); nbuck=int(sys.argv[7]); seed=int(sys.argv[8])
R=Red(v,d); Nsp=sparse_N(R); rng=np.random.default_rng(seed); t0=time.time()
def log(**kw): kw['secs']=round(time.time()-t0,1); print(json.dumps(kw),flush=True)
g=np.unique(np.concatenate([np.load(f) for f in files]))
sgn=np.where((g[:,None]>>np.arange(R.r,dtype=np.uint64))&np.uint64(1),1,-1) if len(g)<0 else None
# row sums from bits (vectorised in blocks)
rs=np.empty(len(g),np.int64)
for s in range(0,len(g),1<<22):
    b=g[s:s+(1<<22)]; ones=np.zeros(len(b),np.int64)
    for q in range(1,R.r): ones+=((b>>np.uint64(q))&np.uint64(1)).astype(np.int64)
    y0=np.where(b&np.uint64(1),1,-1); rs[s:s+len(b)]=y0+d*(2*ones-(R.r-1))
pools={}
for a in set(prof):
    x=g[rs==a]
    if len(x)>Lmax: x=x[rng.choice(len(x),Lmax,replace=False)]
    pools[a]=x
del g,rs
log(stage='pools',sizes={int(a):len(x) for a,x in pools.items()},log2={int(a):round(float(np.log2(len(x))),2) for a,x in pools.items()})
classes=rng.permutation(R.L); chunk=classes[:k]; cc=classes[k]; rest=classes[k+1:]
W=rng.integers(1,2**63,size=len(rest),dtype=np.uint64)|np.uint64(1); G=np.uint64(0x9E3779B97F4A7C15)
D={}
for a,x in pools.items():
    ch=np.empty((len(x),k),np.int8); c=np.empty(len(x),np.int16); h=np.empty(len(x),np.uint64)
    for s in range(0,len(x),1<<21):
        P=paf_from_bits(x[s:s+(1<<21)],*Nsp,R.L,R.r)
        ch[s:s+len(P)]=P[:,chunk]; c[s:s+len(P)]=P[:,cc]
        h[s:s+len(P)]=(P[:,rest].astype(np.int64).astype(np.uint64)*W[None,:]).sum(1)
    D[a]=dict(ch=ch,c=c,h=h)
log(stage='prep')
def keyof(ch16):
    key=np.zeros(ch16.shape[0],np.uint64)
    for i in range(k): key|=(ch16[:,i].astype(np.int8).astype(np.uint8).astype(np.uint64))<<np.uint64(8*i)
    return key
# own keys for B, D slots
sB=keyof(D[prof[1]]['ch'].astype(np.int16)); oB=np.argsort(sB,kind='stable'); sB=sB[oB]
sD=keyof(D[prof[3]]['ch'].astype(np.int16)); oD=np.argsort(sD,kind='stable'); sD=sD[oD]
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
    nA=keyof(tau[None,:]-D[prof[0]]['ch'].astype(np.int16)); nC=keyof(-tau[None,:]-D[prof[2]]['ch'].astype(np.int16))
    A,B,C,Dd=[D[a] for a in prof]
    z=np.zeros(1,np.uint64); zi=np.zeros(1,np.int64)
    nAB=stage1(nA,A['c'],A['h'],sB,oB,B['c'],B['h'],-10**6,10**6,1,G,z,zi,zi,True)
    nCD=stage1(nC,C['c'],C['h'],sD,oD,Dd['c'],Dd['h'],-10**6,10**6,-1,G,z,zi,zi,True)
    if nAB*16+nCD*16>9e9: log(stage='skip',bucket=b,AB=nAB,CD=nCD); continue
    KA=np.empty(nAB,np.uint64); IA=np.empty(nAB,np.int64); JA=np.empty(nAB,np.int64)
    stage1(nA,A['c'],A['h'],sB,oB,B['c'],B['h'],-10**6,10**6,1,G,KA,IA,JA,False)
    KC=np.empty(nCD,np.uint64); IC=np.empty(nCD,np.int64); JC=np.empty(nCD,np.int64)
    stage1(nC,C['c'],C['h'],sD,oD,Dd['c'],Dd['h'],-10**6,10**6,-1,G,KC,IC,JC,False)
    oc=np.argsort(KC); KC=KC[oc]
    pos=np.searchsorted(KC,KA); pos[pos>=len(KC)]=len(KC)-1
    m=np.where(KC[pos]==KA)[0] if len(KC) else []
    found=0
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
    exp=nAB*nCD*2.0**(-3.3*(R.L-1-k))
    exp_tot+=exp
    log(stage='bucket',b=b,tau=tau.tolist(),log2_prob=round(float(-score(ix)),1) if False else round(float(score(ix)),2),AB=int(nAB),CD=int(nCD),cand=int(len(m)),solutions=found,total=nsol,exp_hits_est=round(exp,3),exp_cum=round(exp_tot,2))
    del KA,IA,JA,KC,IC,JC,oc,pos
log(stage='done',solutions=nsol)
