"""Generic probabilistic k-tree (Wagner) 4-list join:  find (a,b,c,d) with PAF0[a]+PAF1[b]+PAF2[c]+PAF3[d] = 0 on all
classes (list 0 already carries -target).  Buckets tau on k 'chunk' classes (AB chunk sum = tau, CD = -tau) are
processed best-first (empirical pair-sum probabilities); within a bucket one more class (cc) splits memory passes;
the remaining classes are matched through a 64-bit linear hash (bloom_match); candidates are verified exactly.
lib: kt_join4(PAF, verify, k=6, nbuckets=10**6, max_solutions=10, log=print, seed=0, budget=1e8)"""
import sys, json, time, heapq, numpy as np, numba as nb
sys.path.insert(0,'scripts')
from kt_join import merge_thist, stage1_merge_t, bloom_match

@nb.njit(parallel=True,cache=True)
def keyof(ch,tau,sgn):
    N=ch.shape[0]; K=ch.shape[1]; out=np.empty(N,np.uint64)
    for n in nb.prange(N):
        w=np.uint64(0)
        for i in range(K):
            val=np.int64(tau[i])+sgn*np.int64(ch[n,i])
            w|=np.uint64(val&255)<<np.uint64(8*i)
        out[n]=w
    return out

def kt_join4(PAF,verify,k=6,nbuckets=10**6,max_solutions=10,log=print,seed=0,budget=1e8,tmax=1e9,target_pairs=1e7):
    rng=np.random.default_rng(seed); L=PAF[0].shape[1]
    spread=[np.ptp(PAF[1][:,c]) for c in range(L)]
    cls=np.argsort(spread)[::-1]
    # adaptive k: per-class probability of the modal AB pair-sum value, from samples
    ns0=200000; i0=rng.integers(0,len(PAF[0]),ns0); j0=rng.integers(0,len(PAF[1]),ns0)
    pmode=[]
    for c in cls:
        vals,cnt=np.unique(PAF[0][i0,c]+PAF[1][j0,c],return_counts=True); pmode.append(cnt.max()/ns0)
    NAB=float(len(PAF[0]))*len(PAF[1]); kk=0; load=NAB
    while kk<min(8,L-2) and load*pmode[kk]>=target_pairs: load*=pmode[kk]; kk+=1
    k=max(1,kk)
    log(json.dumps(dict(stage='kt_setup',L=int(L),k=int(k),best_bucket_AB_est=float(load),pmode=[round(float(x),4) for x in pmode[:10]])))
    chunk=cls[:k]; cc=cls[k]; rest=cls[k+1:]
    W=rng.integers(1,2**63,size=max(len(rest),1),dtype=np.uint64)|np.uint64(1); G=np.uint64(0x9E3779B97F4A7C15)
    D=[]
    for P in PAF:
        h=(P[:,rest].astype(np.int64).astype(np.uint64)*W[:len(rest)]).sum(1) if len(rest) else np.zeros(len(P),np.uint64)
        D.append(dict(ch=np.ascontiguousarray(P[:,chunk]).astype(np.int16),c=np.ascontiguousarray(P[:,cc]).astype(np.int16),h=h))
    z=np.zeros(k,np.int64)
    sB=keyof(D[1]['ch'],z,1); oB=np.argsort(sB,kind='stable'); sB=sB[oB]
    sD=keyof(D[3]['ch'],z,1); oD=np.argsort(sD,kind='stable'); sD=sD[oD]
    ns=2_000_000
    def pairdist(a1,a2):
        i=rng.integers(0,len(PAF[a1]),ns); j=rng.integers(0,len(PAF[a2]),ns)
        S=D[a1]['ch'][i].astype(np.int32)+D[a2]['ch'][j].astype(np.int32)
        return [dict(zip(*np.unique(S[:,q],return_counts=True))) for q in range(k)]
    dAB=pairdist(0,1); dCD=pairdist(2,3); opts=[]
    for q in range(k):
        vals=sorted(set(dAB[q])&{-x for x in dCD[q]},key=lambda t:-(dAB[q][t]/ns)*(dCD[q].get(-t,0)/ns))
        opts.append([(int(t),np.log2(dAB[q][t]/ns)+np.log2(dCD[q][-t]/ns)) for t in vals[:24]])
    if any(len(o)==0 for o in opts): log(json.dumps(dict(stage='kt',note='no common chunk values'))); return []
    score=lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
    best=tuple(0 for _ in range(k)); heap=[(-score(best),best)]; seen={best}
    T=int(max(np.abs(D[0]['c']).max()+np.abs(D[1]['c']).max(),np.abs(D[2]['c']).max()+np.abs(D[3]['c']).max()))+2
    found=[]; t0=time.time(); nb_=0; ncand=0
    while heap and nb_<nbuckets and len(found)<max_solutions and time.time()-t0<tmax:
        s_,ix=heapq.heappop(heap); nb_+=1
        for q in range(k):
            if ix[q]+1<len(opts[q]):
                nx=ix[:q]+(ix[q]+1,)+ix[q+1:]
                if nx not in seen: seen.add(nx); heapq.heappush(heap,(-score(nx),nx))
        tau=np.array([opts[q][ix[q]][0] for q in range(k)],np.int64)
        nA=keyof(D[0]['ch'],tau,-1); oA=np.argsort(nA); nA=nA[oA]
        nC=keyof(D[2]['ch'],-tau,-1); oC=np.argsort(nC); nC=nC[oC]
        hAB=merge_thist(nA,oA,D[0]['c'],sB,oB,D[1]['c'],64,T); hCD=merge_thist(nC,oC,D[2]['c'],sD,oD,D[3]['c'],64,T)
        both=hAB*hCD[::-1]
        passes=[]; cur=[]; load=0
        for t in np.arange(-T,T+1)[np.argsort(-both)]:
            if both[t+T]==0: continue
            w_=hAB[t+T]+hCD[T-t]
            if load+w_>budget and cur: passes.append(cur); cur=[]; load=0
            cur.append(int(t)); load+=w_
        if cur: passes.append(cur)
        for pl in passes:
            mA=np.zeros(2*T+1,np.bool_); mC=np.zeros(2*T+1,np.bool_)
            for t in pl: mA[t+T]=True; mC[-t+T]=True
            KA,IA,JA=stage1_merge_t(nA,oA,D[0]['c'],D[0]['h'],sB,oB,D[1]['c'],D[1]['h'],1,G,64,mA,T)
            KC,IC,JC=stage1_merge_t(nC,oC,D[2]['c'],D[2]['h'],sD,oD,D[3]['c'],D[3]['h'],-1,G,64,mC,T)
            if len(KA)==0 or len(KC)==0: continue
            ia,ic,_=bloom_match(KA,KC,int(min(30,max(20,np.ceil(np.log2(len(KC)))+2))))
            ncand+=len(ia)
            for q,qc in zip(ia,ic):
                idx=(int(IA[q]),int(JA[q]),int(IC[qc]),int(JC[qc]))
                if verify(idx): found.append(idx)
                if len(found)>=max_solutions: break
            if len(found)>=max_solutions: break
        if nb_%20==0 or nb_<=3: log(json.dumps(dict(stage='kt_progress',buckets=nb_,AB=int(hAB.sum()),CD=int(hCD.sum()),candidates=ncand,solutions=len(found),secs=round(time.time()-t0,1))))
    log(json.dumps(dict(stage='kt_done',buckets=nb_,candidates=ncand,solutions=len(found),secs=round(time.time()-t0,1))))
    return found
