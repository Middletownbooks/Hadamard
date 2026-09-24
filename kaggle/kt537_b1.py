# =====================================================================================
# !!! DO NOT RUN AS CONFIGURED: v=537, d=22, r=35 was later PROVEN INFEASIBLE (no solutions exist):
# !!! the Z8 compression cannot satisfy sum_k P_k(w) = 4m+4 at w = 134, 268 (logs/b1_compression_537.log).
# !!! The code is kept as a generic b=1 join; screen any new configuration with b1red.compression_feasible first.
# b=1 k-tree join for v=537 (GS core length m=536, row sums (2,0,0,0), sum_k PAF_k(s) = -4, s != 0),
# reduced by an order-22 multiplier subgroup U of (Z/536)* with r=35 orbits (L=18 lag classes).
# Self-contained (numpy + numba).  Kaggle-ready.
#
# POOLS ARE COMPLETE: every U-invariant +-1 block with the right row sum and P(w) <= 4m+4 at all w
#   is enumerated (meet-in-the-middle, ~1 min; 2^24.19 blocks per row sum).  Every solution in this
#   configuration is a 4-tuple from these pools, so buckets partition the whole search space.
# Model estimate [E]: ~2^17 solution 4-tuples in the space (partial-join validated through 8 classes).
# Validated in the source repo: recovers the real v=419 core when planted (same code path, d=15 r=42).
#
# CHUNKING: several instances with the same CONFIG and INSTANCE=0..NINST-1 split the buckets
#   (bucket b -> instance b % NINST).  A different seed = a different choice of chunk classes
#   (an independent partition of the same space).  The source repo runs seed=0 locally, so use seed>=1 here.
# RESUME: add the previous output as an input dataset and set RESUME_FROM to its kt537 directory.
# OUTPUT: OUT/solutions_v537.txt (JSON lines, each verified) and OUT/log_inst{INSTANCE}.jsonl
# =====================================================================================
import os, sys, json, time, glob, heapq, shutil
from math import gcd, log as log_, log2 as log2_, pi
import numpy as np, numba as nb

CONFIG = dict(
    v=537, d=22, r=35, which=1,   # which: index of the (d,r) subgroup (1 = the configuration costed at log2E~17)
    k=6,                          # chunk classes in the bucket key
    mem_budget_bytes=10.0e9,      # stage-1 pass budget (24 B per stored pair); fine on 30 GB
    logm=30,                      # bloom filter size 2^logm bytes
    INSTANCE=0, NINST=1,
    max_hours=11.4,
    nbuckets=2000000,
    OUT='/kaggle/working/kt537', RESUME_FROM=None,
    seed=1,
)
if not os.path.exists('/kaggle') and os.environ.get('KT_OUT'): CONFIG['OUT']=os.environ['KT_OUT']
for key in list(CONFIG):
    ev=os.environ.get('KT_'+key.upper())
    if ev is not None: CONFIG[key]=type(CONFIG[key])(ev) if CONFIG[key] is not None else ev
C=CONFIG; T_START=time.time()

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

def units(m): return [a for a in range(1,m) if gcd(a,m)==1]

def closure(gens,m):
    S={1}; frontier=[1]
    while frontier:
        nf=[]
        for x in frontier:
            for g in gens:
                y=x*g%m
                if y not in S: S.add(y); nf.append(y)
        frontier=nf
    return frozenset(S)

def subgroups(m):
    us=units(m); cyc={}
    for a in us: cyc.setdefault(closure([a],m),a)
    subs=set(cyc)
    cl=list(cyc)
    for i in range(len(cl)):
        for j in range(i+1,len(cl)):
            if not cl[i]<=cl[j] and not cl[j]<=cl[i]:
                subs.add(closure(list(cl[i])[:1]+[cyc[cl[i]],cyc[cl[j]]],m))
    return subs

def orbits(U,m):
    lab=-np.ones(m,int); orbs=[]
    for x in range(m):
        if lab[x]<0:
            o=sorted({x*u%m for u in U})
            for y in o: lab[y]=len(orbs)
            orbs.append(o)
    return lab,orbs

def lagclasses(U,m):
    seen=set(); reps=[]; sizes=[]
    for s in range(1,m):
        if s in seen: continue
        c={s*u%m for u in U}|{(-s*u)%m for u in U}; seen|=c; reps.append(s); sizes.append(len(c))
    return np.array(reps),np.array(sizes)

def count_rowsum(sz,a):
    from collections import defaultdict
    D={0:1}
    for s_ in sz:
        E=defaultdict(int)
        for x,c in D.items(): E[x+s_]+=c; E[x-s_]+=c
        D=E
    return D.get(a,0)

def find_U(m,d,r,which=0):
    c=[sorted(U) for U in subgroups(m) if len(U)==d and len(orbits(U,m)[1])==r]
    c.sort(); return c[which], len(c)

class RedC:
    def __init__(s,m,U):
        s.m=m; s.U=sorted(U); s.d=len(U)
        lab,orbs=orbits(U,m); s.orb=lab; s.orbits=orbs; s.r=len(orbs)
        s.size=np.array([len(o) for o in orbs])
        s.lags,s.lagsize=lagclasses(U,m); s.L=len(s.lags)
        # frequency classes: same as lag classes, plus omega=0
        s.freqs=np.concatenate([[0],s.lags])
        x=np.arange(m); E=np.exp(2j*np.pi*np.outer(s.freqs,x)/m)          # (F,m)
        s.S=np.stack([E[:,o].sum(1) for o in orbs],1)                       # (F,r)
    def expand(s,Y): return Y[...,s.orb]
    def paf(s,X):
        F=np.fft.fft(X.astype(float),axis=-1); R=np.fft.ifft(np.abs(F)**2,axis=-1).real
        return np.rint(R).astype(np.int64)
    def count_rowsum(s,a): return count_rowsum(list(s.size),a)

def _tables(R,lo,hi):
    nb_=hi-lo; n=1<<nb_
    bits=((np.arange(n)[:,None]>>np.arange(nb_)[None,:])&1)*2-1                  # (n, nb_) signs
    rs=bits@R.size[lo:hi]
    F=bits@R.S[1:,lo:hi].T                                                       # skip omega=0
    return rs.astype(np.int64), np.ascontiguousarray(F.real), np.ascontiguousarray(F.imag)

@nb.njit(parallel=True,cache=True)
def _scan(rsH,FrH,FiH,rsLs,oL,FrL,FiL,starts,ends,a,bound,shift,write,cnt,out,off):
    NH=rsH.shape[0]; F=FrH.shape[1]
    for h in nb.prange(NH):
        t=a-rsH[h]
        # binary search for t in sorted rsLs
        lo=np.searchsorted(rsLs,t); hi=np.searchsorted(rsLs,t+1)
        c=0; n=off[h] if write else 0
        for p in range(lo,hi):
            l=oL[p]; ok=True
            for f in range(F):
                xr=FrH[h,f]+FrL[l,f]; xi=FiH[h,f]+FiL[l,f]
                if xr*xr+xi*xi>bound: ok=False; break
            if ok:
                if write: out[n]=(np.uint64(h)<<np.uint64(shift))|np.uint64(l); n+=1
                else: c+=1
        if not write: cnt[h]=c
    return 0

def enum_pool(R,a,bound,hkeep=None):
    """hkeep: optional sorted array of high-half indices to scan (random subsample of the full pool)."""
    h=R.r//2
    rsL,FrL,FiL=_tables(R,0,h); rsH,FrH,FiH=_tables(R,h,R.r)
    if hkeep is not None: rsH,FrH,FiH=rsH[hkeep],FrH[hkeep],FiH[hkeep]
    oL=np.argsort(rsL,kind='stable'); rsLs=rsL[oL]
    cnt=np.zeros(len(rsH),np.int64); dummy=np.zeros(1,np.uint64); off=np.zeros(len(rsH),np.int64)
    _scan(rsH,FrH,FiH,rsLs,oL,FrL,FiL,0,0,a,float(bound),h,False,cnt,dummy,off)
    off=np.concatenate([[0],np.cumsum(cnt)[:-1]]); out=np.empty(int(cnt.sum()),np.uint64)
    _scan(rsH,FrH,FiH,rsLs,oL,FrL,FiL,0,0,a,float(bound),h,True,cnt,out,off)
    if hkeep is not None:
        hi=np.asarray(hkeep,np.uint64)[(out>>np.uint64(h)).astype(np.int64)]
        out=(hi<<np.uint64(h))|(out&np.uint64((1<<h)-1))
    out.sort(); return out

def bits_to_Y(bits,r):
    return (((bits[:,None]>>np.arange(r,dtype=np.uint64)[None,:])&np.uint64(1)).astype(np.int64))*2-1

def sparse_N(R):
    rows=[];cols=[];cls=[];cnt=[]
    orb=R.orb; m=R.m
    for c,s in enumerate(R.lags):
        i=np.arange(m); j=(i+s)%m
        Z=np.zeros((R.r,R.r),np.int64); np.add.at(Z,(orb[i],orb[j]),1)
        q,p=np.nonzero(Z)
        for a,b in zip(q,p): cls.append(c); rows.append(a); cols.append(b); cnt.append(Z[a,b])
    return np.array(cls),np.array(rows),np.array(cols),np.array(cnt)

@nb.njit(parallel=True,cache=True)
def paf_from_bits(bits,cls,rows,cols,cnt,L,r):
    N=bits.shape[0]; out=np.zeros((N,L),np.int16)
    for n in nb.prange(N):
        w=bits[n]
        for e in range(cls.shape[0]):
            yq=1 if (w>>np.uint64(rows[e]))&np.uint64(1) else -1
            yp=1 if (w>>np.uint64(cols[e]))&np.uint64(1) else -1
            out[n,cls[e]]+=cnt[e]*yq*yp
    return out

@nb.njit(parallel=True,cache=True)
def merge_thist(sX,oX,cX,sY,oY,cY,nchunk,T):
    """histogram over t=cX+cY (offset T) of all merge matches"""
    NX=sX.shape[0]; NY=sY.shape[0]
    b0=np.linspace(0,NX,nchunk+1).astype(np.int64)
    for b in range(1,nchunk):
        x=max(b0[b],b0[b-1])
        while x<NX and x>0 and sX[x]==sX[x-1]: x+=1
        b0[b]=x
    H=np.zeros((nchunk,2*T+1),np.int64)
    for b in nb.prange(nchunk):
        lo=b0[b]; hi=b0[b+1]
        if lo>=hi: continue
        j=np.searchsorted(sY,sX[lo]); i=lo
        while i<hi and j<NY:
            a=sX[i]; bb=sY[j]
            if a<bb: i+=1
            elif a>bb: j+=1
            else:
                j2=j
                while j2<NY and sY[j2]==a: j2+=1
                i2=i
                while i2<hi and sX[i2]==a: i2+=1
                for ii in range(i,i2):
                    cx=cX[oX[ii]]
                    for jj in range(j,j2):
                        H[b,cx+cY[oY[jj]]+T]+=1
                i=i2; j=j2
    return H.sum(0)

@nb.njit(parallel=True,cache=True)
def stage1_merge_t(sX,oX,cX,hX,sY,oY,cY,hY,sign,G,nchunk,mask,T):
    """as stage1_merge but only emits pairs with tlo<=t<thi"""
    NX=sX.shape[0]; NY=sY.shape[0]
    b0=np.linspace(0,NX,nchunk+1).astype(np.int64)
    for b in range(1,nchunk):
        x=max(b0[b],b0[b-1])
        while x<NX and x>0 and sX[x]==sX[x-1]: x+=1
        b0[b]=x
    cnt=np.zeros(nchunk,np.int64)
    for pas in range(2):
        if pas==1:
            off=np.zeros(nchunk+1,np.int64)
            for b in range(nchunk): off[b+1]=off[b]+cnt[b]
            K=np.empty(off[nchunk],np.uint64); I=np.empty(off[nchunk],np.int64); J=np.empty(off[nchunk],np.int64)
        else:
            off=np.zeros(nchunk+1,np.int64); K=np.empty(0,np.uint64); I=np.empty(0,np.int64); J=np.empty(0,np.int64)
        for b in nb.prange(nchunk):
            lo=b0[b]; hi=b0[b+1]
            if lo>=hi: continue
            n=off[b]; c=0
            j=np.searchsorted(sY,sX[lo]); i=lo
            while i<hi and j<NY:
                a=sX[i]; bb=sY[j]
                if a<bb: i+=1
                elif a>bb: j+=1
                else:
                    j2=j
                    while j2<NY and sY[j2]==a: j2+=1
                    i2=i
                    while i2<hi and sX[i2]==a: i2+=1
                    for ii in range(i,i2):
                        xi=oX[ii]; cx=cX[xi]
                        for jj in range(j,j2):
                            yj=oY[jj]; t=cx+cY[yj]
                            if mask[t+T]:
                                if pas==0: c+=1
                                else:
                                    kk=hX[xi]+hY[yj]+np.uint64(np.int64(t))*G
                                    if sign<0: kk=np.uint64(0)-kk
                                    K[n]=kk; I[n]=xi; J[n]=yj; n+=1
                    i=i2; j=j2
            if pas==0: cnt[b]=c
    return K,I,J

@nb.njit(parallel=True,cache=True)
def bloom_match(KA,KC,logm):
    """return (ia, ic) index pairs with KA[ia]==KC[ic]; parallel byte-bloom on KC, filter KA, exact scan."""
    m=np.int64(1)<<logm; mask=np.uint64(m-1)
    bl=np.zeros(m,np.uint8)
    NC=KC.shape[0]; NA=KA.shape[0]
    for i in nb.prange(NC):
        bl[np.int64((KC[i]*np.uint64(0xD6E8FEB86659FD93))>>np.uint64(64-logm))]=1
    nch=256
    ba=np.linspace(0,NA,nch+1).astype(np.int64); cnt=np.zeros(nch,np.int64)
    for b in nb.prange(nch):
        c=0
        for i in range(ba[b],ba[b+1]):
            if bl[np.int64((KA[i]*np.uint64(0xD6E8FEB86659FD93))>>np.uint64(64-logm))]: c+=1
        cnt[b]=c
    off=np.zeros(nch+1,np.int64)
    for b in range(nch): off[b+1]=off[b]+cnt[b]
    surv=np.empty(off[nch],np.int64)
    for b in nb.prange(nch):
        n=off[b]
        for i in range(ba[b],ba[b+1]):
            if bl[np.int64((KA[i]*np.uint64(0xD6E8FEB86659FD93))>>np.uint64(64-logm))]: surv[n]=i; n+=1
    sk=KA[surv]; o=np.argsort(sk); sk=sk[o]; surv=surv[o]
    bc=np.linspace(0,NC,nch+1).astype(np.int64); cnt2=np.zeros(nch,np.int64)
    NS=sk.shape[0]
    for b in nb.prange(nch):
        c=0
        for i in range(bc[b],bc[b+1]):
            p=np.searchsorted(sk,KC[i])
            while p<NS and sk[p]==KC[i]: c+=1; p+=1
        cnt2[b]=c
    off2=np.zeros(nch+1,np.int64)
    for b in range(nch): off2[b+1]=off2[b]+cnt2[b]
    ia=np.empty(off2[nch],np.int64); ic=np.empty(off2[nch],np.int64)
    for b in nb.prange(nch):
        n=off2[b]
        for i in range(bc[b],bc[b+1]):
            p=np.searchsorted(sk,KC[i])
            while p<NS and sk[p]==KC[i]:
                ia[n]=surv[p]; ic[n]=i; n+=1; p+=1
    return ia,ic,NS


# ---------------------------------------------------------------- main
def main():
    v,d,r,which,k=C['v'],C['d'],C['r'],C['which'],C['k']; m=v-1
    OUT=C['OUT']; os.makedirs(OUT,exist_ok=True)
    if C['RESUME_FROM']:
        for f in glob.glob(os.path.join(C['RESUME_FROM'],'*.jsonl'))+glob.glob(os.path.join(C['RESUME_FROM'],'solutions_*.txt')):
            dst=os.path.join(OUT,os.path.basename(f))
            if not os.path.exists(dst): shutil.copy(f,dst)
    logf=open(os.path.join(OUT,f'log_inst{C["INSTANCE"]}.jsonl'),'a')
    def log(**kw):
        kw['secs']=round(time.time()-T_START,1); s=json.dumps(kw); print(s,flush=True); logf.write(s+'\n'); logf.flush()
    U,nU=find_U(m,d,r,which); R=RedC(m,U); bound=4*m+4
    log(stage='config',config={a:(list(b) if isinstance(b,tuple) else b) for a,b in C.items()},L=R.L,U=U[:8],n_subgroups=nU)
    pools={a:enum_pool(R,a,bound) for a in (2,0)}
    log(stage='pools',sizes={a:len(x) for a,x in pools.items()},log2={a:round(float(np.log2(len(x))),2) for a,x in pools.items()})
    Nsp=sparse_N(R)
    rng=np.random.default_rng(C['seed'])
    cls=rng.permutation(R.L); chunk=cls[:k]; cc=cls[k]; rest=cls[k+1:]
    W=rng.integers(1,2**63,size=len(rest),dtype=np.uint64)|np.uint64(1); G=np.uint64(0x9E3779B97F4A7C15)
    Dd={}
    for a,x in pools.items():
        P=paf_from_bits(x,*Nsp,R.L,R.r)
        if a==2: P+=4                                   # target -4 -> 0
        h=np.empty(len(x),np.uint64)
        for s_ in range(0,len(x),1<<20): h[s_:s_+(1<<20)]=(P[s_:s_+(1<<20)][:,rest].astype(np.int64).astype(np.uint64)*W[None,:]).sum(1)
        Dd[a]=dict(ch=np.ascontiguousarray(P[:,chunk]),c=np.ascontiguousarray(P[:,cc]),h=h,
                   P=P[np.random.default_rng(1).integers(0,len(x),300000)].astype(np.int64))
        del P
    def nbits(idx):
        mu=Dd[2]['P'][:,idx].mean(0)+3*Dd[0]['P'][:,idx].mean(0)
        Cv=np.atleast_2d(np.cov(Dd[2]['P'][:,idx].T.astype(float)))+3*np.atleast_2d(np.cov(Dd[0]['P'][:,idx].T.astype(float)))
        w,V=np.linalg.eigh(Cv); keep=w>1e-6*w.max(); pr=V.T@mu
        lp=keep.sum()*log_(4)-0.5*np.sum(np.log(2*pi*w[keep]))-0.5*np.sum(pr[keep]**2/w[keep])
        return -lp/log_(2)
    n_all=nbits(np.arange(R.L)); n_rest=n_all-nbits(chunk)
    log(stage='model',n_all=round(n_all,1),n_rest=round(n_rest,1),log2E=round(float(np.log2(len(pools[2]))+3*np.log2(len(pools[0]))-n_all),1),chunk=chunk.tolist(),cc=int(cc))
    for a in Dd: del Dd[a]['P']
    c0=m%4; z=np.zeros(k,np.int64)
    sO=keyof(Dd[0]['ch'],z,1,c0); oO=np.argsort(sO,kind='stable'); sO=sO[oO]
    ns=2_000_000
    def pairdist(a1,a2):
        i=rng.integers(0,len(pools[a1]),ns); j=rng.integers(0,len(pools[a2]),ns)
        S=Dd[a1]['ch'][i].astype(np.int32)+Dd[a2]['ch'][j].astype(np.int32)
        return [dict(zip(*np.unique(S[:,q],return_counts=True))) for q in range(k)]
    dAB=pairdist(2,0); dCD=pairdist(0,0); opts=[]
    for q in range(k):
        vals=sorted(set(dAB[q])&{-x for x in dCD[q]},key=lambda t:-(dAB[q][t]/ns)*(dCD[q].get(-t,0)/ns))
        opts.append([(int(t),log2_(dAB[q][t]/ns)+log2_(dCD[q][-t]/ns)) for t in vals[:24]])
    score=lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
    best=tuple(0 for _ in range(k)); heap=[(-score(best),best)]; seen={best}; buckets=[]
    while heap and len(buckets)<C['nbuckets']:
        s,ix=heapq.heappop(heap); buckets.append(ix)
        for q in range(k):
            if ix[q]+1<len(opts[q]):
                nx=ix[:q]+(ix[q]+1,)+ix[q+1:]
                if nx not in seen: seen.add(nx); heapq.heappush(heap,(-score(nx),nx))
    done=set()
    for f in glob.glob(os.path.join(OUT,'log_inst*.jsonl')):
        for l in open(f):
            try: j=json.loads(l)
            except Exception: continue
            if j.get('stage')=='bucket' and j.get('seed',C['seed'])==C['seed']: done.add(j['b'])
    log(stage='resume',done=len(done))
    solf=os.path.join(OUT,f'solutions_v{v}.txt'); T=2*m+16; A=Dd[2]; O=Dd[0]
    nsol=0; exp_tot=0.0
    for b,ix in enumerate(buckets):
        if b%C['NINST']!=C['INSTANCE'] or b in done: continue
        if time.time()-T_START>C['max_hours']*3600: log(stage='time_limit'); break
        tau=np.array([opts[q][ix[q]][0] for q in range(k)],np.int64)
        nA=keyof(A['ch'],tau,-1,c0); oA=np.argsort(nA); nA=nA[oA]
        nC=keyof(O['ch'],-tau,-1,c0); oC=np.argsort(nC); nC=nC[oC]
        hAB=merge_thist(nA,oA,A['c'],sO,oO,O['c'],64,T); hCD=merge_thist(nC,oC,O['c'],sO,oO,O['c'],64,T)
        both=hAB*hCD[::-1]; cAB=int(hAB.sum()); cCD=int(hCD.sum()); expb=float(cAB)*float(cCD)*2.0**(-n_rest)
        budget=C['mem_budget_bytes']/24; passes=[]; cur=[]; load=0
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
            ia,ic,_=bloom_match(KA,KC,C['logm']); ncand+=len(ia)
            for q,qc in zip(ia,ic):
                bits=np.array([pools[2][IA[q]],pools[0][JA[q]],pools[0][IC[qc]],pools[0][JC[qc]]],np.uint64)
                X=R.expand(bits_to_Y(bits,R.r)); tot=R.paf(X).sum(0)
                if np.all(tot[1:]==-4):
                    found+=1; nsol+=1
                    with open(solf,'a') as fo:
                        fo.write(json.dumps(dict(v=v,d=d,r=r,U=U,seed=C['seed'],bucket=b,rowsums=X.sum(1).tolist(),
                            blocks=[''.join('+' if t_>0 else '-' for t_ in x) for x in X]))+'\n')
                    log(stage='SOLUTION',b=b)
            del KA,IA,JA,KC,IC,JC
        exp_tot+=expb
        log(stage='bucket',b=b,seed=C['seed'],tau=tau.tolist(),AB=cAB,CD=cCD,passes=len(passes),cand=int(ncand),
            solutions=found,total=nsol,exp_hits=round(expb,4),exp_cum=round(exp_tot,3))
    log(stage='done',solutions=nsol)
    if os.path.exists(solf): verify_file(solf)

def verify_file(path):
    """independent check: GS b=1 condition sum_k PAF_k(s) = -4 (s != 0), row sums (+-2,0,0,0)"""
    ok=0; bad=0
    for l in open(path):
        j=json.loads(l); X=np.array([[1 if c=='+' else -1 for c in b] for b in j['blocks']])
        F=np.fft.fft(X,axis=1); P=np.rint(np.fft.ifft(np.abs(F)**2,axis=1).real).sum(0)
        good=bool(np.all(P[1:]==-4) and sorted(np.abs(X.sum(1)).tolist())==[0,0,0,2])
        ok+=good; bad+=not good
    print(f'verify {path}: {ok} valid, {bad} invalid',flush=True)

if __name__=='__main__': main()
