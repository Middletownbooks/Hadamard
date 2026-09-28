# =====================================================================================
# v=523 profile A (35,17,17,17) -- version 2, paste the WHOLE file into ONE notebook cell and run.
#
# The ONLY line you normally edit is  INSTANCE = ...  just below the imports.  Use a different INSTANCE in
# each notebook (0, 1, 2, ...).  Every INSTANCE works through ALL buckets on its own independently seeded pools,
# so any INSTANCE value is valid (the old bug where INSTANCE >= 1 processed zero buckets cannot happen: see the
# assignment check printed at start, which aborts BEFORE the hours of generation if an instance would get nothing).
#
# Output handling (the problem with the earlier notebooks): by default the multi-GB pool/prep caches now live in
# /tmp and are NOT saved, so /kaggle/working only holds a small log and the solutions file, and the saved output
# is tiny.  Every solution is ALSO printed in full to the cell output (lines starting 'SOLUTION ') and verified
# (H H^T = 4v I) at the end of the run, so nothing depends on downloading output files.
# Set SAVE_CACHE=True only if you want to resume a chain across sessions (then output is ~10 GB, see README).
#
# Flat-pool k-tree (Wagner) join for b=0 Goethals-Seidel cores, reduced by the order-d
# multiplier subgroup U of (Z/v)*.  Self-contained (numpy + numba).  Kaggle-ready.
#
# Validated in the source repo on v=487, d=9, profile (37,17,17,1): recovers the planted
# Alpoge core, and found 2 new inequivalent solutions in ~40 buckets (2.8 min/bucket, 4 cores).
#
# CHUNKING: run several instances with the same CONFIG except INSTANCE=0..NINST-1.
#   - each instance generates its OWN pools (seeded by INSTANCE) -> independent coverage
#   - bucket b is processed by instance (b % NINST); buckets are ordered most-probable first
# RESUME: everything is written under OUT (default /kaggle/working/kt).  To continue in a new
#   session, add the previous session's output as an input dataset and set RESUME_FROM to its
#   kt directory: pools/prep are reused and already-finished buckets are skipped.
# OUTPUT: OUT/solutions_v{v}.txt  (JSON lines, each verified: sum PAF = 0 and H H^T = 4v I)
#         OUT/log_inst{INSTANCE}.jsonl (per-bucket stage-1 sizes, candidates, expected hits)
# =====================================================================================
import os, sys, json, time, glob, heapq, shutil
import numpy as np, numba as nb
from multiprocessing import Pool

INSTANCE = 0          # <<< EDIT THIS: 0, 1, 2, ... (a different value in each notebook)
SAVE_CACHE = False    # False: caches in /tmp (small saved output).  True: caches in /kaggle/working (resumable, ~10 GB)

CONFIG = dict(
    v=523, d=9, profile=(35, 17, 17, 17),   # b=0 profile (row sums, positive); A,B,C,D slots
    rho=0.39,                 # flatness threshold: PAF-energy ratio sqrt(F/(m L)) <= rho
    k=9,                      # chunk classes (first 8 go in the 64-bit key)
    Lmax=2**27,               # max pool size per row sum (raise with more RAM: ~27 B/block + join buffers)
    gen_minutes=90.0,         # pool generation wall time (all cores)
    walk=300000,              # tabu steps per walk
    mem_budget_bytes=10.0e9,   # stage-1 pass budget (24 B per stored pair); ~10e9 is fine on 30 GB
    INSTANCE=INSTANCE, NINST=1,   # NINST=1: each INSTANCE = independent run over all buckets (recommended)
    max_hours=11.4,           # stop before Kaggle's 12 h limit
    nbuckets=100000,
    OUT='/kaggle/working/kt_A_35-17-17-17', RESUME_FROM=None,
    GEN_DIR='/tmp/kt_gen',    # raw generator chunks (not persisted; deleted after pools are built to save /kaggle/working space)
    seed=0,
    SAVE_CACHE=SAVE_CACHE, CACHE_ROOT='/tmp/kt_cache',
)
if os.path.exists('/kaggle') is False and os.environ.get('KT_OUT'): CONFIG['OUT']=os.environ['KT_OUT']
for key in list(CONFIG):                      # allow env overrides, e.g. KT_INSTANCE=3
    ev=os.environ.get('KT_'+key.upper())
    if ev is not None:
        CONFIG[key]=type(CONFIG[key])(ev) if CONFIG[key] is not None and not isinstance(CONFIG[key],tuple) else (tuple(int(x) for x in ev.split(',')) if isinstance(CONFIG[key],tuple) else ev)
C=CONFIG; T_START=time.time()

# ---------------------------------------------------------------- reduced setting
def primroot(p):
    fs=[q for q in range(2,p) if (p-1)%q==0 and all(q%r for r in range(2,int(q**.5)+1))]
    for g in range(2,p):
        if all(pow(g,(p-1)//q,p)!=1 for q in fs): return g
class Red:
    def __init__(s,m,d):
        assert (m-1)%d==0
        s.m=m; s.d=d; g=primroot(m); h=pow(g,(m-1)//d,m); s.U=sorted({pow(h,k,m) for k in range(d)})
        s.orb=np.zeros(m,int); s.orbits=[[0]]
        for j in range((m-1)//d):
            o=sorted({pow(g,j,m)*u%m for u in s.U})
            for i in o: s.orb[i]=len(s.orbits)
            s.orbits.append(o)
        s.r=len(s.orbits); seen=set(); s.lags=[]
        for j in range((m-1)//d):
            x=pow(g,j,m)
            if x in seen: continue
            seen|={x*u%m for u in s.U}|{(-x*u)%m for u in s.U}; s.lags.append(x)
        s.lags=np.array(s.lags); s.L=len(s.lags)
    def expand(s,Y): return Y[...,s.orb]
    def paf(s,X):
        F=np.fft.fft(X.astype(float),axis=-1); return np.rint(np.fft.ifft(np.abs(F)**2,axis=-1).real).astype(np.int64)
    def psd(s,X): return np.abs(np.fft.fft(X.astype(float),axis=-1))**2

def build_M(R):
    M=np.zeros((R.L,R.r,R.r),np.int32)
    for c,s in enumerate(R.lags):
        i=np.arange(R.m); np.add.at(M[c],(R.orb[i],R.orb[(i+s)%R.m]),1)
    return M+M.transpose(0,2,1)
def sparse_N(R):
    cls=[];rows=[];cols=[];cnt=[]
    for c,s in enumerate(R.lags):
        i=np.arange(R.m); Z=np.zeros((R.r,R.r),np.int64); np.add.at(Z,(R.orb[i],R.orb[(i+s)%R.m]),1)
        q,p=np.nonzero(Z)
        for a,b in zip(q,p): cls.append(c); rows.append(a); cols.append(b); cnt.append(Z[a,b])
    return np.array(cls),np.array(rows),np.array(cols),np.array(cnt)

# ---------------------------------------------------------------- flat-pool generation
@nb.njit(cache=True)
def _init(y,M):
    r=y.shape[0]; L=M.shape[0]; G=np.zeros((L,r),np.int64); P=np.zeros(L,np.int64)
    for c in range(L):
        for q in range(r):
            s=0
            for p in range(r): s+=M[c,q,p]*y[p]
            G[c,q]=s
        t=0
        for q in range(r): t+=y[q]*G[c,q]
        P[c]=t//2
    return G,P
@nb.njit(cache=True)
def walk(y,M,F0,nsteps,tenure,seed,d,targets,out,nout):
    np.random.seed(seed); r=y.shape[0]; L=M.shape[0]; G,P=_init(y,M)
    F=0
    for c in range(L): F+=P[c]*P[c]
    rs=y[0]
    for q in range(1,r): rs+=d*y[q]
    tabu=np.zeros(r,np.int64); n=nout
    for step in range(nsteps):
        best=1<<62; bq=-1; nb_=0
        for q in range(r):
            if tabu[q]>step: continue
            yq=y[q]; dF=0
            for c in range(L):
                dd=-2*yq*G[c,q]+2*M[c,q,q]; dF+=2*P[c]*dd+dd*dd
            if dF<best or (dF==best and np.random.random()<1.0/(nb_+1)):
                if dF<best: nb_=1
                else: nb_+=1
                best=dF; bq=q
        yq=y[bq]
        for c in range(L):
            dd=-2*yq*G[c,bq]+2*M[c,bq,bq]; F+=2*P[c]*dd+dd*dd; P[c]+=dd
        for c in range(L):
            for p in range(r): G[c,p]-=2*yq*M[c,p,bq]
        y[bq]=-yq; tabu[bq]=step+1+tenure+np.random.randint(tenure//2+1)
        rs+=-2*yq*(1 if bq==0 else d)
        if F<=F0:
            a=abs(rs); ok=False
            for t in targets:
                if a==t: ok=True
            if ok:
                sgn=1 if rs>0 else -1; w=np.uint64(0)
                for q in range(r):
                    if y[q]*sgn>0: w|=np.uint64(1)<<np.uint64(q)
                out[n]=w; n+=1
                if n>=out.shape[0]: return n
    return n
def gen_worker(args):
    wid,secs,outdir=args
    R=Red(C['v'],C['d']); M=build_M(R); F0=int(C['rho']**2*R.m*R.L)
    targets=np.array(sorted(set(C['profile'])),np.int64); rng=np.random.default_rng(1000003*C['INSTANCE']+wid+17*C['seed'])
    buf=np.zeros(1<<24,np.uint64); n=0; ci=0; t0=time.time(); it=0
    def flush(b,ci):
        b=np.unique(b); ones=np.zeros(len(b),np.int64)
        for q in range(1,R.r): ones+=((b>>np.uint64(q))&np.uint64(1)).astype(np.int64)
        rs=np.where(b&np.uint64(1),1,-1)+R.d*(2*ones-(R.r-1))
        for a in targets: np.save(f'{outdir}/gen_w{wid}_rs{a}_c{ci:04d}.npy',b[rs==a])
    while time.time()-t0<secs:
        y=rng.choice([-1,1],size=R.r).astype(np.int64)
        n=walk(y,M,F0,C['walk'],max(3,R.r//10),int(rng.integers(1<<31)),R.d,targets,buf,n); it+=1
        if n>len(buf)-(1<<21): flush(buf[:n],ci); ci+=1; n=0
    flush(buf[:n],ci); return it

# ---------------------------------------------------------------- join kernels
@nb.njit(parallel=True,cache=True)
def paf_from_bits(bits,cls,rows,cols,cnt,L):
    N=bits.shape[0]; out=np.zeros((N,L),np.int16)
    for n in nb.prange(N):
        w=bits[n]
        for e in range(cls.shape[0]):
            yq=1 if (w>>np.uint64(rows[e]))&np.uint64(1) else -1
            yp=1 if (w>>np.uint64(cols[e]))&np.uint64(1) else -1
            out[n,cls[e]]+=cnt[e]*yq*yp
    return out
@nb.njit(parallel=True,cache=True)
def keyof_nb(ch,tau,sgn):
    N=ch.shape[0]; K=ch.shape[1]; out=np.empty(N,np.uint64)
    for n in nb.prange(N):
        w=np.uint64(0)
        for i in range(K):
            val=np.int64(tau[i])+sgn*np.int64(ch[n,i]); w|=np.uint64(val&255)<<np.uint64(8*i)
        out[n]=w
    return out
@nb.njit(cache=True)
def _bounds(sX,nchunk):
    NX=sX.shape[0]; b0=np.linspace(0,NX,nchunk+1).astype(np.int64)
    for b in range(1,nchunk):
        x=max(b0[b],b0[b-1])
        while x<NX and x>0 and sX[x]==sX[x-1]: x+=1
        b0[b]=x
    return b0
@nb.njit(parallel=True,cache=True)
def merge_thist(sX,oX,cX,sY,oY,cY,nchunk,T):
    NY=sY.shape[0]; b0=_bounds(sX,nchunk); H=np.zeros((nchunk,2*T+1),np.int64)
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
                    for jj in range(j,j2): H[b,cx+cY[oY[jj]]+T]+=1
                i=i2; j=j2
    return H.sum(0)
@nb.njit(parallel=True,cache=True)
def stage1(sX,oX,cX,hX,sY,oY,cY,hY,sign,G,nchunk,mask,T):
    NY=sY.shape[0]; b0=_bounds(sX,nchunk); cnt=np.zeros(nchunk,np.int64)
    off=np.zeros(nchunk+1,np.int64); K=np.empty(0,np.uint64); I=np.empty(0,np.int64); J=np.empty(0,np.int64)
    for pas in range(2):
        if pas==1:
            for b in range(nchunk): off[b+1]=off[b]+cnt[b]
            K=np.empty(off[nchunk],np.uint64); I=np.empty(off[nchunk],np.int64); J=np.empty(off[nchunk],np.int64)
        for b in nb.prange(nchunk):
            lo=b0[b]; hi=b0[b+1]
            if lo>=hi: continue
            n=off[b]; c=0; j=np.searchsorted(sY,sX[lo]); i=lo
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
    m=np.int64(1)<<logm; bl=np.zeros(m,np.uint8); NC=KC.shape[0]; NA=KA.shape[0]; H=np.uint64(0xD6E8FEB86659FD93)
    for i in nb.prange(NC): bl[np.int64((KC[i]*H)>>np.uint64(64-logm))]=1
    nch=256; ba=np.linspace(0,NA,nch+1).astype(np.int64); cnt=np.zeros(nch,np.int64)
    for b in nb.prange(nch):
        c=0
        for i in range(ba[b],ba[b+1]):
            if bl[np.int64((KA[i]*H)>>np.uint64(64-logm))]: c+=1
        cnt[b]=c
    off=np.zeros(nch+1,np.int64)
    for b in range(nch): off[b+1]=off[b]+cnt[b]
    surv=np.empty(off[nch],np.int64)
    for b in nb.prange(nch):
        n=off[b]
        for i in range(ba[b],ba[b+1]):
            if bl[np.int64((KA[i]*H)>>np.uint64(64-logm))]: surv[n]=i; n+=1
    sk=KA[surv]; o=np.argsort(sk); sk=sk[o]; surv=surv[o]; NS=sk.shape[0]
    bc=np.linspace(0,NC,nch+1).astype(np.int64); cnt2=np.zeros(nch,np.int64)
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
            while p<NS and sk[p]==KC[i]: ia[n]=surv[p]; ic[n]=i; n+=1; p+=1
    return ia,ic

# ---------------------------------------------------------------- driver
def main():
    v,d,prof,k=C['v'],C['d'],tuple(C['profile']),C['k']
    OUT=C['OUT']; os.makedirs(OUT,exist_ok=True)
    if C['NINST']<1: raise SystemExit('NINST must be >= 1')
    if C['INSTANCE']>=C['NINST']:
        print(f"NOTE: INSTANCE={C['INSTANCE']} >= NINST={C['NINST']}: this instance takes buckets b % NINST == {C['INSTANCE']%C['NINST']} "
              f"(with NINST=1 that is ALL buckets, on this instance's own independent pools).",flush=True)
    # bucket-assignment check BEFORE any expensive work
    mine=sum(1 for b in range(min(C['nbuckets'],1000)) if b%C['NINST']==C['INSTANCE']%C['NINST'])
    print(f"ASSIGNMENT CHECK: INSTANCE={C['INSTANCE']} NINST={C['NINST']} -> this run takes {mine} of the first "
          f"{min(C['nbuckets'],1000)} buckets (pools seeded by INSTANCE)",flush=True)
    if mine==0: raise SystemExit('ASSIGNMENT CHECK FAILED: this instance would process no buckets; fix INSTANCE/NINST')
    gdir=os.path.join(C['GEN_DIR'],f'gen_inst{C["INSTANCE"]}')
    cache=os.path.join(OUT if C['SAVE_CACHE'] else C['CACHE_ROOT'],f'cache_inst{C["INSTANCE"]}')
    if C['RESUME_FROM']:
        src=os.path.join(C['RESUME_FROM'],f'cache_inst{C["INSTANCE"]}')
        if os.path.isdir(src) and glob.glob(os.path.join(src,'h_*.npy')):
            cache=src                     # read the finished pools/prep in place (read-only input dataset): no copy, no regeneration
        elif os.path.isdir(src) and not os.path.isdir(cache): shutil.copytree(src,cache)
        for f in glob.glob(os.path.join(C['RESUME_FROM'],'*.jsonl'))+glob.glob(os.path.join(C['RESUME_FROM'],'solutions_*.txt')):
            if not os.path.exists(os.path.join(OUT,os.path.basename(f))): shutil.copy(f,OUT)
    os.makedirs(gdir,exist_ok=True)
    if not cache.startswith('/kaggle/input'): os.makedirs(cache,exist_ok=True)
    logf=open(os.path.join(OUT,f'log_inst{C["INSTANCE"]}.jsonl'),'a')
    def log(**kw):
        kw['t']=round(time.time()-T_START,1); s=json.dumps(kw); print(s,flush=True); logf.write(s+'\n'); logf.flush()
    log(stage='config',config={a:(list(b) if isinstance(b,tuple) else b) for a,b in C.items()})
    R=Red(v,d); Nsp=sparse_N(R)
    rng=np.random.default_rng(C['seed']+7919*C['INSTANCE'])      # pool subsampling only
    rng_cls=np.random.default_rng(12345+C['seed'])                 # chunk classes + hash weights (fixed; must match cache)
    rng_pd=np.random.default_rng(777+C['seed'])                    # bucket-order estimation
    cp=lambda nm: os.path.join(cache,nm+'.npy')
    # 1. pools
    if not all(os.path.exists(cp(f'pool_{a}')) for a in set(prof)):
        if not glob.glob(os.path.join(gdir,'gen_w*_rs*_c*.npy')):
            ncpu=os.cpu_count() or 4; log(stage='gen_start',workers=ncpu,minutes=C['gen_minutes'])
            with Pool(ncpu) as p: p.map(gen_worker,[(w,C['gen_minutes']*60,gdir) for w in range(ncpu)])
        pools={}
        for a in set(prof):
            x=np.unique(np.concatenate([np.load(f) for f in glob.glob(os.path.join(gdir,f'gen_w*_rs{a}_c*.npy'))]+[np.zeros(0,np.uint64)]))
            if len(x)>C['Lmax']: x=x[rng.choice(len(x),C['Lmax'],replace=False)]
            np.save(cp(f'pool_{a}'),x); log(stage='pool',rowsum=a,size=len(x),log2=round(float(np.log2(max(len(x),1))),2)); del x
        shutil.rmtree(gdir,ignore_errors=True)
    pools={a:np.load(cp(f'pool_{a}'),mmap_mode='r').view(np.ndarray) for a in set(prof)}
    # 2. prep
    classes=rng_cls.permutation(R.L); chunk=classes[:k]; cc=classes[k]; rest=classes[k+1:]
    W=rng_cls.integers(1,2**63,size=len(rest),dtype=np.uint64)|np.uint64(1); G=np.uint64(0x9E3779B97F4A7C15)
    if not all(os.path.exists(cp(f'h_{a}')) for a in set(prof)):
        for a,x in pools.items():
            ch=np.empty((len(x),k),np.int8); c_=np.empty(len(x),np.int16); h=np.empty(len(x),np.uint64)
            for s in range(0,len(x),1<<21):
                P=paf_from_bits(np.ascontiguousarray(x[s:s+(1<<21)]),*Nsp,R.L)
                ch[s:s+len(P)]=P[:,chunk]; c_[s:s+len(P)]=P[:,cc]
                h[s:s+len(P)]=(P[:,rest].astype(np.int64).astype(np.uint64)*W[None,:]).sum(1)
            np.save(cp(f'ch_{a}'),ch); np.save(cp(f'c_{a}'),c_); np.save(cp(f'h_{a}'),h); del ch,c_,h
            log(stage='prep',rowsum=a)
    D={a:dict(ch=np.load(cp(f'ch_{a}'),mmap_mode='r').view(np.ndarray),c=np.load(cp(f'c_{a}')),h=np.load(cp(f'h_{a}'),mmap_mode='r').view(np.ndarray)) for a in set(prof)}
    KE=min(k,8)
    own=lambda ch: keyof_nb(np.ascontiguousarray(ch[:,:KE]),np.zeros(KE,np.int64),1)
    partner=lambda ch,tau: keyof_nb(np.ascontiguousarray(ch[:,:KE]),np.asarray(tau,np.int64)[:KE],-1)
    sB=own(D[prof[1]]['ch']); oB=np.argsort(sB,kind='stable'); sB=sB[oB]
    sD=own(D[prof[3]]['ch']); oD=np.argsort(sD,kind='stable'); sD=sD[oD]
    # 3. bucket order (best-first by empirical per-class pair-sum probabilities)
    ns=2_000_000
    def pairdist(a1,a2):
        i=rng_pd.integers(0,len(pools[a1]),ns); j=rng_pd.integers(0,len(pools[a2]),ns)
        S=D[a1]['ch'][i].astype(np.int16)+D[a2]['ch'][j].astype(np.int16)
        return [dict(zip(*np.unique(S[:,q],return_counts=True))) for q in range(k)]
    dAB=pairdist(prof[0],prof[1]); dCD=pairdist(prof[2],prof[3]); opts=[]
    for q in range(k):
        vals=sorted(set(dAB[q])&{-x for x in dCD[q]},key=lambda t:-(dAB[q][t]*dCD[q].get(-t,0)))
        opts.append([(t,np.log2(dAB[q][t]/ns)+np.log2(dCD[q][-t]/ns)) for t in vals[:12]])
    score=lambda ix: sum(opts[q][ix[q]][1] for q in range(k))
    start=tuple(0 for _ in range(k)); heap=[(-score(start),start)]; seen={start}; buckets=[]
    while heap and len(buckets)<C['nbuckets']:
        s,ix=heapq.heappop(heap); buckets.append(ix)
        for q in range(k):
            if ix[q]+1<len(opts[q]):
                nx=ix[:q]+(ix[q]+1,)+ix[q+1:]
                if nx not in seen: seen.add(nx); heapq.heappush(heap,(-score(nx),nx))
    done=set()
    for f in glob.glob(os.path.join(OUT,f'log_inst{C["INSTANCE"]}.jsonl')):
        for l in open(f):
            try:
                z=json.loads(l)
                if z.get('stage')=='bucket': done.add(z['b'])
            except Exception: pass
    log(stage='ready',buckets=len(buckets),done_before=len(done))
    solf=open(os.path.join(OUT,f'solutions_v{v}.txt'),'a'); nsol=0; exp_cum=0.0; T=300
    A,B,Cc,Dd=[D[a] for a in prof]
    for b,ix in enumerate(buckets):
        if b%C['NINST']!=C['INSTANCE']%C['NINST'] or b in done: continue
        if time.time()-T_START>C['max_hours']*3600: log(stage='time_limit'); break
        if b in (0,1,2) or b%25==0: print(f"progress: bucket {b}, solutions so far {nsol}, hours {round((time.time()-T_START)/3600,2)}",flush=True)
        tau=np.array([opts[q][ix[q]][0] for q in range(k)],np.int64)
        nA=partner(A['ch'],tau); oA=np.argsort(nA); nA=nA[oA]
        nC=partner(Cc['ch'],-tau); oC=np.argsort(nC); nC=nC[oC]
        hAB=merge_thist(nA,oA,A['c'],sB,oB,B['c'],64,T); hCD=merge_thist(nC,oC,Cc['c'],sD,oD,Dd['c'],64,T)
        both=hAB*hCD[::-1]; exp_b=float(both.sum()*2.0**(-3.42*(R.L-1-KE-1)))
        budget=C['mem_budget_bytes']/24; passes=[]; cur=[]; load=0
        for t in np.arange(-T,T+1)[np.argsort(-both)]:
            if both[t+T]==0: continue
            w=hAB[t+T]+hCD[T-t]
            if load+w>budget and cur: passes.append(cur); cur=[]; load=0
            cur.append(int(t)); load+=w
        if cur: passes.append(cur)
        found=0; ncand=0
        for pl in passes:
            mA=np.zeros(2*T+1,np.bool_); mC=np.zeros(2*T+1,np.bool_)
            for t in pl: mA[t+T]=True; mC[T-t]=True
            KA,IA,JA=stage1(nA,oA,A['c'],A['h'],sB,oB,B['c'],B['h'],1,G,64,mA,T)
            KC,IC,JC=stage1(nC,oC,Cc['c'],Cc['h'],sD,oD,Dd['c'],Dd['h'],-1,G,64,mC,T)
            ia,ic=bloom_match(KA,KC,30); ncand+=len(ia)
            for q,qc in zip(ia,ic):
                bits=[pools[prof[0]][IA[q]],pools[prof[1]][JA[q]],pools[prof[2]][IC[qc]],pools[prof[3]][JC[qc]]]
                Y=np.array([[1 if (int(bb)>>qq)&1 else -1 for qq in range(R.r)] for bb in bits]); X=R.expand(Y)
                if np.all(R.paf(X).sum(0)[1:]==0):
                    found+=1; nsol+=1
                    rec=dict(v=v,d=d,profile=list(prof),instance=C['INSTANCE'],bucket=b,
                             blocks={kk:''.join('+' if t>0 else '-' for t in x) for kk,x in zip('ABCD',X)})
                    solf.write(json.dumps(rec)+'\n'); solf.flush(); log(stage='SOLUTION',b=b)
                    print('SOLUTION '+json.dumps(rec),flush=True)
            del KA,IA,JA,KC,IC,JC
        exp_cum+=exp_b
        log(stage='bucket',b=b,tau=tau.tolist(),AB=int(hAB.sum()),CD=int(hCD.sum()),passes=len(passes),
            cand=int(ncand),solutions=found,total=nsol,exp_hits=round(exp_b,5),exp_cum=round(exp_cum,4))
    log(stage='end',solutions=nsol)
    solf.close()
    sp=os.path.join(OUT,f'solutions_v{v}.txt')
    if nsol and os.path.exists(sp):
        print('Verifying all solutions in',sp,flush=True); verify_file(sp)
    nb_done=0
    for l in open(os.path.join(OUT,'log_inst%d.jsonl'%C['INSTANCE'])):
        if '"stage": "bucket"' in l: nb_done+=1
    print('RUN SUMMARY: instance=%d buckets_done=%d solutions=%d'%(C['INSTANCE'],nb_done,nsol),flush=True)

def verify_file(path):
    """independent check of every solution: full GS array H H^T = 4v I"""
    for l in open(path):
        s=json.loads(l); X=np.array([[1 if c=='+' else -1 for c in s['blocks'][kk]] for kk in 'ABCD']); m=X.shape[1]
        circ=lambda x: np.array([np.roll(x,i) for i in range(m)]); A,B,Cm,Dm=map(circ,X); Rm=np.fliplr(np.eye(m,dtype=int))
        H=np.block([[A,B@Rm,Cm@Rm,Dm@Rm],[-B@Rm,A,Dm.T@Rm,-Cm.T@Rm],[-Cm@Rm,-Dm.T@Rm,A,B.T@Rm],[-Dm@Rm,Cm.T@Rm,-B.T@Rm,A]])
        print('bucket',s['bucket'],'Hadamard of order',4*m,':',np.array_equal(H@H.T,4*m*np.eye(4*m,dtype=int)))

if __name__=='__main__' or 'ipykernel' in sys.modules:
    main()
