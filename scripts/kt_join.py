"""k-tree (Wagner) 4-list join on flat pools, b=0 reduced (v,d).
Slots A,B,C,D with row sums a1..a4 (pools may coincide).  Lag classes split into:
  chunk K (k<=8 classes): stage 1 requires A+B = 0 and C+D = 0 on K
  bucket class c:          stage 1 emits t = A_c+B_c ; stage 2 requires t_CD = -t
  rest:                    linear 64-bit hash; stage 2 requires h_AB + h_CD = 0
Candidates are verified exactly from the packed orbit signs.
usage: kt_join.py v d profile(csv) k pools_dir [npasses]"""
import sys, os, time, json, numpy as np, numba as nb
sys.path.insert(0,'scripts')
from redlib import Red

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

def prep(bits,Nsp,R,chunk,cc,rest,W,tau):
    P=paf_from_bits(bits,*Nsp,R.L,R.r)
    assert np.abs(P).max()<127
    ch=P[:,chunk].astype(np.int8)
    key=np.zeros(len(bits),np.uint64)
    for i in range(len(chunk)): key|=(ch[:,i].astype(np.uint8).astype(np.uint64))<<np.uint64(8*i)
    nkey=np.zeros(len(bits),np.uint64); nch=(tau[None,:]-ch.astype(np.int16)).astype(np.int8)   # partner key: tau - x
    for i in range(len(chunk)): nkey|=(nch[:,i].astype(np.uint8).astype(np.uint64))<<np.uint64(8*i)
    h=(P[:,rest].astype(np.int64).astype(np.uint64)*W[None,:]).sum(1)
    return dict(key=key,nkey=nkey,c=P[:,cc].astype(np.int16),h=h)

@nb.njit(cache=True)
def stage1(nkeyX,cX,hX,skeyY,ordY,cY,hY,tlo,thi,sign,G,outK,outI,outJ,count_only):
    """for each x: y with key(y)=nkey(x); emit if tlo<=t<thi. returns count"""
    n=0; NY=skeyY.shape[0]
    for i in range(nkeyX.shape[0]):
        k=nkeyX[i]
        lo=np.searchsorted(skeyY,k); 
        while lo<NY and skeyY[lo]==k:
            j=ordY[lo]; t=cX[i]+cY[j]
            if tlo<=t and t<thi:
                if not count_only:
                    kk=hX[i]+hY[j]+np.uint64(np.int64(t))*G
                    if sign<0: kk=np.uint64(0)-kk
                    outK[n]=kk; outI[n]=i; outJ[n]=j
                n+=1
            lo+=1
    return n

def run(v,d,prof,k,pools,npasses=1,seed=0,log=print,tau=None):
    R=Red(v,d); Nsp=sparse_N(R); rng=np.random.default_rng(seed)
    classes=rng.permutation(R.L); chunk=classes[:k]; cc=classes[k]; rest=classes[k+1:]
    W=rng.integers(1,2**63,size=len(rest),dtype=np.uint64)|np.uint64(1)
    G=np.uint64(0x9E3779B97F4A7C15)
    t0=time.time(); D={}
    tau=np.full(k,2,np.int16) if tau is None else np.asarray(tau,np.int16)
    # slots A,C need partner keys (tau - x) and (-tau - x); B,D need own keys
    for s_,a in enumerate(prof):
        tt=tau if s_==0 else -tau
        D[(s_,a)]=prep(pools[a],Nsp,R,chunk,cc,rest,W,tt)
    log(dict(stage='prep',secs=round(time.time()-t0,1),sizes={int(a):len(pools[a]) for a in set(prop) } if False else {int(a):len(pools[a]) for a in set(prof)}))
    A,B,C,Dd=[D[(s_,a)] for s_,a in enumerate(prof)]
    def sorted_side(Y):
        o=np.argsort(Y['key'],kind='stable'); return Y['key'][o],o
    sB,oB=sorted_side(B); sD,oD=sorted_side(Dd)
    # t range
    tmax=int(max(np.abs(X['c']).max() for X in D.values()))*2+1
    run.chunk=chunk; run.cc=cc
    edges=np.linspace(-tmax,tmax+1,npasses+1).astype(int)
    hits=[]; tot1=[0,0]
    for p in range(npasses):
        lo,hi=edges[p],edges[p+1]
        nAB=stage1(A['nkey'],A['c'],A['h'],sB,oB,B['c'],B['h'],lo,hi,1,G,np.zeros(1,np.uint64),np.zeros(1,np.int64),np.zeros(1,np.int64),True)
        nCD=stage1(C['nkey'],C['c'],C['h'],sD,oD,Dd['c'],Dd['h'],-hi+1,-lo+1,-1,G,np.zeros(1,np.uint64),np.zeros(1,np.int64),np.zeros(1,np.int64),True)
        KA=np.zeros(nAB,np.uint64); IA=np.zeros(nAB,np.int64); JA=np.zeros(nAB,np.int64)
        KC=np.zeros(nCD,np.uint64); IC=np.zeros(nCD,np.int64); JC=np.zeros(nCD,np.int64)
        stage1(A['nkey'],A['c'],A['h'],sB,oB,B['c'],B['h'],lo,hi,1,G,KA,IA,JA,False)
        stage1(C['nkey'],C['c'],C['h'],sD,oD,Dd['c'],Dd['h'],-hi+1,-lo+1,-1,G,KC,IC,JC,False)
        tot1[0]+=nAB; tot1[1]+=nCD
        oc=np.argsort(KC); KCs=KC[oc]
        pos=np.searchsorted(KCs,KA); pos[pos>=len(KCs)]=max(len(KCs)-1,0)
        m=np.where((len(KCs)>0)&(KCs[pos]==KA))[0] if len(KCs) else np.array([],int)
        for q in m:
            pp=pos[q]
            while pp<len(KCs) and KCs[pp]==KA[q]:
                hits.append((IA[q],JA[q],IC[oc[pp]],JC[oc[pp]])); pp+=1
        log(dict(stage='pass',p=p,t=[int(lo),int(hi)],AB=int(nAB),CD=int(nCD),cand=len(hits),secs=round(time.time()-t0,1)))
        del KA,IA,JA,KC,IC,JC
    # verify
    sols=[]
    for (i,j,kk,l) in hits:
        bits=[pools[prof[0]][i],pools[prof[1]][j],pools[prof[2]][kk],pools[prof[3]][l]]
        Y=np.array([[1 if (int(b)>>q)&1 else -1 for q in range(R.r)] for b in bits])
        X=R.expand(Y); tot=R.paf(X).sum(0)
        if np.all(tot[1:]==0): sols.append([''.join('+' if t>0 else '-' for t in y) for y in Y])
    log(dict(stage='done',stage1=tot1,candidates=len(hits),solutions=len(sols),secs=round(time.time()-t0,1)))
    return sols,tot1

@nb.njit(parallel=True,cache=True)
def stage1_par(nkeyX,cX,hX,skeyY,ordY,cY,hY,sign,G,nchunk):
    NX=nkeyX.shape[0]; NY=skeyY.shape[0]
    bounds=np.linspace(0,NX,nchunk+1).astype(np.int64)
    cnt=np.zeros(nchunk,np.int64)
    for b in nb.prange(nchunk):
        n=0
        for i in range(bounds[b],bounds[b+1]):
            k=nkeyX[i]; lo=np.searchsorted(skeyY,k)
            while lo<NY and skeyY[lo]==k:
                n+=1; lo+=1
        cnt[b]=n
    off=np.zeros(nchunk+1,np.int64)
    for b in range(nchunk): off[b+1]=off[b]+cnt[b]
    K=np.empty(off[nchunk],np.uint64); I=np.empty(off[nchunk],np.int64); J=np.empty(off[nchunk],np.int64)
    for b in nb.prange(nchunk):
        n=off[b]
        for i in range(bounds[b],bounds[b+1]):
            k=nkeyX[i]; lo=np.searchsorted(skeyY,k)
            while lo<NY and skeyY[lo]==k:
                j=ordY[lo]; t=cX[i]+cY[j]
                kk=hX[i]+hY[j]+np.uint64(np.int64(t))*G
                if sign<0: kk=np.uint64(0)-kk
                K[n]=kk; I[n]=i; J[n]=j; n+=1; lo+=1
    return K,I,J

@nb.njit(cache=True)
def _merge_count(sX,sY,lo,hi):
    # count matches for X-range [lo,hi) against sorted Y (both sorted)
    NY=sY.shape[0]; j=np.searchsorted(sY,sX[lo]); n=0; i=lo
    while i<hi and j<NY:
        a=sX[i]; b=sY[j]
        if a<b: i+=1
        elif a>b: j+=1
        else:
            j2=j
            while j2<NY and sY[j2]==a: j2+=1
            i2=i
            while i2<hi and sX[i2]==a: i2+=1
            n+=(i2-i)*(j2-j); i=i2; j=j2
    return n

@nb.njit(parallel=True,cache=True)
def stage1_merge(sX,oX,cX,hX,sY,oY,cY,hY,sign,G,nchunk):
    """sX: sorted partner keys of X (oX: original idx); sY: sorted own keys of Y. Sequential merge join."""
    NX=sX.shape[0]; NY=sY.shape[0]
    b0=np.linspace(0,NX,nchunk+1).astype(np.int64)
    # align chunk boundaries to key changes
    for b in range(1,nchunk):
        x=max(b0[b],b0[b-1])
        while x<NX and x>0 and sX[x]==sX[x-1]: x+=1
        b0[b]=x
    cnt=np.zeros(nchunk,np.int64)
    for b in nb.prange(nchunk):
        if b0[b]<b0[b+1]: cnt[b]=_merge_count(sX,sY,b0[b],b0[b+1])
    off=np.zeros(nchunk+1,np.int64)
    for b in range(nchunk): off[b+1]=off[b]+cnt[b]
    K=np.empty(off[nchunk],np.uint64); I=np.empty(off[nchunk],np.int64); J=np.empty(off[nchunk],np.int64)
    for b in nb.prange(nchunk):
        lo=b0[b]; hi=b0[b+1]
        if lo>=hi: continue
        n=off[b]; j=np.searchsorted(sY,sX[lo]); i=lo
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
                    xi=oX[ii]
                    for jj in range(j,j2):
                        yj=oY[jj]; t=cX[xi]+cY[yj]
                        kk=hX[xi]+hY[yj]+np.uint64(np.int64(t))*G
                        if sign<0: kk=np.uint64(0)-kk
                        K[n]=kk; I[n]=xi; J[n]=yj; n+=1
                i=i2; j=j2
    return K,I,J

@nb.njit(parallel=True,cache=True)
def merge_count_total(sX,sY,nchunk):
    NX=sX.shape[0]
    b0=np.linspace(0,NX,nchunk+1).astype(np.int64)
    for b in range(1,nchunk):
        x=max(b0[b],b0[b-1])
        while x<NX and x>0 and sX[x]==sX[x-1]: x+=1
        b0[b]=x
    cnt=np.zeros(nchunk,np.int64)
    for b in nb.prange(nchunk):
        if b0[b]<b0[b+1]: cnt[b]=_merge_count(sX,sY,b0[b],b0[b+1])
    return cnt.sum()

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
