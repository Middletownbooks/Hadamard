"""Exact stage B (lift) for the c=6 architecture, as a 4-list join.
Given a compressed solution xc (level t = h/p) and symmetry u:
 1. per block: ALL +-1 S-invariant blocks with K-fibre sums xc[k] (connected components of the fibre-orbit graph,
    brute force per component, Cartesian product), filtered by P(chi) <= 4v on characters nontrivial on H
 2. PAF at the lag-class reps (orbits of <S,-1>) and PSD at character-class reps for every survivor
 3. U = (A - T) + C pairs passing the pair PSD bound -> 64-bit linear hash of the PAF vector, sorted
    V = B + D pairs passing the bound -> look up -hash(V); exact verification of every candidate
usage (lib): run(h,p,u,xc) ; __main__ validates on the real 347 compression (h=57,p=3,u=11)."""
import sys, time, json, numpy as np, numba as nb
sys.path.insert(0,'scripts')

def setup(h,p,u):
    src=open('scripts/bge3_lift.py').read(); cut=src.index("t0=time.time(); stats=dict(")
    g={}; sys_argv=sys.argv; sys.argv=['bge3_lift.py',str(h),str(p),str(u),'0','1']
    exec(src[:cut].replace('cache=True','cache=False'),g); sys.argv=sys_argv
    return g

def components(A):
    """A (fibres x orbits) footprint: connected components over orbits (orbits linked by sharing a fibre)."""
    nf,no=A.shape; par=list(range(no))
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for fi in range(nf):
        qs=np.nonzero(A[fi])[0]
        for q in qs[1:]:
            a,b=f(qs[0]),f(q)
            if a!=b: par[a]=b
    comp={}
    for q in range(no): comp.setdefault(f(q),[]).append(q)
    return list(comp.values())

def enum_lifts(A,target):
    comps=components(A); parts=[]
    for cq in comps:
        fis=np.nonzero(A[:,cq].sum(1))[0]; sub=A[np.ix_(fis,cq)]
        n=len(cq); assert n<=24, 'component too large'
        Y=((np.arange(1<<n)[:,None]>>np.arange(n)[None,:])&1)*2-1
        ok=(Y@sub.T==target[fis][None,:]).all(1)
        parts.append((cq,Y[ok].astype(np.int8)))
        if not ok.any(): return np.zeros((0,A.shape[1]),np.int8)
    total=int(np.prod([len(Y) for _,Y in parts]))
    out=np.empty((total,A.shape[1]),np.int8); rep=total
    tile=1
    for cq,Y in parts:
        rep//=len(Y)
        blockv=np.repeat(Y,rep,axis=0)                    # each option repeated rep times
        out[:,cq]=np.tile(blockv,(tile,1)); tile*=len(Y)
    return out

@nb.njit(parallel=True,cache=True)
def pair_filter_count(PA,PB,bound):
    n=PA.shape[0]; cnt=np.zeros(n,np.int64)
    for i in nb.prange(n):
        c=0
        for j in range(PB.shape[0]):
            ok=True
            for f in range(PA.shape[1]):
                if PA[i,f]+PB[j,f]>bound: ok=False; break
            if ok: c+=1
        cnt[i]=c
    return cnt

@nb.njit(parallel=True,cache=True)
def pair_emit(PA,PB,HA,HB,bound,off,K,I,J):
    n=PA.shape[0]
    for i in nb.prange(n):
        c=off[i]
        for j in range(PB.shape[0]):
            ok=True
            for f in range(PA.shape[1]):
                if PA[i,f]+PB[j,f]>bound: ok=False; break
            if ok: K[c]=HA[i]+HB[j]; I[c]=i; J[c]=j; c+=1
    return 0

@nb.njit(parallel=True,cache=True)
def probe(PA,PB,HA,HB,bound,Ks,nch):
    """count/emit V pairs (i,j) with -(HA[i]+HB[j]) in sorted Ks"""
    n=PA.shape[0]; hits=np.zeros(n,np.int64)
    for i in nb.prange(n):
        c=0
        for j in range(PB.shape[0]):
            ok=True
            for f in range(PA.shape[1]):
                if PA[i,f]+PB[j,f]>bound: ok=False; break
            if ok:
                key=np.uint64(0)-(HA[i]+HB[j]); p=np.searchsorted(Ks,key)
                if p<Ks.shape[0] and Ks[p]==key: c+=1
        hits[i]=c
    return hits

def run(h,p,u,xc,log=print,maxpairs=3e8):
    g=setup(h,p,u); co=g['co']; m=g['m']; v=g['v']; lab=g['lab']; orbs=g['orbs']; ADD=g['ADD']; T_full=g['T_full']
    A=g['footB'].T; t0=time.time()
    # lag classes (orbits of <S,-1> on G\0) and character classes
    e=int(np.lcm.reduce([2,3,h])); perm=g['idx'](co*(u%e)); neg=g['NEG']
    cl=-np.ones(m,int); reps=[]
    for x in range(1,m):
        if cl[x]<0:
            S=set(); y=x
            while y not in S: S.add(y); S.add(int(neg[y])); y=int(perm[y])
            for z in S: cl[z]=len(reps)
            reps.append(x)
    reps=np.array(reps); L=len(reps)
    mods=np.array([2,3,h]); chars=co[co[:,2]!=0]
    # character classes: same orbit structure (multiplication by u, negation) on the dual
    cperm=g['idx'](chars*(u%e)); ccl={}; creps=[]
    for ch in chars:
        key=g['idx'](ch)
        if key in ccl: continue
        S=set(); y=int(key)
        while y not in S: S.add(y); S.add(int(neg[y])); y=int(perm[y])
        for z in S: ccl[z]=len(creps)
        creps.append(int(key))
    creps=np.array(creps)
    ph=np.exp(2j*np.pi*((co[:,None,:]*co[creps][None,:,:])/mods).sum(-1))
    lists=[]
    for k in range(4):
        Y=enum_lifts(A,xc[k]); F=Y[:,lab].astype(np.float32)
        PS=np.empty((len(Y),len(creps)),np.float32)
        for s in range(0,len(Y),50000): PS[s:s+50000]=(np.abs(F[s:s+50000]@ph)**2).astype(np.float32)
        ok=(PS<=4*v+1e-3).all(1); Y=Y[ok]; PS=PS[ok]
        Fk=Y[:,lab].astype(np.int64)
        PAF=np.empty((len(Y),L),np.int64)
        for s in range(0,len(Y),20000):
            Fs=np.fft.fftn(Fk[s:s+20000].reshape(-1,2,3,h).astype(float),axes=(1,2,3))
            R=np.rint(np.fft.ifftn(np.abs(Fs)**2,axes=(1,2,3)).real).reshape(len(Fs),-1)
            PAF[s:s+20000]=R[:,reps].astype(np.int64)
        lists.append(dict(Y=Y,PS=PS,PAF=PAF))
        log(json.dumps(dict(stage='list',block=k,lifts=int(ok.size),psd_pass=int(ok.sum()),secs=round(time.time()-t0,1))))
    # linear hash; A carries -T
    rng=np.random.default_rng(1); W=rng.integers(1,2**63,size=L,dtype=np.uint64)|np.uint64(1)
    tgt=T_full[reps]
    H=[((lists[k]['PAF']-(tgt if k==0 else 0)).astype(np.int64).astype(np.uint64)*W).sum(1) for k in range(4)]
    bound=np.float32(4*v+1e-3)
    # U = A + C, V = B + D
    PS=[lists[k]['PS'] for k in range(4)]
    # store the pairing with the smaller raw product; stream the other one through it
    if len(PS[1])*len(PS[3]) <= len(PS[0])*len(PS[2]): st,pr=(1,3),(0,2)
    else: st,pr=(0,2),(1,3)
    cnt=pair_filter_count(PS[st[0]],PS[st[1]],bound); nU=int(cnt.sum())
    log(json.dumps(dict(stage='stored_pairs',blocks=st,n=nU,of=int(len(PS[st[0]]))*int(len(PS[st[1]])),secs=round(time.time()-t0,1))))
    if nU>maxpairs: log(json.dumps(dict(stage='abort',reason='too many stored pairs'))); return None
    off=np.concatenate([[0],np.cumsum(cnt)[:-1]]).astype(np.int64)
    K=np.empty(nU,np.uint64); I=np.empty(nU,np.int32); J=np.empty(nU,np.int32)
    pair_emit(PS[st[0]],PS[st[1]],H[st[0]],H[st[1]],bound,off,K,I,J)
    o=np.argsort(K); K=K[o]; I=I[o]; J=J[o]
    hits=probe(PS[pr[0]],PS[pr[1]],H[pr[0]],H[pr[1]],bound,K,64)
    log(json.dumps(dict(stage='probe',candidates=int(hits.sum()),secs=round(time.time()-t0,1))))
    sols=[]
    for i in np.nonzero(hits)[0]:
        for j in range(len(PS[pr[1]])):
            if (PS[pr[0]][i]+PS[pr[1]][j]>bound).any(): continue
            key=np.uint64(0)-(H[pr[0]][i]+H[pr[1]][j]); pz=np.searchsorted(K,key)
            while pz<len(K) and K[pz]==key:
                idxs={st[0]:I[pz],st[1]:J[pz],pr[0]:i,pr[1]:j}
                Fq=np.array([lists[k]['Y'][idxs[k]][lab] for k in range(4)]).astype(np.int64)
                P=sum(np.array([int(np.dot(Fq[k],Fq[k][ADD[:,gg]])) for gg in range(m)]) for k in range(4))
                if np.array_equal(P[1:],T_full[1:]): sols.append(Fq)
                pz+=1
    log(json.dumps(dict(stage='done',solutions=len(sols),secs=round(time.time()-t0,1))))
    return sols


@nb.njit(parallel=True,cache=True)
def _emit_by_value(P0,P1,H0,H1,v0,v1s,o1,bound,taus,sign,cnt,off,K,I,J,write):
    """pairs (i in list0, j in list1 via sorted values v1s/o1) with v0[i]+v1[j] = sign*tau for tau in taus,
    passing the PSD bound; key = H0+H1."""
    n=P0.shape[0]
    for i in nb.prange(n):
        c=0; w=off[i] if write else 0
        for tt in range(taus.shape[0]):
            want=sign*taus[tt]-v0[i]
            lo=np.searchsorted(v1s,want); hi=np.searchsorted(v1s,want,side='right')
            for pz in range(lo,hi):
                j=o1[pz]; ok=True
                for f in range(P0.shape[1]):
                    if P0[i,f]+P1[j,f]>bound: ok=False; break
                if ok:
                    if write:
                        K[w]=H0[i]+H1[j]; I[w]=i; J[w]=j; w+=1
                    else: c+=1
        if not write: cnt[i]=c
    return 0

def join_batched(PS,PAF,H,cls0,bound,verify,log=print,batch=1e8):
    """4-list join sum_k PAF_k = 0 (list 0 carries -target), batched on the value of class cls0."""
    V=[np.ascontiguousarray(P[:,cls0]).astype(np.int64) for P in PAF]
    st,pr=((1,3),(0,2)) if len(PS[1])*len(PS[3])<=len(PS[0])*len(PS[2]) else ((0,2),(1,3))
    o1=np.argsort(V[st[1]]).astype(np.int64); v1s=V[st[1]][o1]; o3=np.argsort(V[pr[1]]).astype(np.int64); v3s=V[pr[1]][o3]
    from collections import Counter
    va,ca=np.unique(V[st[0]],return_counts=True); vb,cb=np.unique(V[st[1]],return_counts=True)
    hist=Counter()
    for a,x in zip(va,ca):
        for b_,y in zip(vb,cb): hist[int(a+b_)]+=int(x)*int(y)
    pa=np.unique(V[pr[0]]); pb=np.unique(V[pr[1]]); psums=set(int(a+b_) for a in pa for b_ in pb)
    taus=[tv for tv in sorted(hist) if -tv in psums]
    batches=[]; cur=[]; load=0
    for tv in taus:
        if load+hist[tv]>batch and cur: batches.append(cur); cur=[]; load=0
        cur.append(tv); load+=hist[tv]
    if cur: batches.append(cur)
    found=[]; tot_st=0; tot_c=0
    z=np.zeros(1,np.int64); zu=np.zeros(1,np.uint64); zi=np.zeros(1,np.int32)
    for bt in batches:
        tb=np.array(bt,np.int64)
        cnt=np.zeros(len(PS[st[0]]),np.int64)
        _emit_by_value(PS[st[0]],PS[st[1]],H[st[0]],H[st[1]],V[st[0]],v1s,o1,bound,tb,1,cnt,z,zu,zi,zi,False)
        n=int(cnt.sum()); tot_st+=n
        if n==0: continue
        off=np.concatenate([[0],np.cumsum(cnt)[:-1]]).astype(np.int64)
        K=np.empty(n,np.uint64); I=np.empty(n,np.int32); J=np.empty(n,np.int32)
        _emit_by_value(PS[st[0]],PS[st[1]],H[st[0]],H[st[1]],V[st[0]],v1s,o1,bound,tb,1,cnt,off,K,I,J,True)
        o=np.argsort(K); K=K[o]; I=I[o]; J=J[o]
        cntp=np.zeros(len(PS[pr[0]]),np.int64)
        _emit_by_value(PS[pr[0]],PS[pr[1]],H[pr[0]],H[pr[1]],V[pr[0]],v3s,o3,bound,tb,-1,cntp,z,zu,zi,zi,False)
        if int(cntp.sum())==0: continue
        # probe side in row chunks (bounded memory)
        cum=np.cumsum(cntp); r0=0
        while r0<len(cntp):
            base=cum[r0-1] if r0>0 else 0
            r1=int(np.searchsorted(cum,base+5e7,side='right')); r1=max(r1,r0+1)
            sub=np.arange(r0,min(r1,len(cntp)))
            mm=int(cntp[sub].sum())
            if mm:
                subcnt=cntp[sub].copy(); suboff=np.concatenate([[0],np.cumsum(subcnt)[:-1]]).astype(np.int64)
                Kp=np.empty(mm,np.uint64); Ip=np.empty(mm,np.int32); Jp=np.empty(mm,np.int32)
                Pm=np.ascontiguousarray(PS[pr[0]][sub]); Hm=np.ascontiguousarray(H[pr[0]][sub]); Vm=np.ascontiguousarray(V[pr[0]][sub])
                _emit_by_value(Pm,PS[pr[1]],Hm,H[pr[1]],Vm,v3s,o3,bound,tb,-1,subcnt,suboff,Kp,Ip,Jp,True)
                need=np.uint64(0)-Kp
                pos=np.searchsorted(K,need); hitm=(pos<len(K))&(K[np.minimum(pos,len(K)-1)]==need)
                for q in np.nonzero(hitm)[0]:
                    pz=pos[q]
                    while pz<len(K) and K[pz]==need[q]:
                        idx={st[0]:int(I[pz]),st[1]:int(J[pz]),pr[0]:int(sub[Ip[q]]),pr[1]:int(Jp[q])}
                        tot_c+=1
                        if verify(idx): found.append(idx)
                        pz+=1
            r0=min(r1,len(cntp))
    log(json.dumps(dict(stage='join_batched',batches=len(batches),stored_pairs=int(tot_st),candidates=int(tot_c),solutions=len(found))))
    return found

if __name__=='__main__':
    h,p,u=57,3,11
    g=setup(h,p,u); T=np.load('data/template347.npz'); Fr=np.array(T['F']).astype(np.int64)
    xc=np.zeros((4,g['mc']),np.int64)
    for k in range(4): np.add.at(xc[k],g['pi'],Fr[k])
    sols=run(h,p,u,xc,log=lambda s: print(s,flush=True))
    if sols:
        np.save('data/stageB_347_solutions.npy',np.array(sols))
        same=[any(np.array_equal(s[k],Fr[k]) for k in range(4)) for s in sols]
        print(json.dumps(dict(found=len(sols),matches_real_blocks=same)),flush=True)
