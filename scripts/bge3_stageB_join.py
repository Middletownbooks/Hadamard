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
