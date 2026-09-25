"""Exact stage A (compressed problem) for the c=6 architecture, as enumeration + 4-list join.
Compressed blocks xc on Gc = Z2 x Z3 x Z_t (t = h/p): values in {-p..p} (parity p), constant on S-orbits
(S = <u> acting by scalar multiplication), Q-coset sums = 347 pattern + shift, P(chi) <= 4v on characters nontrivial
on the Z_t part.  Join: sum_k PAFc_k = Tc at every compressed lag class (norm at lag 0 included).
lib: solve_A(h,p,u,log) -> list of (shift, xc[4 x mc]) compressed solutions."""
import sys, json, time, itertools, numpy as np, numba as nb
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup

def half_enum(vals,n):
    """all value assignments for n orbits (each from vals): (len(vals)^n, n) int8"""
    if n==0: return np.zeros((1,0),np.int8)
    g=np.array(np.meshgrid(*[vals]*n,indexing='ij')).reshape(n,-1).T
    return g.astype(np.int8)


@nb.njit(parallel=True,cache=True)
def _match(need,k2s,o2s,F1r,F1i,F2r,F2i,bound,cnt,I,J,off,write):
    n=need.shape[0]; nc=F1r.shape[1]
    for i in nb.prange(n):
        lo=np.searchsorted(k2s,need[i]); hi=np.searchsorted(k2s,need[i],side='right')
        c=0; w=off[i] if write else 0
        for pz in range(lo,hi):
            j=o2s[pz]; ok=True
            for f in range(nc):
                a=F1r[i,f]+F2r[j,f]; b=F1i[i,f]+F2i[j,f]
                if a*a+b*b>bound: ok=False; break
            if ok:
                if write: I[w]=i; J[w]=j; w+=1
                else: c+=1
        if not write: cnt[i]=c
    return 0

def enum_block(g,cs_k,v,p,maxout=5*10**7):
    h_=g['h']; rc=g['rc']; footA=g['footA']; mc=g['mc']; labc=g['labc']; coc=g['coc']; t=g['t']
    mods=np.array([2,3,t]); chars=coc[coc[:,2]!=0]
    ph=np.exp(2j*np.pi*((coc[:,None,:]*chars[None,:,:])/mods).sum(-1))       # (element, char)
    OP=np.zeros((rc,len(chars)),complex)                                        # orbit -> Fourier contribution
    for x in range(mc): OP[labc[x]]+=ph[x]
    vals=np.arange(-p,p+1,2)
    order=np.argsort(-np.abs(OP).sum(1)); n1=rc//2; o1=order[:n1]; o2=order[n1:]
    Y1=half_enum(vals,len(o1)); Y2=half_enum(vals,len(o2))
    S1=Y1.astype(np.int64)@footA[o1]; S2=Y2.astype(np.int64)@footA[o2]
    F1=Y1.astype(float)@OP[o1]; F2=Y2.astype(float)@OP[o2]
    # match coset sums S1 + S2 = cs_k with numba (keys = coset-sum vectors encoded in int64)
    base=4*h_+1
    enc=lambda S: ((S+2*h_)*(base**np.arange(6))).sum(1).astype(np.int64)
    k1=enc(S1); k2=enc(S2); o2s=np.argsort(k2); k2s=k2[o2s]
    need=enc(np.array(cs_k)[None,:]-S1)
    F1r=np.ascontiguousarray(F1.real); F1i=np.ascontiguousarray(F1.imag); F2r=np.ascontiguousarray(F2.real); F2i=np.ascontiguousarray(F2.imag)
    bound=4*v+1e-6
    cnt=np.zeros(len(S1),np.int64); dummyI=np.zeros(1,np.int64); dummyJ=np.zeros(1,np.int64)
    _match(need,k2s,o2s,F1r,F1i,F2r,F2i,bound,cnt,dummyI,dummyJ,np.zeros(len(S1),np.int64),False)
    tot=int(cnt.sum()); trunc=tot>maxout
    if tot==0: return np.zeros((0,rc),np.int8), False
    if trunc: tot=0; cnt[:]=0
    off=np.concatenate([[0],np.cumsum(cnt)[:-1]]).astype(np.int64)
    I=np.empty(max(tot,1),np.int64); J=np.empty(max(tot,1),np.int64)
    if not trunc: _match(need,k2s,o2s,F1r,F1i,F2r,F2i,bound,cnt,I,J,off,True)
    Y=np.empty((tot,rc),np.int8); Y[:,o1]=Y1[I[:tot]]; Y[:,o2]=Y2[J[:tot]]
    return Y, trunc

def solve_A(h,p,u,log=print,max_solutions=200,shift_limit=None):
    g=setup(h,p,u); v=g['v']; mc=g['mc']; labc=g['labc']; coc=g['coc']; t=g['t']
    ADDc=g['ADDc']; NEGc=g['NEGc']; Tc=g['Tc']; shifts=g['shifts']; e=int(np.lcm.reduce([2,3,t]))
    # compressed lag classes (orbits of <S,-1> on Gc, including 0)
    permc=g['idxc'](coc*(u%e)); cl=-np.ones(mc,int); reps=[]
    for x in range(mc):
        if cl[x]<0:
            S=set(); y=x
            while y not in S: S.add(y); S.add(int(NEGc[y])); y=int(permc[y])
            for z in S: cl[z]=len(reps)
            reps.append(x)
    reps=np.array(reps)
    mods=np.array([2,3,t]); chars=coc[coc[:,2]!=0]
    ph=np.exp(2j*np.pi*((coc[:,None,:]*chars[None,:,:])/mods).sum(-1))
    sols=[]; t0=time.time()
    shifts=sorted(shifts,key=lambda x: sum(abs(y) for y in x[0]))
    for si,(s,cs) in enumerate(shifts[:shift_limit] if shift_limit else shifts):
        lists=[]
        for k in range(4):
            Y,trunc=enum_block(g,cs[k],v,p)
            X=Y[:,labc].astype(np.int64)                                       # expanded compressed blocks
            PS=(np.abs(X.astype(float)@ph)**2).astype(np.float32)
            PAF=np.empty((len(X),len(reps)),np.int64)
            for a in range(0,len(X),20000):
                Fs=np.fft.fftn(X[a:a+20000].reshape(-1,2,3,t).astype(float),axes=(1,2,3))
                R=np.rint(np.fft.ifftn(np.abs(Fs)**2,axes=(1,2,3)).real).reshape(len(Fs),-1)
                PAF[a:a+20000]=R[:,reps].astype(np.int64)
            lists.append(dict(X=X,PS=PS,PAF=PAF,trunc=trunc))
        log(json.dumps(dict(stage='A_lists',shift=[int(x) for x in s],sizes=[len(L['X']) for L in lists],truncated=[L['trunc'] for L in lists],secs=round(time.time()-t0,1))))
        if any(len(L['X'])==0 for L in lists): continue
        # join: (A - Tc) + B + C + D = 0 on all compressed lag classes, pair PSD bound
        rng=np.random.default_rng(3); W=rng.integers(1,2**63,size=len(reps),dtype=np.uint64)|np.uint64(1)
        tg=Tc[reps]
        H=[((lists[k]['PAF']-(tg if k==0 else 0)).astype(np.int64).astype(np.uint64)*W).sum(1) for k in range(4)]
        from kt_generic import kt_join4
        PAFl=[lists[k]['PAF']-(tg if k==0 else 0) for k in range(4)]
        n0=len(sols)
        def verify(idx):
            Xq=np.array([lists[k]['X'][idx[k]] for k in range(4)]).astype(np.int64)
            P=sum(np.array([int(np.dot(Xq[k],Xq[k][ADDc[:,gg]])) for gg in range(mc)]) for k in range(4))
            if np.array_equal(P,Tc): sols.append((s,Xq)); return True
            return False
        kt_join4(PAFl,verify,k=6,nbuckets=int(__import__('os').environ.get('A_BUCKETS','400')),max_solutions=min(max_solutions-len(sols),int(__import__('os').environ.get('PER_SHIFT','1000000'))),log=log,seed=si)
        found=len(sols)-n0; n=0; hits=np.zeros(1)
        log(json.dumps(dict(stage='A_join',shift=[int(x) for x in s],stored_pairs=n,candidates=int(hits.sum()),solutions=found,total=len(sols),secs=round(time.time()-t0,1))))
        if len(sols)>=max_solutions: break
    return sols

if __name__=='__main__':
    h,p,u=[int(x) for x in sys.argv[1:4]]
    sols=solve_A(h,p,u,log=lambda s: print(s,flush=True),max_solutions=int(sys.argv[4]) if len(sys.argv)>4 else 200)
    np.save(f'data/stageA_h{h}_p{p}_u{u}.npy',np.array([x for _,x in sols]) if sols else np.zeros((0,)))
    if h==57:   # validation: is the real 347 compression among the solutions?
        g=setup(h,p,u); Fr=np.array(np.load('data/template347.npz')['F']).astype(np.int64)
        xr=np.zeros((4,g['mc']),np.int64)
        for k in range(4): np.add.at(xr[k],g['pi'],Fr[k])
        print(json.dumps(dict(n_solutions=len(sols),contains_real=any(np.array_equal(x,xr) for _,x in sols))),flush=True)
