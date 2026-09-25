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

def enum_block(g,cs_k,v,p,maxout=5*10**7):
    rc=g['rc']; footA=g['footA']; mc=g['mc']; labc=g['labc']; coc=g['coc']; t=g['t']
    mods=np.array([2,3,t]); chars=coc[coc[:,2]!=0]
    ph=np.exp(2j*np.pi*((coc[:,None,:]*chars[None,:,:])/mods).sum(-1))       # (element, char)
    OP=np.zeros((rc,len(chars)),complex)                                        # orbit -> Fourier contribution
    for x in range(mc): OP[labc[x]]+=ph[x]
    vals=np.arange(-p,p+1,2)
    order=np.argsort(-np.abs(OP).sum(1)); n1=rc//2; o1=order[:n1]; o2=order[n1:]
    Y1=half_enum(vals,len(o1)); Y2=half_enum(vals,len(o2))
    S1=Y1.astype(np.int64)@footA[o1]; S2=Y2.astype(np.int64)@footA[o2]
    F1=Y1.astype(float)@OP[o1]; F2=Y2.astype(float)@OP[o2]
    # match coset sums: S1 + S2 = cs_k
    key=lambda S: [tuple(r) for r in S]
    from collections import defaultdict
    idx2=defaultdict(list)
    for j,s in enumerate(key(S2)): idx2[s].append(j)
    out=[]; nout=0; bound=4*v+1e-6
    for i,s in enumerate(key(S1)):
        need=tuple(np.array(cs_k)-np.array(s))
        js=idx2.get(need)
        if not js: continue
        js=np.array(js); P=np.abs(F1[i][None,:]+F2[js])**2
        ok=js[(P<=bound).all(1)]
        if len(ok):
            Y=np.empty((len(ok),rc),np.int8); Y[:,o1]=Y1[i]; Y[:,o2]=Y2[ok]
            out.append(Y); nout+=len(ok)
            if nout>maxout: break
    return (np.concatenate(out) if out else np.zeros((0,rc),np.int8)), nout>maxout

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
        from bge3_stageB_join import pair_filter_count, pair_emit, probe
        bound=np.float32(4*v+1e-3); PS=[L['PS'] for L in lists]
        st,pr=((1,3),(0,2)) if len(PS[1])*len(PS[3])<=len(PS[0])*len(PS[2]) else ((0,2),(1,3))
        cnt=pair_filter_count(PS[st[0]],PS[st[1]],bound); n=int(cnt.sum())
        if n>3e8: log(json.dumps(dict(stage='A_abort',shift=[int(x) for x in s],pairs=n))); continue
        off=np.concatenate([[0],np.cumsum(cnt)[:-1]]).astype(np.int64)
        K=np.empty(n,np.uint64); I=np.empty(n,np.int32); J=np.empty(n,np.int32)
        pair_emit(PS[st[0]],PS[st[1]],H[st[0]],H[st[1]],bound,off,K,I,J)
        o=np.argsort(K); K=K[o]; I=I[o]; J=J[o]
        hits=probe(PS[pr[0]],PS[pr[1]],H[pr[0]],H[pr[1]],bound,K,64)
        found=0
        for i in np.nonzero(hits)[0]:
            for j in range(len(PS[pr[1]])):
                if (PS[pr[0]][i]+PS[pr[1]][j]>bound).any(): continue
                key=np.uint64(0)-(H[pr[0]][i]+H[pr[1]][j]); pz=np.searchsorted(K,key)
                while pz<len(K) and K[pz]==key:
                    idxs={st[0]:I[pz],st[1]:J[pz],pr[0]:i,pr[1]:j}
                    Xq=np.array([lists[k]['X'][idxs[k]] for k in range(4)]).astype(np.int64)
                    P=sum(np.array([int(np.dot(Xq[k],Xq[k][ADDc[:,gg]])) for gg in range(mc)]) for k in range(4))
                    if np.array_equal(P,Tc): sols.append((s,Xq)); found+=1
                    pz+=1
                    if len(sols)>=max_solutions: break
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
