"""Two-sided compression-and-lift for c=6 when h = p*q (e.g. 515: h=85=5*17).
A lift F (+-1 on G = Z2 x Z3 x Z_h, constant on <u>-orbits) has two compressions: xc onto Z_t (t=h/p, from stage A at p)
and yc onto Z_p (fibres z mod p).  Given xc, the lift splits into independent components (fibre-orbit graph); each
component's admissible sign patterns contribute a vector to yc.  Step 1 (this file, lib): the exact distribution of yc
over all lifts of xc (convolution over components), i.e. the lift list is partitioned by yc.
lib: yc_setup(g,p2) ; yc_distribution(g,comps,xck,Ymap) -> dict yc_tuple -> log2 count"""
import numpy as np, math

def yc_setup(g,p2):
    """Z_{p2} fibre map of the full group and per-orbit footprint on it (orbit -> counts per yc cell)."""
    co=g['co']; h=g['h'] if 'h' in g else int(co[:,2].max())+1
    cell=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2)
    A=g['footB'].T                                   # (fibres_xc, orbits)
    no=A.shape[1]; foot=np.zeros((no,6*p2),np.int64)
    for x in range(len(co)): foot[g['lab'][x],cell[x]]+=1
    return foot

def comp_options(A,cq,xck):
    fis=np.nonzero(A[:,cq].sum(1))[0]; sub=A[np.ix_(fis,cq)]; n=len(cq)
    Y=((np.arange(1<<n)[:,None]>>np.arange(n)[None,:])&1)*2-1
    return Y[(Y@sub.T==xck[fis][None,:]).all(1)]

def yc_distribution(g,comps,xck,foot,maxstates=5_000_000):
    """exact distribution {yc (tuple): count} of the Z_p2 compression over all lifts of xck; None if too many states."""
    A=g['footB'].T; dist={(0,)*foot.shape[1]:1}
    # process components with the fewest distinct contributions first
    items=[]
    for cq in comps:
        Y=comp_options(A,cq,xck)
        if len(Y)==0: return {}
        contrib=Y@foot[cq]; u,c=np.unique(contrib,axis=0,return_counts=True)
        items.append([(tuple(int(x) for x in r),int(n)) for r,n in zip(u,c)])
    items.sort(key=len)
    for it in items:
        nd={}
        for k,c in dist.items():
            for r,n in it:
                kk=tuple(a+b for a,b in zip(k,r)); nd[kk]=nd.get(kk,0)+c*n
        dist=nd
        if len(dist)>maxstates: return None
    return dist

def yc_lists(g,comps,xc,p2,S,rng,log=print,chunk=200000):
    """sample S lifts per block, keep PSD-passing, return per block: distinct yc (cells over Z2xZ3xZ_p2) and hit counts."""
    from bge3_residue import lift_sampler
    co=g['co']; v=g['v']; h=int(co[:,2].max())+1; mods=np.array([2,3,h]); chars=co[co[:,2]!=0]
    ph=np.exp(2j*np.pi*((co[:,None,:]*chars[None,:,:])/mods).sum(-1))
    cell=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2); out=[]
    for k in range(4):
        acc={}
        for a in range(0,S,chunk):
            F=lift_sampler(g,comps,xc[k],min(chunk,S-a),rng).astype(float)
            F=F[((np.abs(F@ph)**2)<=4*v+1e-6).all(1)]
            if len(F)==0: continue
            Yc=np.zeros((len(F),6*p2)); np.add.at(Yc.T,cell,F.T); uq,c=np.unique(Yc.astype(np.int64),axis=0,return_counts=True)
            for r,n in zip(map(tuple,uq.tolist()),c.tolist()): acc[r]=acc.get(r,0)+n
        Y=np.array(list(acc.keys()),np.int64).reshape(-1,6*p2); N=np.array(list(acc.values()),np.int64)
        log(dict(stage='yc_list',block=k,distinct=int(len(Y)),hits=int(N.sum()))); out.append((Y,N))
    return out

def yc_join(g,ycl,p2,u,log=print,nbuckets=2000,max_solutions=100):
    """4-list join on the Z_p2-compressed PAF: find (yA,yB,yC,yD) from the lists with sum_k PAF(y_k) = T (compressed)."""
    import sys; sys.path.insert(0,'scripts'); from kt_generic import kt_join4
    co=g['co']; cell=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2)
    co2=np.array([(a,b,z) for a in range(2) for b in range(3) for z in range(p2)]); idx2=lambda X:((X[...,0]%2)*3+(X[...,1]%3))*p2+(X[...,2]%p2)
    T2=np.zeros(6*p2,np.int64); np.add.at(T2,cell,g['T_full'])
    seen=set(); reps=[]
    for x in range(6*p2):
        if x in seen: continue
        y=x; orb=set()
        while y not in orb: orb.add(y); orb.add(int(idx2(-co2[y]))); y=int(idx2(co2[y]*u))
        seen|=orb; reps.append(x)
    reps=np.array(reps)
    def paf(Y):
        F=Y.reshape(-1,2,3,p2).astype(float)
        return np.rint(np.fft.ifftn(np.abs(np.fft.fftn(F,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(Y),-1).astype(np.int64)
    PF=[paf(Y) for Y,_ in ycl]; PAFl=[PF[k][:,reps]-(T2[reps] if k==0 else 0) for k in range(4)]
    sols=[]
    def verify(ix):
        tot=sum(PF[k][ix[k]] for k in range(4))
        if np.array_equal(tot,T2): sols.append(tuple(int(i) for i in ix)); return True
        return False
    kt_join4(PAFl,verify,nbuckets=nbuckets,max_solutions=max_solutions,log=lambda s: log(s) if ('kt_done' in s or 'kt_setup' in s or 'kt_buckets' in s) else None,target_pairs=1e6)
    return sols

def _half_product(parts,no):
    """all combinations over a list of (cq, options Y, contribution C): returns (signs (n,no) int8, contrib (n,d) int64)."""
    total=int(np.prod([len(Y) for _,Y,_ in parts])) if parts else 1
    d=parts[0][2].shape[1] if parts else 0
    S=np.zeros((total,no),np.int8); Cs=np.zeros((total,d),np.int64); rep=total; tile=1
    for cq,Y,C in parts:
        rep//=len(Y); idx=np.tile(np.repeat(np.arange(len(Y)),rep),tile); tile*=len(Y)
        S[:,cq]=Y[idx]; Cs+=C[idx]
    return S,Cs

def enum_lifts_yc(g,comps,xck,yck,foot,maxhalf=2**24):
    """all lifts (orbit signs) of block compression xck whose z-mod-p2 compression equals yck (meet in the middle over
    components, matched on the yc contribution)."""
    A=g['footB'].T; no=A.shape[1]; parts=[]
    for cq in comps:
        Y=comp_options(A,cq,xck)
        if len(Y)==0: return np.zeros((0,no),np.int8)
        parts.append((cq,Y.astype(np.int8),(Y@foot[cq]).astype(np.int64)))
    # balance halves by log size
    parts.sort(key=lambda t:-len(t[1])); H1=[];H2=[]; l1=l2=0.0
    for pt in parts:
        if l1<=l2: H1.append(pt); l1+=np.log2(len(pt[1]))
        else: H2.append(pt); l2+=np.log2(len(pt[1]))
    if max(l1,l2)>np.log2(maxhalf): raise MemoryError(f'half too large: 2^{max(l1,l2):.1f}')
    S1,C1=_half_product(H1,no); S2,C2=_half_product(H2,no)
    rng=np.random.default_rng(7); W=rng.integers(1,2**62,size=C1.shape[1],dtype=np.int64)
    k1=(C1*W).sum(1); k2=((np.asarray(yck)[None,:]-C2)*W).sum(1)
    o=np.argsort(k1); k1s=k1[o]; lo=np.searchsorted(k1s,k2); hi=np.searchsorted(k1s,k2,side='right')
    n=int((hi-lo).sum()); out=np.empty((n,no),np.int8); w=0
    for j in np.nonzero(hi>lo)[0]:
        ii=o[lo[j]:hi[j]]; out[w:w+len(ii)]=S1[ii]+S2[j]; w+=len(ii)
    out=out[:w]
    # exact check of the yc contribution (hash collisions)
    ok=((out.astype(np.int64)@foot)==np.asarray(yck)[None,:]).all(1)
    return out[ok]
