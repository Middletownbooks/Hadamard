"""b>=3 cores (347, 359, 479): anomaly tests with symmetry-aware nulls.
Aut(G) is enumerated in full (images of the standard generators).  For each block: its stabiliser in the
affine group Aut(G) x G (up to sign) and the resulting affine orbits.  Null blocks are random +-1 blocks
constant on the SAME affine orbits with the same row sum, so any known symmetry is built into the null.
Tests
  P. partial symmetries: (sigma, t) with x(sigma g + t) = s x(g) exactly on some H-cosets but not globally
     (count of (sigma,t,coset) triples; and the largest number of cosets covered by a non-global map)
  X. block-pair relations: max over (sigma, t) of |sum_g x_i(sigma g + t) x_j(g)| / m
  Q. which blocks carry the quotient deviation (compressions are forced near-constant: sum of power at
     non-trivial G/H characters = 4b), reported descriptively.
usage: bge2_forensics2.py NNULL"""
import sys, json, itertools, numpy as np
sys.path.insert(0,'scripts')
from hadlib import all_cores
NN=int(sys.argv[1]) if len(sys.argv)>1 else 40
rng=np.random.default_rng(0)

def automorphisms(c):
    mods=np.array(c.mods); co=c.coords(); n=len(mods)
    # candidates for image of e_i: elements x with mods[i]*x == 0
    cand=[np.where((((co*mods[i])%mods)==0).all(1))[0] for i in range(n)]
    auts=[]
    for imgs in itertools.product(*cand):
        I=co[list(imgs)]                                   # (n, n) image coords of generators
        img=(co@I)%mods                                    # image of every element
        perm=c.idx(img)
        if len(np.unique(perm))==c.m: auts.append(perm)
    return auts

def corr_all(y,z,c):
    """A(t) = sum_g y(g+t) z(g) for all t (group FFT)"""
    Y=np.fft.fftn(c.arr(y).astype(float)); Z=np.fft.fftn(c.arr(z).astype(float))
    return np.rint(np.fft.ifftn(Y*np.conj(Z)).real.reshape(-1)).astype(np.int64)

def affine_orbits(x,c,auts,AT):
    parent=np.arange(c.m)
    def f(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    nst=0
    for P in auts:
        A=corr_all(x[P],x,c)                               # sum_g x(sigma(g+t))... see below
        # A(t) = sum_g x(P[g+t]) x(g); exact symmetry g -> P[g+t] when A(t) = m
        for t in np.where(A==c.m)[0]:
            nst+=1; img=P[AT[:,t]]
            for g in range(c.m):
                a,b=f(g),f(int(img[g]))
                if a!=b: parent[a]=b
    lab=np.array([f(g) for g in range(c.m)]); _,lab=np.unique(lab,return_inverse=True)
    return lab,nst

def random_like(x,lab,rng):
    sizes=np.bincount(lab); a=int(x.sum())
    while True:
        s=rng.choice([-1,1],size=len(sizes))
        if s@sizes==a: return s[lab]

def partial_stats(x,c,auts,cos,nq,Hs):
    """count (sigma,t,q) exact coset agreements that are not part of a global (sign) symmetry"""
    Z=[x*(cos==q) for q in range(nq)]
    cnt=0; best=0
    for P in auts:
        y=x[P]
        Aq=np.array([corr_all(y,Z[q],c) for q in range(nq)])      # (nq, m)
        full=(np.abs(Aq)==Hs)
        glob=(Aq==Hs).all(0)|(Aq==-Hs).all(0)
        k=full[:,~glob].sum(0)
        cnt+=int(k.sum()); best=max(best,int(k.max()) if k.size else 0)
    return cnt,best

def pair_stat(xi,xj,c,auts):
    return max(int(np.abs(corr_all(xi[P],xj,c)).max()) for P in auts)

for c in all_cores():
    if c.b<3: continue
    m=c.m; co=c.coords(); X=c.X
    AT=c.idx(co[:,None,:]+co[None,:,:])
    S=c.paf().sum(0); H=np.array(sorted({0}|set(np.where(S==-4*c.b)[0].tolist())))
    cos=-np.ones(m,int); nq=0
    for g in range(m):
        if cos[g]<0: cos[AT[g,H]]=nq; nq+=1
    auts=automorphisms(c)
    labs=[]; nst=[]
    for k in range(4):
        lab,ns=affine_orbits(X[k],c,auts,AT); labs.append(lab); nst.append(ns)
    rec=dict(v=c.v,b=c.b,group=c.mods,m=m,H=len(H),n_aut=len(auts),stab_sizes=nst,orbit_counts=[int(l.max()+1) for l in labs])
    # P: partial symmetries
    real=[partial_stats(X[k],c,auts,cos,nq,len(H)) for k in range(4)]
    nulls=[[partial_stats(random_like(X[k],labs[k],rng),c,auts,cos,nq,len(H)) for _ in range(NN)] for k in range(4)]
    rec['partial']=[dict(block=k,count=real[k][0],max_cosets=real[k][1],
                         null_mean=round(float(np.mean([z[0] for z in nulls[k]])),1),
                         p_count=round(float(np.mean([z[0]>=real[k][0] for z in nulls[k]])),3),
                         p_maxcos=round(float(np.mean([z[1]>=real[k][1] for z in nulls[k]])),3)) for k in range(4)]
    print(json.dumps(rec),flush=True)
    # X: pair relations
    pr=[]
    for i in range(4):
        for j in range(i+1,4):
            a=pair_stat(X[i],X[j],c,auts)
            nl=[pair_stat(random_like(X[i],labs[i],rng),random_like(X[j],labs[j],rng),c,auts) for _ in range(NN)]
            pr.append(dict(pair=[i,j],max_agree=a,null_mean=round(float(np.mean(nl)),1),null_max=int(max(nl)),p=round(float(np.mean([z>=a for z in nl])),3)))
    rec['pairs']=pr
    comp=np.array([[int(X[k][cos==q].sum()) for q in range(nq)] for k in range(4)])
    rec['compressions']=comp.tolist(); rec['quot_dev_by_block']=[float(((comp[k]-comp[k].mean())**2).sum()) for k in range(4)]
    print(json.dumps(rec),flush=True)
