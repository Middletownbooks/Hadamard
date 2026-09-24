"""Exact slice relations in the b>=3 cores.  A slice = restriction of a block to one H-coset.
Relation: x_j(g) = s * chi(g) * x_i(sigma(g) + t) for ALL g in coset q of block j, where sigma in Aut(G),
t in G, chi a real (+-1) character of G, s = +-1.  (Maps coset q of block j onto coset sigma(q)+t of block i.)
Counts, per core: distinct (block_i, coset, block_j, coset) pairs related exactly, split into
  within-block not explained by a global symmetry of that block, and between-block;
null = random blocks constant on each block's own affine orbits (same row sums).
usage: bge2_slices.py NNULL [v ...]"""
import sys, json, itertools, numpy as np
sys.path.insert(0,'scripts')
from hadlib import all_cores
src=open('scripts/bge2_forensics2.py').read().split('for c in all_cores():')[0]; exec(src)
NN=int(sys.argv[1]); VS=[int(x) for x in sys.argv[2:]]
rng=np.random.default_rng(5)

def real_chars(c):
    co=c.coords(); mods=np.array(c.mods); out=[]
    even=[i for i,n in enumerate(mods) if n%2==0]
    for bits in itertools.product([0,1],repeat=len(even)):
        e=np.zeros(c.m,int)
        for b,i in zip(bits,even):
            if b: e+=co[:,i]*(1 if mods[i]==2 else 1)      # (-1)^{x_i} is a character iff n_i even
        out.append(np.where(e%2==0,1,-1))
    return out

def slice_relations(Xs,c,auts,chis,cos,nq,Hs):
    """set of (i,qi,j,qj) with an exact relation (any sigma, chi, t, sign)"""
    rel=set(); G=c.m
    FZ=[[np.conj(np.fft.fftn(c.arr((Xs[j]*(cos==q)).astype(float)))) for q in range(nq)] for j in range(4)]
    for P in auts:
        for chi in chis:
            for i in range(4):
                Y=np.fft.fftn(c.arr((Xs[i][P]*chi).astype(float)))   # y(g) = x_i(sigma g) chi(g)
                for j in range(4):
                    for q in range(nq):
                        A=np.rint(np.fft.ifftn(Y*FZ[j][q]).real.reshape(-1))   # sum_{g in q} y(g+t) x_j(g)
                        for t in np.where(np.abs(A)==Hs)[0]:
                            # coset of block i hit: sigma(q + t) ... record the image coset label
                            g0=np.where(cos==q)[0][0]; qi=int(cos[P[AT_[g0,t]]])
                            if (i,qi)!=(j,q): rel.add((i,qi,j,q))
    return rel

out=[]
for c in all_cores():
    if c.b<3 or (VS and c.v not in VS): continue
    m=c.m; co=c.coords(); X=c.X
    AT_=c.idx(co[:,None,:]+co[None,:,:])
    S=c.paf().sum(0); H=np.array(sorted({0}|set(np.where(S==-4*c.b)[0].tolist())))
    cos=-np.ones(m,int); nq=0
    for g in range(m):
        if cos[g]<0: cos[AT_[g,H]]=nq; nq+=1
    auts=automorphisms(c); chis=real_chars(c)
    labs=[affine_orbits(X[k],c,auts,AT_)[0] for k in range(4)]
    def summarize(rel):
        w=sorted({(i,qi,j,q) for (i,qi,j,q) in rel if i==j}); b=sorted({(i,qi,j,q) for (i,qi,j,q) in rel if i!=j})
        return w,b
    rw,rb=summarize(slice_relations(X,c,auts,chis,cos,nq,len(H)))
    nw=[];nb=[]
    for _ in range(NN):
        Y=np.array([random_like(X[k],labs[k],rng) for k in range(4)])
        w,b=summarize(slice_relations(Y,c,auts,chis,cos,nq,len(H))); nw.append(len(w)); nb.append(len(b))
    rec=dict(v=c.v,b=c.b,group=c.mods,nq=nq,n_aut=len(auts),n_real_chars=len(chis),
             within=len(rw),within_null=nw,between=len(rb),between_null=nb,
             between_pairs=rb[:40],within_pairs=rw[:40])
    print(json.dumps(rec),flush=True)
