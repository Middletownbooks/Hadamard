"""Forensics on the published b>=3 cores (479, 347) with respect to the CRT split H = Z_h1 x Z_h2.
(1) compressions onto Q x Z_h1 and Q x Z_h2: entry variance vs a permuted-block null;
(2) rank test: top-singular-value energy share of each (block, coset) slice as an h1 x h2 matrix, against
    (a) i.i.d. signs and (b) the symmetry-aware null (random blocks on the same affine orbits, same row sum).
Result [M]: (1) within 0.8-1.7 sigma of null; (2) 15-16 sigma vs i.i.d., but z = -0.2 / +0.8 vs the symmetry-aware null:
the apparent low-rank structure is entirely the multiplier symmetry.
usage: bge3_forensics3.py [NNULL]"""
import sys, numpy as np
sys.path.insert(0,'scripts')
src=open('scripts/bge2_forensics2.py').read().split('for c in all_cores():')[0]; exec(src)
from hadlib import Core
NN=int(sys.argv[1]) if len(sys.argv)>1 else 200
def crt(X,Qd,h1,h2):
    nq=int(np.prod(Qd)); Xr=X.reshape(4,nq,h1*h2); M=np.zeros((4,nq,h1,h2))
    for z in range(h1*h2): M[:,:,z%h1,z%h2]=Xr[:,:,z]
    return M
def share(M):
    s=np.linalg.svd(M,compute_uv=False); return float((s[...,0]**2/(s**2).sum(-1)).mean())
def compvar(M,axis):
    Z=M.sum(axis).reshape(4,-1); return float(((Z-Z.mean(1,keepdims=True))**2).mean())
rng=np.random.default_rng(3)
for path,Qd,h1,h2 in [('data/core_v479_b3.txt',(2,2),7,17),('data/core_v347_b5.txt',(2,3),3,19)]:
    c=Core(path); co=c.coords(); AT=c.idx(co[:,None,:]+co[None,:,:]); auts=automorphisms(c)
    labs=[affine_orbits(c.X[k],c,auts,AT)[0] for k in range(4)]
    M=crt(c.X,Qd,h1,h2)
    real=(compvar(M,-1),compvar(M,-2),share(M))
    perm=np.array([[compvar(Mr,-1),compvar(Mr,-2),share(Mr)] for Mr in
                   (crt(np.array([rng.permutation(c.X[k]) for k in range(4)]),Qd,h1,h2) for _ in range(NN))])
    sym=np.array([[compvar(Mr,-1),compvar(Mr,-2),share(Mr)] for Mr in
                  (crt(np.array([random_like(c.X[k],labs[k],rng) for k in range(4)]),Qd,h1,h2) for _ in range(NN))])
    for i,nm in enumerate(['var A (onto Q x Z%d)'%h1,'var B (onto Q x Z%d)'%h2,'top-SV share of h1 x h2 slices']):
        print('%s %-34s real %.3f | permuted null %.3f+-%.3f (z=%.1f) | symmetry-aware null %.3f+-%.3f (z=%.1f)'%(
            path.split('/')[-1],nm,real[i],perm[:,i].mean(),perm[:,i].std(),(real[i]-perm[:,i].mean())/perm[:,i].std(),
            sym[:,i].mean(),sym[:,i].std(),(real[i]-sym[:,i].mean())/sym[:,i].std()),flush=True)
