"""Z2^2-component decomposition of the b=3 core v=479 (G = Z2^2 x Z119, H = Z119).
xhat_{k,chi}(z) = sum_{a,b} chi(a,b) x_k(a,b,z); for each chi the 4 blocks' components are complementary on Z119:
sum_k |FT xhat_{k,chi}(w)|^2 = 4v (w != 0).  Test exact relations xhat_{i,chi}(mu z + t) = s xhat_{j,chi'}(z)
between components (mu in (Z/119)*, t in Z119, s = +-1), within and between blocks; null = random blocks with
each block's own affine-orbit structure (same row sums)."""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from hadlib import all_cores
src=open('scripts/bge2_forensics2.py').read().split('for c in all_cores():')[0]; exec(src)
from math import gcd
NN=int(sys.argv[1]) if len(sys.argv)>1 else 20
c=[x for x in all_cores() if x.v==479][0]; co=c.coords(); n=119
chis=[(0,0),(0,1),(1,0),(1,1)]
def comps(x):
    out=np.zeros((4,n),int)
    for ci,(p,q) in enumerate(chis):
        s=(-1)**((p*co[:,0]+q*co[:,1])%2); np.add.at(out[ci],co[:,2],s*x)
    return out
units=[u for u in range(1,n) if gcd(u,n)==1]
def rel(C):
    """list of (i,ci,j,cj) with exact affine relation (excluding identical index)"""
    R=[]; FT={}
    for j in range(4):
        for cj in range(4): FT[(j,cj)]=np.conj(np.fft.fft(C[j][cj].astype(float)))
    for i in range(4):
        for ci in range(4):
            for u in units:
                y=C[i][ci][(np.arange(n)*u)%n].astype(float); Y=np.fft.fft(y); nrm=(y*y).sum()
                for j in range(4):
                    for cj in range(4):
                        if (i,ci)==(j,cj): continue
                        if (C[j][cj]**2).sum()!=nrm: continue
                        A=np.rint(np.fft.ifft(Y*FT[(j,cj)]).real)
                        if np.abs(A).max()==nrm: R.append((i,ci,j,cj))
    return sorted(set(R))
C=[comps(x) for x in c.X]
P=[[np.abs(np.fft.fft(C[k][ci]))**2 for k in range(4)] for ci in range(4)]
print('complementarity check: per chi, min/max of sum_k P over w!=0:',[(round(float(sum(P[ci])[1:].min()),3),round(float(sum(P[ci])[1:].max()),3)) for ci in range(4)])
print('component norms',[[int((C[k][ci]**2).sum()) for ci in range(4)] for k in range(4)])
print('component value sets',[[sorted(set(C[k][ci].tolist())) for ci in range(4)] for k in range(4)])
R=rel(C); w=[r for r in R if r[0]==r[2]]; b=[r for r in R if r[0]!=r[2]]
print('real: within-block',len(w),w[:12],'| between-block',len(b),b[:12],flush=True)
auts=automorphisms(c); AT=c.idx(co[:,None,:]+co[None,:,:])
labs=[affine_orbits(c.X[k],c,auts,AT)[0] for k in range(4)]
rng=np.random.default_rng(7); nw=[];nb=[]
for _ in range(NN):
    Y=[comps(random_like(c.X[k],labs[k],rng)) for k in range(4)]
    R=rel(Y); nw.append(sum(1 for r in R if r[0]==r[2])); nb.append(sum(1 for r in R if r[0]!=r[2]))
print('null within',nw,' between',nb)
