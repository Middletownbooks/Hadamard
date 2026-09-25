"""Extract a reusable template from the (sign-corrected) v=347 matrix: corner C (20x20), border patterns K, Kc per
(core block, coset (a,b)), and the layout of the 16 core sub-blocks: M_IJ[g,h] = rho(g) gam(h) s F_k(sigma(g,h)),
sigma in {h-g, g+h, g-h, -g-h}, with F_k the row-0 functions of row-block 0.  rho, gam are checked to be real
characters of Z2 x Z3 x Z57 (so the template transfers to Z2 x Z3 x Z_h).  Saves data/template347.npz."""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from extract_group import load, G, solve
N=1388; b=5; mods=(2,3,57); off=4*b; v=N//4; m=v-b; g=G(mods); co=g.coords
exec(open('scripts/bge3_skeleton.py').read().split('def completions')[0].split("cores={")[0])
def block_signs(M,off,m,g):
    E=[]
    for c in range(4):
        sub=M[off:off+m,off+c*m:off+(c+1)*m]; f=solve(sub,g,c>0)[0]; I=np.arange(m)
        F=f[g.idx(g.coords[I][:,None,:]+g.coords[I][None,:,:])] if c>0 else f[g.idx(g.coords[I][None,:,:]-g.coords[I][:,None,:])]
        T=sub*F; E.append(T[0,:]*T[0,0])
    return np.concatenate(E)
M=load(N); E=block_signs(M,off,m,g); D=block_signs(M.T.copy(),off,m,g)
Mt=M.copy(); Mt[:,off:]*=E[None,:]; Mt[off:,:]*=D[:,None]
assert np.array_equal(Mt@Mt.T,N*np.eye(N,dtype=Mt.dtype))
core=Mt[off:,off:]; F=[core[0,J*m:(J+1)*m].copy() for J in range(4)]
I_=np.arange(m); Gc=co
sig={'h-g':g.idx(Gc[None,:,:]-Gc[:,None,:]),'g+h':g.idx(Gc[:,None,:]+Gc[None,:,:]),
     'g-h':g.idx(Gc[:,None,:]-Gc[None,:,:]),'-g-h':g.idx(-Gc[:,None,:]-Gc[None,:,:])}
chars={'1':np.ones(m,int),'a':np.where(Gc[:,0]%2==0,1,-1)}
layout={}
for Ib in range(4):
    for Jb in range(4):
        S=core[Ib*m:(Ib+1)*m,Jb*m:(Jb+1)*m]; found=None
        for k in range(4):
            for sn,ix in sig.items():
                R=S*F[k][ix]
                rho=R[:,0]*R[0,0]; gam=R[0,:]
                if np.array_equal(R,np.outer(rho,gam)):
                    cr=[n for n,ch in chars.items() if np.array_equal(rho,ch) or np.array_equal(rho,-ch)]
                    cg=[n for n,ch in chars.items() if np.array_equal(gam*R[0,0],ch) or np.array_equal(gam,ch) or np.array_equal(gam,-ch)]
                    if cr and cg:
                        s=int(rho[0]*gam[0]) if False else 1
                        rs=int(np.sign(rho[0])); gs=int(np.sign(gam[0]))
                        found=dict(k=k,sigma=sn,rho=cr[0],rho_sign=int(rho[0]*chars[cr[0]][0]),gam=cg[0],gam_sign=int(gam[0]*chars[cg[0]][0]))
                        break
            if found: break
        layout[f'{Ib}{Jb}']=found
print(json.dumps(layout))
assert all(layout.values()), 'some sub-block did not fit'
# border patterns per class (block, coset a,b) and corner
cos=(co[:,0]*3+co[:,1])                                   # coset index 0..5 = (a,b)
K=np.array([[Mt[p,off+J*m:off+(J+1)*m][cos==q][0] for J in range(4) for q in range(6)] for p in range(off)])
Kc=np.array([[Mt[off+J*m:off+(J+1)*m,p][cos==q][0] for J in range(4) for q in range(6)] for p in range(off)])
chk=all(len(set(Mt[p,off+J*m:off+(J+1)*m][cos==q].tolist()))==1 for p in range(off) for J in range(4) for q in range(6))
chk2=all(len(set(Mt[off+J*m:off+(J+1)*m,p][cos==q].tolist()))==1 for p in range(off) for J in range(4) for q in range(6))
comp=np.array([[int(F[J][cos==q].sum()) for q in range(6)] for J in range(4)])
print('border constant on cosets:',chk,chk2,' compressions of row-block-0 functions:',comp.tolist(),' row sums',comp.sum(1).tolist())
np.savez('data/template347.npz',C=Mt[:off,:off],K=K,Kc=Kc,comp=comp,layout=json.dumps(layout),F=np.array(F))
