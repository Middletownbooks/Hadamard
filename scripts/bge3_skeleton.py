"""Skeleton of the b>=3 matrices: corner Hadamard matrix of order 4b, border pattern K = 4b rows of a Hadamard
matrix of order 4c (c=b+1), and how the core's coset sums complete it.
 - completions: all +-1 vectors orthogonal to K's rows (the 4 missing rows of an H(4c))
 - quotient matrix Q (order 4b+4c): border/corner entries, and coset-sums of the sign-corrected core
 - is the border symmetric (column pattern = row pattern)?  corner type (Paley?) by invariants"""
import sys, json, itertools, numpy as np
sys.path.insert(0,'scripts')
from extract_group import load, G, solve
from hadlib import Core
cores={1916:(3,(2,2,119),'data/core_v479_b3.txt'),1388:(5,(2,3,57),'data/core_v347_b5.txt'),1436:(7,(2,2,8,11),'data/core_v359_b7.txt')}
def block_signs(M,off,m,g):
    E=[]
    for c in range(4):
        sub=M[off:off+m,off+c*m:off+(c+1)*m]; f=solve(sub,g,c>0)[0]; I=np.arange(m)
        F=f[g.idx(g.coords[I][:,None,:]+g.coords[I][None,:,:])] if c>0 else f[g.idx(g.coords[I][None,:,:]-g.coords[I][:,None,:])]
        T=sub*F; E.append(T[0,:]*T[0,0])
    return np.concatenate(E)
def completions(K):
    n=K.shape[1]; U,s,Vt=np.linalg.svd(K.astype(float)); B=Vt[np.sum(s>1e-8):].T      # n x 4 basis of complement
    out=set()
    for rows in itertools.combinations(range(n),B.shape[1]):
        Bs=B[list(rows)]
        if abs(np.linalg.det(Bs))<1e-6: continue
        for sg in itertools.product([1,-1],repeat=B.shape[1]):
            a=np.linalg.solve(Bs,np.array(sg,float)); x=B@a
            if np.allclose(np.abs(x),1,atol=1e-6): out.add(tuple(np.rint(x).astype(int)))
        break
    return sorted(out)
def hinv(Hm):
    """cheap invariant profile of a Hadamard matrix: distribution of |4-row profiles| (Hall/Kimura style)"""
    n=Hm.shape[0]; cnt={}
    for a,b,c_ in itertools.combinations(range(min(n,16)),3):
        k=abs(int((Hm[0]*Hm[a]*Hm[b]*Hm[c_]).sum())) if 0 not in (a,b,c_) else None
        if k is not None: cnt[k]=cnt.get(k,0)+1
    return dict(sorted(cnt.items()))
for N,(b,mods,cf) in cores.items():
    M=load(N); v=N//4; off=4*b; m=v-b; g=G(mods); c=b+1; h=m//c
    core=Core(cf); Sg=core.paf().sum(0); H=[int(x) for x in np.where(Sg==-4*b)[0]]+[0]
    cos=-np.ones(m,int); nq=0
    for x in range(m):
        if cos[x]<0:
            for y in H: cos[g.add(x,y)]=nq
            nq+=1
    E=block_signs(M,off,m,g); D=block_signs(M.T.copy(),off,m,g)
    Mt=M.copy(); Mt[:,off:]*=E[None,:]; Mt[off:,:]*=D[:,None]
    ok=np.array_equal(Mt@Mt.T,N*np.eye(N,dtype=Mt.dtype))
    corner=Mt[:off,:off]
    cls=np.concatenate([k*nq+cos for k in range(4)])          # coset-class of each core index
    K=np.array([[Mt[p,off:][cls==j][0] for j in range(4*nq)] for p in range(off)])
    Kc=np.array([[Mt[off:,p][cls==j][0] for j in range(4*nq)] for p in range(off)])
    comp=completions(K)
    # quotient: coset sums of the core
    Z=np.zeros((4*nq,4*nq),int)
    for i in range(4*nq):
        r=off+np.where(cls==i)[0][0]
        for j in range(4*nq): Z[i,j]=int(Mt[r,off:][cls==j].sum())
    # relation: does Z (coset sums, divided by h? no: entries are sums over h columns) span the completion space?
    Ccomp=np.array(comp) if comp else np.zeros((0,4*nq))
    proj=None
    if len(comp):
        B=np.array(comp,float).T; P=B@np.linalg.pinv(B)
        proj=float(np.abs(Z.astype(float)-(P@Z.T.astype(float)).T).max())
    rec=dict(v=v,b=b,c=c,h=h,G=mods,matrix_ok=ok,corner_is_hadamard=bool(np.array_equal(corner@corner.T,off*np.eye(off,dtype=int))),
             K_rows_orthogonal=bool(np.array_equal(K@K.T,4*c*np.eye(off,dtype=int))),border_symmetric=[bool(np.array_equal(K,Kc)),bool(np.array_equal(K,-Kc))],
             n_completion_vectors=len(comp),Z_minus_projection_on_completion_space=proj,
             Z_value_set=sorted(set(Z.ravel().tolist())))
    print(json.dumps(rec),flush=True)
    if len(comp)>=4:
        Hfull=np.vstack([K,np.array(comp[:])]) 
        # choose 4 mutually orthogonal completion vectors
        best=None
        for q in itertools.combinations(range(len(comp)),4):
            Cq=np.array([comp[i] for i in q])
            if np.array_equal(Cq@Cq.T,4*c*np.eye(4,dtype=int)): best=Cq; break
        if best is not None:
            W=np.vstack([K,best]); print('  completed H(%d) orthogonal:'%(4*c),bool(np.array_equal(W@W.T,4*c*np.eye(4*c,dtype=int))))
            print('  completion rows:'); [print('   ',''.join('+' if t>0 else '-' for t in r)) for r in best]
    print('  core coset-sum matrix Z (first row of each core row-block; entries = sums over h=%d cosets):'%h)
    for R in range(4): print('   ',Z[R*nq].tolist())
    np.save(f'data/skeleton_{N}.npy',dict(K=K,Kc=Kc,corner=corner,Z=Z,comp=np.array(comp)),allow_pickle=True)
