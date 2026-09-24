"""Border structure of the b>=3 matrices.  Undo the normalisation on the core columns (and, via the transpose,
the core rows) using the rank-1 sign pattern of each group-developed block, then read every border row
(resp. column) across the 4 core blocks: is it constant on H-cosets?  If so the border is a 4b x 4c sign matrix
(c = b+1 cosets per block).  Also reports the corner 4b x 4b block's Gram structure."""
import sys, json, itertools, numpy as np
sys.path.insert(0,'scripts')
from extract_group import load, G, solve
from hadlib import Core
cores={1388:(5,(2,3,57),'data/core_v347_b5.txt'),1436:(7,(2,2,8,11),'data/core_v359_b7.txt'),1916:(3,(2,2,119),'data/core_v479_b3.txt')}

def block_signs(M,off,m,g):
    """per core column block c: column sign vector E_c (from row block 0), up to a global sign"""
    E=[]
    for c in range(4):
        sub=M[off:off+m,off+c*m:off+(c+1)*m]; fs=solve(sub,g,c>0); f=fs[0]
        I=np.arange(m)
        F=f[g.idx(g.coords[I][:,None,:]+g.coords[I][None,:,:])] if c>0 else f[g.idx(g.coords[I][None,:,:]-g.coords[I][:,None,:])]
        T=sub*F; E.append(T[0,:]*T[0,0])
    return E

def border_pattern(M,b,mods,H):
    m=int(np.prod(mods)); off=4*b; g=G(mods); E=block_signs(M,off,m,g)
    cos=-np.ones(m,int); nq=0
    for x in range(m):
        if cos[x]<0:
            for h in H: cos[g.add(x,h)]=nq
            nq+=1
    S=np.zeros((off,4*nq),int); const=True
    for p in range(off):
        for c in range(4):
            w=M[p,off+c*m:off+(c+1)*m]*E[c]
            for q in range(nq):
                vals=set(w[cos==q].tolist())
                if len(vals)!=1: const=False; S[p,c*nq+q]=0
                else: S[p,c*nq+q]=vals.pop()
    return S,const,nq,cos

out={}
for N,(b,mods,cf) in cores.items():
    M=load(N); core=Core(cf); Sg=core.paf().sum(0); H=[int(h) for h in np.where(Sg==-4*b)[0]]+[0]
    Sr,cr,nq,cos=border_pattern(M,b,mods,H)          # border rows x core columns
    Sc,cc,_,_=border_pattern(M.T.copy(),b,mods,H)     # border columns x core rows
    corner=M[:4*b,:4*b]
    rec=dict(N=N,v=N//4,b=b,group=mods,c=nq,rows_constant_on_cosets=cr,cols_constant_on_cosets=cc)
    # Gram of border rows over core part, and over the corner
    Grow=Sr@Sr.T; Gcol=Sc@Sc.T
    rec['row_pattern_gram_values']=sorted(set(np.round(Grow[np.triu_indices(4*b,1)]).astype(int).tolist()))
    rec['corner_row_gram_values']=sorted(set((corner@corner.T)[np.triu_indices(4*b,1)].tolist()))
    # per core block: the 4b x c sign pattern; its rank and column-sum structure
    blocks=[Sr[:,k*nq:(k+1)*nq] for k in range(4)]
    rec['row_pattern_rank_per_block']=[int(np.linalg.matrix_rank(B_)) for B_ in blocks]
    rec['row_pattern_rank_all']=int(np.linalg.matrix_rank(Sr))
    rec['row_pattern_distinct_columns']=int(len({tuple(col) for col in Sr.T}))
    rec['col_pattern_distinct_columns']=int(len({tuple(col) for col in Sc.T}))
    print(json.dumps(rec),flush=True)
    np.save(f'data/border_rows_{N}.npy',Sr); np.save(f'data/border_cols_{N}.npy',Sc); np.save(f'data/border_corner_{N}.npy',corner)
    print('border-row pattern (rows = border rows, cols = [block A cosets | B | C | D]):')
    for row in Sr: print('  ',''.join('+' if t>0 else ('-' if t<0 else '?') for t in row))
