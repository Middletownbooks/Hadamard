"""Canonicalise the border (remove the real-character twist u so the completion rows are H4 (x) 1_c) and test
whether K = H4 (x) L (Kronecker: each K row, reshaped 4 x c, has rank 1) and whether Kc ~ K."""
import sys, json, itertools, numpy as np
for N,v,c in [(1916,479,4),(1388,347,6),(1436,359,8)]:
    S=np.load(f'data/skeleton_{N}.npy',allow_pickle=True).item(); K=S['K']; Kc=S['Kc']; comp=S['comp']; Z=S['Z']
    # u = coset pattern of the first completion vector's first block (up to sign)
    best=None
    for q in itertools.combinations(range(len(comp)),4):
        Cq=comp[list(q)]
        if np.array_equal(Cq@Cq.T,4*c*np.eye(4,dtype=int)): best=Cq; break
    u=best[0][:c]*best[0][0]
    U=np.tile(u,4)
    Kp=K*U[None,:]; Kcp=Kc*U[None,:]; Cp=best*U[None,:]
    blk_const=all(len(set(r[j*c:(j+1)*c]))==1 for r in Cp for j in range(4))
    ranks=[int(np.linalg.matrix_rank(r.reshape(4,c).astype(float))) for r in Kp]
    bal=[int(np.abs(r.reshape(4,c).sum(1)).max()) for r in Kp]
    # row-space comparison K vs Kc
    same_space=np.linalg.matrix_rank(np.vstack([Kp,Kcp]).astype(float))
    print(json.dumps(dict(v=v,c=c,u=u.tolist(),completion_block_constant_after_untwist=blk_const,
        K_row_ranks_as_4xc=sorted(set(ranks)),count_rank1=ranks.count(1),max_block_sum=max(bal),
        rank_K=int(np.linalg.matrix_rank(Kp.astype(float))),rank_K_plus_Kc=int(same_space))))
    print('  K (untwisted), rows as [A|B|C|D] cosets:')
    for r in Kp: print('   ',' '.join(''.join('+' if t>0 else '-' for t in r[j*c:(j+1)*c]) for j in range(4)))
    Zp=Z*U[None,:]*U[:,None]
    print('  core coset sums (untwisted), first row of each row-block:')
    for R in range(4): print('   ',[Zp[R*c][j*c:(j+1)*c].tolist() for j in range(4)])
