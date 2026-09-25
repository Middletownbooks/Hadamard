"""Exact stage B for the c=6 architecture: given a compressed solution xc (level t = h/p) and symmetry u, enumerate
ALL +-1 S-invariant blocks whose K-fibre sums equal xc[k] (a product of independent choices), filter by the PSD
bound (P(w) <= 4v for w nontrivial on H; the H-trivial part is fixed by the compression), and save per-block lists
of lag-class PAF vectors for a 4-list join.   lib: enum_lifts(A, target) -> array of orbit-sign vectors."""
import numpy as np
def enum_lifts(A,target,limit=5*10**7):
    """A: (fibres, orbits) footprint matrix; all y in {+-1}^orbits with A y = target.  Backtracking with
    capacity pruning; returns int8 array (n, orbits)."""
    order=np.argsort(-A.sum(0)); Ao=A[:,order]
    rem=np.concatenate([np.cumsum(Ao[:,::-1],axis=1)[:,::-1],np.zeros((A.shape[0],1),int)],axis=1)
    out=[]; stack=[(0,np.zeros(A.shape[0],int),[])]
    while stack:
        i,cur,assign=stack.pop()
        if i==Ao.shape[1]:
            if (cur==target).all(): out.append(assign)
            if len(out)>=limit: break
            continue
        for s in (1,-1):
            new=cur+s*Ao[:,i]
            if (np.abs(target-new)<=rem[:,i+1]).all() and ((target-new-rem[:,i+1])%2==0).all():
                stack.append((i+1,new,assign+[s]))
    Y=np.array(out,np.int8).reshape(-1,A.shape[1]); Yo=np.empty_like(Y); Yo[:,order]=Y
    return Yo
