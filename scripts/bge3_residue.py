"""Forced-residue (lift-level) feasibility filter for compressed solutions of the c=6 architecture.
For a compressed solution xc, the PAF of every lift of block k at certain lag classes is fixed modulo a power of 2
(orbit repetition x fibre-sum parity).  Sample lifts per block, find per class the largest 2^j (<= 2^JMAX) at which the
residue is constant, and require sum_k r_k = T (mod min_k 2^j_k) at every class.  Necessary condition for a lift."""
import numpy as np
def lift_sampler(g,comps,xck,ns,rng):
    A=g['footB'].T; opts=[]
    for cq in comps:
        fis=np.nonzero(A[:,cq].sum(1))[0]; sub=A[np.ix_(fis,cq)]; n=len(cq)
        Y=((np.arange(1<<n)[:,None]>>np.arange(n)[None,:])&1)*2-1
        ok=Y[(Y@sub.T==xck[fis][None,:]).all(1)]
        if len(ok)==0: return None
        opts.append((cq,ok))
    Ys=np.empty((ns,A.shape[1]),int)
    for cq,ok in opts: Ys[:,cq]=ok[rng.integers(0,len(ok),ns)]
    return Ys[:,g['lab']]
def residue_ok(g,comps,reps,xc,ns=48,jmax=6,seed=0):
    rng=np.random.default_rng(seed); ADD=g['ADD']; T=g['T_full']; R=[]; M=[]
    for k in range(4):
        F=lift_sampler(g,comps,xc[k],ns,rng)
        if F is None: return False,'no lift'
        P=np.array([[int(np.dot(f,f[ADD[:,r]])) for r in reps] for f in F])
        mk=np.ones(len(reps),np.int64); rk=np.zeros(len(reps),np.int64)
        for ci in range(len(reps)):
            for j in range(jmax,1,-1):
                if len(set((P[:,ci]%(2**j)).tolist()))==1: mk[ci]=2**j; rk[ci]=P[0,ci]%(2**j); break
        R.append(rk); M.append(mk)
    Mmin=np.minimum.reduce(M); tot=sum(R)
    bad=[int(reps[ci]) for ci in range(len(reps)) if Mmin[ci]>1 and (tot[ci]-T[reps[ci]])%Mmin[ci]!=0]
    return len(bad)==0, bad
