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

def class_reps(g,u,h):
    """lag-class representatives under x -> u*x and x -> -x (same classes as stage B)."""
    co=g['co']; m=g['m']; NEG=g['NEG']; e=int(np.lcm.reduce([2,3,h])); perm=g['idx'](co*(u%e))
    cl=-np.ones(m,int); reps=[]
    for x in range(1,m):
        if cl[x]<0:
            S=set(); y=x
            while y not in S: S.add(y); S.add(int(NEG[y])); y=int(perm[y])
            for z in S: cl[z]=len(reps)
            reps.append(x)
    return np.array(reps)

def joint_residue_ok(g,comps,reps,xc,h,ns=20000,M=64,maxres=16,seed=0):
    """Joint (multi-class) version of residue_ok.  Per block, the vector of lift PAFs at the lag reps, reduced mod M,
    takes few values on the 'forced' classes (<= maxres residues per block); require some choice of one sampled residue
    vector per block with r0+r1+r2+r3 = T (mod M) on all those classes simultaneously.  Necessary condition, up to
    sampling of the residue vectors (sizes reported)."""
    rng=np.random.default_rng(seed); T=g['T_full'][reps]; R=[]
    for k in range(4):
        F=lift_sampler(g,comps,xc[k],ns,rng)
        if F is None: return False,dict(reason='no lift')
        Fo=F.astype(float).reshape(ns,2,3,h)
        P=np.rint(np.fft.ifftn(np.abs(np.fft.fftn(Fo,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(ns,-1)[:,reps].astype(np.int64)
        R.append(P%M)
    C=[c for c in range(len(reps)) if all(len(np.unique(R[k][:,c]))<=maxres for k in range(4))]
    # full residue vectors are near-unique per sample, so check every PAIR of classes jointly (catches coupled residues)
    import itertools
    def sumset(X,Y): return np.unique(((X[:,None,:]+Y[None,:,:])%M).reshape(-1,X.shape[1]),axis=0)
    bad=[]
    for c1,c2 in itertools.combinations(C,2):
        V=[np.unique(R[k][:,[c1,c2]],axis=0) for k in range(4)]
        AB=sumset(V[0],V[1]); CD=sumset((T[[c1,c2]][None,:]-V[2])%M,-V[3]%M)
        if not (set(map(tuple,AB.tolist()))&set(map(tuple,CD.tolist()))): bad.append((int(reps[c1]),int(reps[c2])))
    info=dict(classes=len(C),pairs=len(C)*(len(C)-1)//2,bad_pairs=len(bad),first_bad=bad[:4])
    return len(bad)==0,info
