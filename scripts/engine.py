"""W2 compiled local-search engine for b=0 reduced (U-orbit) GS cores.
State: y[4,r] orbit signs.  Orbit-level PAF on L lag classes: PAF_k(c)=1/2 y^T M_c y, M_c=N_c+N_c^T.
Flip q in block k: dPAF_c = -2 y_q (M_c y)_q + 2 M_c[q,q].
G[k,c,q] = (M_c y_k)_q maintained incrementally (O(L r) per flip).
Modes
  S3 joint:   E = sum_c R_c^2 + lam * sum_k sum_c PAF_k(c)^2        (R = sum_k PAF_k)
  S2 shared:  block k (in turn) minimises sum_c (PAF_k(c) - t_k(c))^2 + mu*sum_c PAF_k(c)^2,
              t_k = PAF_k - R/4 frozen at the start of block k's turn; `inner` flips per turn.
"""
import numpy as np, numba as nb, time, sys
sys.path.insert(0,'scripts') if 'scripts' not in sys.path else None
from redlib import Red

def build(R):
    m,r,L=R.m,R.r,R.L
    M=np.zeros((L,r,r),np.int32)
    orb=R.orb
    for c,s in enumerate(R.lags):
        i=np.arange(m); j=(i+s)%m
        np.add.at(M[c],(orb[i],orb[j]),1)
        M[c]+=M[c].T.copy() if False else 0
    M=M+M.transpose(0,2,1)
    return M

@nb.njit(cache=True)
def init_state(y,M):
    K,r=y.shape; L=M.shape[0]
    G=np.zeros((K,L,r),np.int64); P=np.zeros((K,L),np.int64)
    for k in range(K):
        for c in range(L):
            for q in range(r):
                s=0
                for p in range(r): s+=M[c,q,p]*y[k,p]
                G[k,c,q]=s
            t=0
            for q in range(r): t+=y[k,q]*G[k,c,q]
            P[k,c]=t//2
    return G,P

@nb.njit(cache=True)
def flip(y,G,P,M,k,q):
    L=M.shape[0]; r=y.shape[1]; yq=y[k,q]
    for c in range(L):
        P[k,c]+= -2*yq*G[k,c,q]+2*M[c,q,q]
    for c in range(L):
        for p in range(r):
            G[k,c,p]-=2*yq*M[c,p,q]
    y[k,q]=-yq

@nb.njit(cache=True)
def run_joint(y,M,lam,maxsteps,tenure,seed,noise):
    """tabu search on E=sum R^2+lam*sum PAF^2 ; returns (steps, solved)"""
    np.random.seed(seed)
    K,r=y.shape; L=M.shape[0]
    G,P=init_state(y,M)
    R=np.zeros(L,np.int64)
    for c in range(L):
        for k in range(K): R[c]+=P[k,c]
    tabu=np.zeros((K,r),np.int64)
    d=np.zeros(L,np.int64)
    for step in range(maxsteps):
        rr=0
        for c in range(L): rr+=R[c]*R[c]
        if rr==0: return step,True
        best=1<<62; bk=-1; bq=-1; nb_=0
        for k in range(K):
            for q in range(r):
                if tabu[k,q]>step: continue
                yq=y[k,q]; dE=0
                for c in range(L):
                    dd=-2*yq*G[k,c,q]+2*M[c,q,q]
                    dE+=2*R[c]*dd+dd*dd+lam*(2*P[k,c]*dd+dd*dd)
                if dE<best or (dE==best and np.random.random()<1.0/(nb_+1)):
                    if dE<best: nb_=1
                    else: nb_+=1
                    best=dE; bk=k; bq=q
        if noise>0 and np.random.random()<noise:
            bk=np.random.randint(K); bq=np.random.randint(r)
        yq=y[bk,bq]
        for c in range(L):
            R[c]+=-2*yq*G[bk,c,bq]+2*M[c,bq,bq]
        flip(y,G,P,M,bk,bq)
        tabu[bk,bq]=step+1+tenure+np.random.randint(max(1,tenure//2+1))
    return maxsteps,False

@nb.njit(cache=True)
def run_shared(y,M,mu,maxsteps,inner,tenure,seed,gauss_seidel):
    """S2: residual-sharing block coordinate descent (tabu inside each block's turn)."""
    np.random.seed(seed)
    K,r=y.shape; L=M.shape[0]
    G,P=init_state(y,M)
    R=np.zeros(L,np.int64)
    for c in range(L):
        for k in range(K): R[c]+=P[k,c]
    tabu=np.zeros((K,r),np.int64)
    T=np.zeros(L,np.float64); R0=np.zeros(L,np.int64)
    step=0
    while step<maxsteps:
        for c in range(L): R0[c]=R[c]          # Jacobi: residual frozen for the sweep
        for kk in range(K):
            k=kk
            rr=0
            for c in range(L): rr+=R[c]*R[c]
            if rr==0: return step,True
            for c in range(L):
                T[c]=P[k,c]-(R[c] if gauss_seidel else R0[c])/4.0
            for it in range(inner):
                best=1e300; bq=-1; nb_=0
                for q in range(r):
                    if tabu[k,q]>step: continue
                    yq=y[k,q]; dE=0.0
                    for c in range(L):
                        dd=-2*yq*G[k,c,q]+2*M[c,q,q]
                        e0=P[k,c]-T[c]
                        dE+=2*e0*dd+dd*dd+mu*(2*P[k,c]*dd+dd*dd)
                    if dE<best or (dE==best and np.random.random()<1.0/(nb_+1)):
                        if dE<best: nb_=1
                        else: nb_+=1
                        best=dE; bq=q
                if bq<0: break
                if best>0 and it>0: break      # block has reached a local optimum of its sub-objective
                yq=y[k,bq]
                for c in range(L): R[c]+=-2*yq*G[k,c,bq]+2*M[c,bq,bq]
                flip(y,G,P,M,k,bq)
                tabu[k,bq]=step+1+tenure; step+=1
                rr=0
                for c in range(L): rr+=R[c]*R[c]
                if rr==0: return step,True
    return step,False

def random_y(R,rng): return rng.choice([-1,1],size=(4,R.r)).astype(np.int64)

def fingerprint(R,y):
    X=R.expand(y); P=R.psd(X); S=P[:,1:]/(4*R.m)
    pa=R.paf(X)[:,1:]
    return float(np.mean((S-0.25)**2)), float(np.mean(pa.std(1))/np.sqrt(R.m))

@nb.njit(cache=True)
def run_joint_sched(y,M,lam0,tau,lam_min,maxsteps,tenure,seed):
    """tabu on E=sum R^2 + lam(t) sum PAF^2, lam(t)=max(lam_min, lam0*exp(-t/tau))"""
    np.random.seed(seed)
    K,r=y.shape; L=M.shape[0]
    G,P=init_state(y,M)
    R=np.zeros(L,np.int64)
    for c in range(L):
        for k in range(K): R[c]+=P[k,c]
    tabu=np.zeros((K,r),np.int64)
    for step in range(maxsteps):
        rr=0
        for c in range(L): rr+=R[c]*R[c]
        if rr==0: return step,True
        lam=max(lam_min,lam0*np.exp(-step/tau))
        best=1e300; bk=-1; bq=-1; nb_=0
        for k in range(K):
            for q in range(r):
                if tabu[k,q]>step: continue
                yq=y[k,q]; dE=0.0
                for c in range(L):
                    dd=-2*yq*G[k,c,q]+2*M[c,q,q]
                    dE+=2*R[c]*dd+dd*dd+lam*(2*P[k,c]*dd+dd*dd)
                if dE<best or (dE==best and np.random.random()<1.0/(nb_+1)):
                    if dE<best: nb_=1
                    else: nb_+=1
                    best=dE; bk=k; bq=q
        yq=y[bk,bq]
        for c in range(L):
            R[c]+=-2*yq*G[bk,c,bq]+2*M[c,bq,bq]
        flip(y,G,P,M,bk,bq)
        tabu[bk,bq]=step+1+tenure+np.random.randint(max(1,tenure//2+1))
    return maxsteps,False

@nb.njit(cache=True)
def run_hinge(y,M,lam,F0,maxsteps,tenure,seed):
    """tabu on E = sum_c R_c^2 + lam * sum_k max(0, F_k - F0), F_k = sum_c PAF_k(c)^2.
    Zero-energy states = solutions whose blocks all satisfy F_k <= F0."""
    np.random.seed(seed)
    K,r=y.shape; L=M.shape[0]
    G,P=init_state(y,M)
    R=np.zeros(L,np.int64); F=np.zeros(K,np.int64)
    for c in range(L):
        for k in range(K): R[c]+=P[k,c]; F[k]+=P[k,c]*P[k,c]
    tabu=np.zeros((K,r),np.int64)
    bestE=1e300
    for step in range(maxsteps):
        rr=0
        for c in range(L): rr+=R[c]*R[c]
        if rr==0:
            ok=True
            for k in range(K):
                if F[k]>F0: ok=False
            if ok: return step,True
        best=1e300; bk=-1; bq=-1; nb_=0
        for k in range(K):
            pen0=max(0,F[k]-F0)
            for q in range(r):
                if tabu[k,q]>step: continue
                yq=y[k,q]; dR=0; dF=0
                for c in range(L):
                    dd=-2*yq*G[k,c,q]+2*M[c,q,q]
                    dR+=2*R[c]*dd+dd*dd; dF+=2*P[k,c]*dd+dd*dd
                dE=dR+lam*(max(0,F[k]+dF-F0)-pen0)
                if dE<best or (dE==best and np.random.random()<1.0/(nb_+1)):
                    if dE<best: nb_=1
                    else: nb_+=1
                    best=dE; bk=k; bq=q
        yq=y[bk,bq]
        for c in range(L):
            dd=-2*yq*G[bk,c,bq]+2*M[c,bq,bq]
            R[c]+=dd; F[bk]+=2*P[bk,c]*dd+dd*dd
        flip(y,G,P,M,bk,bq)
        tabu[bk,bq]=step+1+tenure+np.random.randint(max(1,tenure//2+1))
    return maxsteps,False

@nb.njit(cache=True)
def flat_tabu(y,M,maxsteps,tenure,seed,rowsum_fixed):
    """single block: minimise F = sum_c PAF(c)^2 by tabu; flips (rowsum_fixed=0) or swaps (1).
    returns best F and best y (y overwritten with best)."""
    np.random.seed(seed)
    r=y.shape[1]; L=M.shape[0]
    G,P=init_state(y,M)
    F=0
    for c in range(L): F+=P[0,c]*P[0,c]
    bestF=F; besty=y.copy()
    tabu=np.zeros(r,np.int64)
    dd1=np.zeros((r,L),np.int64)
    for step in range(maxsteps):
        if rowsum_fixed==0:
            best=1<<62; bq=-1; nb_=0
            for q in range(r):
                if tabu[q]>step: continue
                yq=y[0,q]; dF=0
                for c in range(L):
                    dd=-2*yq*G[0,c,q]+2*M[c,q,q]; dF+=2*P[0,c]*dd+dd*dd
                if F+dF<bestF: best=dF; bq=q; break      # aspiration
                if dF<best or (dF==best and np.random.random()<1.0/(nb_+1)):
                    if dF<best: nb_=1
                    else: nb_+=1
                    best=dF; bq=q
            yq=y[0,bq]
            for c in range(L):
                dd=-2*yq*G[0,c,bq]+2*M[c,bq,bq]; F+=2*P[0,c]*dd+dd*dd
            flip(y,G,P,M,0,bq); tabu[bq]=step+1+tenure+np.random.randint(tenure//2+1)
        else:
            # swap: flip a (+1 orbit) and a (-1 orbit) of equal size (orbit 0 excluded)
            for q in range(r):
                yq=y[0,q]
                for c in range(L): dd1[q,c]=-2*yq*G[0,c,q]+2*M[c,q,q]
            best=1<<62; ba=-1; bb=-1; nb_=0
            for a in range(1,r):
                if y[0,a]!=1 or tabu[a]>step: continue
                for b in range(1,r):
                    if y[0,b]!=-1 or tabu[b]>step: continue
                    dF=0
                    for c in range(L):
                        dd=dd1[a,c]+dd1[b,c]+4*y[0,a]*y[0,b]*M[c,a,b]
                        dF+=2*P[0,c]*dd+dd*dd
                    if dF<best or (dF==best and np.random.random()<1.0/(nb_+1)):
                        if dF<best: nb_=1
                        else: nb_+=1
                        best=dF; ba=a; bb=b
            if ba<0: continue
            for q in (ba,bb):
                yq=y[0,q]
                for c in range(L):
                    dd=-2*yq*G[0,c,q]+2*M[c,q,q]; F+=2*P[0,c]*dd+dd*dd
                flip(y,G,P,M,0,q); tabu[q]=step+1+tenure+np.random.randint(tenure//2+1)
        if F<bestF:
            bestF=F; besty[:]=y[:]
    y[:]=besty[:]
    return bestF

@nb.njit(cache=True)
def flat_pool(y,M,F0,nsteps,tenure,seed,stride,out,maxout):
    """tabu walk minimising F=sum_c PAF(c)^2; every `stride` steps, if F<=F0, store y into out.
    Restarts with random perturbation when stuck.  returns number stored."""
    np.random.seed(seed)
    r=y.shape[1]; L=M.shape[0]
    G,P=init_state(y,M)
    F=0
    for c in range(L): F+=P[0,c]*P[0,c]
    tabu=np.zeros(r,np.int64); n=0; last=-stride
    for step in range(nsteps):
        best=1<<62; bq=-1; nb_=0
        for q in range(r):
            if tabu[q]>step: continue
            yq=y[0,q]; dF=0
            for c in range(L):
                dd=-2*yq*G[0,c,q]+2*M[c,q,q]; dF+=2*P[0,c]*dd+dd*dd
            if dF<best or (dF==best and np.random.random()<1.0/(nb_+1)):
                if dF<best: nb_=1
                else: nb_+=1
                best=dF; bq=q
        yq=y[0,bq]
        for c in range(L):
            dd=-2*yq*G[0,c,bq]+2*M[c,bq,bq]; F+=2*P[0,c]*dd+dd*dd
        flip(y,G,P,M,0,bq); tabu[bq]=step+1+tenure+np.random.randint(tenure//2+1)
        if F<=F0 and step-last>=stride:
            out[n,:]=y[0,:]; n+=1; last=step
            if n>=maxout: return n
    return n

@nb.njit(cache=True)
def run_breakout(y,M,lam,maxsteps,tenure,seed,winc):
    """tabu + breakout weights: E = sum_c w_c R_c^2 + lam sum_k sum_c PAF_k(c)^2.
    At a local minimum (no improving non-tabu move) weights of unsatisfied classes grow by winc."""
    np.random.seed(seed)
    K,r=y.shape; L=M.shape[0]
    G,P=init_state(y,M)
    R=np.zeros(L,np.int64)
    for c in range(L):
        for k in range(K): R[c]+=P[k,c]
    w=np.ones(L); tabu=np.zeros((K,r),np.int64)
    for step in range(maxsteps):
        rr=0
        for c in range(L): rr+=R[c]*R[c]
        if rr==0: return step,True
        best=1e300; bk=-1; bq=-1; nb_=0
        for k in range(K):
            for q in range(r):
                if tabu[k,q]>step: continue
                yq=y[k,q]; dE=0.0
                for c in range(L):
                    dd=-2*yq*G[k,c,q]+2*M[c,q,q]
                    dE+=w[c]*(2*R[c]*dd+dd*dd)+lam*(2*P[k,c]*dd+dd*dd)
                if dE<best or (dE==best and np.random.random()<1.0/(nb_+1)):
                    if dE<best: nb_=1
                    else: nb_+=1
                    best=dE; bk=k; bq=q
        if best>=0:
            for c in range(L):
                if R[c]!=0: w[c]+=winc
        yq=y[bk,bq]
        for c in range(L):
            R[c]+=-2*yq*G[bk,c,bq]+2*M[c,bq,bq]
        flip(y,G,P,M,bk,bq)
        tabu[bk,bq]=step+1+tenure+np.random.randint(max(1,tenure//2+1))
    return maxsteps,False

@nb.njit(cache=True)
def run_shared2(y,M,mu,share,maxsteps,inner,tenure,seed,gauss_seidel):
    """S2 (fixed): each block's turn = `inner` tabu steps on
         sum_c (PAF_k(c) - t_k(c))^2 + mu * sum_c PAF_k(c)^2,   t_k = PAF_k - share*R  (frozen at turn start)
       share=0.25 -> block absorbs a quarter of the residual; share=1 -> full absorption (block LNS).
       Returns as soon as the global residual hits 0."""
    np.random.seed(seed)
    K,r=y.shape; L=M.shape[0]
    G,P=init_state(y,M)
    R=np.zeros(L,np.int64)
    for c in range(L):
        for k in range(K): R[c]+=P[k,c]
    tabu=np.zeros((K,r),np.int64)
    T=np.zeros(L); R0=np.zeros(L,np.int64)
    step=0
    while step<maxsteps:
        for c in range(L): R0[c]=R[c]
        order=np.random.permutation(K)
        for kk in range(K):
            k=order[kk]
            for c in range(L):
                T[c]=P[k,c]-share*(R[c] if gauss_seidel else R0[c])
            for it in range(inner):
                best=1e300; bq=-1; nb_=0
                for q in range(r):
                    if tabu[k,q]>step: continue
                    yq=y[k,q]; dE=0.0
                    for c in range(L):
                        dd=-2*yq*G[k,c,q]+2*M[c,q,q]
                        e0=P[k,c]-T[c]
                        dE+=2*e0*dd+dd*dd+mu*(2*P[k,c]*dd+dd*dd)
                    if dE<best or (dE==best and np.random.random()<1.0/(nb_+1)):
                        if dE<best: nb_=1
                        else: nb_+=1
                        best=dE; bq=q
                yq=y[k,bq]
                for c in range(L): R[c]+=-2*yq*G[k,c,bq]+2*M[c,bq,bq]
                flip(y,G,P,M,k,bq)
                tabu[k,bq]=step+1+tenure+np.random.randint(tenure//2+1); step+=1
                rr=0
                for c in range(L): rr+=R[c]*R[c]
                if rr==0: return step,True
    return step,False
