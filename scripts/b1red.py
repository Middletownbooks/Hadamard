"""Reduced (U-orbit) setting for b=1 GS cores of composite length m=v-1.
RedC(m,U): orbits of U on Z_m (orbit 0 = {0}), lag classes = orbits of <U,-1> on Z_m\\{0},
frequency classes (same sets), orbit-sum matrix S for PSD filtering.
find_U(m,d,r): a subgroup of order d with r orbits.  canon(m,U,blocks): translate each real block to U-invariance."""
import sys, numpy as np
sys.path.insert(0,'scripts')
from b1_multiplier_enum import subgroups, orbits, lagclasses, count_rowsum

def find_U(m,d,r,which=0):
    c=[sorted(U) for U in subgroups(m) if len(U)==d and len(orbits(U,m)[1])==r]
    c.sort(); return c[which], len(c)

class RedC:
    def __init__(s,m,U):
        s.m=m; s.U=sorted(U); s.d=len(U)
        lab,orbs=orbits(U,m); s.orb=lab; s.orbits=orbs; s.r=len(orbs)
        s.size=np.array([len(o) for o in orbs])
        s.lags,s.lagsize=lagclasses(U,m); s.L=len(s.lags)
        # frequency classes: same as lag classes, plus omega=0
        s.freqs=np.concatenate([[0],s.lags])
        x=np.arange(m); E=np.exp(2j*np.pi*np.outer(s.freqs,x)/m)          # (F,m)
        s.S=np.stack([E[:,o].sum(1) for o in orbs],1)                       # (F,r)
    def expand(s,Y): return Y[...,s.orb]
    def paf(s,X):
        F=np.fft.fft(X.astype(float),axis=-1); R=np.fft.ifft(np.abs(F)**2,axis=-1).real
        return np.rint(R).astype(np.int64)
    def count_rowsum(s,a): return count_rowsum(list(s.size),a)

def canon(R,blocks):
    """blocks: (4,m) +-1. Return (4,r) orbit signs after translating each block to U-invariance (None if impossible)."""
    out=[]
    for x in blocks:
        ok=None
        for t in range(R.m):
            y=np.roll(x,-t)
            if all(np.array_equal(y,y[(np.arange(R.m)*u)%R.m]) for u in R.U):
                ok=y; break
        if ok is None: return None
        out.append(np.array([ok[o[0]] for o in R.orbits]))
    return np.array(out)

def read_core(path):
    L={l.split('=')[0]:l.split('=')[1].strip() for l in open(path) if '=' in l and l[0] in 'ABCD'}
    return np.array([[1 if c=='+' else -1 for c in L[k]] for k in 'ABCD'])

def psd_feasible(R,maxrel=22):
    """Exact obstruction test. For each frequency class w != 0, let Z = orbits with S_q(w) != 0.
    If |Z| <= maxrel, enumerate all sign patterns on Z (keeping only those whose partial row sum can be
    completed to the target by the other orbits) -> the exact set of achievable P(w) for row sums 2 and 0.
    Infeasible if no P_A + P_B + P_C + P_D = 4m+4.  Returns (feasible, list of (w, |Z|, ok))."""
    from b1_multiplier_enum import signed_sums
    bad=[]; info=[]
    for f in range(1,len(R.freqs)):
        Sq=R.S[f]; Z=np.where(np.abs(Sq)>1e-7)[0]
        if len(Z)>maxrel: continue
        others=[int(R.size[q]) for q in range(R.r) if q not in set(Z.tolist())]
        ss=signed_sums(others)
        Y=((np.arange(1<<len(Z))[:,None]>>np.arange(len(Z))[None,:])&1)*2-1
        ps=Y@R.size[Z]; Pw=np.rint(np.abs(Y@Sq[Z])**2).astype(np.int64)
        V={a:set(Pw[[(a-p) in ss for p in ps]].tolist()) for a in (2,0)}
        S3={x+y for x in V[0] for y in V[0]}; S3={x+y for x in S3 for y in V[0]}
        ok=any((4*R.m+4-pa) in S3 for pa in V[2])
        info.append((int(R.freqs[f]),len(Z),ok))
        if not ok: bad.append(int(R.freqs[f]))
    return len(bad)==0, info

def joint_feasible(R,maxz=26,nsamp=20000):
    """Exact necessary condition, joint over the integer-valued frequency classes W (P(w) depends only on the
    orbits Z_w with S_q(w) != 0).  Enumerate sign patterns on Z = union Z_w (row sum completable by the other
    orbits, P(w) <= 4m+4 for w in W) -> exact achievable vectors (P(w))_{w in W} for row sums 2 and 0;
    feasible iff some A + B + C + D vector equals (4m+4,...).  Returns (verdict, W, |Z|); verdict None if |Z|>maxz."""
    from b1_multiplier_enum import signed_sums
    Y=np.random.default_rng(0).choice([-1,1],size=(nsamp,R.r)); P=np.abs(Y@R.S[1:].T)**2
    W=[f+1 for f in range(P.shape[1]) if np.all(np.abs(P[:,f]-np.rint(P[:,f]))<1e-6)]
    if not W: return True,[],0
    Z=sorted(set(np.where(np.abs(R.S[W])>1e-7)[1].tolist()))
    if len(Z)>maxz: return None,[int(R.freqs[w]) for w in W],len(Z)
    others=[int(R.size[q]) for q in range(R.r) if q not in set(Z)]; ss=signed_sums(others)
    B=4*R.m+4; vecs={}
    for s in range(0,1<<len(Z),1<<20):
        idx=np.arange(s,min(1<<len(Z),s+(1<<20)))
        Yz=((idx[:,None]>>np.arange(len(Z))[None,:])&1)*2-1
        ps=Yz@R.size[Z]; Pw=np.rint(np.abs(Yz@R.S[W][:,Z].T)**2).astype(np.int64)
        okb=(Pw<=B).all(1)
        for a in (2,0):
            sel=okb&np.array([(a-p) in ss for p in ps])
            vecs.setdefault(a,set()).update(map(tuple,Pw[sel]))
    V0=np.array(sorted(vecs[0])); V2=np.array(sorted(vecs[2]))
    if len(V0)==0 or len(V2)==0: return False,[int(R.freqs[w]) for w in W],len(Z)
    tgt=np.full(len(W),B)
    S2=np.unique((V0[:,None,:]+V0[None,:,:]).reshape(-1,len(W)),axis=0)
    S2=S2[(S2<=B).all(1)]
    need={tuple(x) for x in (tgt-V2)}          # need B+C+D = tgt - A
    S2set={tuple(x) for x in S2}
    for z in V0:
        for nv in need:
            if tuple(np.array(nv)-z) in S2set: return True,[int(R.freqs[w]) for w in W],len(Z)
    return False,[int(R.freqs[w]) for w in W],len(Z)

def compression_feasible(R,t,maxstates=6e7):
    """Exact necessary condition via the Z_t compression xbar = M y (M[c,q] = #elements of orbit q = c mod t).
    DP over orbits of the achievable xbar (row sum = sum xbar); P at w = (m/t) j depends only on xbar.
    Requires sum_k P_k(w) = 4m+4 jointly at all j != 0 with every P_k(w) <= 4m+4."""
    m=R.m; assert m%t==0
    M=np.zeros((t,R.r),np.int64)
    for q,o in enumerate(R.orbits):
        for x in o: M[x%t,q]+=1
    cur=np.zeros((1,t),np.int16)
    for q in range(R.r):
        col=M[:,q].astype(np.int16); cur=np.ascontiguousarray(np.concatenate([cur+col,cur-col]))
        key=cur.view(np.dtype((np.void,2*t))).ravel()
        _,i=np.unique(key,return_index=True); cur=cur[i]
        if len(cur)>maxstates: return None,len(cur)
    cur=cur.astype(np.int64)
    rs=cur.sum(1); j=np.arange(1,t); E=np.exp(-2j*np.pi*np.outer(j,np.arange(t))/t)
    B=4*m+4; out={}
    for a in (2,0):
        X=cur[rs==a]; P=np.rint(np.abs(X@E.T)**2).astype(np.int64); P=P[(P<=B).all(1)]
        out[a]=np.unique(P,axis=0)
    V0,V2=out[0],out[2]
    if len(V0)==0 or len(V2)==0: return False,(len(V2),len(V0))
    S2=np.unique((V0[:,None,:]+V0[None,:,:]).reshape(-1,t-1),axis=0); S2=S2[(S2<=B).all(1)]
    S2set={tuple(x) for x in S2}; tgt=np.full(t-1,B)
    for a_ in V2:
        rem=tgt-a_
        for z in V0:
            if tuple(rem-z) in S2set: return True,(len(V2),len(V0))
    return False,(len(V2),len(V0))
