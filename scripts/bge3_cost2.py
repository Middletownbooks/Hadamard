"""Refined cost model for b>=3 cores: exact Q-compression probability (DP over orbits) + Gaussian model for the
remaining constraints CONDITIONED on the compression (samples drawn from the exact compression fibre).

Why: bge3_cost.py fits one Gaussian to the lag-class PAF vector of random S-invariant blocks. The compression
constraints (sum_k |xbar^_k(psi)|^2 = 4b, tiny) are a large deviation of quadratic functionals (coset sums of PAF =
PAF of the compression), where a Gaussian fit can be off by many bits. Here the compression is handled exactly
(it is LINEAR in the orbit signs: xbar = M y), and the Gaussian is used only for what is left.

Output per (config, symmetry group):
  log2W      = log2 sum over admissible compression quadruples of prod_k N_k(xbar_k)   [exact; N_k = #blocks]
  n_comp     = sum_k log2 N_k(a_k) - log2 W restricted to the profile    [bits paid to the compression]
  n_lift     = -log2 P(remaining PAF constraints | compressions)          [Gaussian, conditioned samples]
  log2E      = log2W - n_lift                                              [compare with bge3_cost log2E]
usage: bge3_cost2.py v c Qmods d r L [group_index]     env NS (samples/block, default 4000), TOPQ (quads to lift-cost)
"""
import sys, os, json, itertools, numpy as np
from math import log2, log, pi
from collections import defaultdict
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; d0=int(sys.argv[4]); r0=int(sys.argv[5]); L0=int(sys.argv[6])
gsel=int(sys.argv[7]) if len(sys.argv)>7 else 0
NS=int(os.environ.get('NS','4000')); TOPQ=int(os.environ.get('TOPQ','6')); MAXST=float(os.environ.get('MAXST','4e6'))
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(r0),str(r0)]
exec(src[:cut])
cands=[(key,val) for key,val in seen.items() if val[0]==d0 and val[1]==r0 and len(lagclasses(val[3]))==L0]
assert cands, 'no symmetry group with that (d,r,L)'
key,(d,r,lab,grp,abel)=cands[gsel]
sz=np.bincount(lab); reps=lagclasses(grp); L=len(reps)
tgt=np.array([-4*b if int(x) in Hset else 4 for x in reps],float)
nQ=len(Q)
# ---- compression structure: classes = cosets of H = Q-part of each element ----
cls=np.ravel_multi_index(tuple(np.moveaxis(co[:,:nQ],-1,0)),Q) if nQ>1 else co[:,0].copy()
Mfull=np.zeros((c,r),np.int64); np.add.at(Mfull,(cls,lab),1)
# class-orbits under S (xbar is constant on them)
par=np.arange(c)
def fnd(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
rep_el=np.array([np.where(cls==k)[0][0] for k in range(c)])
for g in grp:
    img=cls[g[rep_el]]
    for k in range(c):
        x_,y_=fnd(k),fnd(int(img[k]))
        if x_!=y_: par[x_]=y_
corb=np.array([fnd(k) for k in range(c)]); _,corb=np.unique(corb,return_inverse=True); nco=int(corb.max())+1
orbrep=np.array([np.where(corb==o)[0][0] for o in range(nco)])
M=Mfull[orbrep]                       # (nco, r): xbar on class-orbit o = sum_j M[o,j] y_j
csize=np.bincount(corb)               # classes per class-orbit
# characters of Q, as functions on class-orbits
chiv=np.array(np.unravel_index(np.arange(c),Q)).T; Qm=np.array(Q)
Efull=np.exp(-2j*np.pi*((chiv[:,None,:]*chiv[None,:,:])/Qm).sum(-1))      # (class, psi)
E=np.stack([Efull[corb==o].sum(0) for o in range(nco)])                    # (class-orbit, psi)
# ---- feasible row sums (profiles) ----
profs=profiles(sz); Aok=sorted({s*a for p in profs for a in p for s in (1,-1)})
# ---- exact DP over orbits: N(xbar) for xbar with |xbar^(psi)|^2 <= 4b (psi != 1) and row sum in Aok ----
FM=(M.T.astype(float))@E                       # (r, psi) contribution of orbit j (sign +1)
order=np.argsort(-np.abs(FM).sum(1)); Mo=M[:,order]; FMo=FM[order]
rem=np.vstack([np.abs(FMo[q:]).sum(0) for q in range(r)]+[np.zeros(c)])
bound=np.full(c,np.sqrt(4*b)); bound[0]=np.sqrt(4*v)
st={tuple([0]*nco):1.0}; overflow=False
for q in range(r):
    col=Mo[:,q]; ns=defaultdict(float)
    for s_,w in st.items():
        s=np.array(s_)
        for sg in (1,-1):
            ns[tuple(s+sg*col)]+=w
    keys=np.array(list(ns.keys())); vals=np.array(list(ns.values()))
    F=keys.astype(float)@E
    keep=(np.abs(F)-rem[q+1][None,:] <= bound[None,:]+1e-6).all(1)
    keep&=(np.abs(F[:,0].real[:,None]-np.array(Aok,float)[None,:]) <= rem[q+1][0]+1e-6).any(1)
    st={tuple(k):w for k,w in zip(keys[keep],vals[keep])}
    if len(st)>MAXST: overflow=True; break
if overflow:
    print(json.dumps(dict(v=v,c=c,Q=Q,h=h,d=d,r=r,L=L,note='DP state overflow',states=len(st)))); sys.exit()
X=np.array(list(st.keys()),np.int64); Nx=np.array(list(st.values()))
Fx=X.astype(float)@E; Px=np.rint(np.abs(Fx)**2*1e6)/1e6
ax=np.rint(Fx[:,0].real).astype(np.int64)
# ---- admissible quadruples: sum_k Px[:,1:] = 4b, sum a^2 = 4v ----
W=np.concatenate([ax[:,None]**2,Px[:,1:]],axis=1); T=np.concatenate([[4*v],np.full(c-1,4*b)])
keyf=lambda A: [tuple(np.round(t,4)) for t in A]
pairs=defaultdict(list)
for i in range(len(W)):
    S=W[i]+W[i:]
    ok=(S<=T+1e-6).all(1)
    for jj in np.where(ok)[0]:
        pairs[tuple(np.round(S[jj],4))].append((i,i+int(jj)))
quads=[]
for k2,lst in pairs.items():
    need=tuple(np.round(T-np.array(k2),4))
    if need in pairs and k2<=need:
        for (i,j) in lst:
            for (k,l) in pairs[need]:
                quads.append((i,j,k,l))
# weight of each unordered quadruple, counting orderings of the 4 blocks
def mult(q):
    from collections import Counter
    from math import factorial
    cnt=Counter(q); m_=factorial(4)
    for x in cnt.values(): m_//=factorial(x)
    return m_
qw=[(np.log2(Nx[list(q)]).sum()+log2(mult(q)),q) for q in quads]
qw.sort(reverse=True)
log2W=float(np.log2(sum(2.0**(w-qw[0][0]) for w,_ in qw))+qw[0][0]) if qw else float('-inf')
# block counts per row sum (for comparison with bge3_cost)
from b1_multiplier_enum import count_rowsum
rec=dict(v=v,b=b,c=c,Q=Q,h=h,d=d,r=r,L=L,class_orbits=nco,n_compressions=len(X),n_quads=len(quads),log2W=round(log2W,2))
if not quads:
    rec['note']='no admissible compression quadruple'; print(json.dumps(rec)); sys.exit()
# ---- conditional sampling within a compression fibre: start from DP backtrack, then swap moves ----
rng=np.random.default_rng(1)
types=defaultdict(list)
for j in range(r): types[tuple(M[:,j])].append(j)
def fibre_start(xbar):
    """random y in {+-1}^r with M y = xbar (randomised greedy + repair)"""
    for attempt in range(2000):
        y=rng.choice([-1,1],r); cur=M@y; diff=xbar-cur
        for it in range(20*r):
            if not diff.any(): return y
            j=int(rng.integers(r)); dj=-2*y[j]*M[:,j]
            if np.abs(diff-dj).sum()<np.abs(diff).sum() or (rng.random()<0.05 and np.abs(diff-dj).sum()==np.abs(diff).sum()):
                y[j]=-y[j]; diff=diff-dj
        if not diff.any(): return y
    return None
# compression-preserving moves: (a) swap two orbits of the same signature; (b) change the +count of 2-3 signature
# groups along a small integer kernel relation sum_g D_g M_g = 0 (needed: count vectors are not unique per xbar)
glist_=list(types.items()); Mg_=np.array([mc for mc,_ in glist_]); gsz_=[len(o) for _,o in glist_]
rels=[]
for kk_ in (2,3):
    for gs in itertools.combinations(range(len(glist_)),kk_):
        for D in itertools.product(*[range(-min(3,gsz_[g]),min(3,gsz_[g])+1) for g in gs]):
            if any(x==0 for x in D): continue
            if not (sum(D[i]*Mg_[g] for i,g in enumerate(gs))!=0).any(): rels.append((gs,D))
        if len(rels)>400: break
def fibre_samples(xbar,n):
    y=fibre_start(np.asarray(xbar))
    if y is None: return None
    out=np.empty((n,r),np.int64); groups=[g for g in types.values() if len(g)>1]
    for s in range(n):
        for _ in range(3*r):
            if rels and rng.random()<0.3:
                gs,D=rels[int(rng.integers(len(rels)))]; ok=True; flips=[]
                for g,dg in zip(gs,D):
                    orbs=np.array(glist_[g][1]); want=-1 if dg>0 else 1          # flip orbits currently at 'want'
                    cand=orbs[y[orbs]==want]
                    if len(cand)<abs(dg): ok=False; break
                    flips.extend(rng.choice(cand,abs(dg),replace=False).tolist())
                if ok:
                    for j in flips: y[j]=-y[j]
            elif groups:
                g=groups[int(rng.integers(len(groups)))]
                i,j=rng.choice(g,2,replace=False)
                if y[i]!=y[j]: y[i],y[j]=y[j],y[i]
        out[s]=y
    return out
def paf_reps(Y):
    Xs=Y[:,lab].reshape((-1,)+tuple(mods)).astype(float)
    F=np.fft.fftn(Xs,axes=tuple(range(1,len(mods)+1)))
    P=np.fft.ifftn(np.abs(F)**2,axes=tuple(range(1,len(mods)+1))).real
    return np.rint(P.reshape(len(Y),-1))[:,reps]
best=None; lifts=[]
cache={}
for w,q in qw[:TOPQ]:
    S=[]
    for i in q:
        if i not in cache:
            Ys=fibre_samples(X[i],NS)
            cache[i]=None if Ys is None else (lambda P:(P.mean(0),np.cov(P.T)))(paf_reps(Ys))
        S.append(cache[i])
    if any(s is None for s in S): lifts.append(dict(w=round(w,1),note='no fibre start')); continue
    mu=sum(s[0] for s in S)-tgt; C=sum(s[1] for s in S)
    wv,V=np.linalg.eigh(C); keep=wv>1e-6*wv.max(); pr=V.T@mu
    nullbad=bool((~keep).any() and np.abs(pr[~keep]).max()>1e-3)
    lp=keep.sum()*log(4)-0.5*np.sum(np.log(2*pi*wv[keep]))-0.5*np.sum(pr[keep]**2/wv[keep])
    n_lift=-lp/log(2)
    lifts.append(dict(log2Nquad=round(w,1),rowsums=[int(ax[i]) for i in q],n_lift=round(n_lift,1),null_dirs=int((~keep).sum()),null_mismatch=nullbad))
    if not nullbad and (best is None or w-n_lift>best[0]): best=(w-n_lift,n_lift,w)
rec['top_quads']=lifts
if best:
    # log2E: weight-sum over all quads using the best quad's n_lift as representative
    rec['n_lift']=round(best[1],1); rec['log2E']=round(log2W-best[1],1)
    lN=sum(log2(max(count_rowsum(list(sz),a),1)) for a in [ax[i] for i in qw[0][1]])
    rec['log2N_rowsums_topquad']=round(lN,1); rec['n_comp_topquad']=round(lN-qw[0][0],1)
print(json.dumps(rec,default=int))
