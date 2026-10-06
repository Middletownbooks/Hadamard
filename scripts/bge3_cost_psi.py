"""[generalized-Psi variant: target sum_k PAF_k(q,z) = -4 sum_{psi in Psi} psi(q) off 0; b = |Psi|]\nCost model for b>=3 (bordered, two-valued) GS cores over G = Q x Z_h, H = {0} x Z_h, c = |Q| = b+1.
Symmetry: cyclic subgroups S = <(alpha, u)> of Aut(Q) x (Z/h)* (alpha an automorphism of Q, u a unit mod h).
Blocks constant on S-orbits; lag classes = orbits of <S, -1> on G\\{0}; target -4b on H\\{0}, +4 off H.
Gaussian estimate (as for b=1): n = -log2 P(4 random S-invariant blocks with the row-sum profile hit the target),
log2E = sum_k log2 N_k(a_k) - n, over row-sum profiles sum a_k^2 = 4v.
usage: bge3_cost_psi.py v Qmods(csv) all|basis [rmin rmax]      env NS (samples per row sum, default 3000)"""
import sys, os, json, itertools, numpy as np
from math import gcd, log, log2, pi
sys.path.insert(0,'scripts')
from b1_multiplier_enum import signed_sums, count_rowsum
NS=int(os.environ.get('NS','3000'))
v=int(sys.argv[1]); Q=[int(x) for x in sys.argv[2].split(',')]; c=int(np.prod(Q))
PSI=sys.argv[3]                     # 'all' (original architecture) or 'basis' (Q elementary 2-group, Psi = coordinate characters)
rmin=int(sys.argv[4]) if len(sys.argv)>4 else 20; rmax=int(sys.argv[5]) if len(sys.argv)>5 else 200
Qco_=np.array(np.unravel_index(np.arange(c),Q)).T
if PSI=='all':
    psis=[np.array(p) for p in np.array(np.unravel_index(np.arange(1,c),Q)).T]
elif PSI=='basis':
    assert all(x==2 for x in Q); psis=[np.eye(len(Q),dtype=int)[i] for i in range(len(Q))]
else: raise SystemExit('PSI must be all|basis')
b=len(psis); assert (v-b)%c==0; h=(v-b)//c; m=c*h
# sum_{psi in Psi} psi(q), real for conjugation-closed Psi
psisum=np.array([sum(np.cos(2*np.pi*((q*p)/np.array(Q)).sum()) for p in psis) for q in Qco_])
assert np.allclose(psisum,np.rint(psisum)); psisum=np.rint(psisum).astype(int)
mods=Q+[h]; co=np.array(np.unravel_index(np.arange(m),mods)).T
def idx(X): return np.ravel_multi_index(tuple(np.moveaxis(np.asarray(X)%np.array(mods),-1,0)),mods)
# Aut(Q) by generator images
Qm=np.array(Q); Qco=np.array(np.unravel_index(np.arange(c),Q)).T
cand=[np.where((((Qco*Qm[i])%Qm)==0).all(1))[0] for i in range(len(Q))]
autQ=[]
for imgs in itertools.product(*cand):
    I=Qco[list(imgs)]; img=(Qco@I)%Qm
    if len({tuple(x) for x in img})==c: autQ.append(I)
units=[u for u in range(1,h) if gcd(u,h)==1] if h>1 else [1]
def perm_of(I,u):
    X=np.concatenate([(co[:,:len(Q)]@I)%Qm, (co[:,-1:]*u)%h],axis=1); return idx(X)
# symmetry groups S <= Aut(Q) x (Z/h)*, enumerated in compact pair form (pQ: permutation of Q's c elements, u: unit)
qidx=lambda X: np.ravel_multi_index(tuple(np.moveaxis(np.asarray(X)%Qm,-1,0)),Q)
autP=[tuple(qidx((Qco@I)%Qm).tolist()) for I in autQ]; Ibyp={p_:I for p_,I in zip(autP,autQ)}
idQ=tuple(range(c))
def mul(x,y): return (tuple(x[0][i] for i in y[0]), (x[1]*y[1])%h)
def closure(gens_,cap):
    e=(idQ,1%h if h>1 else 0); els={e}; frontier=[e]
    while frontier:
        nf=[]
        for x in frontier:
            for g_ in gens_:
                y=mul(x,g_)
                if y not in els:
                    els.add(y); nf.append(y)
                    if len(els)>cap: return None
        frontier=nf
    return frozenset(els)
def n_orbits(S):   # Burnside
    tot=0
    for p_,u in S: tot+=sum(1 for i in range(c) if p_[i]==i)*gcd((u-1)%h,h) if h>1 else 1
    return tot//len(S)
GMAX=int(os.environ.get('GMAX','48')); PAIRMAX=int(os.environ.get('PAIRMAX','12'))
groups=set()
for p_ in autP:
    for u in units:
        S=closure([(p_,u)],GMAX)
        if S is not None and len(S)>1: groups.add(S)
cyc=[S for S in groups if len(S)<=PAIRMAX]
gen_of={S:max(S,key=lambda x:(len(closure([x],GMAX) or []))) for S in cyc}
cl=list(cyc)
for i in range(len(cl)):
    for j in range(i+1,len(cl)):
        if cl[i]<=cl[j] or cl[j]<=cl[i]: continue
        S=closure([gen_of[cl[i]],gen_of[cl[j]]],GMAX)
        if S is not None: groups.add(S)
seen={}
for S in groups:
    r=n_orbits(S)
    if not (rmin<=r<=rmax): continue
    if any(psisum[p_[i]]!=psisum[i] for p_,u in S for i in range(c)): continue
    grp=[perm_of(Ibyp[p_],u) for p_,u in S]
    lab=-np.ones(m,int); rr=0
    for x in range(m):
        if lab[x]<0:
            for g in grp: lab[g[x]]=rr
            rr+=1
    abel=all(mul(x,y)==mul(y,x) for x in S for y in S)
    seen[S]=(len(grp),rr,lab,grp,abel)
print(json.dumps(dict(n_symmetry_groups=len(groups),in_r_range=len(seen))),flush=True)
neg=idx(-co)
Hset=set(idx(np.concatenate([np.zeros((h,len(Q)),int),np.arange(h)[:,None]],axis=1)).tolist())
def lagclasses(grp):
    cl=-np.ones(m,int); reps=[]
    for x in range(1,m):
        if cl[x]<0:
            S={int(g[x]) for g in grp}|{int(neg[g[x]]) for g in grp}
            for y in S: cl[y]=len(reps)
            reps.append(x)
    return np.array(reps)
def profiles(sz):
    ss=signed_sums(list(sz)); A=sorted(a for a in ss if a>=0 and (a-m)%2==0 and a*a<=4*v)
    out=[]
    for p in itertools.combinations_with_replacement(A,4):
        if sum(x*x for x in p)==4*v: out.append(p)
    return out
rng=np.random.default_rng(0)
res=[]
for key,(d,r,lab,grp,abel) in sorted(seen.items(),key=lambda kv:(kv[1][1],-kv[1][0])):
    if not (rmin<=r<=rmax): continue
    sz=np.bincount(lab); reps=lagclasses(grp); L=len(reps)
    qidx_=np.ravel_multi_index(tuple(co[reps][:,:len(Q)].T),Q)
    tgt=np.array([-4*psisum[q] for q in qidx_],float)
    profs=profiles(sz)
    profs=sorted(profs,key=lambda p:-sum(log2(max(count_rowsum(list(sz),a),1)) for a in p))[:int(os.environ.get('PROFMAX','6'))]
    if not profs: res.append(dict(d=d,r=r,L=L,note='no row-sum profile')); continue
    stats={}
    def st(a):
        if a in stats: return stats[a]
        # vectorised exact sampling: P[i][s] = fraction of sign patterns on orbits i.. with partial sum s
        S0=int(sz.sum()); off=S0; P=np.zeros((r+1,2*S0+1)); P[r][off]=1.0
        for i in range(r-1,-1,-1):
            z=int(sz[i]); P[i][z:]+=0.5*P[i+1][:2*S0+1-z]; P[i][:2*S0+1-z]+=0.5*P[i+1][z:]
        if abs(a)>S0 or P[0][a+off]==0: stats[a]=None; return None
        rem=np.full(NS,a); Y=np.empty((NS,r),int)
        for i in range(r):
            z=int(sz[i]); i1=rem-z+off; i0=rem+z+off
            p1=np.where((i1>=0)&(i1<=2*S0),P[i+1][np.clip(i1,0,2*S0)],0.0)
            p0=np.where((i0>=0)&(i0<=2*S0),P[i+1][np.clip(i0,0,2*S0)],0.0)
            y=np.where(rng.random(NS)*(p1+p0)<p1,1,-1); Y[:,i]=y; rem=rem-y*z
        out=Y
        X=np.array(out)[:,lab].reshape((-1,)+tuple(mods)).astype(float)
        F=np.fft.fftn(X,axes=tuple(range(1,len(mods)+1))); Pf=np.fft.ifftn(np.abs(F)**2,axes=tuple(range(1,len(mods)+1))).real
        P=np.rint(Pf.reshape(len(out),-1))[:,reps]
        stats[a]=(P.mean(0),np.cov(P.T)); return stats[a]
    print(json.dumps(dict(progress=True,d=d,r=r,L=L,n_profiles=len(profs))),flush=True)
    best=None
    for p in profs:
        S=[st(a) for a in p]
        if any(s is None for s in S): continue
        mu=sum(s[0] for s in S)-tgt; C=sum(s[1] for s in S)
        w,V=np.linalg.eigh(C); keep=w>1e-6*w.max(); pr=V.T@mu
        if (~keep).any() and np.abs(pr[~keep]).max()>1e-6: continue          # exact infeasibility in a null direction
        lp=keep.sum()*log(4)-0.5*np.sum(np.log(2*pi*w[keep]))-0.5*np.sum(pr[keep]**2/w[keep])
        n=-lp/log(2); lN=sum(log2(count_rowsum(list(sz),a)) for a in p); lE=lN-n
        if best is None or lE>best['log2E']: best=dict(profile=[int(x) for x in p],n=round(n,1),log2N=round(lN,1),log2E=round(lE,1))
    res.append(dict(d=d,r=r,L=L,abelian=abel,**(best or dict(note='no feasible profile'))))
    print(json.dumps(dict(v=v,b=b,c=c,Q=Q,h=h,**res[-1]),default=int),flush=True)
    if sum(1 for x in res if x.get('log2E',-1)>0)>=int(os.environ.get('STOPAFTER','999')): break
