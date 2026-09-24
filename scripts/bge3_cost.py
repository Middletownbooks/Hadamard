"""Cost model for b>=3 (bordered, two-valued) GS cores over G = Q x Z_h, H = {0} x Z_h, c = |Q| = b+1.
Symmetry: cyclic subgroups S = <(alpha, u)> of Aut(Q) x (Z/h)* (alpha an automorphism of Q, u a unit mod h).
Blocks constant on S-orbits; lag classes = orbits of <S, -1> on G\\{0}; target -4b on H\\{0}, +4 off H.
Gaussian estimate (as for b=1): n = -log2 P(4 random S-invariant blocks with the row-sum profile hit the target),
log2E = sum_k log2 N_k(a_k) - n, over row-sum profiles sum a_k^2 = 4v.
usage: bge3_cost.py v c Qmods(csv) [rmin rmax]      env NS (samples per row sum, default 3000)"""
import sys, os, json, itertools, numpy as np
from math import gcd, log, log2, pi
sys.path.insert(0,'scripts')
from b1_multiplier_enum import signed_sums, count_rowsum
NS=int(os.environ.get('NS','3000'))
v=int(sys.argv[1]); c=int(sys.argv[2]); Q=[int(x) for x in sys.argv[3].split(',')]
rmin=int(sys.argv[4]) if len(sys.argv)>4 else 20; rmax=int(sys.argv[5]) if len(sys.argv)>5 else 200
b=c-1; h=(v+1)//c-1; m=c*h; assert int(np.prod(Q))==c and m==v-b
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
    ss=signed_sums(list(sz)); A=sorted(a for a in ss if a>=0 and (a-m)%2==0)
    out=[]
    for p in itertools.combinations_with_replacement(A,4):
        if sum(x*x for x in p)==4*v: out.append(p)
    return out
rng=np.random.default_rng(0)
res=[]
for key,(d,r,lab,grp,abel) in sorted(seen.items(),key=lambda kv:(kv[1][1],-kv[1][0])):
    if not (rmin<=r<=rmax): continue
    sz=np.bincount(lab); reps=lagclasses(grp); L=len(reps)
    tgt=np.array([-4*b if int(x) in Hset else 4 for x in reps],float)
    profs=profiles(sz)
    if not profs: res.append(dict(d=d,r=r,L=L,note='no row-sum profile')); continue
    stats={}
    def st(a):
        if a in stats: return stats[a]
        from collections import defaultdict
        W=[defaultdict(int) for _ in range(r+1)]; W[r][0]=1
        for i in range(r-1,-1,-1):
            for s_,cnt in W[i+1].items(): W[i][s_+int(sz[i])]+=cnt; W[i][s_-int(sz[i])]+=cnt
        if W[0].get(a,0)==0: stats[a]=None; return None
        Y=np.empty((NS,r),int)
        for t in range(NS):
            rem=a
            for i in range(r):
                p1=W[i+1].get(rem-int(sz[i]),0); p0=W[i+1].get(rem+int(sz[i]),0)
                if rng.random()*(p1+p0)<p1: Y[t,i]=1; rem-=int(sz[i])
                else: Y[t,i]=-1; rem+=int(sz[i])
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
