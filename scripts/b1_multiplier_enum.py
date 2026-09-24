"""b=1 multiplier configurations for GS cores of length m=v-1 (even):
subgroups U of (Z/m)* (1- and 2-generated), orbits of U on Z_m, feasibility of row sums (+-2,0,0,0)
with U-invariant +-1 blocks (signed orbit-size sums), and a Gaussian estimate of constraint bits
n = -log2 P(4 random U-invariant blocks with row sums (2,0,0,0) satisfy sum_k PAF_k(s) = -4 for all s != 0).
usage: b1_multiplier_enum.py v [v ...]   (env TOPK: how many smallest-r feasible configs to cost per v)"""
import sys, os, json, numpy as np
from math import gcd, log2, pi, log
TOPK=int(os.environ.get('TOPK','4')); NS=int(os.environ.get('NS','4000'))
def units(m): return [a for a in range(1,m) if gcd(a,m)==1]
def closure(gens,m):
    S={1}; frontier=[1]
    while frontier:
        nf=[]
        for x in frontier:
            for g in gens:
                y=x*g%m
                if y not in S: S.add(y); nf.append(y)
        frontier=nf
    return frozenset(S)
def subgroups(m):
    us=units(m); cyc={}
    for a in us: cyc.setdefault(closure([a],m),a)
    subs=set(cyc)
    cl=list(cyc)
    for i in range(len(cl)):
        for j in range(i+1,len(cl)):
            if not cl[i]<=cl[j] and not cl[j]<=cl[i]:
                subs.add(closure(list(cl[i])[:1]+[cyc[cl[i]],cyc[cl[j]]],m))
    return subs
def orbits(U,m):
    lab=-np.ones(m,int); orbs=[]
    for x in range(m):
        if lab[x]<0:
            o=sorted({x*u%m for u in U})
            for y in o: lab[y]=len(orbs)
            orbs.append(o)
    return lab,orbs
def lagclasses(U,m):
    seen=set(); reps=[]; sizes=[]
    for s in range(1,m):
        if s in seen: continue
        c={s*u%m for u in U}|{(-s*u)%m for u in U}; seen|=c; reps.append(s); sizes.append(len(c))
    return np.array(reps),np.array(sizes)
def signed_sums(sizes):
    """set of achievable sum_q eps_q*size_q"""
    S={0}
    for s in sizes: S={x+s for x in S}|{x-s for x in S}
    return S
def count_rowsum(sz,a):
    from collections import defaultdict
    D={0:1}
    for s_ in sz:
        E=defaultdict(int)
        for x,c in D.items(): E[x+s_]+=c; E[x-s_]+=c
        D=E
    return D.get(a,0)
def cost(U,m,rng):
    lab,orbs=orbits(U,m); sz=np.array([len(o) for o in orbs]); r=len(orbs)
    reps,_=lagclasses(U,m); L=len(reps)
    def sample(a,n):
        out=[]; tries=0
        while len(out)<n and tries<400:
            Y=rng.choice([-1,1],size=(20000,r)); rs=Y@sz; Y=Y[rs==a]; out.extend(Y[:n-len(out)]); tries+=1
        return np.array(out)
    stats={}
    for a in (0,2):
        Y=sample(a,NS)
        if len(Y)<200: return None
        X=Y[:,lab]; F=np.fft.fft(X.astype(float),axis=1); P=np.rint(np.fft.ifft(np.abs(F)**2,axis=1).real)[:,reps]
        stats[a]=(P.mean(0),np.cov(P.T))
    mu=3*stats[0][0]+stats[2][0]+4.0     # target -4 -> shift
    C=3*stats[0][1]+stats[2][1]
    w,V=np.linalg.eigh(C); keep=w>1e-6*w.max()
    proj=V.T@mu; D=int(keep.sum())-1                      # one direction implied by the row sums
    idx=np.argsort(-w)[:D]
    lp=D*log(4)-0.5*np.sum(np.log(2*pi*w[idx]))-0.5*np.sum(proj[idx]**2/w[idx])
    null_mis=float(np.abs(proj[~keep]).max()) if (~keep).any() else 0.0
    N0=count_rowsum(sz,0); N2=count_rowsum(sz,2)
    n=-lp/log(2); lE=(log2(N2)+3*log2(N0)-n) if N0 and N2 else None
    return dict(r=r,L=L,dims=D,n_cond=round(float(n),1),log2E=None if lE is None else round(float(lE),1),null_mismatch=round(null_mis,2))
rng=np.random.default_rng(0)
for v in ([int(x) for x in sys.argv[1:]] if __name__=='__main__' else []):
    m=v-1; res=[]
    for U in subgroups(m):
        if len(U)==1: continue
        lab,orbs=orbits(U,m); sz=[len(o) for o in orbs]
        if len(orbs)>400: continue
        ss=signed_sums(sz)
        feas=(0 in ss) and (2 in ss)
        res.append(dict(d=len(U),r=len(orbs),feasible=feas,U=U))
    feas=sorted([x for x in res if x['feasible'] and x['r']<=int(os.environ.get('RMAX','160'))],key=lambda x:(x['r'],-x['d']))
    out=dict(v=v,m=m,n_subgroups=len(res),configs=[])
    for x in feas:
        c=cost(sorted(x['U']),m,rng)
        if c is None: continue
        if c['null_mismatch']>1e-6: continue          # exactly infeasible (degenerate classes miss the target)
        out['configs'].append(dict(d=x['d'],**c))
    print(json.dumps(out),flush=True)
