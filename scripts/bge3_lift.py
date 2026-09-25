"""Compression-and-lift solver for the c=6 (b=5) architecture of v=347:  v = 6(h+1)-1, G = Z2 x Z3 x Z_h,
H = {(0,0,z)}, target sum_k PAF_k = -4b on H\\{0}, +4 off H.  h = p*t, K = subgroup of Z_h of order p.
Stage A: compressed sequences on Gc = Z2 x Z3 x Z_t (entries in {-p..p}, parity p), constant on orbits of the
         symmetry S (scalar multiplier u), Q-coset sums fixed to 347's pattern + shift; compressed PAF target
         Tc(g') = sum over the K-coset of T (norm at g'=0: 4m - 4b(p-1)).  Annealing with +-2 transfers
         between orbits of equal Q-coset footprint.
Stage B: lift to +-1 on G, constant on S-orbits, K-fibre sums equal to the stage-A solution; annealing with
         opposite-sign orbit swaps of equal fibre footprint.  Full energy.
On success: assemble with 347's corner/border (data/template347.npz), verify H H^T = 4v I.
usage: bge3_lift.py h p u [seed] [max_seconds]     env ITERS_A, ITERS_B, LIFT_TRIES"""
import sys, os, json, time, numpy as np, numba as nb
from math import gcd
h=int(sys.argv[1]); p=int(sys.argv[2]); u=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
tmax=float(sys.argv[5]) if len(sys.argv)>5 else 600
assert h%p==0; t=h//p; b=5; c=6; v=c*(h+1)-1; m=c*h; mc=c*t; N=4*v
ITA=int(float(os.environ.get('ITERS_A','2e6'))); ITB=int(float(os.environ.get('ITERS_B','5e6'))); LT=int(os.environ.get('LIFT_TRIES','4'))
T=np.load('data/template347.npz'); comp=T['comp']; layout=json.loads(str(T['layout']))
rng=np.random.default_rng(seed)
def grp(n):
    co=np.array([(a,bb,z) for a in range(2) for bb in range(3) for z in range(n)])
    idx=lambda X: ((np.asarray(X)[...,0]%2)*3+(np.asarray(X)[...,1]%3))*n+(np.asarray(X)[...,2]%n)
    return co,idx
co,idx=grp(h); coc,idxc=grp(t)
ADD=idx(co[:,None,:]+co[None,:,:]).astype(np.int32); NEG=idx(-co).astype(np.int32)
ADDc=idxc(coc[:,None,:]+coc[None,:,:]).astype(np.int32); NEGc=idxc(-coc).astype(np.int32)
Qc=(co[:,0]*3+co[:,1]).astype(np.int32); Qcc=(coc[:,0]*3+coc[:,1]).astype(np.int32)
pi=idxc(np.stack([co[:,0],co[:,1],co[:,2]%t],1)).astype(np.int32)          # compression map G -> Gc
T_full=np.where((co[:,0]==0)&(co[:,1]==0),-4*b,4).astype(np.int64); T_full[0]=4*m
Tc=np.zeros(mc,np.int64); np.add.at(Tc,pi,T_full)
e=int(np.lcm.reduce([2,3,h])); assert gcd(u,e)==1
def orbits(n,co_,idx_,mod):
    perm=idx_(co_*(u%mod)); lab=-np.ones(len(co_),int); orbs=[]
    for x in range(len(co_)):
        if lab[x]<0:
            o=[]; y=x
            while lab[y]<0: lab[y]=len(orbs); o.append(y); y=perm[y]
            orbs.append(o)
    return lab,orbs
lab,orbs=orbits(h,co,idx,e); labc,orbsc=orbits(t,coc,idxc,int(np.lcm.reduce([2,3,t])))
r=len(orbs); rc=len(orbsc)
def pack(orbs_):
    L=max(len(o) for o in orbs_); A=np.full((len(orbs_),L),-1,np.int32)
    for i,o in enumerate(orbs_): A[i,:len(o)]=o
    return A
OM=pack(orbs); OMc=pack(orbsc)
# footprints: stage A (Q-coset counts per compressed orbit), stage B (fibre counts per full orbit)
footA=np.zeros((rc,6),np.int64)
for i,o in enumerate(orbsc):
    for x in o: footA[i,Qcc[x]]+=1
fkeyA=np.unique(footA,axis=0,return_inverse=True)[1].astype(np.int32)
footB=np.zeros((r,mc),np.int64)
for i,o in enumerate(orbs):
    for x in o: footB[i,pi[x]]+=1
fkeyB=np.unique(footB,axis=0,return_inverse=True)[1].astype(np.int32)
# shifts (Q-coset sums = comp + even shift, sum a^2 = 4v)
R_=2*((h+7)//2); a0=comp.sum(1); allowed=[[s for s in range(-R_,R_+1,2) if np.abs(comp[k]+s).max()<=h] for k in range(4)]
shifts=[]
for s0 in allowed[0]:
    for s1 in allowed[1]:
        for s2 in allowed[2]:
            rem=4*v-(a0[0]+6*s0)**2-(a0[1]+6*s1)**2-(a0[2]+6*s2)**2
            if rem<0: continue
            rt=int(round(rem**0.5))
            if rt*rt!=rem: continue
            for a3 in {rt,-rt}:
                if (a3-a0[3])%6==0 and (a3-a0[3])//6 in allowed[3]:
                    s=(s0,s1,s2,(a3-a0[3])//6); shifts.append((s,comp+np.array(s)[:,None]))
@nb.njit(cache=True)
def pafsum(F,ADD):
    K,n=F.shape; P=np.zeros(n,np.int64)
    for k in range(K):
        for g in range(n):
            s=0
            for x in range(n): s+=F[k,x]*F[k,ADD[x,g]]
            P[g]+=s
    return P
@nb.njit(cache=True)
def anneal_A(Y,F,OM,fkey,ADD,NEG,T,bound,seed,iters,T0,T1):
    """integer orbit values Y (K x rc), F expanded; moves: +2 on o1, -2 on o2 (same footprint), |value|<=bound"""
    np.random.seed(seed); K,rc=Y.shape; n=F.shape[1]; P=pafsum(F,ADD)
    E=0
    for g in range(n): d=P[g]-T[g]; E+=d*d
    D=np.zeros(n,np.int64); Temp=T0; best=E
    for it in range(iters):
        k=np.random.randint(K); o1=np.random.randint(rc); o2=-1
        for tries in range(30):
            q=np.random.randint(rc)
            if q!=o1 and fkey[q]==fkey[o1]: o2=q; break
        if o2<0 or Y[k,o1]+2>bound or Y[k,o2]-2<-bound: continue
        for g in range(n): D[g]=0
        for (oo,dl) in ((o1,2),(o2,-2)):
            for j in range(OM.shape[1]):
                x=OM[oo,j]
                if x<0: break
                fx=F[k,x]
                D[0]+=2*dl*fx+dl*dl
                for g in range(1,n): D[g]+=dl*(F[k,ADD[x,g]]+F[k,ADD[x,NEG[g]]])
                F[k,x]=fx+dl
        e=0
        for g in range(n): d=P[g]+D[g]-T[g]; e+=d*d
        if e<=E or np.random.random()<np.exp((E-e)/Temp):
            Y[k,o1]+=2; Y[k,o2]-=2
            for g in range(n): P[g]+=D[g]
            E=e
            if E<best: best=E
            if E==0: return 0,best
        else:
            for (oo,dl) in ((o2,-2),(o1,2)):
                for j in range(OM.shape[1]):
                    x=OM[oo,j]
                    if x<0: break
                    F[k,x]-=dl
        if it%1000==0: Temp=T0*(T1/T0)**(it/iters)
    return E,best
@nb.njit(cache=True)
def anneal_B(Y,F,OM,fkey,ADD,NEG,T,seed,iters,T0,T1):
    np.random.seed(seed); K,r=Y.shape; n=F.shape[1]; P=pafsum(F,ADD)
    E=0
    for g in range(1,n): d=P[g]-T[g]; E+=d*d
    D=np.zeros(n,np.int64); Temp=T0; best=E
    for it in range(iters):
        k=np.random.randint(K); q1=np.random.randint(r); q2=-1
        for tries in range(30):
            qq=np.random.randint(r)
            if qq!=q1 and fkey[qq]==fkey[q1] and Y[k,qq]==-Y[k,q1]: q2=qq; break
        if q2<0: continue
        for g in range(n): D[g]=0
        for qq in (q1,q2):
            for j in range(OM.shape[1]):
                x=OM[qq,j]
                if x<0: break
                fx=F[k,x]
                for g in range(1,n): D[g]+=-2*fx*(F[k,ADD[x,g]]+F[k,ADD[x,NEG[g]]])
                F[k,x]=-fx
        e=0
        for g in range(1,n): d=P[g]+D[g]-T[g]; e+=d*d
        if e<=E or np.random.random()<np.exp((E-e)/Temp):
            Y[k,q1]=-Y[k,q1]; Y[k,q2]=-Y[k,q2]
            for g in range(1,n): P[g]+=D[g]
            E=e
            if E<best: best=E
            if E==0: return 0,best
        else:
            for qq in (q2,q1):
                for j in range(OM.shape[1]):
                    x=OM[qq,j]
                    if x<0: break
                    F[k,x]=-F[k,x]
        if it%1000==0: Temp=T0*(T1/T0)**(it/iters)
    return E,best
def expand(Y,OM,n):
    F=np.zeros((Y.shape[0],n),np.int64)
    for k in range(Y.shape[0]):
        for i in range(OM.shape[0]):
            for x in OM[i]:
                if x<0: break
                F[k,x]=Y[k,i]
    return F
def init_A(cs):
    """integer orbit values (parity p, |.|<=p) with Q-coset sums cs[k] (random greedy repair)"""
    Y=np.empty((4,rc),np.int64)
    for k in range(4):
        for attempt in range(500):
            y=rng.choice(np.arange(-p,p+1,2),size=rc); cur=footA.T@y
            for it in range(40*rc):
                diff=cs[k]-cur
                if not diff.any(): break
                i=rng.integers(rc); step=2*np.sign(diff@footA[i]) if (diff@footA[i])!=0 else 0
                if step==0 or abs(y[i]+step)>p: continue
                new=cur+step*footA[i]
                if np.abs(cs[k]-new).sum()<np.abs(diff).sum(): y[i]+=step; cur=new
            if (cur==cs[k]).all(): Y[k]=y; break
        else: return None
    return Y
def init_B(xc):
    """+-1 orbit signs with fibre sums xc (per Gc element), greedy repair"""
    Y=np.empty((4,r),np.int64)
    for k in range(4):
        for attempt in range(500):
            y=rng.choice([-1,1],size=r); cur=footB.T@y
            for it in range(40*r):
                diff=xc[k]-cur
                if not diff.any(): break
                i=rng.integers(r); dd=diff@footB[i]
                if dd==0 or np.sign(dd)==y[i]: continue
                new=cur-2*y[i]*footB[i]
                if np.abs(xc[k]-new).sum()<np.abs(diff).sum(): y[i]=-y[i]; cur=new
            if (cur==xc[k]).all(): Y[k]=y; break
        else: return None
    return Y
t0=time.time(); stats=dict(A_runs=0,A_solved=0,lift_runs=0,A_best=None,B_best=None); solved=None; tries=0
while time.time()-t0<tmax and solved is None:
    s,cs=shifts[tries%len(shifts)]; tries+=1
    YA=init_A(cs)
    if YA is None: continue
    FA=expand(YA,OMc,mc); stats['A_runs']+=1
    EA,bA=anneal_A(YA,FA,OMc,fkeyA,ADDc,NEGc,Tc,p,seed*7919+tries,ITA,float(4*mc*p),1.0)
    stats['A_best']=bA if stats['A_best'] is None else min(stats['A_best'],bA)
    if EA!=0: continue
    stats['A_solved']+=1
    for lt in range(LT):
        YB=init_B(FA)
        if YB is None: break
        FB=expand(YB,OM,m); stats['lift_runs']+=1
        EB,bB=anneal_B(YB,FB,OM,fkeyB,ADD,NEG,T_full,seed*104729+tries*31+lt,ITB,float(2*m),1.0)
        stats['B_best']=bB if stats['B_best'] is None else min(stats['B_best'],bB)
        if EB==0:
            P=pafsum(FB,ADD)
            if np.array_equal(P[1:],T_full[1:]): solved=FB; break
    if stats['A_runs']%10==0: print(json.dumps(dict(progress=True,secs=round(time.time()-t0),**stats)),flush=True)
res=dict(h=h,p=p,t=t,v=v,u=u,r=r,rc=rc,n_shifts=len(shifts),secs=round(time.time()-t0,1),solved=solved is not None,**stats)
if solved is not None:
    F=solved; off=20; Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=T['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=T['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=T['Kc'][:,J*6+Qc].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    ok=bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64))); res['hadamard_verified']=ok
    with open(f'data/lift_c6_v{v}.txt','w') as fo:
        fo.write(f'v={v} N={N} b={b} m={m} group=2x3x{h} (mixed radix, last index fastest) u={u} p={p} hadamard_verified={ok}\n')
        for k,nm in enumerate('ABCD'): fo.write(nm+'='+''.join('+' if x>0 else '-' for x in F[k])+'\n')
print(json.dumps(res),flush=True)
