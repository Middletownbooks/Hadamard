"""Symmetry-reduced ladder solver for the c=6 (b=5) architecture of v=347 (G = Z2 x Z3 x Z_h, H = {(0,0,z)}).
Blocks are invariant under S = <u> acting by scalar multiplication (u a unit mod lcm(6,h), u = 2 mod 3 so that
347's coset-sum pattern is S-invariant).  Variables = signs of S-orbits.  Coset sums are fixed to 347's pattern plus
even per-block shifts; moves swap the signs of two orbits with equal coset footprint and opposite signs.
Energy = sum_{g != 0} (sum_k PAF_k(g) - T(g))^2.  On success the full matrix is assembled with 347's corner/border
and H H^T = 4v I is verified.
usage: bge3_rung_sym.py h d [seed] [max_seconds]    env ITERS, T0, T1"""
import sys, os, json, time, itertools, numpy as np, numba as nb
from math import gcd
h=int(sys.argv[1]); dS=int(sys.argv[2]); seed=int(sys.argv[3]) if len(sys.argv)>3 else 0; tmax=float(sys.argv[4]) if len(sys.argv)>4 else 600
b=5; c=6; v=c*(h+1)-1; m=c*h; N=4*v
T=np.load('data/template347.npz'); comp=T['comp']; layout=json.loads(str(T['layout']))
co=np.array([(a,bb,z) for a in range(2) for bb in range(3) for z in range(h)])
def idx(X): X=np.asarray(X); return ((X[...,0]%2)*3+(X[...,1]%3))*h+(X[...,2]%h)
ADD=idx(co[:,None,:]+co[None,:,:]).astype(np.int32); NEG=idx(-co).astype(np.int32)
coset=(co[:,0]*3+co[:,1]).astype(np.int32)
tgt=np.where((co[:,0]==0)&(co[:,1]==0),-4*b,4).astype(np.int64); tgt[0]=0
# symmetry: u of multiplicative order dS on G, u = 2 mod 3 (odd)
e=int(np.lcm.reduce([2,3,h])); cand=[]
for u in range(1,e):
    if gcd(u,e)!=1 or u%3!=2: continue
    k=1; x=u
    while x!=1: x=x*u%e; k+=1
    if k==dS: cand.append(u)
if not cand: print(json.dumps(dict(h=h,d=dS,note='no multiplier of that order with u=2 mod 3')),flush=True); sys.exit()
u=int(os.environ['U']) if os.environ.get('U') else cand[seed%len(cand)]
if u not in cand: print(json.dumps(dict(h=h,d=dS,u=u,note='U not a valid multiplier of that order')),flush=True); sys.exit()
perm=idx(co*u)
lab=-np.ones(m,int); orbs=[]
for x in range(m):
    if lab[x]<0:
        o=[]; y=x
        while lab[y]<0: lab[y]=len(orbs); o.append(y); y=perm[y]
        orbs.append(o)
r=len(orbs)
foot=np.zeros((r,6),np.int64)
for q_,o in enumerate(orbs):
    for x in o: foot[q_,coset[x]]+=1
# admissible shifts (target coset sums must be reachable by orbit signs)
from collections import defaultdict
def reachable(target):   # can orbit signs realise the coset-sum vector 'target'?  DP over footprint groups
    groups=defaultdict(int)
    for q_ in range(r): groups[tuple(foot[q_])]+=1
    S={tuple([0]*6)}
    for fp,cnt in groups.items():
        fp=np.array(fp); S={tuple(np.array(s)+(2*k-cnt)*fp) for s in S for k in range(cnt+1)}
        if len(S)>2e6: return True
    return tuple(target) in S
shifts=[]
R_=2*((h+7)//2)
allowed=[[s_ for s_ in range(-R_,R_+1,2) if np.abs(comp[k]+s_).max()<=h] for k in range(4)]
_rc={}
def reach_ok(k,s_):
    if (k,s_) not in _rc: _rc[(k,s_)]=reachable(comp[k]+s_)
    return _rc[(k,s_)]
a0=comp.sum(1)
for s0 in allowed[0]:
    A0=(a0[0]+6*s0)**2
    for s1 in allowed[1]:
        A1=A0+(a0[1]+6*s1)**2
        if A1>4*v: continue
        for s2 in allowed[2]:
            A2=A1+(a0[2]+6*s2)**2; rem=4*v-A2
            if rem<0: continue
            rt=int(round(rem**0.5))
            if rt*rt!=rem: continue
            for a3 in {rt,-rt}:
                if (a3-a0[3])%6: continue
                s3=(a3-a0[3])//6
                if s3 in allowed[3]:
                    s_=(s0,s1,s2,s3); shifts.append((s_,comp+np.array(s_)[:,None]))
if not shifts: print(json.dumps(dict(h=h,d=dS,u=u,r=r,note='no admissible shift reachable under S')),flush=True); sys.exit()
ITERS=int(float(os.environ.get('ITERS','2e7'))); T0_=float(os.environ.get('T0',str(2*m))); T1_=float(os.environ.get('T1','1.0'))
omem=np.full((r,2*dS+2),-1,np.int32)
for q_,o in enumerate(orbs): omem[q_,:len(o)]=o
fkey=np.unique(foot,axis=0,return_inverse=True)[1].astype(np.int32)     # footprint class of each orbit
@nb.njit(cache=True)
def setup(Y,omem,ADD,m):
    K,r=Y.shape; F=np.zeros((K,m),np.int64)
    for k in range(K):
        for q in range(r):
            for j in range(omem.shape[1]):
                x=omem[q,j]
                if x<0: break
                F[k,x]=Y[k,q]
    P=np.zeros(m,np.int64)
    for k in range(K):
        for g in range(m):
            t=0
            for x in range(m): t+=F[k,x]*F[k,ADD[x,g]]
            P[g]+=t
    return F,P
@nb.njit(cache=True)
def chunk(Y,F,P,E,omem,fkey,ADD,NEG,tgt,m,seed,it0,n,total,T0,T1):
    """n annealing steps continuing a run of 'total' steps (geometric cooling); mutates Y,F,P; returns (E, bestE, solved)"""
    np.random.seed(seed); K,r=Y.shape; D=np.zeros(m,np.int64); best=E
    Temp=T0*(T1/T0)**(it0/total)
    for it in range(it0,it0+n):
        k=np.random.randint(K); q1=np.random.randint(r); q2=-1
        for tries in range(30):
            qq=np.random.randint(r)
            if qq!=q1 and fkey[qq]==fkey[q1] and Y[k,qq]==-Y[k,q1]: q2=qq; break
        if q2<0: continue
        for g in range(m): D[g]=0
        for qq in (q1,q2):
            for j in range(omem.shape[1]):
                x=omem[qq,j]
                if x<0: break
                fx=F[k,x]
                for g in range(1,m): D[g]+=-2*fx*(F[k,ADD[x,g]]+F[k,ADD[x,NEG[g]]])
                F[k,x]=-fx
        e=0
        for g in range(1,m): d=P[g]+D[g]-tgt[g]; e+=d*d
        if e<=E or np.random.random()<np.exp((E-e)/Temp):
            Y[k,q1]=-Y[k,q1]; Y[k,q2]=-Y[k,q2]
            for g in range(1,m): P[g]+=D[g]
            E=e
            if E<best: best=E
            if E==0: return E,best,True
        else:
            for qq in (q2,q1):
                for j in range(omem.shape[1]):
                    x=omem[qq,j]
                    if x<0: break
                    F[k,x]=-F[k,x]
        if it%1000==0: Temp=T0*(T1/T0)**(it/total)
    return E,best,False
rng=np.random.default_rng(seed)
def init(cs):
    # random orbit signs hitting the coset-sum target: greedy by footprint class then repair via DP-free random search
    for attempt in range(300):
        Y=np.empty((4,r),np.int64)
        ok=True
        for k in range(4):
            y=rng.choice([-1,1],size=r); tgtk=cs[k]
            cur=foot.T@y
            for it in range(4*r):
                diff=tgtk-cur
                if not diff.any(): break
                q=rng.integers(r)
                y2=-y[q]; new=cur+2*y2*foot[q]
                if np.abs(tgtk-new).sum()<np.abs(diff).sum(): y[q]=y2; cur=new
            if (cur!=tgtk).any(): ok=False; break
            Y[k]=y
        if ok: return Y
    return None
t0=time.time(); tries=0; solved=None; its=0; bestE=None; CH=int(float(os.environ.get('CHUNK','2e5')))
while time.time()-t0<tmax and solved is None:
    s,cs=shifts[tries%len(shifts)]; Y=init(cs); tries+=1
    if Y is None: continue
    F,P=setup(Y,omem,ADD,m); E=int(sum((P[g]-tgt[g])**2 for g in range(1,m))); runbest=E
    for it0 in range(0,ITERS,CH):
        E,b2,ok=chunk(Y,F,P,E,omem,fkey,ADD,NEG,tgt,m,seed*100003+tries*1009+it0//CH,it0,min(CH,ITERS-it0),ITERS,T0_,T1_)
        its+=CH; runbest=min(runbest,b2)
        if (it0//CH)%50==0: print(json.dumps(dict(progress=True,restart=tries,it=it0,E=int(E),run_best=int(runbest),secs=round(time.time()-t0))),flush=True)
        if ok: solved=F; break
        if time.time()-t0>tmax: break
    bestE=runbest if bestE is None else min(bestE,runbest)
res=dict(h=h,v=v,m=m,d=dS,u=u,r=r,n_shifts=len(shifts),iters_per_run=ITERS,restarts=tries,iters=int(its),best_E=bestE,secs=round(time.time()-t0,1),solved=solved is not None)
if solved is not None:
    F=solved; off=20; Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=T['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=T['K'][:,J*6+coset]; Hm[off+J*m:off+(J+1)*m,:off]=T['Kc'][:,J*6+coset].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    ok=bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64))); res['hadamard_verified']=ok
    with open(f'data/rung_c6_v{v}_sym.txt','w') as fo:
        fo.write(f'v={v} N={N} b={b} m={m} group=2x3x{h} (mixed radix, last index fastest) symmetry u={u} d={dS} hadamard_verified={ok}\n')
        for k,nm in enumerate('ABCD'): fo.write(nm+'='+''.join('+' if t>0 else '-' for t in F[k])+'\n')
print(json.dumps(res),flush=True)
