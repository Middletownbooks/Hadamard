"""Ladder rung for the c=6 (b=5) architecture of v=347: v = 6(h+1)-1, G = Z2 x Z3 x Z_h, H = {(0,0,z)}.
1. choose even per-block shifts s so that coset sums = comp347 + s and sum_k a_k^2 = 4v (|coset sum| <= h)
2. simulated annealing on the 4 block functions (swap moves inside a coset keep coset sums fixed),
   energy = sum_{g != 0} (sum_k PAF_k(g) - T(g))^2, T = -4b on H\\{0}, +4 off H
3. assemble the full matrix with 347's corner/border (data/template347.npz) and verify H H^T = 4v I.
usage: bge3_rung.py h [seed] [max_seconds]      output: data/rung_c6_v{v}.txt (on success), JSON log line"""
import sys, json, time, itertools, numpy as np, numba as nb
h=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 0; tmax=float(sys.argv[3]) if len(sys.argv)>3 else 600
b=5; c=6; v=c*(h+1)-1; m=c*h; N=4*v
T=np.load('data/template347.npz'); comp=T['comp']; layout=json.loads(str(T['layout']))
co=np.array([(a,bb,z) for a in range(2) for bb in range(3) for z in range(h)])
def idx(X): X=np.asarray(X); return ((X[...,0]%2)*3+(X[...,1]%3))*h+(X[...,2]%h)
ADD=idx(co[:,None,:]+co[None,:,:]).astype(np.int32); NEG=idx(-co).astype(np.int32)
coset=(co[:,0]*3+co[:,1]).astype(np.int32)
tgt=np.where((co[:,0]==0)&(co[:,1]==0),-4*b,4).astype(np.int64); tgt[0]=0
# 1. shifts
import os
ITERS=int(float(os.environ.get('ITERS','2e7'))); T0_=float(os.environ.get('T0',str(2*m))); T1_=float(os.environ.get('T1','1.0'))
shifts=[]
R_=2*((h+7)//2)
for s_ in itertools.product(range(-R_,R_+1,2),repeat=4):
    cs_=comp+np.array(s_)[:,None]
    if np.abs(cs_).max()>h: continue
    a=cs_.sum(1)
    if (a*a).sum()==4*v: shifts.append((s_,cs_))
if not shifts: print(json.dumps(dict(h=h,v=v,note='no admissible shift')),flush=True); sys.exit()
s,cs=shifts[0]
@nb.njit(cache=True)
def paf(f,ADD):
    m=f.shape[0]; P=np.zeros(m,np.int64)
    for g in range(m):
        t=0
        for x in range(m): t+=f[x]*f[ADD[x,g]]
        P[g]=t
    return P
@nb.njit(cache=True)
def anneal(F,cos_,ADD,NEG,tgt,seed,tmax_iter,T0,T1):
    np.random.seed(seed); K,m=F.shape
    P=np.zeros(m,np.int64)
    for k in range(K):
        for g in range(m):
            t=0
            for x in range(m): t+=F[k,x]*F[k,ADD[x,g]]
            P[g]+=t
    def energy(P):
        e=0
        for g in range(1,m): d=P[g]-tgt[g]; e+=d*d
        return e
    E=energy(P); best=E; Temp=T0
    D=np.zeros(m,np.int64)
    for it in range(tmax_iter):
        k=np.random.randint(K); x=np.random.randint(m)
        # partner y in same coset with opposite sign
        y=-1
        for tries in range(20):
            yy=np.random.randint(m)
            if cos_[yy]==cos_[x] and F[k,yy]==-F[k,x]: y=yy; break
        if y<0: continue
        # delta of flipping x then y (exact, sequential)
        for g in range(m): D[g]=0
        fx=F[k,x]
        for g in range(1,m):
            D[g]+=-2*fx*(F[k,ADD[x,g]]+F[k,ADD[x,NEG[g]]]) if ADD[x,g]!=x else 0
        # apply x flip temporarily
        F[k,x]=-fx
        fy=F[k,y]
        for g in range(1,m):
            D[g]+=-2*fy*(F[k,ADD[y,g]]+F[k,ADD[y,NEG[g]]])
        # correction for 2g=0 double counting: when ADD[x,g]==ADD[x,NEG[g]] the formula counts the pair twice -> exact for PAF sum convention
        e=0
        for g in range(1,m): d=P[g]+D[g]-tgt[g]; e+=d*d
        if e<=E or np.random.random()<np.exp((E-e)/Temp):
            F[k,y]=-fy
            for g in range(1,m): P[g]+=D[g]
            E=e
            if E<best: best=E
            if E==0: return it,0
        else:
            F[k,x]=fx
        if it%1000==0: Temp=T0*(T1/T0)**(it/tmax_iter)
    return tmax_iter,E
rng=np.random.default_rng(seed)
def init():
    F=np.empty((4,m),np.int64)
    for k in range(4):
        for q in range(6):
            ix=np.where(coset==q)[0]; npos=(h+cs[k,q])//2
            vals=np.array([1]*npos+[-1]*(h-npos)); rng.shuffle(vals); F[k,ix]=vals
    return F
t0=time.time(); tries=0; solved=None; its=0
while time.time()-t0<tmax:
    s,cs=shifts[tries%len(shifts)]; F=init(); tries+=1
    it,E=anneal(F,coset,ADD,NEG,tgt,seed*1000+tries,ITERS,T0_,T1_)
    its+=it
    # verify with exact PAF
    P=sum(paf(F[k],ADD) for k in range(4))
    if np.array_equal(P[1:],tgt[1:]): solved=F; break
res=dict(h=h,v=v,m=m,n_shifts=len(shifts),iters_per_run=ITERS,shift=list(s),rowsums=cs.sum(1).tolist(),restarts=tries,iters=int(its),secs=round(time.time()-t0,1),solved=solved is not None)
if solved is not None:
    F=solved; off=20; Hm=np.zeros((N,N),np.int64)
    Hm[:off,:off]=T['C']
    cls=lambda J: J*6+coset
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=T['K'][:,cls(J)]
        Hm[off+J*m:off+(J+1)*m,:off]=T['Kc'][:,cls(J)].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Y,X=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Y,X)]
    ok=bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64)))
    res['hadamard_verified']=ok
    with open(f'data/rung_c6_v{v}.txt','w') as fo:
        fo.write(f'v={v} N={N} b={b} m={m} group=2x3x{h} (mixed radix, last index fastest) hadamard_verified={ok}\n')
        for k,nm in enumerate('ABCD'): fo.write(nm+'='+''.join('+' if t>0 else '-' for t in F[k])+'\n')
    if ok: np.save(f'data/rung_c6_v{v}_H.npy',Hm.astype(np.int8))
print(json.dumps(res),flush=True)
