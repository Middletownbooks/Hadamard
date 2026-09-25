"""Sampled stage B for large lift lists (e.g. v=515, p=5): for a compressed solution xc, draw S random lifts per block
from the product structure (components of the fibre-orbit graph), dedupe, PSD-filter (P <= 4v on characters nontrivial
on H), compute PAF at the lag-class reps, and run the bucketed k-tree join (kt_generic.kt_join4).  Any hit is verified
on the full group, assembled with 347's corner/border (data/template347.npz) and checked H H^T = 4v I.
usage: bge3_stageB_sampled.py h p u sol_index S nbuckets [seed]"""
import sys, os, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from kt_generic import kt_join4
h,p,u,si,S,nbk=[int(float(x)) for x in sys.argv[1:7]]; seed=int(sys.argv[7]) if len(sys.argv)>7 else 0
g=setup(h,p,u); A=g['footB'].T; co=g['co']; m=g['m']; v=g['v']; lab=g['lab']; ADD=g['ADD']; NEG=g['NEG']; T_full=g['T_full']
XC=np.load(os.environ.get('XCFILE',f'data/stageA_h{h}_p{p}_u{u}.npy')); xc=XC[si]
rng=np.random.default_rng(seed); t0=time.time()
def log(**kw): kw['secs']=round(time.time()-t0,1); print(json.dumps(kw),flush=True)
e=int(np.lcm.reduce([2,3,h])); perm=g['idx'](co*(u%e))
cl=-np.ones(m,int); reps=[]
for x in range(1,m):
    if cl[x]<0:
        Sx=set(); y=x
        while y not in Sx: Sx.add(y); Sx.add(int(NEG[y])); y=int(perm[y])
        for z in Sx: cl[z]=len(reps)
        reps.append(x)
reps=np.array(reps)
mods=np.array([2,3,h]); chars=co[co[:,2]!=0]
ccl={}; creps=[]
for ch in chars:
    key=int(g['idx'](ch))
    if key in ccl: continue
    Sx=set(); y=key
    while y not in Sx: Sx.add(y); Sx.add(int(NEG[y])); y=int(perm[y])
    for z in Sx: ccl[z]=len(creps)
    creps.append(key)
ph=np.exp(2j*np.pi*((co[:,None,:]*co[np.array(creps)][None,:,:])/mods).sum(-1))
comps=components(A)
lists=[]
for k in range(4):
    opts=[]
    for cq in comps:
        fis=np.nonzero(A[:,cq].sum(1))[0]; sub=A[np.ix_(fis,cq)]; n=len(cq)
        Y=((np.arange(1<<n)[:,None]>>np.arange(n)[None,:])&1)*2-1
        opts.append((cq,Y[(Y@sub.T==xc[k][fis][None,:]).all(1)].astype(np.int8)))
    Ys=np.empty((S,A.shape[1]),np.int8)
    for cq,ok in opts: Ys[:,cq]=ok[rng.integers(0,len(ok),S)]
    Ys=np.unique(Ys,axis=0)
    keep=np.zeros(len(Ys),bool); PAF=[]
    for a in range(0,len(Ys),20000):
        F=Ys[a:a+20000][:,lab].astype(float)
        P=np.abs(F@ph)**2; ok=(P<=4*v+1e-6).all(1); keep[a:a+20000]=ok
        if ok.any():
            Fo=F[ok].reshape(-1,2,3,h)
            R=np.rint(np.fft.ifftn(np.abs(np.fft.fftn(Fo,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(Fo),-1)[:,reps]
            PAF.append(R.astype(np.int64))
    Yk=Ys[keep]; PAFk=np.concatenate(PAF) if PAF else np.zeros((0,len(reps)),np.int64)
    lists.append(dict(Y=Yk,PAF=PAFk))
    log(stage='B_list',block=k,sampled=int(S),unique=int(len(Ys)),psd_pass=int(len(Yk)))
tg=T_full[reps]; PAFl=[lists[k]['PAF']-(tg if k==0 else 0) for k in range(4)]
sols=[]
def verify(idx):
    Fq=np.array([lists[k]['Y'][idx[k]][lab] for k in range(4)]).astype(np.int64)
    P=sum(np.array([int(np.dot(Fq[k],Fq[k][ADD[:,gg]])) for gg in range(m)]) for k in range(4))
    if np.array_equal(P[1:],T_full[1:]): sols.append(Fq); return True
    return False
kt_join4(PAFl,verify,nbuckets=nbk,max_solutions=1,log=lambda s: print(s,flush=True),seed=seed)
res=dict(h=h,p=p,u=u,v=v,sol_index=si,S=S,solutions=len(sols))
if sols:
    T=np.load('data/template347.npz'); layout=json.loads(str(T['layout'])); F=sols[0]; N=4*v; off=20
    Qc=(co[:,0]*3+co[:,1]).astype(np.int32)
    Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=T['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=T['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=T['Kc'][:,J*6+Qc].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    ok=bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64))); res['hadamard_verified']=ok
    with open(f'data/solution_v{v}_b5.txt','w') as fo:
        fo.write(f'v={v} N={N} b=5 m={m} group=2x3x{h} (mixed radix, last index fastest) u={u} built with 347 corner/border; hadamard_verified={ok}\n')
        for k,nm in enumerate('ABCD'): fo.write(nm+'='+''.join('+' if x>0 else '-' for x in F[k])+'\n')
print(json.dumps(res),flush=True)
