"""Exact lift for the two-sided route: for each yc 4-tuple of a compressed solution, enumerate ALL lifts of each block
matching both compressions (xc: z mod t from stage A; yc: z mod p2), then the exact PSD-filtered join (bge3_stageB_join.run).
Every hit is assembled with 347's corner/border and checked H H^T = 4v I (object check).  Lift lists are cached per
(block, yc).  usage: bge3_twosided_lift.py h p u p2 XCFILE sol_index"""
import sys, os, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components, run
from bge3_twosided import yc_setup, enum_lifts_yc
h,p,u,p2=[int(x) for x in sys.argv[1:5]]; XC=np.load(sys.argv[5]); si=int(sys.argv[6])
g=setup(h,p,u); comps=components(g['footB'].T); foot=yc_setup(g,p2); xc=XC[si]; v=g['v']; m=g['m']; co=g['co']
Z=np.load(f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}.npz'); T=Z['tuples']; t0=time.time()
log=lambda d: print(json.dumps(d),flush=True)
# reduce tuples modulo per-block translation by multiples of t=h/p (keeps xc, u-invariance when (u-1)t = 0 mod h, and PAF;
# shifts yc cyclically in z mod p2): lifts and solutions correspond bijectively, so one representative per class suffices
tt=h//p
if ((u-1)*tt)%h==0 and os.environ.get('CANON','1')=='1':
    def canon(y):
        Y=np.asarray(y).reshape(6,p2); st=tt%p2
        return min(tuple(np.roll(Y,j*st,axis=1).ravel().tolist()) for j in range(p2))
    cmap=[{j:canon(Z[f'Y{k}'][j]) for j in set(T[:,k].tolist())} for k in range(4)]; seenc=set(); keep=[]
    for i,tu in enumerate(T.tolist()):
        c=tuple(cmap[k][tu[k]] for k in range(4))
        if c not in seenc: seenc.add(c); keep.append(i)
    log(dict(stage='canon',tuples=int(len(T)),classes=len(keep))); T=T[keep]
log(dict(stage='start',tuples=int(len(T)),distinct_per_block=[int(len(set(T[:,k]))) for k in range(4)]))
cache={}
def lifts(k,j):
    if (k,j) not in cache:
        cache[(k,j)]=enum_lifts_yc(g,comps,xc[k],Z[f'Y{k}'][j],foot)
        log(dict(stage='lifts',block=k,yc=int(j),n=int(len(cache[(k,j)])),secs=round(time.time()-t0,1)))
    return cache[(k,j)]
def assemble(F):
    import json as js
    Tm=np.load('data/template347.npz'); layout=js.loads(str(Tm['layout'])); N=4*v; off=20; ADD=g['ADD']; NEG=g['NEG']
    Qc=(co[:,0]*3+co[:,1]).astype(np.int32); Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=Tm['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=Tm['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=Tm['Kc'][:,J*6+Qc].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    return Hm, bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64)))
nsol=0
done_t=set()
if os.environ.get('RESUME_LOG') and os.path.exists(os.environ['RESUME_LOG']):
    for l in open(os.environ['RESUME_LOG']):
        if '"stage": "tuple"' in l: done_t.add(json.loads(l)['tuple'])
    log(dict(stage='resume',tuples_done=sorted(done_t)))
for ti,tu in enumerate(T):
    if ti in done_t: continue
    L=[lifts(k,int(tu[k])) for k in range(4)]
    sols=run(h,p,u,xc,log=lambda s: None,lifts=L) or []
    log(dict(stage='tuple',tuple=ti,lift_sizes=[int(len(x)) for x in L],solutions=len(sols),secs=round(time.time()-t0,1)))
    for F in sols:
        Hm,ok=assemble(F); nsol+=1
        log(dict(stage='SOLUTION',tuple=ti,hadamard_verified=ok))
        with open(f'data/solution_v{v}_twosided_u{u}_sol{si}_{nsol}.txt','w') as fo:
            fo.write(f'v={v} N={4*v} b=5 group=Z2xZ3xZ{h} (mixed radix, last index fastest) u={u} built with 347 corner/border; hadamard_verified={ok}\n')
            for k,nm in enumerate('ABCD'): fo.write(nm+'='+''.join('+' if x>0 else '-' for x in F[k])+'\n')
log(dict(stage='done',solutions=nsol,secs=round(time.time()-t0,1)))
