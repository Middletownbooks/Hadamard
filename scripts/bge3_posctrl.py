"""Blind positive control for the c=6 compression-and-lift pipeline on a small rung (see posctrl/PREREGISTRATION.md).
Per slice (h, p, u): exact stage A (all shifts; exact 4-list join over complete compressed-block lists -> ALL
compressions), exact lift of every compression (bge3_stageB_join.run on complete lift lists), residue filters on every
compression (false-rejection check), fast three-sided lift (wc_lift) on every yc class that holds exact solutions
(completeness check), and an object check of every hit (full 4v x 4v matrix, H H^T = 4v I).
usage: bge3_posctrl.py h p u [maxpairs]"""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components, run
from bge3_stageA import enum_block
from bge3_residue import residue_ok, joint_residue_ok, class_reps
from bge3_twosided import yc_setup, yc_distribution, comp_options
from bge3_wc_lift import wc_lift
h,p,u=[int(x) for x in sys.argv[1:4]]; t0=time.time()
g=setup(h,p,u); v=g['v']; m=g['m']; co=g['co']; comps=components(g['footB'].T); labc=g['labc']; t=g['t']; mc=g['mc']
p2=h//p; foot=yc_setup(g,p2); reps=class_reps(g,u,h); Tc=g['Tc']
log=lambda d: print(json.dumps(d),flush=True)
def cpaf(X):   # compressed PAF, all lags, on Z2 x Z3 x Z_t
    if len(X)==0: return np.zeros((0,mc),np.int64)
    F=X.reshape(-1,2,3,t).astype(float)
    return np.rint(np.fft.ifftn(np.abs(np.fft.fftn(F,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(X),-1).astype(np.int64)
def assemble_ok(F,flip=False):
    Tm=np.load('data/template347.npz'); layout=json.loads(str(Tm['layout'])); N=4*v; off=20; ADD=g['ADD']; NEG=g['NEG']
    F=[f.copy() for f in F]
    if flip: F[0][0]*=-1
    Qc=(co[:,0]*3+co[:,1]).astype(np.int32); Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=Tm['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=Tm['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=Tm['Kc'][:,J*6+Qc].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    return bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64))), Hm
# ---- exact stage A
comp_list=[]; rng=np.random.default_rng(5); W=rng.integers(1,2**62,size=mc,dtype=np.int64).astype(np.uint64)
for si,(s,cs) in enumerate(g['shifts']):
    L=[]
    for k in range(4):
        Y,trunc=enum_block(g,cs[k],v,p)
        if trunc:
            res=dict(stage='slice_done',h=h,v=v,p=p,u=u,status='stageA_too_large',note='a compressed block list exceeds 5e7; exhaustive stage A infeasible')
            log(res); json.dump(dict(summary=res),open(f'posctrl/slice_h{h}_p{p}_u{u}.json','w')); sys.exit(0)
        X=Y[:,labc].astype(np.int64); L.append((X,cpaf(X)))
    sizes=[len(x[0]) for x in L]
    if min(sizes)==0: log(dict(stage='A_shift',shift=[int(x) for x in s],sizes=sizes,compressions=0)); continue
    PL=[L[0][1]-Tc,L[1][1],L[2][1],L[3][1]]; H_=[(P_.astype(np.uint64)*W).sum(1) for P_ in PL]
    # pairing with the smallest stored side; stream the other side in chunks
    pairings=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
    best=min(((min(len(H_[a])*len(H_[b]),len(H_[c])*len(H_[d])),(a,b),(c,d)) for (a,b),(c,d) in pairings))
    nst,(a,b),(c,d)=best
    if len(H_[a])*len(H_[b])>len(H_[c])*len(H_[d]): (a,b),(c,d)=(c,d),(a,b)
    if nst>2e8:
        res=dict(stage='slice_done',h=h,v=v,p=p,u=u,status='stageA_join_too_large',shift=[int(x) for x in s],sizes=sizes,stored_pairs=int(nst))
        log(res); json.dump(dict(summary=res),open(f'posctrl/slice_h{h}_p{p}_u{u}.json','w')); sys.exit(0)
    ST=(H_[a][:,None]+H_[b][None,:]).ravel(); o=np.argsort(ST); STs=ST[o]; del ST; found=0
    step=max(1,int(2e7//max(len(H_[d]),1)))
    for c0 in range(0,len(H_[c]),step):
        PR=(np.uint64(0)-(H_[c][c0:c0+step,None]+H_[d][None,:])).ravel(); lo=np.searchsorted(STs,PR); hi=np.searchsorted(STs,PR,side='right')
        for q in np.nonzero(hi>lo)[0]:
            ic,id_=divmod(int(q),len(H_[d])); ic+=c0
            for r in range(lo[q],hi[q]):
                ia,ib=divmod(int(o[r]),len(H_[b])); idx={a:ia,b:ib,c:ic,d:id_}
                if np.array_equal(sum(PL[k][idx[k]] for k in range(4)),np.zeros(mc,np.int64)):
                    comp_list.append((si,np.array([L[k][0][idx[k]] for k in range(4)]))); found+=1
    log(dict(stage='A_shift',shift=[int(x) for x in s],sizes=sizes,compressions=found,secs=round(time.time()-t0,1)))
ncomp_total=len(comp_list); MAXLIFT=int(__import__('os').environ.get('MAXLIFT','3000'))
if ncomp_total>MAXLIFT:
    pick=np.random.default_rng(2026).choice(ncomp_total,MAXLIFT,replace=False); comp_list=[comp_list[i] for i in pick]   # random order (time cap below)
log(dict(stage='A_done',compressions=ncomp_total,lifted=len(comp_list),sampled=ncomp_total>MAXLIFT))
# ---- exact lift of every compression, filters, fast-lift completeness, object checks
cellq=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2)

def lift_count_log2(k,xck):
    A=g['footB'].T; return float(sum(np.log2(max(len(comp_options(A,cq,xck)),1)) for cq in comps))
def lift_via_classes(xc):
    """production path, exact: exact yc distribution per block (DP over components), exact yc-level join,
    then the fast three-sided lift on every valid yc 4-tuple."""
    D=[yc_distribution(g,comps,xc[k],foot) for k in range(4)]
    if any(d is None for d in D): return None
    if any(len(d)==0 for d in D): return []
    Y=[np.array(list(d.keys()),np.int64) for d in D]
    cell2=np.array([(a,b,z) for a in range(2) for b in range(3) for z in range(p2)])
    def paf2(Yk):
        F=Yk.reshape(-1,2,3,p2).astype(float)
        return np.rint(np.fft.ifftn(np.abs(np.fft.fftn(F,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(Yk),-1).astype(np.int64)
    T2=np.zeros(6*p2,np.int64); np.add.at(T2,cellq,g['T_full'])
    P=[paf2(y) for y in Y]; Wy=np.random.default_rng(9).integers(1,2**62,size=6*p2,dtype=np.int64).astype(np.uint64)
    hA=((P[0]-T2).astype(np.uint64)*Wy).sum(1); hB=(P[1].astype(np.uint64)*Wy).sum(1); hC=(P[2].astype(np.uint64)*Wy).sum(1); hD=(P[3].astype(np.uint64)*Wy).sum(1)
    AB=(hA[:,None]+hB[None,:]).ravel(); CD=(np.uint64(0)-(hC[:,None]+hD[None,:])).ravel(); o=np.argsort(CD); CDs=CD[o]
    lo=np.searchsorted(CDs,AB); hi=np.searchsorted(CDs,AB,side='right'); sols=[]
    for q in np.nonzero(hi>lo)[0]:
        ia,ib=divmod(int(q),len(hB))
        for r in range(lo[q],hi[q]):
            ic,id_=divmod(int(o[r]),len(hD))
            if not np.array_equal(P[0][ia]+P[1][ib]+P[2][ic]+P[3][id_],T2): continue
            fs,_=wc_lift(g,comps,foot,xc,[Y[0][ia],Y[1][ib],Y[2][ic],Y[3][id_]],u); sols+=fs
    return sols
nlift=0; nsol=0; false_rej=0; pc3_bad=0; hadamard_bad=0; neg_ok=None; per=[]
LIFT_TMAX=float(__import__('os').environ.get('LIFT_TMAX','1e12')); t_lift=time.time(); nlifted=0
if ncomp_total<=MAXLIFT: comp_list=[comp_list[i] for i in np.random.default_rng(2027).permutation(len(comp_list))]
for ci,(si,xcX) in enumerate(comp_list):
    if time.time()-t_lift>LIFT_TMAX: log(dict(stage='lift_timecap',lifted=ci)); break
    nlifted=ci+1
    xc=np.zeros((4,mc),np.int64)
    for k in range(4): xc[k]=xcX[k]          # compressed values per element of Gc
    big=max(lift_count_log2(k,xc[k]) for k in range(4))>21
    if big:
        try: sols=lift_via_classes(xc); method='classes'
        except MemoryError: sols=None
        if sols is None: per.append(dict(c=ci,shift_index=si,solutions=None,filters_pass=None,method='lift_too_large')); continue
    else:
        sols=run(h,p,u,xc,log=lambda s: None,maxpairs=3e8) or []; method='exact_join'
    ok1,_=residue_ok(g,comps,reps,xc,seed=ci); ok2=ok1
    if ok1:
        for M in (16,64):
            okj,_=joint_residue_ok(g,comps,reps,xc,h,ns=4000,M=M,seed=ci)
            if not okj: ok2=False; break
    if sols:
        nlift+=1; nsol+=len(sols)
        if not ok2: false_rej+=1
        for F in sols:
            ok,_=assemble_ok(F)
            if not ok: hadamard_bad+=1
            if neg_ok is None: neg_ok=not assemble_ok(F,flip=True)[0]
        # PC3: per yc class among the exact solutions, the fast lift must return the same solutions
        classes={}
        for F in sols:
            yc=np.zeros((4,6*p2),np.int64)
            for k in range(4): np.add.at(yc[k],cellq,F[k])
            classes.setdefault(yc.tobytes(),(yc,[]))[1].append(F)
        for key,(yc,Fs) in (classes.items() if method=='exact_join' else []):
            fs,_=wc_lift(g,comps,foot,xc,yc,u)
            A={tuple(np.concatenate(F).tolist()) for F in Fs}; B={tuple(np.concatenate(F).tolist()) for F in fs}
            if A!=B: pc3_bad+=1
    per.append(dict(c=ci,shift_index=si,solutions=len(sols),filters_pass=bool(ok2),method=method))
    if ci%20==0: log(dict(stage='lift_progress',done=ci+1,of=len(comp_list),liftable=nlift,solutions=nsol,secs=round(time.time()-t0,1)))
res=dict(stage='slice_done',h=h,v=v,p=p,u=u,compressions=ncomp_total,lifted=nlifted,sampled=(ncomp_total>MAXLIFT or nlifted<len(comp_list)),liftable=nlift,lift_rate=(nlift/nlifted if nlifted else None),
         solutions=nsol,filter_false_rejections=false_rej,filter_pass_total=int(sum(x['filters_pass'] for x in per)),
         fastlift_mismatches=pc3_bad,hadamard_failures=hadamard_bad,negative_test_rejected=neg_ok,secs=round(time.time()-t0,1))
log(res); json.dump(dict(summary=res,per_compression=per),open(f'posctrl/slice_h{h}_p{p}_u{u}.json','w'))
