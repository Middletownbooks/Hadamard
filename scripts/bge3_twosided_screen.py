"""Exact pairwise forced-residue screen for two-sided yc classes.  For a class (yc_A..yc_D) the lift lists are complete
(all lifts matching xc and yc), so per block the set of (PAF(c1) mod M, PAF(c2) mod M) pairs over ALL its lifts is exact.
A class can hold a solution only if, for every pair of lag classes (c1,c2), some choice r_A + r_B + r_C + r_D = T (mod M).
Necessary condition, exact (no sampling).  Presence tables are cached per (block, yc).
lib: Screen(g, comps, xc, foot, reps, M).check(Ylists_per_block) ; usage: bge3_twosided_screen.py (validation on 347)"""
import sys, time, itertools, numpy as np
sys.path.insert(0,'scripts')
from bge3_twosided import enum_lifts_yc

class Screen:
    def __init__(s,g,comps,xc,foot,reps,Ms=(16,64)):
        s.g=g; s.comps=comps; s.xc=xc; s.foot=foot; s.reps=np.asarray(reps); s.Ms=Ms; s.cache={}
        s.h=int(g['co'][:,2].max())+1; s.T=g['T_full'][s.reps]; s.pairs=list(itertools.combinations(range(len(s.reps)),2))
    def table(s,k,yc):
        key=(k,tuple(int(x) for x in yc))
        if key in s.cache: return s.cache[key]
        Y=enum_lifts_yc(s.g,s.comps,s.xc[k],yc,s.foot); lab=s.g['lab']; R=np.empty((len(Y),len(s.reps)),np.int64)
        for a in range(0,len(Y),20000):
            F=Y[a:a+20000][:,lab].reshape(-1,2,3,s.h).astype(float)
            R[a:a+20000]=np.rint(np.fft.ifftn(np.abs(np.fft.fftn(F,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(F),-1)[:,s.reps]
        tabs={}
        for M in s.Ms:
            Rm=R%M; pres=np.zeros((len(s.pairs),M*M),bool)
            for pi,(c1,c2) in enumerate(s.pairs): pres[pi,np.unique(Rm[:,c1]*M+Rm[:,c2])]=True
            tabs[M]=pres
        s.cache[key]=(len(Y),tabs); return s.cache[key]
    def check(s,ycs):
        """ycs: 4 yc vectors.  returns (ok, info)"""
        T=[s.table(k,ycs[k]) for k in range(4)]
        if any(t[0]==0 for t in T): return False,dict(reason='empty lift list')
        for M in s.Ms:
            idx=np.arange(M*M); a=idx//M; b=idx%M
            for pi,(c1,c2) in enumerate(s.pairs):
                P=[np.nonzero(T[k][1][M][pi])[0] for k in range(4)]
                AB={(int((x//M+y//M)%M),int((x%M+y%M)%M)) for x in P[0] for y in P[1]}
                need=((int(s.T[c1])%M),(int(s.T[c2])%M))
                ok=any(((need[0]-(x//M+y//M))%M,(need[1]-(x%M+y%M))%M) in AB for x in P[2] for y in P[3])
                if not ok: return False,dict(M=M,pair=(int(s.reps[c1]),int(s.reps[c2])),lift_sizes=[t[0] for t in T])
        return True,dict(lift_sizes=[t[0] for t in T])

if __name__=='__main__':
    import json
    from bge3_stageB_join import setup, components
    from bge3_twosided import yc_setup
    from bge3_residue import class_reps
    g=setup(57,3,11); comps=components(g['footB'].T); foot=yc_setup(g,3); co=g['co']
    Fr=np.array(np.load('data/template347.npz')['F']).astype(np.int64)
    xc=np.zeros((4,g['mc']),np.int64)
    for k in range(4): np.add.at(xc[k],g['pi'],Fr[k])
    cell=((co[:,0]%2)*3+(co[:,1]%3))*3+(co[:,2]%3); yr=np.zeros((4,18),np.int64)
    for k in range(4): np.add.at(yr[k],cell,Fr[k])
    S=Screen(g,comps,xc,foot,class_reps(g,11,57))
    ok,info=S.check(yr); print(json.dumps(dict(check='347 real class must pass',ok=ok,**info)))
    # negative test: perturb one block's yc (move 2 between two cells of the same coset) -> expect a fail or an empty list
    y2=yr.copy(); y2[0,0]+=2; y2[0,1]-=2
    ok2,info2=S.check(y2); print(json.dumps(dict(check='perturbed yc (should fail or be empty)',ok=ok2,**info2)))
