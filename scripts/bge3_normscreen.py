"""Norm-budget screen for b>=3 configurations (exact necessary condition, fast).
For the compression x_k onto classes of Q x Z_t (t | h), Parseval plus the two-valued target give the identity
    sum_k ||x_k||^2 = N_t = (4v(ct - c + 1) + 4b(c - 1)) / (ct)           (verified exactly on the 479 and 347 cores)
Each S-orbit lies in one class-orbit, so a block's compression is a free choice of an achievable value per class-orbit
(subset sums of +-orbit sizes).  DP per block over class-orbits on (row sum a, norm); then require a row-sum profile
(a_1..a_4), sum a_k^2 = 4v, with norms summing to N_t.
usage: bge3_normscreen.py v c Qmods d r L [t ...]"""
import sys, os, json, itertools, numpy as np
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; d0=int(sys.argv[4]); r0=int(sys.argv[5]); L0=int(sys.argv[6])
ts=[int(x) for x in sys.argv[7:]] or [1]
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(r0),str(r0)]
exec(src[:cut])
cands=[(key,val) for key,val in seen.items() if val[0]==d0 and val[1]==r0 and len(lagclasses(val[3]))==L0]
for gi,(key,(d,r,lab,grp,abel)) in enumerate(cands):
    sz=np.bincount(lab); profs=profiles(sz)
    rec=dict(v=v,c=c,Q=Q,h=h,d=d,r=r,L=L0,group=gi,n_profiles=len(profs),tests={})
    for t in ts:
        if h%t: continue
        cls=np.ravel_multi_index(tuple(np.moveaxis(np.concatenate([co[:,:len(Q)],co[:,-1:]%t],axis=1),-1,0)),Q+[t])
        K=c*t
        # class-orbits under S
        par=list(range(K))
        def f(a):
            while par[a]!=a: par[a]=par[par[a]]; a=par[a]
            return a
        for g in grp:
            for x in range(m):
                a_,b_=f(int(cls[x])),f(int(cls[g[x]]))
                if a_!=b_: par[a_]=b_
        corb={}
        for k in range(K): corb.setdefault(f(k),[]).append(k)
        orbs=list(corb.values())
        k2o={k:i for i,ks in enumerate(orbs) for k in ks}
        # each S-orbit (label) contributes |orbit|/|class-orbit| to every class of its class-orbit
        vals=[]
        for i,ks in enumerate(orbs):
            parts=[]
            for o in range(r):
                el=np.where(lab==o)[0]; oo={k2o[int(cls[x])] for x in el}
                assert len(oo)==1
                if oo.pop()==i: parts.append(len(el)//len(ks))
            S={0}
            for p_ in parts: S={s+p_ for s in S}|{s-p_ for s in S}
            vals.append((len(ks),sorted(S)))
        # DP over class-orbits: reachable (a, norm)
        R={(0,0)}
        for n_o,V in vals:
            R={(a+n_o*x,nr+n_o*x*x) for (a,nr) in R for x in V}
        Nt=(4*v*(c*t-c+1)+4*b*(c-1))
        assert Nt%(c*t)==0; Nt//=c*t
        byA={}
        for a,nr in R: byA.setdefault(abs(a),set()).add(nr)
        ok=False
        for p in profs:
            sets=[byA.get(a,set()) for a in p]
            if any(not s for s in sets): continue
            s12={x+y for x in sets[0] for y in sets[1] if x+y<=Nt}
            s34={x+y for x in sets[2] for y in sets[3] if x+y<=Nt}
            if any((Nt-x) in s34 for x in s12): ok=True; break
        rec['tests'][t]=ok
        rec.setdefault('N_t',{})[t]=Nt
    print(json.dumps(rec,default=int),flush=True)
