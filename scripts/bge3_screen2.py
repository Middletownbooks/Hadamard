"""[v2: row sums restricted to feasible profiles] Exact obstruction screen for b>=3 configurations (G = Q x Z_h, H = {0} x Z_h, symmetry S from bge3_cost).
Compress onto Q x Z_t (t | h): xbar = M y with M[class, orbit] = #orbit elements in that class.  P at every character
of Q x Z_t depends only on xbar.  Required: sum_k P_k(chi) = 4v if chi is nontrivial on the Z_t part, 4b if trivial on
Z_t but nontrivial on Q, and sum_k a_k^2 = 4v at chi = 0; every P_k(chi) <= that target.
DP over orbits of achievable xbar per block (capped), then a pairwise meet-in-the-middle over the 4 blocks.
usage: bge3_screen.py v c Qmods d r L [t ...]   (selects the symmetry group with that (d, r, L) from bge3_cost's list)"""
import sys, os, json, itertools, numpy as np
from math import gcd
sys.argv_saved=list(sys.argv)
v=int(sys.argv[1]); c=int(sys.argv[2]); Qs=sys.argv[3]; d0=int(sys.argv[4]); r0=int(sys.argv[5]); L0=int(sys.argv[6])
ts=[int(x) for x in sys.argv[7:]] or [1]
MAXST=float(os.environ.get('MAXST','3e6'))
# reuse bge3_cost's group enumeration by exec'ing its setup portion
src=open('scripts/bge3_cost.py').read(); cut=src.index("rng=np.random.default_rng(0)\nres=[]")
sys.argv=['bge3_cost.py',str(v),str(c),Qs,str(r0),str(r0)]
exec(src[:cut])
cands=[(key,val) for key,val in seen.items() if val[0]==d0 and val[1]==r0 and len(lagclasses(val[3]))==L0]
out=[]
for gi,(key,(d,r,lab,grp,abel)) in enumerate(cands):
    sz=np.bincount(lab); rec=dict(v=v,c=c,Q=Q,h=h,d=d,r=r,L=L0,group=gi,tests={})
    profs=profiles(sz); Aok=sorted({s*a for p in profs for a in p for s in (1,-1)})
    rec['n_profiles']=len(profs)
    if not profs: rec['tests']={t:False for t in ts}; out.append(rec); print(json.dumps(rec,default=int),flush=True); continue
    Aok_arr=np.array(Aok,float)
    for t in ts:
        if h%t: continue
        # class of each element in Q x Z_t
        cls=np.ravel_multi_index(tuple(np.moveaxis(np.concatenate([co[:,:len(Q)],co[:,-1:]%t],axis=1),-1,0)),Q+[t])
        K=c*t; M=np.zeros((K,r),np.int64); np.add.at(M,(cls,lab),1)
        # symmetry reduction: S permutes classes; xbar is constant on class-orbits (rows of M equal there)
        cpar=np.arange(K)
        def cf(a_):
            while cpar[a_]!=a_: cpar[a_]=cpar[cpar[a_]]; a_=cpar[a_]
            return a_
        rep_el=np.array([np.where(cls==k)[0][0] for k in range(K)])
        for g in grp:
            img=cls[g[rep_el]]
            for k in range(K):
                x_,y_=cf(k),cf(int(img[k]))
                if x_!=y_: cpar[x_]=y_
        corb=np.array([cf(k) for k in range(K)]); _,corb=np.unique(corb,return_inverse=True); nco=corb.max()+1
        orbrep=np.array([np.where(corb==o)[0][0] for o in range(nco)])
        assert all(np.array_equal(M[k],M[orbrep[corb[k]]]) for k in range(K))
        Mfull=M; M=M[orbrep]                                                    # (nco, r): one row per class-orbit
        rec.setdefault('class_orbits',{})[t]=int(nco)
        # characters of Q x Z_t and their bounds: |F(chi)| <= sqrt(target) at the end
        chiv=np.array(np.unravel_index(np.arange(K),Q+[t])).T; mods_=np.array(Q+[t])
        Efull=np.exp(-2j*np.pi*((np.array(np.unravel_index(np.arange(K),Q+[t])).T[:,None,:]*chiv[None,:,:])/mods_).sum(-1))   # (class, chi)
        E=np.stack([Efull[corb==o].sum(0) for o in range(nco)])                   # (class-orbit, chi): xbar_orbit value * sum of members
        tg=np.where(chiv[:,-1]!=0,4*v,np.where((chiv[:,:-1]!=0).any(1),4*b,4*v)).astype(float)
        FM=(M.T.astype(float))@E                                              # orbit q -> contribution to F(chi)
        rem=np.vstack([np.abs(FM[q:]).sum(0) for q in range(r)]+[np.zeros(K)])  # max remaining |contribution|
        order=np.argsort(-np.abs(FM).sum(1)); FM=FM[order]; Mo=M[:,order]
        rem=np.vstack([np.abs(FM[q:]).sum(0) for q in range(r)]+[np.zeros(K)])
        cur=np.zeros((1,nco),np.int32); ok=True
        for q in range(r):
            col=Mo[:,q].astype(np.int32); cur=np.ascontiguousarray(np.concatenate([cur+col,cur-col]))
            _,i=np.unique(cur.view(np.dtype((np.void,4*nco))).ravel(),return_index=True); cur=cur[i]
            keep=np.empty(len(cur),bool)
            for c0 in range(0,len(cur),500000):                     # chunked: Fc for all states at once can OOM
                Fc=cur[c0:c0+500000].astype(float)@E
                kk=(np.abs(Fc)-rem[q+1][None,:] <= np.sqrt(tg)[None,:]+1e-6).all(1)
                kk&=(np.abs(Fc[:,0].real[:,None]-Aok_arr[None,:]) <= rem[q+1][0]+1e-6).any(1)
                keep[c0:c0+500000]=kk
            cur=cur[keep]
            if len(cur)>MAXST: ok=False; break
            if len(cur)==0: break
        if not ok: rec['tests'][t]=None; continue
        if len(cur)==0: rec['tests'][t]=False; continue            # no achievable block compression
        cur=cur[:,corb]                                                          # expand to all classes
        X=cur.reshape((-1,)+tuple(Q+[t])).astype(float)
        F=np.fft.fftn(X,axes=tuple(range(1,len(Q)+2))); P=np.rint(np.abs(F)**2).reshape(len(cur),-1).astype(np.int64)
        a=cur.sum(1).astype(np.int64)
        # targets per character index (same raveled order as P)
        chi=np.array(np.unravel_index(np.arange(K),Q+[t])).T
        tgt=np.where(chi[:,-1]!=0,4*v,np.where((chi[:,:-1]!=0).any(1),4*b,4*v))
        good=(P[:,1:]<=tgt[None,1:]).all(1)&np.isin(a,Aok)
        P=P[good]; a=a[good]
        V=np.unique(np.concatenate([a[:,None],P[:,1:]],axis=1),axis=0)      # (a, P at nonzero characters)
        # need 4 rows with sum_k a_k^2 = 4v and sum P = tgt[1:]
        W=np.concatenate([V[:,:1]**2,V[:,1:]],axis=1)
        T=np.concatenate([[4*v],tgt[1:]])
        # chunked pair sums, pruned componentwise by the target (all entries are >= 0)
        parts=[]; tot=0; over=False
        for i0 in range(0,len(W),512):
            S=(W[i0:i0+512,None,:]+W[None,:,:]).reshape(-1,W.shape[1]); S=S[(S<=T[None,:]).all(1)]
            if len(S): S=np.unique(S,axis=0); parts.append(S); tot+=len(S)
            if tot>4e7: over=True; break
        if over: rec['tests'][t]=('too many pair sums',int(tot)); continue
        S2=np.unique(np.concatenate(parts),axis=0) if parts else np.zeros((0,W.shape[1]),np.int64)
        key=lambda A: A.view(np.dtype((np.void,A.dtype.itemsize*A.shape[1]))).ravel()
        S2c=np.ascontiguousarray(S2); Nd=np.ascontiguousarray(T[None,:]-S2)
        feas=bool(np.isin(key(Nd),key(S2c)).any()) if len(S2) else False
        rec['tests'][t]=bool(feas)
    out.append(rec); print(json.dumps(rec,default=int),flush=True)
