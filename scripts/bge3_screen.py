"""Exact obstruction screen for b>=3 configurations (G = Q x Z_h, H = {0} x Z_h, symmetry S from bge3_cost).
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
    for t in ts:
        if h%t: continue
        # class of each element in Q x Z_t
        cls=np.ravel_multi_index(tuple(np.moveaxis(np.concatenate([co[:,:len(Q)],co[:,-1:]%t],axis=1),-1,0)),Q+[t])
        K=c*t; M=np.zeros((K,r),np.int64); np.add.at(M,(cls,lab),1)
        # characters of Q x Z_t and their bounds: |F(chi)| <= sqrt(target) at the end
        chiv=np.array(np.unravel_index(np.arange(K),Q+[t])).T; mods_=np.array(Q+[t])
        E=np.exp(-2j*np.pi*((np.array(np.unravel_index(np.arange(K),Q+[t])).T[:,None,:]*chiv[None,:,:])/mods_).sum(-1))   # (class, chi)
        tg=np.where(chiv[:,-1]!=0,4*v,np.where((chiv[:,:-1]!=0).any(1),4*b,4*v)).astype(float)
        FM=(M.T.astype(float))@E                                              # orbit q -> contribution to F(chi)
        rem=np.vstack([np.abs(FM[q:]).sum(0) for q in range(r)]+[np.zeros(K)])  # max remaining |contribution|
        order=np.argsort(-np.abs(FM).sum(1)); FM=FM[order]; Mo=M[:,order]
        rem=np.vstack([np.abs(FM[q:]).sum(0) for q in range(r)]+[np.zeros(K)])
        cur=np.zeros((1,K),np.int32); ok=True
        for q in range(r):
            col=Mo[:,q].astype(np.int32); cur=np.ascontiguousarray(np.concatenate([cur+col,cur-col]))
            _,i=np.unique(cur.view(np.dtype((np.void,4*K))).ravel(),return_index=True); cur=cur[i]
            Fc=cur.astype(float)@E
            keep=(np.abs(Fc)-rem[q+1][None,:] <= np.sqrt(tg)[None,:]+1e-6).all(1)
            cur=cur[keep]
            if len(cur)>MAXST: ok=False; break
            if len(cur)==0: break
        if not ok: rec['tests'][t]=None; continue
        X=cur.reshape((-1,)+tuple(Q+[t])).astype(float)
        F=np.fft.fftn(X,axes=tuple(range(1,len(Q)+2))); P=np.rint(np.abs(F)**2).reshape(len(cur),-1).astype(np.int64)
        a=cur.sum(1).astype(np.int64)
        # targets per character index (same raveled order as P)
        chi=np.array(np.unravel_index(np.arange(K),Q+[t])).T
        tgt=np.where(chi[:,-1]!=0,4*v,np.where((chi[:,:-1]!=0).any(1),4*b,4*v))
        good=(P[:,1:]<=tgt[None,1:]).all(1)&(a*a<=4*v)
        P=P[good]; a=a[good]
        V=np.unique(np.concatenate([a[:,None],P[:,1:]],axis=1),axis=0)      # (a, P at nonzero characters)
        # need 4 rows with sum_k a_k^2 = 4v and sum P = tgt[1:]
        W=np.concatenate([V[:,:1]**2,V[:,1:]],axis=1)
        T=np.concatenate([[4*v],tgt[1:]])
        if len(W)>4000: rec['tests'][t]=('too many vectors',len(W)); continue
        S2=(W[:,None,:]+W[None,:,:]).reshape(-1,W.shape[1]); S2=S2[(S2<=T[None,:]).all(1)]
        S2=np.unique(S2,axis=0); need={tuple(x) for x in (T[None,:]-S2)}
        feas=any(tuple(x) in need for x in S2)
        rec['tests'][t]=bool(feas)
    out.append(rec); print(json.dumps(rec,default=int),flush=True)
