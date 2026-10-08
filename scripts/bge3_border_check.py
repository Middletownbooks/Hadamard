"""Border (corner) compatibility of a b>=3 core's Q-compression.

Fact (verified exactly on the 347 template): with compressed core Cb (the 4c x 4c GS array of the coset sums xbar_k
over Q, same layout as the core), border rows K (4b x 4c) and border columns Kc (4b x 4c, stored like the template),
the corner is forced:   C = -K Cb^T Kc^T / (4c).
K, Kc must be +-1, rows mutually orthogonal and balanced inside each block (they span R^4 (x) 1_c^perp); then C is
automatically sqrt(4b)-orthogonal, so the ONLY extra condition is that every entry of K Cb^T Kc^T equals +-4c.
It depends only on the deviations d_k = xbar_k - mean, and on HOW the compression sits in the GS layout.

Modes
  c4          exhaustive over all c=4 (Q=Z2^2) compressions (3 single +-2 bumps + 1 constant block): every one has a border
  template347 verify the corner formula on data/template347.npz
  c6 N        c=6: (a) 347's compression; (b) 150 core-symmetry transforms of it, all with 347's frame; (c) sampled
              admissible compressions of each profile type with N frames derived from 347's K (random within-block
              perms/signs); Kc solved exactly by enumerating eps in {+-1}^20 (y' = Z^T eps / 20 must be +-1)
usage: bge3_border_check.py {c4|template347|c6 N}"""
import sys, json, itertools, numpy as np
T=np.load('data/template347.npz'); LAYOUT=json.loads(str(T['layout']))
def qgroup(mods):
    mods=np.array(mods); co=np.array(list(itertools.product(*[range(m) for m in mods])))
    def idx(X): X=np.asarray(X)%mods; return np.ravel_multi_index(tuple(np.moveaxis(X,-1,0)),tuple(mods))
    return co,idx
def compressed_core(xbar,mods):
    co,idx=qgroup(mods); c=len(co)
    sig={'h-g':lambda g,x: idx(x-g),'g+h':lambda g,x: idx(g+x),'-g-h':lambda g,x: idx(-g-x),'g-h':lambda g,x: idx(g-x)}
    G,X=np.meshgrid(np.arange(c),np.arange(c),indexing='ij'); M=np.zeros((4*c,4*c),np.int64)
    for key,L in LAYOUT.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        M[I_*c:(I_+1)*c,J_*c:(J_+1)*c]=sg*np.asarray(xbar[L['k']])[sig[L['sigma']](co[G],co[X])]
    return M
def clique(O,cands,k):
    if k==0: return []
    for t,i in enumerate(cands):
        nc=[j for j in cands[t+1:] if O[i,j]]
        if len(nc)>=k-1:
            r=clique(O,nc,k-1)
            if r is not None: return [i]+r
    return None
def border_c4(xbar):
    co,_=qgroup((2,2))
    chars=[np.array([(-1)**(a*x+bb*y) for (x,y) in co]) for (a,bb) in [(0,1),(1,0),(1,1)]]
    bv=[s*ch for ch in chars for s in (1,-1)]
    Y=np.array([np.concatenate(t) for t in itertools.product(bv,repeat=4)]); Y=Y[[i for i in range(len(Y)) if Y[i,0]==1]]
    O=(Y@Y.T==0); A=compressed_core(xbar,(2,2)).T; F=np.abs(Y@A@Y.T)==16
    sol=[None]
    def dfs(Kr,Kc_,cols):
        if sol[0] is not None: return
        if len(Kr)==12:
            cc=clique(O,cols,12)
            if cc is not None: sol[0]=(Kr[:],cc)
            return
        if clique(O,cols,12) is None: return
        for t,i in enumerate(Kc_):
            nk=[j for j in Kc_[t+1:] if O[i,j]]
            if len(nk)<12-len(Kr)-1: continue
            nc=[j for j in cols if F[i,j]]
            if len(nc)<12: continue
            Kr.append(i); dfs(Kr,nk,nc); Kr.pop()
    dfs([],list(range(len(Y))),list(range(len(Y))))
    if sol[0] is None: return None
    K=Y[sol[0][0]]; Kc=Y[sol[0][1]]; C=-(K@A@Kc.T)//16
    assert (np.abs(C)==1).all() and np.array_equal(C@C.T,12*np.eye(12,dtype=int))
    return K,Kc,C
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='template347':
        C=T['C']; K=T['K']; Kc=T['Kc']; Cb=compressed_core(T['comp'],(2,3))
        print(json.dumps(dict(max_abs_err=int(np.abs(-K@Cb.T@Kc.T/24-C).max()))))
    elif mode=='c4':
        ok=0; tot=0
        for const in range(4):
            for qs in itertools.product(range(4),repeat=3):
                for sg in itertools.product((1,-1),repeat=3):
                    xbar=np.zeros((4,4),np.int64); bi=0
                    for k in range(4):
                        if k==const: continue
                        xbar[k,qs[bi]]=2*sg[bi]; bi+=1
                    tot+=1; ok+=border_c4(xbar) is not None
        print(json.dumps(dict(c=4,compressions_tested=tot,with_border=ok)))
    elif mode=='c6':
        # (a) 347's own compression with 347's frame; (b) 150 core-symmetry transforms of 347's compression with 347's frame;
        # (c) sampled admissible compressions of every profile type, N random frames each (frames = 347's K with random
        #     within-block permutations and signs; Kc solved exactly).
        N=int(sys.argv[2]) if len(sys.argv)>2 else 20
        K0=T['K'].astype(np.int64); i=np.arange(2**19,dtype=np.int64)
        EPS=np.concatenate([np.ones((len(i),1),np.int64),1-2*((i[:,None]>>np.arange(19))&1)],axis=1)
        def solve(K,A):
            Yv=EPS@(K@A); ok=(np.abs(Yv)==20).all(1); V=Yv[ok]//20
            if len(V)<20: return False
            O=(V@V.T==0); return clique(O,list(range(len(V))),20) is not None
        rng=np.random.default_rng(0)
        def frame():
            K=K0.copy()
            for J in range(4):
                p=rng.permutation(6); s_=rng.choice([-1,1]); K[:,J*6:(J+1)*6]=s_*K[:,J*6:(J+1)*6][:,p]
            return K
        comp=T['comp']; mods=(2,3)
        print(json.dumps(dict(test='347 own compression, 347 frame',border=solve(K0,compressed_core(comp,mods).T))),flush=True)
        co=[(a,b_) for a in range(2) for b_ in range(3)]
        def tr(x,t,rev):
            s_=-1 if rev else 1
            return np.array([x[co.index(((s_*a+t[0])%2,(s_*b_+t[1])%3))] for (a,b_) in co])
        ok=0
        for it in range(150):
            perm=rng.permutation(4)
            x=np.array([int(rng.choice([-1,1]))*tr(comp[perm[k]],(int(rng.integers(2)),int(rng.integers(3))),bool(rng.integers(2))) for k in range(4)])
            ok+=solve(K0,compressed_core(x,mods).T)
        print(json.dumps(dict(test='core-symmetry transforms of 347 compression, 347 frame',with_border=int(ok),of=150)),flush=True)
        # admissible c=6 deviation patterns E (xbar = 1 + 2E), entries 0..2
        Qco=np.array(co); Fm=np.exp(-2j*np.pi*((Qco[:,None,:]*Qco[None,:,:])/np.array(mods)).sum(-1))
        cands=[]
        for E in itertools.product(range(3),repeat=6):
            E=np.array(E)
            if E.min()!=0: continue
            P=np.abs(E@Fm)**2
            if (P[1:]<=5+1e-9).all(): cands.append((E,np.rint(P[1:]*1e4)/1e4))
        from collections import defaultdict
        pair=defaultdict(list)
        for a in range(len(cands)):
            for b2 in range(a,len(cands)):
                sm=cands[a][1]+cands[b2][1]
                if (sm<=5+1e-6).all(): pair[tuple(np.round(sm,3))].append((a,b2))
        sols=set()
        for k2,lst in pair.items():
            need=tuple(np.round(5-np.array(k2),3))
            for (a,b2) in lst:
                for (c2,d2) in pair.get(need,[]): sols.add(tuple(sorted((a,b2,c2,d2))))
        byprof=defaultdict(list)
        for s4 in sols:
            byprof[tuple(sorted(int(round(6*((cands[t][0]-cands[t][0].mean())**2).sum())) for t in s4))].append(s4)
        for prof,lst in sorted(byprof.items()):
            found=0; tested=0
            for pi in rng.choice(len(lst),min(6,len(lst)),replace=False):
                xbar=np.array([cands[t][0] for t in lst[pi]])*2+1
                A=compressed_core(xbar,mods).T; tested+=1
                if any(solve(K0 if t==0 else frame(),A) for t in range(N)): found+=1
            print(json.dumps(dict(profile=prof,admissible_multisets=len(lst),tested=tested,frames_each=N,with_border=found)),flush=True)
