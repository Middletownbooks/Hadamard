"""Method-indicating anomaly tests for the b>=2 cores (347 b=5, 359 b=7, 479 b=3).
1. H = {g : sum_k PAF_k(g) = -4b} u {0}; G/H structure; compressions xbar_k(c) = sum_{g in c} x_k(g).
2. PSD shares at the characters trivial on H (= compression spectrum) vs Dirichlet(1,1,1,1) null.
3. Coset-restriction relations: for each block, max over (scalar u, t in H, sign) of the agreement between
   its restrictions to two H-cosets; null = shuffle within cosets (keeps compressions).
4. Block-pair relations: max over (scalar u, t in G, sign) of agreement between blocks i and j;
   null = random blocks with the same row sums.
5. Per-coset constancy / balance counts.
Writes JSON lines."""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from hadlib import all_cores
rng=np.random.default_rng(0); NNULL=int(sys.argv[1]) if len(sys.argv)>1 else 200

def elem_order(c,g):
    co=c.coords()[g]; mods=np.array(c.mods); o=1
    for x,n in zip(co,mods): o=np.lcm(o,n//np.gcd(x,n))
    return int(o)

def add_table(c):
    co=c.coords(); return c.idx(co[:,None,:]+co[None,:,:])          # (m,m) g+h

def best_agree(x,y,maps,T):
    """max over maps u (perm arrays) and translations t in index set T of |sum x(u g + t)*y(g)| over g in dom"""
    return max(np.abs(np.array([np.dot(x[Pm[:,t]],y) for t in range(Pm.shape[1])])).max() for Pm in maps)

for c in all_cores():
    if c.b<2: continue
    m=c.m; S=c.paf().sum(0); H=np.array(sorted({0}|set(np.where(S==-4*c.b)[0].tolist())))
    AT=add_table(c); co=c.coords()
    isgrp=all(AT[h1,h2] in set(H.tolist()) for h1 in H for h2 in H)
    # cosets
    lab=-np.ones(m,int); reps=[]
    for g in range(m):
        if lab[g]<0: lab[AT[g,H]]=len(reps); reps.append(g)
    nq=len(reps)
    qord=[elem_order(c,g) for g in reps]                       # orders in G (upper bounds); quotient orders:
    def qorder(g):
        k=1; x=g
        while lab[x]!=lab[0]: x=AT[x,g]; k+=1
        return k
    qo=sorted(qorder(g) for g in reps)
    horders=sorted({elem_order(c,h) for h in H})
    X=c.X; comp=np.array([[int(X[k][lab==q].sum()) for q in range(nq)] for k in range(4)])
    rec=dict(v=c.v,b=c.b,group=c.mods,m=m,H_size=len(H),H_is_subgroup=isgrp,quotient_orders=qo,H_elem_orders=horders,
             rowsums=X.sum(1).tolist(),compressions=comp.tolist(),comp_norm=int((comp**2).sum()))
    # 2. PSD shares at characters trivial on H
    P=c.psd(); mods=np.array(c.mods)
    fr=co  # frequency coords
    triv=np.array([all((np.sum(fr[w]*co[h]*(np.lcm.reduce(mods)//mods))%np.lcm.reduce(mods)==0) for h in H) for w in range(m)])
    Pq=P[:,triv]; tot=Pq.sum(0); nz=tot>1e-9
    sh=(Pq[:,nz]/tot[nz]); rec['n_quot_freqs']=int(triv.sum()); rec['quot_share_var']=round(float(sh.var(0).mean()),4)
    rec['quot_shares']=np.round(sh.T,3).tolist()
    # 5. constancy: cosets where a block is constant, and exact-balance counts
    Hs=len(H); rec['coset_abs_sums_min']=int(np.abs(comp).min()); rec['coset_constant']=int((np.abs(comp)==Hs).sum())
    # 3. coset-restriction relations. Order each coset as rep+h for h in H (h in fixed order).
    units=c.units(); Hidx={h:i for i,h in enumerate(H)}
    Hmaps=[]
    for u in units:
        uh=c.idx(co[H]*u)
        if not all(x in Hidx for x in uh): continue
        # Pm[:,t] = index (within H order) of u*h + t_h
        Pm=np.array([[Hidx[AT[uh[i],H[t]]] for t in range(Hs)] for i in range(Hs)])
        Hmaps.append(Pm)
    # dedupe maps
    seen=set(); Hm=[]
    for Pm in Hmaps:
        key=Pm[:,0].tobytes()
        if key not in seen: seen.add(key); Hm.append(Pm)
    rec['n_H_multipliers']=len(Hm)
    def restr(x): return np.array([x[AT[reps[q],H]] for q in range(nq)])
    def coset_stat(x):
        Rr=restr(x); best=[]
        for q1 in range(nq):
            for q2 in range(q1+1,nq): best.append(best_agree(Rr[q1],Rr[q2],Hm,None))
        return max(best), best
    real=[coset_stat(X[k]) for k in range(4)]
    null=[]
    for _ in range(NNULL):
        k=rng.integers(4); y=X[k].copy()
        for q in range(nq):
            ii=np.where(lab==q)[0]; y[ii]=rng.permutation(y[ii])
        null.append(coset_stat(y)[0])
    null=np.array(null)
    rec['coset_rel']=[dict(block=k,max_agree=int(real[k][0]),p=round(float((null>=real[k][0]).mean()),3)) for k in range(4)]
    rec['coset_rel_null_mean']=round(float(null.mean()),1); rec['coset_rel_null_max']=int(null.max())
    # 4. block pair relations under scalar multipliers + G translations
    Gmaps=[]; seenG=set()
    for u in units:
        ug=c.idx(co*u); key=ug.tobytes()
        if key in seenG: continue
        seenG.add(key); Gmaps.append(AT[ug,:])   # Pm[g,t] = u g + t
    def pair_stat(x,y): return best_agree(x,y,Gmaps,None)
    pr=[]
    for i in range(4):
        for j in range(i+1,4): pr.append(((i,j),int(pair_stat(X[i],X[j]))))
    nullp=[]
    for _ in range(max(20,NNULL//10)):
        i,j=rng.choice(4,2,replace=False); xi=rng.permutation(X[i]); xj=rng.permutation(X[j]); nullp.append(pair_stat(xi,xj))
    nullp=np.array(nullp)
    rec['pair_rel']=[dict(pair=p,max_agree=a,p=round(float((nullp>=a).mean()),3)) for p,a in pr]
    rec['pair_rel_null_mean']=round(float(nullp.mean()),1); rec['n_G_multipliers']=len(Gmaps)
    print(json.dumps(rec),flush=True)
