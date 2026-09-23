"""W0 fingerprint census over all extracted cores.  Writes logs/census.json and prints a table."""
import sys, json, numpy as np, itertools, time
from math import lgamma, log
sys.path.insert(0,'scripts')
from hadlib import *
from dirichlet_null import cdf

def log10_binom_tail(k,n):
    # P(X>=k), X~Bin(n,1/2)
    from math import comb
    import math
    return math.log10(sum(comb(n,j) for j in range(k,n+1)))-n*math.log10(2)

def common_U(c):
    st=[dict(stabiliser(c,k)) for k in range(4)]
    U=sorted(set(st[0])&set(st[1])&set(st[2])&set(st[3]))
    return U,st

def normalise(c,U,st):
    """choose translation per block so each block is U-invariant about 0 where possible"""
    if len(U)<=1: return c.X.copy(),[0]*4
    co=c.coords(); Y=[];shifts=[]
    for k in range(4):
        best=None
        for t in range(c.m):
            y=c.X[k][c.idx(co+co[t])]
            if all(np.array_equal(y[c.idx(co*u)],y) for u in U): best=t;break
        shifts.append(best); Y.append(c.X[k] if best is None else c.X[k][c.idx(co+co[best])])
    return np.array(Y),shifts

def char_classes(c,U):
    """orbits of <U,-1> acting on nonzero characters (character coords scaled by u)"""
    co=c.coords(); seen=np.zeros(c.m,bool); seen[0]=True; reps=[]
    mults=set(U)|{1}; mults|={-u for u in mults}
    for i in range(1,c.m):
        if seen[i]: continue
        orb={int(c.idx(co[i]*u)) for u in mults}
        for j in orb: seen[j]=True
        real = int(c.idx(-co[i]))==i
        reps.append((i,len(orb),real))
    return reps

def neighbours_raw(c,Xn):
    """count solutions at Hamming distance 1 and 2 (arbitrary entry flips) preserving the PAF-sum target"""
    m=c.m; co=c.coords()
    T=c.idx(co[:,None,:]+co[None,:,:])   # T[i,s]=i+s
    Tm=c.idx(co[:,None,:]-co[None,:,:])  # i-s
    D=[]
    for k in range(4):
        x=Xn[k]
        d=-2*x[:,None]*(x[T]+x[Tm]); d[:,0]=0
        D.append(d)
    D=np.concatenate(D)   # (4m, m)
    n1=int(np.sum(~D[:,1:].any(1)))
    # different blocks: D_a = -D_b
    h={}
    for a in range(4*m): h.setdefault(D[a].tobytes(),[]).append(a)
    n2x=0
    for a in range(4*m):
        for bb in h.get((-D[a]).tobytes(),[]):
            if bb//m != a//m and bb>a: n2x+=1
    # same block
    n2s=0
    for k in range(4):
        x=Xn[k]; Dk=D[k*m:(k+1)*m]
        for i in range(m):
            S=Dk[i][None,:]+Dk[i+1:]          # pairs (i,j>i)
            nz=(S[:,1:]!=0).sum(1)
            for jj in np.where(nz<=2)[0]:
                j=i+1+jj; s=S[jj].copy()
                s[int(c.idx(co[j]-co[i]))]+=4*x[i]*x[j]; s[int(c.idx(co[i]-co[j]))]+=4*x[i]*x[j]
                if not s[1:].any(): n2s+=1
    return n1,n2x,n2s

def neighbours_orbit(c,Y,U):
    """solutions at distance 1,2 in the U-orbit space (flip whole orbits of the normalised blocks)"""
    co=c.coords(); orb=-np.ones(c.m,int); reps=[]
    for i in range(c.m):
        if orb[i]<0:
            o={int(c.idx(co[i]*u)) for u in U}
            for j in o: orb[j]=len(reps)
            reps.append(sorted(o))
    r=len(reps); M=np.zeros((r,c.m)); 
    for q,o in enumerate(reps): M[q,o]=1
    base=c.paf(Y); tgt=base.sum(0)
    def pafs(Z): return c.paf(Z)
    # single orbit flips
    D=[]
    for k in range(4):
        Z=np.repeat(Y[k][None,:],r,0)*(1-2*M).astype(int)
        D.append(pafs(Z)-base[k][None,:])
    n1=sum(int(np.sum(~d[:,1:].any(1))) for d in D)
    h={}
    for k in range(4):
        for q in range(r): h.setdefault(D[k][q].tobytes(),[]).append((k,q))
    n2x=sum(1 for k in range(4) for q in range(r) for (kk,qq) in h.get((-D[k][q]).tobytes(),[]) if kk>k)
    n2s=0
    for k in range(4):
        for q in range(r):
            Z=np.repeat((Y[k]*(1-2*M[q]).astype(int))[None,:],r-q-1,0)*(1-2*M[q+1:]).astype(int)
            if len(Z)==0: continue
            d=pafs(Z)-base[k][None,:]
            n2s+=int(np.sum(~d[:,1:].any(1)))
    return r,n1,n2x,n2s

def half_shift(c,Xn):
    if c.m%2 or len(c.mods)>1: return None
    h=c.m//2; out=[]
    for k in range(4):
        for l in range(4):
            agree=int(np.sum(Xn[k]==np.roll(Xn[l],-h)))
            a=max(agree,c.m-agree); out.append(('ABCD'[k]+'ABCD'[l], agree, round(float(log10_binom_tail(a,c.m)+np.log10(2)),1)))
    return out

rows=[]
for c in all_cores():
    t0=time.time()
    U,st=common_U(c); d=len(U)
    Y,shifts=normalise(c,U,st)
    R=c.paf(Y); P=c.psd(Y)
    tot=P.sum(0); S=P[:,1:]/tot[1:]            # shares per nonzero character
    share_var=float(np.mean((S-0.25)**2))
    pair={}
    for a,b in itertools.combinations(range(4),2):
        pair['ABCD'[a]+'ABCD'[b]]=round(float(np.var(S[a]+S[b])),4)
    reps=char_classes(c,U)
    lp=0.0; nreal=0
    for (i,sz,real) in reps:
        s=P[:,i]/P[:,i].sum(); t=float(((s-0.25)**2).sum())
        lp+=float(np.log10(cdf(t,0.5 if real else 1.0)[0])); nreal+=real
    from scipy.stats import chi2
    fisher=float(chi2.logsf(-2*lp*np.log(10),2*len(reps))/np.log(10))
    blocks=[]
    for k in range(4):
        p=R[k][1:]
        blocks.append(dict(paf_std=round(float(p.std()),2), ratio=round(float(p.std()/np.sqrt(c.m)),3),
                           maxabs=int(np.abs(p).max()), psd_range=[round(float(P[k,1:].min()/c.m),2),round(float(P[k,1:].max()/c.m),2)]))
    rs=np.abs(Y.sum(1)); special=bool(len(set(rs.tolist()))<=2 and max(np.bincount(rs))>=3)
    nb=neighbours_raw(c,Y)
    nbo=neighbours_orbit(c,Y,U) if d>1 else None
    row=dict(core=c.name,v=c.v,b=c.b,m=c.m,group=c.mods,U=U,d=d,shifts=shifts,rowsums=Y.sum(1).tolist(),special=special,
             blocks=blocks,mean_ratio=round(float(np.mean([bk['ratio'] for bk in blocks])),3),
             share_var=round(share_var,4),pair_var=pair,classes=len(reps),real_classes=nreal,log10_rarity=round(lp,1),log10_null_expect=round(-0.4343*len(reps),1),log10_fisher_p=round(fisher,1),
             nb_raw=nb,nb_orbit=nbo,half_shift=half_shift(c,Y),secs=round(time.time()-t0,1))
    rows.append(row); print(json.dumps(row),flush=True)
json.dump(rows,open('logs/census.json','w'),indent=1)
