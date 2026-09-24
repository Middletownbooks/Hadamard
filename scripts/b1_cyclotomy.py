"""Cyclotomy test for the b=1 half-shift structure of v=167, 179 (h = m/2 prime).
For block k: S_k={j in Z_h: x(j)=x(j+h)}, u_k = x on S_k, w_k = x on the complement.
(1) support purity: for each subgroup H<(Z/h)*, |H|>=2, and translation t, purity(t,H) = fraction of points of Z_h\{t}
    lying in H-cosets (translated by t) that are entirely inside or entirely outside S_k. 1.0 = exact union of classes.
    Compared with the same statistic for random supports of the same size (null max over t).
(2) sign/character: for each multiplicative character chi of (Z/h)* of order e in {2,...} and t, agreement of sign(u_k(j))
    with eps*chi(j-t) on S_k (real characters only: e=2 Legendre), and same for w_k. Null: max over t,eps for random signs."""
import sys, numpy as np, json
sys.path.insert(0,'scripts')
from hadlib import all_cores
from redlib import primroot
rng=np.random.default_rng(0)
def subgroups(h):
    g=primroot(h); out=[]
    for dd in range(2,h):
        if (h-1)%dd==0:
            gen=pow(g,(h-1)//dd,h); out.append((dd,sorted({pow(gen,i,h) for i in range(dd)})))
    return out
def cosets(h,H):
    seen=set(); res=[]
    for a in range(1,h):
        if a in seen: continue
        c={a*x%h for x in H}; seen|=c; res.append(sorted(c))
    return res
def purity(Sset,h,cos,t):
    ins=np.zeros(h,bool); ins[[s for s in Sset]]=True
    good=0
    for c in cos:
        pts=[(x+t)%h for x in c]; v=ins[pts]
        if v.all() or (~v).all(): good+=len(c)
    return good/(h-1)
def legendre(a,h): return 0 if a%h==0 else (1 if pow(a,(h-1)//2,h)==1 else -1)
out=[]
for c in all_cores():
    if c.b!=1 or c.m not in (166,178): continue
    h=c.m//2; X=c.X
    subs=[(dd,H) for dd,H in subgroups(h) if dd<=h//4]
    for k in range(4):
        x=X[k]; S=[j for j in range(h) if x[j]==x[j+h]]
        res=dict(v=c.v,h=h,block='ABCD'[k],size=len(S))
        for dd,H in subs:
            cos=cosets(h,H)
            best=max(purity(S,h,cos,t) for t in range(h))
            nullv=[max(purity(list(rng.choice(h,len(S),replace=False)),h,cos,t) for t in range(h)) for _ in range(60)]
            res[f'purity_|H|={dd}']=(round(best,3),'null95',round(float(np.quantile(nullv,0.95)),3),'nullmax',round(float(max(nullv)),3))
        # Legendre sign test on u (support) and w (complement)
        chi=np.array([legendre(j,h) for j in range(h)])
        for name,idx in (('u',[j for j in range(h) if x[j]==x[j+h]]),('w',[j for j in range(h) if x[j]!=x[j+h]])):
            idx=np.array(idx); sg=x[idx]
            best=max(abs(int(np.sum(sg*chi[(idx-t)%h]))) for t in range(h))
            nl=[max(abs(int(np.sum(rng.choice([-1,1],len(idx))*chi[(idx-t)%h]))) for t in range(h)) for _ in range(200)]
            res[f'legendre_{name}']=(best,'of',len(idx),'null95',int(np.quantile(nl,0.95)))
        print(json.dumps(res),flush=True); out.append(res)
