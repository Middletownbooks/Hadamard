"""v=359 (b=7, G = Z2 x Z2 x Z8 x Z11): stabiliser of the core under the FULL Aut(G) = Aut(Z2^2 x Z8) x (Z/11)*,
up to per-block translation (earlier census searched scalar multipliers only)."""
import sys, itertools, numpy as np; sys.path.insert(0,'scripts')
from hadlib import all_cores, translate_match
c=[x for x in all_cores() if x.v==359][0]; co=c.coords()        # coords (a,b,e,f): a,b in Z2, e in Z8, f in Z11
P=[(a,b,e) for a in range(2) for b in range(2) for e in range(8)]
def add(x,y): return ((x[0]+y[0])%2,(x[1]+y[1])%2,(x[2]+y[2])%8)
def mul(k,x):
    r=(0,0,0)
    for _ in range(k): r=add(r,x)
    return r
tors2=[x for x in P if mul(2,x)==(0,0,0)]
auts=[]
for i1 in tors2:
    for i2 in tors2:
        for i3 in P:
            img={}
            for (a,b,e) in P: img[(a,b,e)]=add(add(mul(a,i1),mul(b,i2)),mul(e,i3))
            if len(set(img.values()))==32: auts.append(img)
print('|Aut(Z2^2 x Z8)| =',len(auts))
units11=list(range(1,11))
stab={k:[] for k in range(4)}
for ai,A in enumerate(auts):
    Pimg=np.array([A[(int(x[0]),int(x[1]),int(x[2]))] for x in co[:,:3]])
    for u in units11:
        img=np.concatenate([Pimg,(co[:,3:4]*u)%11],axis=1); perm=c.idx(img)
        for k in range(4):
            if translate_match(c.X[k],c.X[k][perm],c): stab[k].append((ai,u))
common=set(stab[0])&set(stab[1])&set(stab[2])&set(stab[3])
print('per-block stabiliser sizes',[len(stab[k]) for k in range(4)],' common',len(common))
scal=[(ai,u) for ai,u in common if all(auts[ai][x]==mul(1,x) or True for x in P)]
# which common elements are scalar (x -> s x on the 2-part, u on Z11 with s = u mod 8 consistency not required here)
def is_scalar(A):
    for s in (1,3,5,7):
        if all(A[x]==mul(s,x) for x in P): return s
    return None
print('common elements (aut index, u, scalar-on-2-part?):',[(ai,u,is_scalar(auts[ai])) for ai,u in sorted(common)])
