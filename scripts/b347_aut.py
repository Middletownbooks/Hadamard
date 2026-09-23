"""v=347 (b=5, G=Z2 x Z3 x Z57 = Z2 x Z3^2 x Z19): stabiliser under GL(2,3) x (Z/19)* up to per-block translation; affine orbits."""
import sys,itertools,numpy as np; sys.path.insert(0,'scripts')
from hadlib import all_cores, translate_match
c=[x for x in all_cores() if x.v==347][0]; co=c.coords()   # (a mod2, b mod3, e mod57)
z=co[:,0]; x3=co[:,1]; y3=co[:,2]%3; w=co[:,2]%19
def back(z,x3,y3,w):
    e=(y3*19*pow(19,-1,3)+w*3*pow(3,-1,19))%57
    return c.idx(np.stack([z%2,x3%3,e],1))
assert np.array_equal(back(z,x3,y3,w),np.arange(c.m))
GL=[np.array(M).reshape(2,2) for M in itertools.product(range(3),repeat=4) if (M[0]*M[3]-M[1]*M[2])%3]
units=range(1,19)
res={}
for gi,S in enumerate(GL):
    nx=(S[0,0]*x3+S[0,1]*y3)%3; ny=(S[1,0]*x3+S[1,1]*y3)%3
    for u in units:
        perm=back(z,nx,ny,(w*u)%19)
        ok=[bool(translate_match(c.X[k],c.X[k][perm],c)) for k in range(4)]
        if all(ok): res[(gi,u)]=perm
print('GL(2,3) size',len(GL),' common stabiliser order',len(res))
for k in range(4):
    parent=list(range(c.m))
    def f(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    for key,perm in res.items():
        t=co[translate_match(c.X[k],c.X[k][perm],c)[0]]
        sh=c.idx((co+t)%np.array(c.mods))
        for i in range(c.m):
            a,b=f(int(perm[i])),f(int(sh[i]))
            if a!=b: parent[a]=b
    roots={f(i) for i in range(c.m)}
    ok=all(len({c.X[k][i] for i in range(c.m) if f(i)==rt})==1 for rt in roots)
    print('block',k,'affine orbits',len(roots),'constant',ok)
