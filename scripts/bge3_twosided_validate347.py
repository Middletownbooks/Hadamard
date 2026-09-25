"""Validate the two-sided yc stage on v=347 (h=57=3*19, xc = z mod 19 from the real core, second side z mod 3):
sampled yc lists + yc join must contain the real core's yc tuple."""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_twosided import yc_lists, yc_join
S=int(float(sys.argv[1])) if len(sys.argv)>1 else 2e6
h,p,u,p2=57,3,11,3; g=setup(h,p,u); comps=components(g['footB'].T); co=g['co']
Fr=np.array(np.load('data/template347.npz')['F']).astype(np.int64)
xc=np.zeros((4,g['mc']),np.int64)
for k in range(4): np.add.at(xc[k],g['pi'],Fr[k])
cell=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2)
yreal=np.zeros((4,6*p2),np.int64)
for k in range(4): np.add.at(yreal[k],cell,Fr[k])
log=lambda d: print(json.dumps(d) if isinstance(d,dict) else d,flush=True)
ycl=yc_lists(g,comps,xc,p2,S,np.random.default_rng(0),log=log)
for k,(Y,N) in enumerate(ycl): log(dict(block=k,real_yc_in_list=bool((Y==yreal[k]).all(1).any())))
sols=yc_join(g,ycl,p2,u,log=log)
real=[any((ycl[k][0][s[k]]==yreal[k]).all() for k in range(4)) and all((ycl[k][0][s[k]]==yreal[k]).all() for k in range(4)) for s in sols]
log(dict(yc_tuples=len(sols),contains_real=any(real)))
