"""Two-sided yc stage for a compressed solution: sampled yc lists (z mod p2) + yc join; saves the yc lists and tuples.
usage: bge3_twosided_run.py h p u p2 XCFILE sol_index S [nbuckets]"""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_twosided import yc_lists, yc_join
h,p,u,p2=[int(x) for x in sys.argv[1:5]]; XC=np.load(sys.argv[5]); si=int(sys.argv[6]); S=int(float(sys.argv[7])); nbk=int(sys.argv[8]) if len(sys.argv)>8 else 2000
g=setup(h,p,u); comps=components(g['footB'].T)
log=lambda d: print(json.dumps(d) if isinstance(d,dict) else d,flush=True)
ycl=yc_lists(g,comps,XC[si],p2,S,np.random.default_rng(1000+si),log=log)
sols=yc_join(g,ycl,p2,u,log=log,nbuckets=nbk)
np.savez(f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}.npz',**{f'Y{k}':ycl[k][0] for k in range(4)},**{f'N{k}':ycl[k][1] for k in range(4)},tuples=np.array(sols,np.int64).reshape(-1,4))
log(dict(stage='yc_done',sol=si,yc_tuples=len(sols)))
