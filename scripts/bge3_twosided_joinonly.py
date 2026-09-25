"""Re-run the yc join on saved yc lists (data/twosided_*.npz) with a larger bucket budget.
usage: bge3_twosided_joinonly.py h p u p2 sol_index nbuckets"""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup
from bge3_twosided import yc_join
h,p,u,p2,si,nbk=[int(x) for x in sys.argv[1:7]]
f=f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}.npz'; Z=np.load(f); ycl=[(Z[f'Y{k}'],Z[f'N{k}']) for k in range(4)]
g=setup(h,p,u); log=lambda d: print(json.dumps(d) if isinstance(d,dict) else d,flush=True)
sols=yc_join(g,ycl,p2,u,log=log,nbuckets=nbk,max_solutions=1000)
np.savez(f,**{k:Z[k] for k in Z.files if k!='tuples'},tuples=np.array(sols,np.int64).reshape(-1,4))
log(dict(stage='yc_done',sol=si,yc_tuples=len(sols),nbuckets=nbk))
