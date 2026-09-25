# validate enum_lifts_yc + exact join on 347: fixing the real yc must recover the real core
import sys,json,numpy as np; sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components, run
from bge3_twosided import yc_setup, enum_lifts_yc
g=setup(57,3,11); comps=components(g['footB'].T); foot=yc_setup(g,3); co=g['co']
Fr=np.array(np.load('data/template347.npz')['F']).astype(np.int64)
xc=np.zeros((4,g['mc']),np.int64)
for k in range(4): np.add.at(xc[k],g['pi'],Fr[k])
cell=((co[:,0]%2)*3+(co[:,1]%3))*3+(co[:,2]%3); yr=np.zeros((4,18),np.int64)
for k in range(4): np.add.at(yr[k],cell,Fr[k])
L=[enum_lifts_yc(g,comps,xc[k],yr[k],foot) for k in range(4)]
print('lifts per block with real yc',[len(x) for x in L],flush=True)
sols=run(57,3,11,xc,log=lambda s: print(s,flush=True),lifts=L)
print('solutions',len(sols),'real recovered',any(all(np.array_equal(s[k],Fr[k]) for k in range(4)) for s in sols))
