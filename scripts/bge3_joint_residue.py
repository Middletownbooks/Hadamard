"""Joint forced-residue filter on compressed solutions.  usage: bge3_joint_residue.py 347  |  bge3_joint_residue.py 515 XCFILE
347: validates on the real compression (must pass).  515: tests each compressed solution in XCFILE."""
import sys, json, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_residue import class_reps, joint_residue_ok
if sys.argv[1]=='347':
    h,p,u=57,3,11; g=setup(h,p,u); Fr=np.array(np.load('data/template347.npz')['F']).astype(np.int64)
    xc=np.zeros((4,g['mc']),np.int64)
    for k in range(4): np.add.at(xc[k],g['pi'],Fr[k])
    XC=[xc]
else:
    h,p,u=85,5,53; g=setup(h,p,u); XC=np.load(sys.argv[2])
comps=components(g['footB'].T); reps=class_reps(g,u,h)
for i,xc in enumerate(XC):
    for M in (16,32,64,128):
        ok,info=joint_residue_ok(g,comps,reps,xc,h,M=M,seed=i)
        print(json.dumps(dict(v=sys.argv[1],sol=i,M=M,ok=bool(ok),**info)),flush=True)
