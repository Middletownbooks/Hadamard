"""Stage B on the stage-A compressed solutions for v=515 (h=85, p=5, u=53): list sizes, then join."""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, enum_lifts, components
XC=np.load('data/stageA_h85_p5_u53.npy')
g=setup(85,5,53); A=g['footB'].T
comps=components(A); print(json.dumps(dict(n_components=len(comps),max_component=max(len(c) for c in comps),orbits=g['r'])),flush=True)
for si,xc in enumerate(XC[:int(sys.argv[1]) if len(sys.argv)>1 else 1]):
    sizes=[]
    for k in range(4):
        tot=1
        for cq in comps:
            fis=np.nonzero(A[:,cq].sum(1))[0]; sub=A[np.ix_(fis,cq)]; n=len(cq)
            Y=((np.arange(1<<n)[:,None]>>np.arange(n)[None,:])&1)*2-1
            tot*=int((Y@sub.T==xc[k][fis][None,:]).all(1).sum())
        sizes.append(tot)
    print(json.dumps(dict(solution=si,lift_counts=sizes,log2=[round(float(np.log2(max(x,1))),1) for x in sizes])),flush=True)
