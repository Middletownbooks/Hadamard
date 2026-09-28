"""Probe for a further exact partition of the two-sided lift lists: for a yc class, take each block's complete
PSD-passing lift list and compress it along a Q direction:
  zc: sum over the Z2 coordinate  -> sequence on Z3 x Z_h, values in {-2,0,2}
  wc: sum over the Z3 coordinate  -> sequence on Z2 x Z_h, values in {-3,-1,1,3}
Reports per block the list size, the number of distinct compressed vectors, and the effective number of cells
(sum n)^2/sum n^2.  The compressed PAF equations are implied by the full ones, so the lift join can be partitioned by
these cells exactly; the gain depends on these counts and on how many cell 4-tuples satisfy the compressed equations.
usage: bge3_threesided_probe.py h p u p2 XCFILE sol_index class_index"""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_twosided import yc_setup, enum_lifts_yc
h,p,u,p2=[int(x) for x in sys.argv[1:5]]; XC=np.load(sys.argv[5]); si=int(sys.argv[6]); ci=int(sys.argv[7])
g=setup(h,p,u); comps=components(g['footB'].T); foot=yc_setup(g,p2); co=g['co']; v=g['v']; lab=g['lab']
Z=np.load(f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}.npz'); C=np.load(f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}_classes.npy'); t=C[ci]
mods=np.array([2,3,h]); chars=co[co[:,2]!=0]; ph=np.exp(2j*np.pi*((co[:,None,:]*chars[None,:,:])/mods).sum(-1))
zcell=(co[:,1]%3)*h+co[:,2]; wcell=(co[:,0]%2)*h+co[:,2]
out=dict(sol=si,cls=ci,blocks=[]); t0=time.time()
LZ=[]; LW=[]
for k in range(4):
    Y=enum_lifts_yc(g,comps,XC[si][k],Z[f'Y{k}'][t[k]],foot); keep=[]
    for a in range(0,len(Y),50000):
        F=Y[a:a+50000][:,lab].astype(np.float32); keep.append(((np.abs(F@ph)**2)<=4*v+1e-3).all(1))
    Y=Y[np.concatenate(keep)]; F=Y[:,lab].astype(np.int16)
    zc=np.zeros((len(F),3*h),np.int16); np.add.at(zc.T,zcell,F.T); wc=np.zeros((len(F),2*h),np.int16); np.add.at(wc.T,wcell,F.T)
    res={}
    for nm,M in (('zc',zc),('wc',wc)):
        _,cnt=np.unique(M,axis=0,return_counts=True); res[nm]=dict(distinct=int(len(cnt)),C_eff_log2=round(float(np.log2(cnt.sum()**2/(cnt.astype(float)**2).sum())),1))
    out['blocks'].append(dict(block=k,psd_lifts=int(len(F)),**res)); print(json.dumps(out['blocks'][-1]),flush=True)
print(json.dumps(dict(stage='done',secs=round(time.time()-t0,1))),flush=True)
