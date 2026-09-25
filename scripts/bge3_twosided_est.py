"""Feasibility estimate for the two-sided lift (see bge3_twosided.py).  For one compressed solution xc (Z_t level),
sample lifts per block, keep PSD-passing ones, and compute their Z_p2 compression yc.  Reports per block the effective
number of yc cells C_eff = (sum n)^2 / sum n^2 (mass-weighted), and estimates q = P(four random PSD-passing yc satisfy
the Z_p2-compressed PAF equations) from per-class pair-sum distributions (independent classes).  The partitioned exact
search costs ~ C^2 q x (unpartitioned); worthwhile only if log2(C^2 q) << 0.
usage: bge3_twosided_est.py h p u p2 XCFILE sol_index S"""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_residue import lift_sampler
h,p,u,p2=[int(x) for x in sys.argv[1:5]]; XC=np.load(sys.argv[5]); si=int(sys.argv[6]); S=int(float(sys.argv[7]))
g=setup(h,p,u); comps=components(g['footB'].T); co=g['co']; v=g['v']; xc=XC[si]; rng=np.random.default_rng(si)
mods=np.array([2,3,h]); chars=co[co[:,2]!=0]; ph=np.exp(2j*np.pi*((co[:,None,:]*chars[None,:,:])/mods).sum(-1))
cell=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2)
# Z_p2-compressed PAF target and lag classes (under x -> -x and u)
co2=np.array([(a,b,z) for a in range(2) for b in range(3) for z in range(p2)]); idx2=lambda X:((X[...,0]%2)*3+(X[...,1]%3))*p2+(X[...,2]%p2)
T_full=g['T_full']; T2=np.zeros(6*p2,np.int64); np.add.at(T2,cell,T_full)
seen=set(); reps=[]
for x in range(6*p2):
    if x in seen: continue
    y=x; orb=set()
    while y not in orb: orb.add(y); orb.add(int(idx2(-co2[y]))); y=int(idx2(co2[y]*u))
    seen|=orb; reps.append(x)
reps=np.array(reps)
def paf2(Yc):
    F=Yc.reshape(-1,2,3,p2).astype(float)
    return np.rint(np.fft.ifftn(np.abs(np.fft.fftn(F,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(Yc),-1)[:,reps].astype(np.int64)
t0=time.time(); PA=[]; out=dict(h=h,p=p,u=u,p2=p2,sol=si,S=S,lag_classes=len(reps),blocks=[])
for k in range(4):
    Fs=[]
    for a in range(0,S,200000):
        F=lift_sampler(g,comps,xc[k],min(200000,S-a),rng).astype(float)
        ok=((np.abs(F@ph)**2)<=4*v+1e-6).all(1); Fs.append(F[ok])
    F=np.concatenate(Fs); Yc=np.zeros((len(F),6*p2)); np.add.at(Yc.T,cell,F.T)
    Yc=Yc.astype(np.int64); uq,cnt=np.unique(Yc,axis=0,return_counts=True)
    ceff=float(cnt.sum())**2/float((cnt.astype(float)**2).sum())
    PA.append(paf2(Yc)); out['blocks'].append(dict(psd_pass=int(len(F)),distinct_yc=int(len(uq)),C_eff_log2=round(float(np.log2(ceff)),1),
         singletons=int((cnt==1).sum())))
    print(json.dumps(dict(block=k,**out['blocks'][-1],secs=round(time.time()-t0,1))),flush=True)
tg=T2[reps]; ns=400000; lq=0.0
i=[rng.integers(0,len(P),ns) for P in PA]
AB=PA[0][i[0]]+PA[1][i[1]]-tg[None,:]; CD=-(PA[2][i[2]]+PA[3][i[3]])
for c in range(len(reps)):
    va,ca=np.unique(AB[:,c],return_counts=True); d=dict(zip(*np.unique(CD[:,c],return_counts=True)))
    lq+=np.log2(max(sum(n/ns*d.get(x,0)/ns for x,n in zip(va,ca)),1e-300))
out['log2_q_indep']=round(float(lq),1)
C=[b['C_eff_log2'] for b in out['blocks']]; out['log2_C2q_note']='compare sum over blocks of log2 C (pairing cost) with q'
out['log2_prod_C']=round(float(sum(C)),1)
print(json.dumps(out),flush=True)
