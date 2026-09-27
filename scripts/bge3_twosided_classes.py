"""Enumerate ALL yc 4-tuple classes of a compressed solution (two-sided route): run the yc join over every common
bucket, keeping one representative per class (per-block translation by multiples of h/p, which fixes xc and, when
(u-1)(h/p) = 0 mod h, u-invariance; lifts and solutions correspond bijectively within a class).  No tuple cap.
Writes data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{i}_classes.npy (rows = tuple indices into the saved yc lists).
usage: bge3_twosided_classes.py h p u p2 sol_index [nbuckets]"""
import sys, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup
from kt_generic import kt_join4
h,p,u,p2,si=[int(x) for x in sys.argv[1:6]]; nbk=int(sys.argv[6]) if len(sys.argv)>6 else 10**6
Z=np.load(f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}.npz'); Ys=[Z[f'Y{k}'] for k in range(4)]
g=setup(h,p,u); co=g['co']; tt=h//p; assert ((u-1)*tt)%h==0
cell=((co[:,0]%2)*3+(co[:,1]%3))*p2+(co[:,2]%p2)
co2=np.array([(a,b,z) for a in range(2) for b in range(3) for z in range(p2)]); idx2=lambda X:((X[...,0]%2)*3+(X[...,1]%3))*p2+(X[...,2]%p2)
T2=np.zeros(6*p2,np.int64); np.add.at(T2,cell,g['T_full'])
seen=set(); reps=[]
for x in range(6*p2):
    if x in seen: continue
    y=x; orb=set()
    while y not in orb: orb.add(y); orb.add(int(idx2(-co2[y]))); y=int(idx2(co2[y]*u))
    seen|=orb; reps.append(x)
reps=np.array(reps)
def paf(Y):
    F=Y.reshape(-1,2,3,p2).astype(float)
    return np.rint(np.fft.ifftn(np.abs(np.fft.fftn(F,axes=(1,2,3)))**2,axes=(1,2,3)).real).reshape(len(Y),-1).astype(np.int64)
PF=[paf(Y) for Y in Ys]; PAFl=[PF[k][:,reps]-(T2[reps] if k==0 else 0) for k in range(4)]
st=tt%p2
def canon(y): Y=np.asarray(y).reshape(6,p2); return min(tuple(np.roll(Y,j*st,axis=1).ravel().tolist()) for j in range(p2))
ccache=[{} for _ in range(4)]
def cblock(k,i):
    if i not in ccache[k]: ccache[k][i]=canon(Ys[k][i])
    return ccache[k][i]
classes={}; ntup=[0]; t0=time.time()
def verify(ix):
    tot=sum(PF[k][ix[k]] for k in range(4))
    if not np.array_equal(tot,T2): return False
    ntup[0]+=1; key=tuple(cblock(k,int(ix[k])) for k in range(4))
    if key in classes: return False
    classes[key]=tuple(int(i) for i in ix); return True
log=lambda s: print(s,flush=True) if ('kt_done' in s or 'kt_buckets' in s) else None
kt_join4(PAFl,verify,nbuckets=nbk,max_solutions=10**9,log=log,target_pairs=1e6)
R=np.array(list(classes.values()),np.int64).reshape(-1,4)
np.save(f'data/twosided_h{h}_p{p}_u{u}_q{p2}_sol{si}_classes.npy',R)
print(json.dumps(dict(stage='classes_done',sol=si,tuples=ntup[0],classes=len(R),secs=round(time.time()-t0,1))),flush=True)
