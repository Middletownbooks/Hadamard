import sys,json,numpy as np; sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_stageA import enum_block
from bge3_residue import residue_ok, joint_residue_ok, class_reps
L=open('data/rung_c6_v47_sym.txt').read().split('\n')
F=np.array([[1 if c=='+' else -1 for c in l.split('=')[1]] for l in L[1:5]])
g=setup(7,7,29); comps=components(g['footB'].T); co=g['co']; m=g['m']; v=g['v']; xc=np.zeros((4,g['mc']),np.int64)
for k in range(4): np.add.at(xc[k],g['pi'],F[k])
labc=g['labc']
# 1. stage-A completeness: is the known compression among the enumerated compressed blocks for its shift?
hit=None
for si,(s,cs) in enumerate(g['shifts']):
    if not all(np.array_equal(np.bincount((g['coc'][:,0]*3+g['coc'][:,1]),weights=xc[k],minlength=6).astype(int),np.array(cs[k])) for k in range(4)): continue
    ok=[]
    for k in range(4):
        Y,tr=enum_block(g,cs[k],v,7); X=Y[:,labc]; ok.append(bool((X==xc[k][None,:]).all(1).any()))
    hit=(si,[int(x) for x in s],ok)
print(json.dumps(dict(check='known compression in stage-A block lists',result=hit)))
# 2. filters
reps=class_reps(g,29,7)
print(json.dumps(dict(check='filters on known compression',residue_ok=bool(residue_ok(g,comps,reps,xc)[0]),joint16=bool(joint_residue_ok(g,comps,reps,xc,7,M=16)[0]),joint64=bool(joint_residue_ok(g,comps,reps,xc,7,M=64)[0]))))
# 3. assemble with the pipeline template code
ADD=g['ADD']; NEG=g['NEG']; Tm=np.load('data/template347.npz'); layout=json.loads(str(Tm['layout'])); N=4*v; off=20
Qc=(co[:,0]*3+co[:,1]).astype(np.int32); Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=Tm['C']
for J in range(4):
    Hm[:off,off+J*m:off+(J+1)*m]=Tm['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=Tm['Kc'][:,J*6+Qc].T
sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
for key,Lx in layout.items():
    I_,J_=int(key[0]),int(key[1]); sg=Lx['rho_sign']*Lx['gam_sign']
    Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[Lx['k']][sig[Lx['sigma']](Yg,Xg)]
P=sum(np.array([int(np.dot(F[k],F[k][ADD[:,gg]])) for gg in range(m)]) for k in range(4))
print(json.dumps(dict(check='known solution: PAF equations and pipeline assembly',paf_ok=bool(np.array_equal(P[1:],g['T_full'][1:])),hadamard=bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64))))))
