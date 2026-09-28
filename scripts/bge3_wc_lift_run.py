"""Driver: fast exact (three-sided) lift over the yc classes of one compressed solution (v=515 two-sided route).
Resumable: classes already in the log are skipped.  Hits are assembled with 347's corner/border and checked
H H^T = 4v I; solutions saved to data/solution_v515_wc_u{u}_sol{s}_c{c}_{i}.txt.
usage: bge3_wc_lift_run.py sol_index CLASSES_FILE LOG [class indices csv or 'all']"""
import sys, os, json, time, numpy as np
sys.path.insert(0,'scripts')
from bge3_stageB_join import setup, components
from bge3_twosided import yc_setup
from bge3_wc_lift import wc_lift
si=int(sys.argv[1]); C=np.load(sys.argv[2]); logp=sys.argv[3]; sel=sys.argv[4] if len(sys.argv)>4 else 'all'
h,p,u,p2=85,5,161,5
g=setup(h,p,u); comps=components(g['footB'].T); foot=yc_setup(g,p2); co=g['co']; m=g['m']; v=g['v']
xc=np.load('data/config_screen_h85_p5_u161_pass.npy')[si]; Z=np.load(f'data/twosided_h85_p5_u161_q5_sol{si}.npz')
todo=list(range(len(C))) if sel=='all' else [int(x) for x in sel.split(',')]
done=set()
if os.path.exists(logp):
    for l in open(logp):
        try: done.add(json.loads(l)['class'])
        except Exception: pass
def assemble(F):
    Tm=np.load('data/template347.npz'); layout=json.loads(str(Tm['layout'])); N=4*v; off=20; ADD=g['ADD']; NEG=g['NEG']
    Qc=(co[:,0]*3+co[:,1]).astype(np.int32); Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=Tm['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=Tm['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=Tm['Kc'][:,J*6+Qc].T
    sig={'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    return bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64)))
logf=open(logp,'a')
for c in todo:
    if c in done: continue
    t=C[c]; sols,st=wc_lift(g,comps,foot,xc,[Z[f'Y{k}'][t[k]] for k in range(4)],u)
    rec=dict(sol=si,**{'class':int(c)},**st); logf.write(json.dumps(rec)+'\n'); logf.flush(); print(json.dumps(rec),flush=True)
    for i,F in enumerate(sols):
        ok=assemble(F); print(json.dumps(dict(stage='SOLUTION',sol=si,**{'class':int(c)},hadamard_verified=ok)),flush=True)
        with open(f'data/solution_v515_wc_u{u}_sol{si}_c{c}_{i}.txt','w') as fo:
            fo.write(f'v=515 N=2060 b=5 group=Z2xZ3xZ85 (mixed radix, last index fastest) u={u} compression={si} class={c}; built with 347 corner/border; hadamard_verified={ok}\n')
            for k,nm in enumerate('ABCD'): fo.write(nm+'='+''.join('+' if x>0 else '-' for x in F[k])+'\n')
