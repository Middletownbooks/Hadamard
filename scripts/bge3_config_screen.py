"""Configuration-level lift screen for the c=6 architecture (option 2): for G = Z2 x Z3 x Z_h, compression level p and
multiplier u, generate stage-A compressed solutions across shifts (PER_SHIFT per shift, A_TMAX seconds), then apply the
per-class and pairwise-joint forced-residue filters.  A configuration with compressions but no survivor, all failing
at a common lag signature, is reported dead [E: sampled filter].
usage: bge3_config_screen.py h p u nA     env PER_SHIFT (3), A_TMAX (1800), A_JTMAX (150, per-shift join), A_BUCKETS (200)"""
import sys, os, json, time, collections, numpy as np
sys.path.insert(0,'scripts')
os.environ.setdefault('PER_SHIFT','3'); os.environ.setdefault('A_TMAX','1800'); os.environ.setdefault('A_JTMAX','150'); os.environ.setdefault('A_BUCKETS','200')
from bge3_stageA import solve_A
from bge3_stageB_join import setup, components
from bge3_residue import residue_ok, joint_residue_ok, class_reps
h,p,u,nA=[int(x) for x in sys.argv[1:5]]; t0=time.time()
g=setup(h,p,u); comps=components(g['footB'].T); reps=class_reps(g,u,h)
quiet=lambda s: None if '"kt_progress"' in s or '"kt_setup"' in s or '"kt_buckets"' in s else print(s,flush=True)
sols=solve_A(h,p,u,log=quiet,max_solutions=nA)
res=collections.Counter(); sig=collections.Counter(); good=[]
for i,(s,xc) in enumerate(sols):
    ok,bad=residue_ok(g,comps,reps,xc,seed=i)
    if not ok:
        k='no_lift' if isinstance(bad,str) else 'per_class'; res[k]+=1; sig[(k,str(bad[:2]))]+=1; continue
    for M in (16,64):
        ok,info=joint_residue_ok(g,comps,reps,xc,h,ns=8000,M=M,seed=i)
        if not ok: res[f'joint{M}']+=1; sig[(f'joint{M}',str(info['first_bad'][:2]))]+=1; break
    if ok: res['pass']+=1; good.append(xc)
    print(json.dumps(dict(sol=i,shift=[int(x) for x in s],ok=bool(ok))),flush=True)
out=dict(stage='config_summary',h=h,p=p,u=u,r=int(g['r']),rc=int(g['rc']),n_shifts=len(g['shifts']),compressions=len(sols),
         shifts_with_sols=len({tuple(int(x) for x in s) for s,_ in sols}),outcome=dict(res),top_signatures=[[k[0],k[1],n] for k,n in sig.most_common(4)],
         secs=round(time.time()-t0,1))
print(json.dumps(out),flush=True)
if good: np.save(f'data/config_screen_h{h}_p{p}_u{u}_pass.npy',np.array(good))
