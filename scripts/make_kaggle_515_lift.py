"""Build kaggle/kt515_u161_sol5_lift.py: a self-contained, paste-into-one-cell notebook for the exact two-sided lift
on the 56 yc classes of v=515 (c=6, h=85, p=5, u=161) compression 5.  Embeds the three solver modules verbatim and
the small data files (base64)."""
import base64, io, json, sys, numpy as np
SOL=int(sys.argv[1]) if len(sys.argv)>1 else 5
mods={p:open(p).read() for p in ['scripts/bge3_lift.py','scripts/bge3_stageB_join.py','scripts/bge3_twosided.py']}
Z=np.load(f'data/twosided_h85_p5_u161_q5_sol{SOL}.npz'); C=np.load(f'data/twosided_h85_p5_u161_q5_sol{SOL}_classes.npy')
yc=np.stack([np.stack([Z[f'Y{k}'][C[i,k]] for k in range(4)]) for i in range(len(C))]).astype(np.int16)   # (56,4,30)
xc=np.load('data/config_screen_h85_p5_u161_pass.npy')[SOL].astype(np.int16)
def b64(**arrs):
    buf=io.BytesIO(); np.savez_compressed(buf,**arrs); return base64.b64encode(buf.getvalue()).decode()
T=np.load('data/template347.npz'); tmpl=base64.b64encode(open('data/template347.npz','rb').read()).decode()
job=b64(xc=xc,yc=yc)
out=f'''# =====================================================================================
# v=515 (Hadamard order 2060), c=6 architecture, G = Z2 x Z3 x Z85, multiplier u=161, compression {SOL}:
# EXACT two-sided lift of the {len(yc)} yc classes.  Paste this WHOLE file into ONE notebook cell and run.
#
# What it does, per class: enumerate ALL +-1 lifts of each block matching both compressions (z mod 17 from stage A,
# z mod 5 = the class), PSD-filter, exact 4-list join (batched).  A hit is assembled with 347's corner/border into the
# full 2060 x 2060 matrix and checked H H^T = 2060 I.  A class with 0 solutions is an exact null.
#
# EDIT ONLY THESE (just below the imports):
#   INSTANCE, NINST : run NINST notebooks in parallel with INSTANCE = 0, 1, ..., NINST-1 (same NINST in all).
#                     Instance i takes classes c with c % NINST == i.  0 <= INSTANCE < NINST is checked BEFORE any work.
#   CLASS_LIST      : optional explicit list, e.g. [0, 5, 9]; overrides INSTANCE/NINST.
#   RESUME_FROM     : None, or the path of a previous run's output folder added as an input dataset
#                     (e.g. '/kaggle/input/<dataset>/lift515'); classes already logged there are skipped.
#   SELFTEST        : True runs a 1-minute check on the known v=347 matrix first (must print SELFTEST PASS).
# Expect roughly 3-6 hours per class on 4 cores [estimate]; the run stops starting new classes after MAX_HOURS.
# OUTPUT (small): /kaggle/working/lift515/log_inst<I>.jsonl (one line per finished class) and, if any,
#   solutions_*.txt.  Every class result and every SOLUTION is also printed in the cell output.
# =====================================================================================
import os, sys, json, time, base64, io, warnings, numpy as np
warnings.filterwarnings('ignore', category=RuntimeWarning)   # uint64 hash wraparound is intentional

INSTANCE = 0          # <<< EDIT: 0 .. NINST-1
NINST = {min(4,len(yc))}             # <<< EDIT: number of parallel notebooks (same value in all of them)
CLASS_LIST = None     # <<< optional, e.g. [0, 7, 13]
RESUME_FROM = None    # <<< optional previous output folder
SELFTEST = True
MAX_HOURS = 11.0      # do not START a class after this many hours (Kaggle limit is 12 h)
JOIN_BATCH = 1e9      # stored pairs per batch in the exact join (~16 GB at 1e9; lower if memory is short)

WORK = os.environ.get('LIFT515_WORK','/kaggle/working/lift515_code'); OUT = os.environ.get('LIFT515_OUT','/kaggle/working/lift515_c{SOL}')
MODULES = {json.dumps(mods)}
TEMPLATE_B64 = "{tmpl}"
JOB_B64 = "{job}"

def _prepare():
    os.makedirs(os.path.join(WORK,'scripts'),exist_ok=True); os.makedirs(os.path.join(WORK,'data'),exist_ok=True); os.makedirs(OUT,exist_ok=True)
    for p,src in MODULES.items(): open(os.path.join(WORK,p),'w').write(src)
    open(os.path.join(WORK,'data/template347.npz'),'wb').write(base64.b64decode(TEMPLATE_B64))
    os.chdir(WORK); sys.path.insert(0,'scripts'); os.environ['JOIN_BATCH']=str(JOIN_BATCH)

def _assemble(g,F):
    co=g['co']; m=g['m']; v=g['v']; ADD=g['ADD']; NEG=g['NEG']
    Tm=np.load('data/template347.npz'); layout=json.loads(str(Tm['layout'])); N=4*v; off=20
    Qc=(co[:,0]*3+co[:,1]).astype(np.int32); Hm=np.zeros((N,N),np.int64); Hm[:off,:off]=Tm['C']
    for J in range(4):
        Hm[:off,off+J*m:off+(J+1)*m]=Tm['K'][:,J*6+Qc]; Hm[off+J*m:off+(J+1)*m,:off]=Tm['Kc'][:,J*6+Qc].T
    sig={{'h-g':lambda y,x: ADD[x,NEG[y]],'g+h':lambda y,x: ADD[y,x],'-g-h':lambda y,x: NEG[ADD[y,x]],'g-h':lambda y,x: ADD[y,NEG[x]]}}
    Yg,Xg=np.meshgrid(np.arange(m),np.arange(m),indexing='ij')
    for key,L in layout.items():
        I_,J_=int(key[0]),int(key[1]); sg=L['rho_sign']*L['gam_sign']
        Hm[off+I_*m:off+(I_+1)*m,off+J_*m:off+(J_+1)*m]=sg*F[L['k']][sig[L['sigma']](Yg,Xg)]
    return Hm, bool(np.array_equal(Hm@Hm.T,N*np.eye(N,dtype=np.int64)))

def _selftest():
    from bge3_stageB_join import setup, components, run
    from bge3_twosided import yc_setup, enum_lifts_yc
    g=setup(57,3,11); comps=components(g['footB'].T); foot=yc_setup(g,3); co=g['co']
    Fr=np.array(np.load('data/template347.npz')['F']).astype(np.int64); xc=np.zeros((4,g['mc']),np.int64)
    for k in range(4): np.add.at(xc[k],g['pi'],Fr[k])
    cell=((co[:,0]%2)*3+(co[:,1]%3))*3+(co[:,2]%3); yr=np.zeros((4,18),np.int64)
    for k in range(4): np.add.at(yr[k],cell,Fr[k])
    L=[enum_lifts_yc(g,comps,xc[k],yr[k],foot) for k in range(4)]
    sols=run(57,3,11,xc,log=lambda s: None,lifts=L) or []
    real=any(all(np.array_equal(s[k],Fr[k]) for k in range(4)) for s in sols)
    ok=real and all(_assemble(g,s)[1] for s in sols)
    print(('SELFTEST PASS' if ok else 'SELFTEST FAIL')+f': v=347 lifts {{[len(x) for x in L]}}, solutions {{len(sols)}}, real core recovered {{real}}, all assemble to Hadamard {{ok}}',flush=True)
    if not ok: raise SystemExit('SELFTEST FAILED - do not trust this environment')

def main():
    t0=time.time(); _prepare()
    J=np.load(io.BytesIO(base64.b64decode(JOB_B64))); xc=J['xc'].astype(np.int64); YC=J['yc'].astype(np.int64); ncl=len(YC)
    if CLASS_LIST is not None: mine=[int(c) for c in CLASS_LIST]
    else:
        if not (isinstance(NINST,int) and NINST>=1 and 0<=INSTANCE<NINST): raise SystemExit(f'ASSIGNMENT CHECK FAILED: need 0 <= INSTANCE < NINST (got INSTANCE={{INSTANCE}}, NINST={{NINST}})')
        mine=[c for c in range(ncl) if c%NINST==INSTANCE]
    print(f'ASSIGNMENT CHECK: {{len(mine)}} of {{ncl}} classes for this run: {{mine}}',flush=True)
    if not mine: raise SystemExit('ASSIGNMENT CHECK FAILED: no classes assigned')
    if SELFTEST: _selftest()
    from bge3_stageB_join import setup, components, run
    from bge3_twosided import yc_setup, enum_lifts_yc
    tag=('list' if CLASS_LIST is not None else f'inst{{INSTANCE}}of{{NINST}}'); logp=os.path.join(OUT,f'log_{{tag}}.jsonl')
    done=set()
    for d in [OUT]+([RESUME_FROM] if RESUME_FROM else []):
        if d and os.path.isdir(d):
            for f in os.listdir(d):
                if f.startswith('log_') and f.endswith('.jsonl'):
                    for l in open(os.path.join(d,f)):
                        try: done.add(json.loads(l)['class'])
                        except Exception: pass
    todo=[c for c in mine if c not in done]; print(f'already done (skipped): {{sorted(set(mine)&done)}}; to run: {{todo}}',flush=True)
    g=setup(85,5,161); comps=components(g['footB'].T); foot=yc_setup(g,5); cache={{}}; logf=open(logp,'a')
    for c in todo:
        if time.time()-t0>MAX_HOURS*3600: print('MAX_HOURS reached: not starting further classes',flush=True); break
        tc=time.time(); L=[]
        for k in range(4):
            key=(k,tuple(YC[c,k].tolist()))
            if key not in cache: cache[key]=enum_lifts_yc(g,comps,xc[k],YC[c,k],foot)
            L.append(cache[key])
        print(f'class {{c}}: lift list sizes {{[len(x) for x in L]}} (built in {{time.time()-tc:.0f}} s); exact join running...',flush=True)
        sols=run(85,5,161,xc,log=lambda s: print('   ',s,flush=True) if ('batched' in s or 'stored_pairs' in s) else None,lifts=L) or []
        rec=dict(v=515,u=161,compression={SOL},**{{'class':int(c)}},lift_sizes=[int(len(x)) for x in L],solutions=len(sols),hours=round((time.time()-tc)/3600,2))
        logf.write(json.dumps(rec)+'\\n'); logf.flush(); print('CLASS RESULT '+json.dumps(rec),flush=True)
        for i,F in enumerate(sols):
            Hm,ok=_assemble(g,F); s=dict(v=515,u=161,compression={SOL},**{{'class':int(c)}},hadamard_verified=ok,group='Z2xZ3xZ85 mixed radix, last index fastest',
                  blocks={{nm:''.join('+' if x>0 else '-' for x in F[k]) for k,nm in enumerate('ABCD')}})
            print('SOLUTION '+json.dumps(s),flush=True)
            open(os.path.join(OUT,f'solution_class{{c}}_{{i}}.txt'),'w').write(json.dumps(s)+'\\n')
    print(f'RUN SUMMARY: finished {{sum(1 for _ in open(logp))}} class lines in {{logp}}; hours {{(time.time()-t0)/3600:.2f}}',flush=True)

if __name__=='__main__' or 'ipykernel' in sys.modules:
    main()
'''
open(f'kaggle/kt515_u161_sol{SOL}_lift.py','w').write(out); print('written',len(out)//1024,'KB')
