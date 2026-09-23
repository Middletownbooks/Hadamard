"""Normalise restricted-range WL (ln g on [lo,hi]) to the full-range WL curves by matching log2 counts on overlap bins."""
import json,sys,numpy as np
full=[json.load(open(f)) for f in sys.argv[1].split(',')]
rr=[json.load(open(f)) for f in sys.argv[2].split(',')]
ov_lo,ov_hi=float(sys.argv[3]),float(sys.argv[4])
def full_cum(rho):
    return np.mean([np.interp(rho,np.array(f['ratio_upper']),np.array(f['log2_count_le'])) for f in full])
out={}
for r in rr:
    ru=np.array(r['ratio_upper']); lng=np.array(r['lng']); vis=lng>0
    # log2 counts per bin up to a constant: lng/ln2 ; cumulative (log-sum-exp)
    l2=np.where(vis,lng/np.log(2),-np.inf)
    cum=np.log2(np.cumsum(2.0**(l2-l2[vis].max())))+l2[vis].max()
    sel=(ru>=ov_lo)&(ru<=ov_hi)
    off=np.mean([full_cum(x) for x in ru[sel]]-cum[sel])
    resid=np.std([full_cum(x) for x in ru[sel]]-cum[sel])
    c=cum+off
    print(f"seed {r['seed']}: flips {r['flips']} lnf {r['lnf']:.1e} secs {r['secs']} lowest visited {ru[vis][0]:.3f} overlap-fit resid {resid:.2f} bits")
    for rho in [0.28,0.29,0.30,0.31,0.32,0.34,0.36,0.38,0.40]:
        out.setdefault(rho,[]).append(float(np.interp(rho,ru,c)))
print('rho  log2 #blocks(ratio<=rho), per seed')
for rho,v in out.items(): print(rho,[round(x,1) for x in v])
