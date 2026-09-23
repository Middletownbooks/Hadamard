"""Summarise W1 null runs: per v, share variance / PAF ratio / Fisher distribution of uniform solutions,
KS test of per-class T statistic against the Dirichlet(1,1,1,1) CDF."""
import sys, json, glob, numpy as np
sys.path.insert(0,'scripts')
from dirichlet_null import cdf
from scipy import stats
rows=[]
for f in sorted(glob.glob('logs/w1_v*_d*.jsonl'),key=lambda p:int(p.split('_v')[1].split('_')[0])):
    S=[json.loads(l) for l in open(f)]
    if not S: continue
    v,d=S[0]['v'],S[0]['d']
    sv=np.array([s['sv'] for s in S]); ra=np.array([s['ratio'] for s in S]); fi=np.array([s['fisher'] for s in S])
    T=np.concatenate([s['T'] for s in S]); U=cdf(T)
    ks=stats.kstest(U,'uniform')
    profs={}
    for s in S: profs.setdefault(tuple(s['prof']),[]).append(s['sv'])
    rows.append(dict(v=v,d=d,n=len(S),sv_mean=float(sv.mean()),sv_sd=float(sv.std()),sv_min=float(sv.min()),sv_p01=float(np.quantile(sv,0.01)),
                     ratio_mean=float(ra.mean()),fisher_median=float(np.median(fi)),fisher_min=float(fi.min()),
                     ks_stat=float(ks.statistic),ks_p=float(ks.pvalue),per_profile={str(k):(len(x),round(float(np.mean(x)),4)) for k,x in profs.items()}))
    r=rows[-1]
    print(f"v={v} d={d} n={r['n']} sv mean {r['sv_mean']:.4f} sd {r['sv_sd']:.4f} min {r['sv_min']:.4f} p01 {r['sv_p01']:.4f} | ratio {r['ratio_mean']:.3f} | log10 Fisher p median {r['fisher_median']:.2f} min {r['fisher_min']:.2f} | KS(T vs Dir) D={r['ks_stat']:.3f} p={r['ks_p']:.2g}")
    print('     per profile (n, mean sv):',r['per_profile'])
json.dump(rows,open('logs/w1_summary.json','w'),indent=1)
