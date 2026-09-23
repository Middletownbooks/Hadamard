"""Generate report/tables.md from logs."""
import json, glob, collections, numpy as np, os
out=[]
# W0
R=json.load(open('logs/census.json'))
out.append('## W0 census (all numbers measured)\n')
out.append('| v | b | m | G | d | r | row sums | special | PAF std/√m per block | mean | share var | log10 ∏p | null E[log10 ∏p] | log10 Fisher p | nbrs raw ≤2 | nbrs orbit ≤2 |')
out.append('|'+'---|'*16)
for r in R:
    out.append(f"| {r['v']} | {r['b']} | {r['m']} | {'×'.join('Z'+str(x) for x in r['group'])} | {r['d']} | {r['nb_orbit'][0] if r['nb_orbit'] else r['m']} | {tuple(r['rowsums'])} | {'yes' if r['special'] else 'no'} | {', '.join('%.3f'%b['ratio'] for b in r['blocks'])} | {r['mean_ratio']:.3f} | {r['share_var']:.4f} | {r['log10_rarity']} | {r['log10_null_expect']} | {r['log10_fisher_p']} | {sum(r['nb_raw'])} | {sum(r['nb_orbit'][1:]) if r['nb_orbit'] else '–'} |")
# W1
out.append('\n## W1 null ladder (measured; uniform solutions from pool-MITM, exhaustive at v=43)\n')
out.append('| v | d | n bits (Gaussian) | #solutions sampled | share var mean | sd | min | 0.0375(1−1.45d/m) | PAF ratio mean | P(log10 Fisher ≤ −2) (nominal 0.01) | KS D (T vs Dirichlet) |')
out.append('|'+'---|'*11)
nb={}
for f in ['logs/nbits.jsonl','logs/nbits_ladder.jsonl']:
    if os.path.exists(f):
        for l in open(f): x=json.loads(l); nb[(x['v'],x['d'])]=x['n_bits']
S=json.load(open('logs/w1_summary.json'))
for s in S:
    sols=[json.loads(l) for l in open(f"logs/w1_v{s['v']}_d{s['d']}.jsonl")]
    fi=np.array([z['fisher'] for z in sols])
    out.append(f"| {s['v']} | {s['d']} | {nb.get((s['v'],s['d']),'')} | {s['n']} | {s['sv_mean']:.4f} | {s['sv_sd']:.4f} | {s['sv_min']:.4f} | {0.0375*(1-1.45*s['d']/s['v']):.4f} | {s['ratio_mean']:.3f} | {(fi<=-2).mean():.3f} | {s['ks_stat']:.3f} |")
# W2
out.append('\n## W2 candidate runs (measured)\n')
out.append('| alg | v | d | n | params | solved/runs | median steps (all runs, cap counts as cap) | median secs | median share var of solutions | median PAF ratio |')
out.append('|'+'---|'*10)
rows=collections.defaultdict(list)
for f in glob.glob('logs/grid1.jsonl')+glob.glob('logs/w2_*.jsonl'):
    for l in open(f):
        try: r=json.loads(l)
        except: continue
        if r.get('kind')=='floor': continue
        alg=r.get('alg') or r.get('kind')
        if alg=='s3sched': alg='S3sched'
        keys=[k for k in ('lam','lam0','tau','lam_min','rho','winc','mu','share','inner','beta','eta','upd','tenure') if k in r]
        p=', '.join(f"{k}={r[k]}" for k in keys)
        rows[(alg,r['v'],r['d'],p)].append(r)
for (alg,v,d,p),rs in sorted(rows.items(),key=lambda t:(t[0][1],t[0][0],t[0][3])):
    ok=[r for r in rs if r.get('solved')]
    st=[r.get('steps',r.get('iters',0)) for r in rs]
    out.append(f"| {alg} | {v} | {d} | {nb.get((v,d),'')} | {p} | {len(ok)}/{len(rs)} | {int(np.median(st))} | {np.median([r['secs'] for r in rs]):.1f} | {('%.4f'%np.median([r['sv'] for r in ok])) if ok else '–'} | {('%.3f'%np.median([r['ratio'] for r in ok])) if ok else '–'} |")
# floors
out.append('\n## Single-block flatness floor (tabu, 400k steps, 4 seeds; measured)\n')
fl=collections.defaultdict(list)
for l in open('logs/grid1.jsonl'):
    r=json.loads(l)
    if r['kind']=='floor': fl[(r['v'],r['d'])].append(r['ratio'])
real={(c['v'],c['d']):c['mean_ratio'] for c in R}
out.append('| v | d | best floor ratio | median floor ratio | real core mean ratio | real/best floor |')
out.append('|'+'---|'*6)
for k in sorted(fl):
    x=np.array(fl[k]); rr=real.get(k)
    out.append(f"| {k[0]} | {k[1]} | {x.min():.3f} | {np.median(x):.3f} | {rr if rr else '–'} | {('%.2f'%(rr/x.min())) if rr else '–'} |")
os.makedirs('report',exist_ok=True); open('report/tables.md','w').write('\n'.join(out)+'\n'); print('\n'.join(out[-40:]))
