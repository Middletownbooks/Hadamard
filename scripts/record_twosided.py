"""Append finished two-sided runs for multiplier U to data/twosided_515_u{U}_results.jsonl (idempotent)."""
import sys, json, glob, re, os
U=sys.argv[1]; out=f'data/twosided_515_u{U}_results.jsonl'
have=set()
if os.path.exists(out): have={json.loads(l)['sol'] for l in open(out)}
new=[]
for f in glob.glob(f'logs/twosided/u{U}_sol*.log'):
    m=re.search(r'_sol(\d+)\.log$',f)
    if not m: continue
    L=[json.loads(l) for l in open(f) if l.startswith('{')]
    done=[d for d in L if d.get('stage')=='yc_done']
    if not done or int(m.group(1)) in have: continue
    new.append(dict(sol=int(m.group(1)),yc_samples_per_block=2e7,yc_distinct=[d['distinct'] for d in L if d.get('stage')=='yc_list'],
                    buckets=next((d['n_common'] for d in L if d.get('stage')=='kt_buckets'),0),yc_tuples=done[-1]['yc_tuples']))
with open(out,'a') as fo:
    for d in sorted(new,key=lambda d:d['sol']): fo.write(json.dumps(d)+'\n'); print(json.dumps(d))
