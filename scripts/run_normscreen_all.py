"""Run bge3_normscreen over every feasible (log2E > 0) b>=3 config of the open orders, all proper divisors t of h."""
import json, glob, subprocess, sys
OPEN=[523, 571, 573, 599, 643, 647, 659, 669, 719, 739]
rows={}
for f in glob.glob('logs/bge3_cost_open/v*.log'):
    for l in open(f):
        if '"log2E"' not in l: continue
        r=json.loads(l)
        if r['v'] in OPEN and r['log2E']>0:
            key=(r['v'],r['c'],tuple(r['Q']),r['d'],r['r'],r['L'])
            rows[key]=min(rows.get(key,1e9),r['n'])
for (v,c,Q,d,r,L),n in sorted(rows.items(),key=lambda kv:(kv[0][0],kv[1])):
    h=(v+1)//c-1; ts=[str(t) for t in range(1,h) if h%t==0]
    out=subprocess.run(['timeout','1800','python3','scripts/bge3_normscreen.py',str(v),str(c),','.join(map(str,Q)),str(d),str(r),str(L)]+ts,capture_output=True,text=True)
    recs=[json.loads(x) for x in out.stdout.splitlines() if x.startswith('{"v"')]
    print(json.dumps(dict(v=v,c=c,Q=list(Q),d=d,r=r,L=L,n=n,h=h,groups=[x['tests'] for x in recs],
                          survives=any(all(t.values()) for t in (x['tests'] for x in recs)),
                          err=None if recs else out.stderr[-300:])),flush=True)
