"""For each open order: feasible (log2E>0) b>=3 configurations from logs/bge3_cost_open, cheapest n first;
screen each with the exact quotient test (t=1, plus t=5 when 5 | h) until one survives.  usage: bge3_screen_ladder.py v [v ...]"""
import sys, json, glob, subprocess
for v in [int(x) for x in sys.argv[1:]]:
    rows={}
    for f in glob.glob(f'logs/bge3_cost_open/v{v}_*.log'):
        for l in open(f):
            if '"log2E"' not in l: continue
            r=json.loads(l)
            if r['log2E']<=0: continue
            key=(r['c'],tuple(r['Q']),r['d'],r['r'],r['L'])
            rows[key]=min(rows.get(key,1e9),r['n'])
    tried=0
    for (c,Q,d,r,L),n in sorted(rows.items(),key=lambda kv:kv[1]):
        h=(v+1)//c-1; ts=['1']+(['5'] if h%5==0 else [])
        out=subprocess.run(['timeout','1800','python3','scripts/bge3_screen.py',str(v),str(c),','.join(map(str,Q)),str(d),str(r),str(L)]+ts,capture_output=True,text=True).stdout
        res=[json.loads(l) for l in out.splitlines() if l.startswith('{"v"')]
        verdicts=[x['tests'] for x in res]
        anypass=any(all(val is not False for val in t.values()) for t in verdicts)
        decided=any(all(val is True for val in t.values()) for t in verdicts)
        print(json.dumps(dict(v=v,c=c,Q=Q,d=d,r=r,L=L,n=n,screen=verdicts,survives=anypass,decided_pass=decided)),flush=True)
        tried+=1
        if decided or tried>=8: break
