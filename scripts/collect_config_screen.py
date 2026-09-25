"""Collect config-screen outcomes from logs/config_screen/*.log into data/config_screen_515.jsonl (committed record)."""
import glob, json, re
out=[]
for f in sorted(glob.glob('logs/config_screen/h*_p*_u*.log')):
    h,p,u=map(int,re.findall(r'h(\d+)_p(\d+)_u(\d+)',f)[0]); L=open(f).read().splitlines()
    summ=[json.loads(l) for l in L if '"config_summary"' in l]
    if summ: out.append(summ[-1]); continue
    err=[l for l in L if 'Error' in l]
    if err: out.append(dict(stage='config_error',h=h,p=p,u=u,error=err[-1][:200],note='stage-A half-enumeration intractable at this compression level'))
with open('data/config_screen_515.jsonl','w') as fo:
    for d in out: fo.write(json.dumps(d)+'\n')
print(len(out),'records')
