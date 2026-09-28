"""Collect fast exact lift results (logs/wc/sol*.log) into data/wc_lift_515_u161.jsonl."""
import json, glob
out=[]
for f in sorted(glob.glob('logs/wc/sol*.log')):
    for l in open(f):
        try: out.append(json.loads(l))
        except Exception: pass
with open('data/wc_lift_515_u161.jsonl','w') as fo:
    for d in out: fo.write(json.dumps(d)+'\n')
from collections import Counter
print(len(out),'classes', dict(Counter(d['sol'] for d in out)),'solutions',sum(d['solutions'] for d in out))
