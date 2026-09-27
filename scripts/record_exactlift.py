"""Collect exact two-sided lift results (logs/twosided/lift_*.log) into data/twosided_515_u161_exactlift_all.jsonl."""
import json, glob, os
out=[]
for f in sorted(glob.glob('logs/twosided/lift_s*.log')):
    job=os.path.basename(f)[5:-4]
    for l in open(f):
        if '"stage": "tuple"' in l:
            d=json.loads(l); out.append(dict(u=161,job=job,class_index=d['tuple'],lift_sizes=d['lift_sizes'],solutions=d['solutions']))
with open('data/twosided_515_u161_exactlift_all.jsonl','w') as fo:
    for d in out: fo.write(json.dumps(d)+'\n')
print(len(out),'classes;',sum(d['solutions'] for d in out),'solutions')
