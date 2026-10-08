"""Summarise exhaustive slice sweeps (JSON lines from sweep_b1.py) into a markdown table."""
import json, glob, collections
rows = collections.OrderedDict()
for f in ['logs/lp_all_r30.log', 'logs/b0b1_cyclic_r30.log', 'logs/b1_537_Z2cu_Z67.log', 'logs/b1_669_Z2sq_Z167.log'] + sorted(glob.glob('logs/b1_nc_*.log')):
    try:
        lines = open(f).read().splitlines()
    except FileNotFoundError:
        continue
    for l in lines:
        if not l.startswith('{'): continue
        d = json.loads(l)
        key = (d['kind'], d['G'])
        R = rows.setdefault(key, dict(slices=0, sol=0, undecided=0, rmax=0, surv=0, errors=0))
        if 'error' in d:
            R['errors'] += 1; continue
        R['slices'] += 1; R['sol'] += d['nsol']; R['rmax'] = max(R['rmax'], d['r']); R['surv'] += d['surv']
        if d.get('undecided'): R['undecided'] += 1
print('| kind | group G | slices (r<=30) | max r | blocks surviving PSD+row-sum (total) | solutions | undecided (join cap) |')
print('|---|---|---|---|---|---|---|')
tot = collections.Counter()
for (k, G), R in rows.items():
    print(f"| {k} | {G} | {R['slices']} | {R['rmax']} | {R['surv']} | {R['sol']} | {R['undecided']} |")
    tot[k] += R['slices']
print(); print('totals per kind:', dict(tot))
