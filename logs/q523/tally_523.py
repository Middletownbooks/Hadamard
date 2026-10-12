"""Tally the 523 production logs: S-solutions produced (S-join lines), D-checked (D-check lines, resumed D-checks
including the part done before the interruption, and the progress of live checkpoints)."""
import json, glob, os, math
prod = 0; checked = 0; cores = 0; joins = 0
for f in sorted(glob.glob('logs/q523/prod/p523*_w*.jsonl')):
    for l in open(f):
        d = json.loads(l)
        if d['stage'] == 'S-join': prod += d['s_solutions']; joins += 1
        if d['stage'] == 'D-check': checked += d['checked']
        # a resumed D-check also completes the checks its killed predecessor made (never logged by that process)
        if d['stage'] == 'D-check (resumed)': checked += d['start'] + d['checked']
        if d['stage'] == 'CORE' and d.get('ok'): cores += 1
live = 0
for f in glob.glob('logs/q523/prod/ckpt/p523*_c*.json'):
    live += json.load(open(f))['done']
tot = checked + live
# expected cores: model D_def = 26.14 (section 5) and the D-success measured on 1000 production S-solutions,
# 2^-26.02 (theory/q523_dsuccess_real.py)
print(json.dumps(dict(middle_values=joins, s_solutions_produced=prod, d_checked=tot, d_checked_log2=round(math.log2(max(tot, 1)), 2),
                      expected_cores=round(tot / 2 ** 26.14, 3), expected_cores_measured=round(tot / 2 ** 26.02, 3), cores=cores)))
