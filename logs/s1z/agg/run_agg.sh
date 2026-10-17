#!/bin/bash
# aggregated cheap-level compression screens over the second-tier configurations (restart-safe per multiset)
cd /home/user/Hadamard
python3 - <<'PY'
import json, subprocess, sys
cfgs = json.load(open('logs/s1z/ts/configs.json'))
skip = {'739_c4_q23_s2_H8_U[]', '719_c6_q17_s2_H7_U[[neg,[1]]]'}
for o in cfgs:
    if o['tag'] in skip: continue
    lv = o['levels'][:-1]   # all but the full level
    out = f"logs/s1z/agg/agg_{o['tag']}.jsonl"
    subprocess.run([sys.executable, 'theory/s1z_zqscreen.py', json.dumps(o['cfg']), out, '--levels', json.dumps(lv)],
                   stdout=open(out.replace('.jsonl', '.out'), 'a'), stderr=subprocess.STDOUT)
PY
