#!/bin/bash
# Re-run every core-search record whose latest verdict is UNDECIDED with the fixed code (window filter on the first
# y-quadruple too; streaming join).  Results append to logs/fable_fresh/tcore.jsonl (the latest record per key counts).
cd /home/user/Hadamard
L=logs/fable_fresh
python3 - > $L/tcore_redo_todo.txt <<'PY'
import json
last = {}
for l in open('logs/fable_fresh/tcore.jsonl'):
    try: r = json.loads(l); last[r['key']] = r
    except Exception: pass
for k, r in last.items():
    if r['verdict'].startswith('UNDECIDED'): print(k)
PY
while IFS= read -r KEY; do
  echo "$(date -u +%T) redo $KEY" >> $L/tcore_redo.log
  (ulimit -v 6000000; OMP_NUM_THREADS=1 timeout 5400 nice -n 4 python3 theory/fable_tcore.py "$KEY" --stream --maxlift 1e6 --maxpairs 5e7 >> $L/tcore_redo.log 2>&1)
done < $L/tcore_redo_todo.txt
echo "$(date -u +%T) redo done" >> $L/tcore_redo.log
