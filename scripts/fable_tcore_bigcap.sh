#!/bin/bash
# Final core pass for the slices still UNDECIDED: 5e6 lifts per block, streaming join, small cache, 11 GB cap,
# one at a time, after the main core queue has finished.
cd /home/user/Hadamard
L=logs/fable_fresh
until grep -q "tcore queue done" $L/tcore.log 2>/dev/null; do sleep 300; done
python3 - > $L/tcore_bigcap_todo.txt <<'PY'
import json
last = {}
for l in open('logs/fable_fresh/tcore.jsonl'):
    try: r = json.loads(l); last[r['key']] = r
    except Exception: pass
for k, r in last.items():
    if r['verdict'].startswith('UNDECIDED'): print(k)
PY
while IFS= read -r KEY; do
  echo "$(date -u +%T) bigcap $KEY" >> $L/tcore_bigcap.log
  (ulimit -v 11000000; OMP_NUM_THREADS=2 timeout 10800 nice -n 4 python3 theory/fable_tcore.py "$KEY" --stream --maxlift 5e6 --maxpairs 5e7 --cache 4 >> $L/tcore_bigcap.log 2>&1)
done < $L/tcore_bigcap_todo.txt
echo "$(date -u +%T) bigcap done" >> $L/tcore_bigcap.log
