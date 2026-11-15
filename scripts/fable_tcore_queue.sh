#!/bin/bash
# Exact core search (theory/fable_tcore.py) for every census record that passed G2 (and G3 or capped there).
# Restart-safe: keys already in logs/fable_fresh/tcore.jsonl are skipped.  Re-scans the census every pass.
cd /home/user/Hadamard
L=logs/fable_fresh; touch $L/tcore.jsonl
while true; do
  python3 - > $L/tcore_todo.txt <<'PY'
import json
done = set()
for l in open('logs/fable_fresh/tcore.jsonl'):
    try: done.add(json.loads(l)['key'])
    except Exception: pass
R = {}
for l in open('logs/fable_fresh/tcensus3.jsonl'):
    try: r = json.loads(l); R[r['key']] = r
    except Exception: pass
for k, r in R.items():
    if k in done: continue
    if r.get('alive') is True or (r.get('alive') is None and str(r.get('gate', '')).startswith('G3')):
        if r['r'] <= 48: print(k)
PY
  n=$(wc -l < $L/tcore_todo.txt)
  if [ "$n" -eq 0 ]; then
    if grep -q "tcensus3 done" $L/tcensus3.log; then echo "$(date -u +%T) tcore queue done" >> $L/tcore.log; break; fi
    sleep 300; continue
  fi
  while IFS= read -r KEY; do
    echo "$(date -u +%T) start $KEY" >> $L/tcore.log
    (ulimit -v 6000000; OMP_NUM_THREADS=2 timeout 3600 nice -n 3 python3 theory/fable_tcore.py "$KEY" --planted --maxlift 1e6 >> $L/tcore.log 2>&1)
    rc=$?
    python3 - "$KEY" "$rc" <<'PY2'
import json, sys
key, rc = sys.argv[1], sys.argv[2]
have = set()
for l in open('logs/fable_fresh/tcore.jsonl'):
    try: have.add(json.loads(l)['key'])
    except Exception: pass
if key not in have:
    open('logs/fable_fresh/tcore.jsonl', 'a').write(json.dumps(dict(key=key, verdict=f'UNDECIDED (rc={rc}, no record)')) + '\n')
PY2
  done < $L/tcore_todo.txt
done
