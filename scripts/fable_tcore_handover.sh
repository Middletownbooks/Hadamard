#!/bin/bash
# When the core queue has no candidate left (the a = 5 census pass is very slow and adds few), stop the queue loop,
# mark it done so that the large-cap pass starts, and let the large-cap pass run alone (memory).
cd /home/user/Hadamard
L=logs/fable_fresh
while true; do
  n=$(python3 - <<'PY'
import json
last = set()
for l in open('logs/fable_fresh/tcore.jsonl'):
    try: last.add(json.loads(l)['key'])
    except Exception: pass
R = {}
for l in open('logs/fable_fresh/tcensus3.jsonl'):
    try: r = json.loads(l); R[r['key']] = r
    except Exception: pass
print(sum(1 for k, r in R.items() if k not in last and r['r'] <= 48 and (r.get('alive') is True or (r.get('alive') is None and str(r.get('gate', '')).startswith('G3')))))
PY
)
  if [ "$n" -eq 0 ] && ! ps -eo args | grep -E "^python3 theory/fable_tcore.py" | grep -v grep > /dev/null; then break; fi
  sleep 120
done
for p in $(ps -eo pid,args | grep -E "fable_tcore_queue.sh" | grep -v grep | awk '{print $1}'); do kill $p; done
echo "$(date -u +%T) tcore queue done (all current candidates; the a = 5 census is still adding slowly)" >> $L/tcore.log
