#!/bin/bash
# Second pass over the tier-1/2 keys that the first pass left undecided by MemoryError: one at a time, 12 GB cap,
# 90 minutes.  Waits for the first pass.  Results to logs/fable_fresh/tier12_stage1_pass2.jsonl.
cd /home/user/Hadamard
OUT=logs/fable_fresh/tier12_stage1_pass2.jsonl; LOG=logs/fable_fresh/tier12_stage1_pass2.log; touch $OUT
while ps -eo args | grep -E "fable_tier12_queue.sh" | grep -v grep > /dev/null; do sleep 60; done
python3 - > logs/fable_fresh/tier12_pass2_keys.txt <<'PY'
import json
seen = {}
for l in open('logs/fable_fresh/tier12_stage1.jsonl'):
    r = json.loads(l); seen[json.dumps(r['key'])] = r
for k, r in seen.items():
    if 'MemoryError' in json.dumps(r.get('B')) or 'rc=1 ' in json.dumps(r.get('B')): print(k)
PY
while IFS= read -r KEY; do
  [ -z "$KEY" ] && continue
  if grep -qF "$KEY" $OUT; then continue; fi
  echo "$(date -u +%F' '%T) start $KEY" >> $LOG
  OMP_NUM_THREADS=1 timeout 5400 nice -n 5 python3 theory/open_screen_stage1_exact.py "$KEY" --methods B --strong --planted 1 \
      --aslimit 12 --out $OUT >> $LOG 2>&1
  rc=$?; echo "$(date -u +%F' '%T) end $KEY rc=$rc" >> $LOG
  if ! grep -qF "$KEY" $OUT; then echo "{\"key\": $KEY, \"B\": {\"undecided\": \"rc=$rc (timeout 5400 s or memory 12 GB)\"}, \"phase\": \"capped\"}" >> $OUT; fi
done < logs/fable_fresh/tier12_pass2_keys.txt
echo "$(date -u +%F' '%T) pass2 done" >> $LOG
