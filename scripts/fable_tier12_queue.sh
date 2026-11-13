#!/bin/bash
# Exact stage 1 (method B, --strong, one planted control) of the tier-1 and tier-2 slices of the open-order screen
# that survive the axis-energy gate and have no exact decision yet (keys: logs/fable_fresh/tier12_keys.txt).
# Restart-safe (a key with a line in the out file is skipped).  Waits for the tier-0 queue (memory).
cd /home/user/Hadamard
OUT=logs/fable_fresh/tier12_stage1.jsonl
LOG=logs/fable_fresh/tier12_stage1.log
touch $OUT
while ps -eo args | grep -E "fable_tier0_queue.sh" | grep -v grep > /dev/null; do sleep 30; done
while IFS= read -r KEY; do
  [ -z "$KEY" ] && continue
  if grep -qF "$KEY" $OUT; then continue; fi
  echo "$(date -u +%F' '%T) start $KEY" >> $LOG
  OMP_NUM_THREADS=1 timeout 2700 nice -n 5 python3 theory/open_screen_stage1_exact.py "$KEY" --methods B --strong --planted 1 \
      --aslimit 7 --out $OUT >> $LOG 2>&1
  rc=$?
  echo "$(date -u +%F' '%T) end $KEY rc=$rc" >> $LOG
  if ! grep -qF "$KEY" $OUT; then echo "{\"key\": $KEY, \"B\": {\"undecided\": \"rc=$rc (timeout 2700 s or memory 7 GB)\"}, \"phase\": \"capped\"}" >> $OUT; fi
done < logs/fable_fresh/tier12_keys.txt
echo "$(date -u +%F' '%T) tier12 queue done" >> $LOG
