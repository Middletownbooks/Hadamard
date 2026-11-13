#!/bin/bash
# Second worker for the tier-0 stage-1 decisions: the last four keys in reverse order, one planted control each.
# Writes to the same out file, so scripts/fable_tier0_queue.sh skips what is done here.
cd /home/user/Hadamard
OUT=logs/fable_fresh/tier0_stage1.jsonl
LOG=logs/fable_fresh/tier0_stage1_w2.log
while IFS= read -r KEY; do
  [ -z "$KEY" ] && continue
  if grep -qF "$KEY" $OUT; then echo "$(date -u +%F' '%T) skip $KEY" >> $LOG; continue; fi
  echo "$(date -u +%F' '%T) start $KEY" >> $LOG
  OMP_NUM_THREADS=1 timeout 14400 nice -n 5 python3 theory/open_screen_stage1_exact.py "$KEY" --methods B --strong --planted 1 \
      --aslimit 7 --out $OUT >> $LOG 2>&1
  echo "$(date -u +%F' '%T) end $KEY rc=$?" >> $LOG
done <<'KEYS'
["k4", 1187, [4], [8], 37, 4, [["neg", [0, 1]]]]
["b1", 947, [], [11, 2], 43, 14, [[0, 10]]]
["b1", 2315, [], [13, 2], 89, 8, [[0, 6]]]
["b1", 1423, [], [3, 3, 2], 79, 6, [["neg", [0, 1]]]]
KEYS
echo "$(date -u +%F' '%T) worker 2 done" >> $LOG
