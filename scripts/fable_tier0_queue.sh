#!/bin/bash
# Exact stage-1 decision (method B, --strong) of the eight tier-0 slices of report/open_screen_answer.md that were
# listed "not run".  Restart-safe: a key with a line in the out file is skipped.  One job at a time (memory).
cd /home/user/Hadamard
OUT=logs/fable_fresh/tier0_stage1.jsonl
LOG=logs/fable_fresh/tier0_stage1.log
touch $OUT
while IFS= read -r KEY; do
  [ -z "$KEY" ] && continue
  if grep -qF "$KEY" $OUT; then echo "$(date -u +%F' '%T) skip $KEY" >> $LOG; continue; fi
  echo "$(date -u +%F' '%T) start $KEY" >> $LOG
  timeout 10800 nice -n 5 python3 theory/open_screen_stage1_exact.py "$KEY" --methods B --strong --planted 3 \
      --aslimit 11 --out $OUT >> $LOG 2>&1
  echo "$(date -u +%F' '%T) end $KEY rc=$?" >> $LOG
done <<'KEYS'
["k4", 1039, [4], [7], 37, 4, [["neg", [0, 1]]]]
["b1", 1427, [], [23, 2], 31, 10, [[0, 22]]]
["b1", 1259, [], [17, 2], 37, 12, [[0, 16]]]
["b1", 1893, [], [11, 2, 2], 43, 6, [[0, 10]]]
["b1", 1423, [], [3, 3, 2], 79, 6, [["neg", [0, 1]]]]
["b1", 2315, [], [13, 2], 89, 8, [[0, 6]]]
["b1", 947, [], [11, 2], 43, 14, [[0, 10]]]
["k4", 1187, [4], [8], 37, 4, [["neg", [0, 1]]]]
KEYS
echo "$(date -u +%F' '%T) queue done" >> $LOG
