#!/bin/bash
# CP-SAT level gate (theory/zq_xcheck_run.py, 20 s per level and row pattern) on each order's survivors once its
# exact-level pass is done (logs/zq/gates/run.done).  Restart-safe.
cd /home/user/Hadamard
while true; do
  for tag in $(cut -d' ' -f1 logs/zq/gates/run.done 2>/dev/null); do
    grep -q "^$tag$" logs/zq/xcheck/run.done 2>/dev/null && continue
    python3 theory/zq_xcheck_run.py $tag --secs 20 >> logs/zq/xcheck/$tag.log 2>&1 && echo $tag >> logs/zq/xcheck/run.done
  done
  [ $(sort -u logs/zq/xcheck/run.done 2>/dev/null | wc -l) -ge 21 ] && break
  sleep 120
done
