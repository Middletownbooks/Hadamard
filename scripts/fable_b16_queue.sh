#!/bin/bash
# Borders of the 2396 slice (normal form) invariant under each of the five classes of subgroups of order 3.
cd /home/user/Hadamard
L=logs/fable_fresh
for w in d a r ra rd; do
  out=$L/b16_$w.log
  if grep -q "^summary" $out 2>/dev/null; then continue; fi
  while ps -eo args | grep -E "fable_border16.py" | grep -v grep > /dev/null; do sleep 20; done
  if grep -q "^summary" $out 2>/dev/null; then continue; fi
  OMP_NUM_THREADS=1 nice -n 3 python3 theory/fable_border16.py $w > $out 2>&1
  echo "$(date -u +%T) order-3 class $w: $(tail -1 $out)" >> $L/qborder_queue.log
done
