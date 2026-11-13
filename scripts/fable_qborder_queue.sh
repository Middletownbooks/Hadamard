#!/bin/bash
# After the Q8 run: the three cyclic order-4 ansaetze for the 2396 slice, case 0, one at a time.
cd /home/user/Hadamard
L=logs/fable_fresh
while pgrep -f "fable_border13.py q599" > /dev/null; do sleep 20; done
for sub in z1 z2 z3; do
  out=$L/qborder_${sub}_q599_case0.log
  if grep -q "invariant border for" $out 2>/dev/null; then continue; fi
  nice -n 3 python3 theory/fable_border14.py sub=$sub q599 0 300000000 > $out 2>&1
done
echo "$(date -u +%T) qborder queue done" >> $L/qborder_queue.log
