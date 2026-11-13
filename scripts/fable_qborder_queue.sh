#!/bin/bash
# Symmetric border searches for the 2396 slice (fable_border14.py, with the symmetry reduction of the first row),
# one at a time: Q8 (cross-check of the unreduced run), then the three cyclic groups of order 4.
cd /home/user/Hadamard
L=logs/fable_fresh
CASE=${1:-0}
for sub in q8 z1 z2 z3; do
  out=$L/qborder14_${sub}_q599_case${CASE}.log
  if grep -q "invariant border for" $out 2>/dev/null; then continue; fi
  nice -n 3 python3 theory/fable_border14.py sub=$sub q599 $CASE 400000000 > $out 2>&1
  echo "$(date -u +%T) $sub case $CASE: $(tail -1 $out)" >> $L/qborder_queue.log
done
echo "$(date -u +%T) qborder queue done (case $CASE)" >> $L/qborder_queue.log
