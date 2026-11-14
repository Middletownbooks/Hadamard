#!/bin/bash
# Plan item 2: single-character borders without any symmetry for the 2396 slice (solutions 0 and 2160).
cd /home/user/Hadamard
L=logs/fable_fresh
for c in 0 2160; do
  out=$L/qborder14_id_single_q599_case$c.log
  if grep -q "invariant border for" $out 2>/dev/null; then continue; fi
  OMP_NUM_THREADS=2 nice -n 3 python3 theory/fable_border14.py sub=id mode=single q599 $c 100000000000000 >> $out 2>&1
  echo "$(date -u +%T) id single case $c: $(grep 'invariant border for' $out | tail -1)" >> $L/qborder_queue.log
done
