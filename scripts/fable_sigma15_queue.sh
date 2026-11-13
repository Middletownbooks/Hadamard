#!/bin/bash
# sigma-invariant borders for the 2396 slice by clique search: NFIX = 7 (cross-check of the CP-SAT proof), 1, 4, 10.
cd /home/user/Hadamard
L=logs/fable_fresh
CASE=${1:-0}
for nf in 7 1 4 10; do
  out=$L/sigma15_q599_case${CASE}_f${nf}.log
  if grep -q "sigma-invariant border" $out 2>/dev/null; then continue; fi
  OMP_NUM_THREADS=1 nice -n 3 python3 theory/fable_border15.py q599 $CASE $nf 20000000000 > $out 2>&1
  echo "$(date -u +%T) sigma NFIX $nf case $CASE: $(tail -1 $out)" >> $L/qborder_queue.log
done
