#!/bin/bash
# Atlas, continued: the order-4 cyclic runs for |T| = 7 need more than 11 GB (2.1e6 rows; killed), so they are
# skipped; the remaining runs with a 9 GB address-space cap.
cd /home/user/Hadamard
L=logs/fable_fresh/atlas; mkdir -p $L
echo "$(date -u +%T) t7new_z1/z2/z3: not run (out of memory at 11 GB on the 2.1e6-row problem)" >> $L/atlas.log
run() { name=$1; shift; out=$L/$name.log
  if grep -q "invariant border for" $out 2>/dev/null; then return; fi
  (ulimit -v 9000000; OMP_NUM_THREADS=1 timeout 10800 nice -n 4 python3 theory/fable_border14.py "$@" > $out 2>&1)
  echo "$(date -u +%T) $name: $(grep 'invariant border for' $out | tail -1 | cut -c1-200) $(grep -m1 -o 'MemoryError' $out)" >> $L/atlas.log; }
run t7flat_q8  sub=q8 synthetic 4 0,1,2,3,4,5,6 4,2,2,2 1000000000000
run t11a_q8s   sub=q8 mode=single synthetic 4 0,1,2,3,4,5,6,7,8,9,10 6,2,2,0 1000000000000
run t11b_q8s   sub=q8 mode=single synthetic 4 0,1,2,3,4,5,6,8,9,10,12 6,2,2,0 1000000000000
run t7new_q8s  sub=q8 mode=single synthetic 4 0,1,2,3,4,5,8 4,2,2,2 1000000000000
echo "$(date -u +%T) atlas done" >> $L/atlas.log
