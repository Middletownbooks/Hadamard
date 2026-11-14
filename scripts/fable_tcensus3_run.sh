#!/bin/bash
# Plan item 1: census of small general-T slices (a = 6, 5, 4; r <= 64), restart-safe.
cd /home/user/Hadamard
for a in 6 5 4; do
  OMP_NUM_THREADS=1 nice -n 4 python3 theory/fable_tcensus3.py --a $a --rmax 64 --qorbmax 12 --redo-caps >> logs/fable_fresh/tcensus3.log 2>&1
done
echo "$(date -u +%T) tcensus3 done" >> logs/fable_fresh/tcensus3.log
