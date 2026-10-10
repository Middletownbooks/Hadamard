#!/bin/bash
# batch of 523 stage-1 solutions with the free-composition CP-SAT model (scripts/cyclo_stage1_cpsat.py), 3 workers each
cd /home/user/Hadamard
for seed in 2 3 4 5 6 7 8 9 10 11 12 13; do
  out=logs/q523/y523/y_cpsat_s${seed}.json
  [ -s $out ] && continue
  timeout 1500 python3 scripts/cyclo_stage1_cpsat.py 18 29 4 --time 1400 --workers 3 --seed $seed > $out 2>&1
done
