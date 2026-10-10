#!/bin/bash
# batch of 523 stage-1 solutions (fixed lag-0 compositions, several seeds), one at a time with 3 workers
cd /home/user/Hadamard
for seed in 11 12 13 14 15 16 17 18; do
  for comp in "13:6,15:4,27:0,29:0,1:62" "13:10,15:1,27:0,29:0,1:61" "13:2,15:7,27:0,29:0,1:63" "13:3,15:3,27:1,29:0,1:65"; do
    tag=$(echo $comp | tr ':,' '__')
    out=logs/q523/y523/y_${tag}_s${seed}.json
    [ -s $out ] && continue
    timeout 1300 python3 theory/q523_stage1.py 18 29 4 --comp $comp --time 1200 --workers 3 --seed $seed > $out 2>&1
  done
done
