#!/bin/bash
cd /home/user/Hadamard; mkdir -p logs/q523/y235
for comp in "5:12,7:6,11:2,13:0,1:52" "5:10,7:7,11:2,13:0,1:53" "5:11,7:8,11:0,13:1,1:52" "5:9,7:10,11:1,13:0,1:52" "5:12,7:5,11:1,13:1,1:53"; do
  for seed in 2; do
    tag=$(echo $comp | tr ':,' '__'); out=logs/q523/y235/y_${tag}_s${seed}.json
    [ -s $out ] && continue
    timeout 1300 python3 theory/q523_stage1.py 18 13 4 --comp $comp --time 1200 --workers 1 --seed $seed > $out 2>&1
  done
done
