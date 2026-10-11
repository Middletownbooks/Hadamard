#!/bin/bash
# one batch of 523 stage-1 solutions per lag-0 composition (theory/q523_ymitm.py), for the composition screen
cd /home/user/Hadamard
for comp in "13:10,15:1,1:61" "13:2,15:7,1:63" "13:7,27:1,1:64" "13:3,15:3,27:1,1:65" "13:5,15:1,29:1,1:65" "13:1,15:4,29:1,1:66" "15:2,27:2,1:68" "13:2,27:1,29:1,1:68" "15:1,29:2,1:69"; do
  timeout 400 python3 theory/q523_ymitm.py 18 29 4 --comp $comp --N 20 --mid 3 --cents 3 --secs 300 --seed 3 --maxper 5 --out logs/q523/y523 > logs/q523/y523/log_$(echo $comp | tr ':,' '-_').jsonl 2>&1
done
