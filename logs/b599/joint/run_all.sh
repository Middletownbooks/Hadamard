#!/bin/bash
# exact joint stage 1 of the 598 doubly-staged b = 1 slices (theory/nd_joint_exhaust.py), plants first
cd /home/user/Hadamard
for mu in 49 31 25 45 51; do
  timeout 1800 python3 theory/nd_joint_exhaust.py 598 $mu 13 23 --plant 3 > logs/b599/joint/plant_598_$mu.log 2>&1
  timeout 14400 python3 theory/nd_joint_exhaust.py 598 $mu 13 23 --dir logs/b599/joint --workers 2 > logs/b599/joint/run_598_$mu.log 2>&1
done
