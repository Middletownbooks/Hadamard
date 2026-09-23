#!/bin/bash
cd /home/user/Hadamard
sleep 60
while ps -eo args | grep -q "[k]t_gen.py 523"; do sleep 10; done
sleep 5
KTCACHE=data/kt523/cache python3 -u scripts/kt_run.py 523 9 35,17,17,17 9 6.7e7 "data/kt523/gen_s*.npy" 2000 0 > logs/kt523_join.log 2>&1
