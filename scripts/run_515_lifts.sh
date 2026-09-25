#!/bin/bash
cd /home/user/Hadamard
for i in 0 3 6 9; do echo $i; done | xargs -P2 -I{} bash -c 'XCFILE=data/stageA_515_feasible12.npy python3 -u scripts/bge3_stageB_sampled.py 85 5 53 {} 4e6 400 {} > logs/stageB_515_feas{}.log 2>&1'
echo LIFTS_DONE > logs/stageB_515_lifts_done
