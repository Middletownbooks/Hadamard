#!/bin/bash
# positive control v=95: remaining slices, 2 at a time
cd /home/user/Hadamard
printf "5 17\n5 7\n3 7\n5 11\n3 11\n5 19\n3 19\n5 29\n3 29\n" | xargs -P2 -L1 bash -c 'NUMBA_NUM_THREADS=2 timeout 14400 python3 -u scripts/bge3_posctrl.py 15 $0 $1 > logs/posctrl/h15_p$0_u$1.log 2>&1'
touch logs/posctrl/done95
