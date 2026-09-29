#!/bin/bash
# positive control: v=131 and v=215 high-order slices, serial, 4 threads
cd /home/user/Hadamard
printf "21 3 11\n21 3 19\n21 7 5\n21 7 11\n21 7 19\n35 5 17\n35 5 23\n35 5 37\n35 5 73\n35 7 17\n35 7 23\n35 7 37\n35 7 73\n35 7 11\n35 7 19\n35 7 59\n35 7 79\n35 7 101\n35 7 149\n" | \
  xargs -P1 -L1 bash -c 'NUMBA_NUM_THREADS=2 MAXLIFT=1000000 LIFT_TMAX=2400 timeout 10800 python3 -u scripts/bge3_posctrl.py $0 $1 $2 > logs/posctrl/h$0_p$1_u$2.log 2>&1'
touch logs/posctrl/done_rungs
