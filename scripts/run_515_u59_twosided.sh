#!/bin/bash
# two-sided yc stage on the remaining u=59 joint-feasible compressions (1..11), 2 at a time, full bucket budget
cd /home/user/Hadamard; mkdir -p logs/twosided
seq ${FIRST:-1} ${LAST:-11} | xargs -P2 -I{} bash -c "NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 -u scripts/bge3_twosided_run.py 85 5 59 5 data/config_screen_h85_p5_u59_pass.npy {} 2e7 100000 > logs/twosided/u59_sol{}.log 2>&1"
touch logs/twosided/done_u59
