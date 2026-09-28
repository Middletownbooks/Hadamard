#!/bin/bash
# fast exact lift over all remaining u=161 classes; jobs in logs/wc/jobs.txt ("sol classes_file log class_csv"), 2 at a time
cd /home/user/Hadamard
cat logs/wc/jobs.txt | xargs -P2 -L1 bash -c 'NUMBA_NUM_THREADS=2 python3 -u scripts/bge3_wc_lift_run.py "$0" "$1" "$2" "$3" >> logs/wc/out_sol$0.txt 2>&1'
touch logs/wc/done
