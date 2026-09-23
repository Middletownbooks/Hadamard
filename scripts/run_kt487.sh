#!/bin/bash
cd /home/user/Hadamard
until [ $(ls data/kt487/gen_s*.npy 2>/dev/null | wc -l) -ge 3 ]; do sleep 10; done
sleep 5
python3 -u scripts/kt_run.py 487 9 37,17,17,1 9 6.7e7 "data/kt487/gen_s*.npy" 400 0 > logs/kt487_join.log 2>&1
