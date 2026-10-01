#!/bin/bash
# Follow-on W(65) queue: shards whose Kaggle records were lost. Waits for the first local queue (PID $1) to exit.
cd /home/user/Hadamard
while kill -0 $1 2>/dev/null; do sleep 120; done
xargs -P4 -L1 bash -c 'o=logs/wil65/wil_n65_c$0_s$1of32.json; [ -f $o ] && exit 0; NUMBA_NUM_THREADS=1 WIL_CACHE_GB=2 PYTHONPATH=scripts python3 -u scripts/wil_search.py 65 5 13 --case $0 --shard $1/32 --out $o > logs/wil65/c$0_s$1.log 2>&1; echo "exit $?" >> logs/wil65/c$0_s$1.log' < logs/wil65/queue2.txt
