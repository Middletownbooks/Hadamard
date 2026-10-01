#!/bin/bash
# Local W(65) shard queue (K=32): 4 at a time; a shard is skipped if its JSON exists. Order: logs/wil65/queue.txt
cd /home/user/Hadamard
xargs -P4 -L1 bash -c 'o=logs/wil65/wil_n65_c$0_s$1of32.json; [ -f $o ] && exit 0; NUMBA_NUM_THREADS=1 WIL_CACHE_GB=2 PYTHONPATH=scripts python3 -u scripts/wil_search.py 65 5 13 --case $0 --shard $1/32 --out $o > logs/wil65/c$0_s$1.log 2>&1; echo "exit $?" >> logs/wil65/c$0_s$1.log' < logs/wil65/queue.txt
