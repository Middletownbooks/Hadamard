#!/bin/bash
# re-run the census after the first pass (the a = 6 run was killed for memory at p = 719); restart-safe
cd /home/user/Hadamard
while ps -eo args | grep -E "fable_tcensus3_run.sh" | grep -v grep > /dev/null; do sleep 60; done
sed -i '/tcensus3 done/d' logs/fable_fresh/tcensus3.log
(ulimit -v 7000000; scripts/fable_tcensus3_run.sh)
