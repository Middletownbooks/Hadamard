#!/bin/bash
# a = 6 is complete; a = 5 and a = 4 run as separate passes.  When both have ended, mark the census done
# (the core queue and the large-cap pass wait for this marker).
cd /home/user/Hadamard
while ps -eo args | grep -E "fable_tcensus3.py --a [45]" | grep -v grep > /dev/null; do sleep 120; done
echo "$(date -u +%T) tcensus3 done" >> logs/fable_fresh/tcensus3.log
