#!/bin/bash
# wait for seed 2 of the 361 loop to finish, stop the loop (no further seeds), then run the 599 border job alone
cd /home/user/Hadamard
until [ -f /home/user/skew_work/p361_s2/DONE ] || grep -q "seed 2 exit" logs/fable_fresh/p361_loop.log; do sleep 30; done
pkill -f "loop361.sh" ; sleep 2; pkill -f "skew_kt.py theory/specs/p361_torus_a.json"
echo "$(date -u +%F' '%T) after361: seed 2 finished, loop stopped" >> logs/fable_fresh/p361_loop.log
nice -n 5 python3 theory/fable_border7.py logs/fable_fresh/q599_solutions.json 0 37 7200 1 > logs/fable_fresh/q599_border7_0.log 2>&1
echo "$(date -u +%F' '%T) border7 exit $?" >> logs/fable_fresh/q599_border7_0.log
