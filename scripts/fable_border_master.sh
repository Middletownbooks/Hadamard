#!/bin/bash
# After the case-0 queues: the same symmetric border searches for a representative with a non-constant Walsh
# profile on T (solution index 2160; the 32 constant profiles are equivalent to case 0 by signed block permutations).
cd /home/user/Hadamard
while ps -eo args | grep -E "fable_qborder_queue.sh|fable_border14.py" | grep -v grep > /dev/null; do sleep 30; done
scripts/fable_qborder_queue.sh 2160
while ps -eo args | grep -E "fable_sigma15_queue.sh|fable_border15.py" | grep -v grep > /dev/null; do sleep 30; done
scripts/fable_sigma15_queue.sh 2160
echo "$(date -u +%T) master done" >> logs/fable_fresh/qborder_queue.log
