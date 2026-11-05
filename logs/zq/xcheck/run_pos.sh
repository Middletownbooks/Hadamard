#!/bin/bash
# CP-SAT levels at 600 s on the alive structures with log2 E_Q >= 0
cd /home/user/Hadamard
while read tag keys; do
  python3 theory/zq_xcheck_run.py $tag --secs 600 --keys $keys >> logs/zq/xcheck/$tag.log 2>&1
done < logs/zq/xcheck/positive_keys.txt
echo POSDONE >> logs/zq/xcheck/run.done
