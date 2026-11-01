#!/bin/bash
# border sampler, round 2 (theory/zq_border_fixx.py, GS ansatz with X sampled): q = 7 at 1 h per pattern (seed 2),
# then the first 6 patterns of q = 9 (Z_9, Z_3^2) and q = 11 at 600 s each
cd /home/user/Hadamard
python3 theory/zq_border_fixx.py logs/zq/patterns_q7.json --tmax 3600 --tone 20 --seed 2 --out logs/zq/border_fixx_q7_s2.jsonl >> logs/zq/border_run.log 2>&1
python3 theory/zq_border_fixx.py logs/zq/patterns_q9.json --only 2,3,4,5 --tmax 600 --tone 20 --stop-first --out logs/zq/border_fixx_q9.jsonl >> logs/zq/border_run.log 2>&1
python3 theory/zq_border_fixx.py logs/zq/patterns_q9_33.json --only 0,1,2,3,4,5 --tmax 600 --tone 20 --stop-first --out logs/zq/border_fixx_q9_33.jsonl >> logs/zq/border_run.log 2>&1
python3 theory/zq_border_fixx.py logs/zq/patterns_q11.json --only 0,1,2,3,4,5 --tmax 600 --tone 30 --stop-first --out logs/zq/border_fixx_q11.jsonl >> logs/zq/border_run.log 2>&1
echo BORDERDONE >> logs/zq/border_run.log
