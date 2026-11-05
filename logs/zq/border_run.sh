#!/bin/bash
# border sampler (theory/zq_border_fixx.py): first find one bordered pattern per q, 600 s per pattern
cd /home/user/Hadamard
python3 theory/zq_border_fixx.py logs/zq/patterns_q7.json --tmax 1200 --tone 20 --out logs/zq/border_fixx_q7.jsonl >> logs/zq/border_run.log 2>&1
python3 theory/zq_border_fixx.py logs/zq/patterns_q9.json --tmax 600 --tone 20 --stop-first --out logs/zq/border_fixx_q9.jsonl >> logs/zq/border_run.log 2>&1
python3 theory/zq_border_fixx.py logs/zq/patterns_q9_33.json --tmax 600 --tone 20 --stop-first --out logs/zq/border_fixx_q9_33.jsonl >> logs/zq/border_run.log 2>&1
python3 theory/zq_border_fixx.py logs/zq/patterns_q11.json --tmax 600 --tone 30 --stop-first --out logs/zq/border_fixx_q11.jsonl >> logs/zq/border_run.log 2>&1
echo BORDERDONE >> logs/zq/border_run.log
