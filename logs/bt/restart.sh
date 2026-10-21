#!/bin/sh
# Restart every b3-push worker after a container restart (each runner resumes from its output file).
cd /home/user/Hadamard
CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
nohup python3 theory/bt_census_run.py logs/bt/census/k4_P4_w3.jsonl --ns 955 --rmin 8 --rmax 160 >> logs/bt/census/k4_P4_w3.log 2>&1 &
nohup python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 739,643 --rmin 8 --rmax 160 --cap 16384 >> logs/bt/census/k4_P4_w4.log 2>&1 &
nohup python3 theory/bt_census_run.py logs/bt/census/k4_P4_w5.jsonl --ns 751,599,719 --rmin 8 --rmax 160 >> logs/bt/census/k4_P4_w5.log 2>&1 &
nohup sh logs/bt/census2/queue2.sh 0 2 > /dev/null 2>&1 &
nohup sh logs/bt/census2/queue2.sh 1 2 > /dev/null 2>&1 &
nohup sh -c "python3 theory/bt_gates_run.py logs/bt/cost_next_gate.jsonl logs/bt/gates_run.jsonl $CEN --cap 34; python3 theory/bt_gates_run.py logs/bt/cost_next_gate.jsonl logs/bt/gates_run_extra.jsonl $CEN --cap 38 --keys 1bf02b0d94af1d4a" >> logs/bt/gates_run.log 2>&1 &
nohup nice -n 5 python3 theory/bt_cost.py logs/bt/cost_next_gate.jsonl $CEN --top 25 --min-logE -10 >> logs/bt/cost_next_gate2.log 2>&1 &
nohup nice -n 10 python3 theory/bt_thm_rigid.py logs/bt/census/thm_rigid.jsonl $CEN >> logs/bt/census/thm_rigid.log 2>&1 &
