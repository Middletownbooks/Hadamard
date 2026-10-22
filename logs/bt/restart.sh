#!/bin/sh
# Restart every live b3-push worker after a container restart (each runner resumes from its output file; the census
# runners cache their structure enumerations next to their output, so a restart skips the enumeration).
cd /home/user/Hadamard
CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
# census
nohup python3 theory/bt_census_run.py logs/bt/census/k4_P4_w8.jsonl --ns 643 --rmin 8 --rmax 160 --cap 4096 >> logs/bt/census/k4_P4_w8.log 2>&1 &
nohup sh -c "python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 739 --rmin 8 --rmax 160 --cap 16384 --min-cap 1024; python3 theory/bt_census_run.py logs/bt/census/k4_P4_w4.jsonl --ns 643 --rmin 8 --rmax 160 --cap 16384 --min-cap 4096" >> logs/bt/census/k4_P4_w4.log 2>&1 &
nohup python3 theory/bt_census2.py logs/bt/census2/k8c_839.jsonl --fam k8 --ns 839 --rmin 12 --rmax 110 --cap 16384 >> logs/bt/census2/k8c_839.log 2>&1 &
nohup nice -n 19 python3 theory/bt_census_split.py logs/bt/census/split_k4.jsonl logs/bt/census/k4_P4_w0.jsonl,logs/bt/census/k4_P4_w1.jsonl,logs/bt/census/k4_P4_w2.jsonl,logs/bt/census/k4_P4_w3.jsonl,logs/bt/census/k4_P4_w5.jsonl --rmin 16 --rmax 50 >> logs/bt/census/split_k4.log 2>&1 &
# screens [P] and exact gates [M]
nohup nice -n 10 python3 theory/bt_thm_rigid.py logs/bt/census/thm_rigid.jsonl $CEN >> logs/bt/census/thm_rigid.log 2>&1 &
nohup sh logs/bt/closure.sh >> logs/bt/closure.log 2>&1 &
nohup sh logs/bt/close_loop.sh > /dev/null 2>&1 &
