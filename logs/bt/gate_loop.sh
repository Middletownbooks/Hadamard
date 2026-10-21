#!/bin/sh
# Repeat cost pass + affordable gates (cap 36) until a round adds no gate record (kills promote new structures into
# the top 25 per order and family that bt_cost covers).
cd /home/user/Hadamard
for round in 1 2 3 4 5 6 7 8; do
  CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
  before=$(cat logs/bt/gates_run.jsonl | wc -l)
  nice -n 5 python3 theory/bt_cost.py logs/bt/cost_v2.jsonl $CEN --top 25 --min-logE -10
  python3 theory/bt_gates_run.py logs/bt/cost_v2.jsonl logs/bt/gates_run.jsonl $CEN --cap 36
  after=$(cat logs/bt/gates_run.jsonl | wc -l)
  echo "gate_loop round $round: $before -> $after gate records"
  [ "$after" -eq "$before" ] && break
done
