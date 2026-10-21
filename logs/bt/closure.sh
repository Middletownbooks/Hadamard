#!/bin/sh
# Exact closure pass [M]: every exact level G/N (N invariant inside K, coarsest first, at most 300 s each; the pair
# bound of bt_cost.py overestimates badly, so levels are tried under a time limit instead of a bound) for each alive
# structure of the families judged "effectively closed" [E] (tier 1), then for the other alive structures with
# log2 E_Q < 0 or no finite estimate (tier 2).  The key lists are regenerated every round (the census keeps adding
# structures); a round that adds no record ends the pass.  Restart-safe (bt_gates_levels.py resumes from its output).
cd /home/user/Hadamard
for round in 1 2 3 4 5 6 7 8 9 10; do
  CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
  before=$(cat logs/bt/gates_run_closure.jsonl 2>/dev/null | wc -l)
  python3 logs/bt/closure_keys.py > logs/bt/closure_keys.txt
  for tier in 1 2; do
    sed -n "${tier}p" logs/bt/closure_keys.txt > logs/bt/closure_keys_t$tier.txt
    nice -n 10 python3 theory/bt_gates_levels.py logs/bt/gates_run_closure.jsonl $CEN --levels auto --keys @logs/bt/closure_keys_t$tier.txt --tmax 300
    echo "closure round $round tier $tier done"
  done
  after=$(cat logs/bt/gates_run_closure.jsonl | wc -l)
  [ "$after" -eq "$before" ] && break
done
