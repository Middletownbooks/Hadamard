#!/bin/sh
# Exact closure pass [M]: every exact level G/N (N invariant inside K, coarsest first, at most 300 s each; the pair
# bound of bt_cost.py overestimates badly, so levels are tried under a time limit instead of a bound) for each alive
# structure of the families judged "effectively closed" [E] (tier 1), then for the other alive structures with
# log2 E_Q < 0 or no finite estimate (tier 2).  Restart-safe (bt_gates_levels.py resumes from its output file).
cd /home/user/Hadamard
CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
[ -s logs/bt/closure_keys.txt ] || python3 logs/bt/closure_keys.py > logs/bt/closure_keys.txt
for tier in 1 2; do
  sed -n "${tier}p" logs/bt/closure_keys.txt > logs/bt/closure_keys_t$tier.txt
  nice -n 10 python3 theory/bt_gates_levels.py logs/bt/gates_run_closure.jsonl $CEN --levels auto --keys @logs/bt/closure_keys_t$tier.txt --tmax 300
  echo "closure tier $tier done"
done
