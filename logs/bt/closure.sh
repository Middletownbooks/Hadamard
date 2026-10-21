#!/bin/sh
# Exact closure pass [M]: next-gate costs and every exact gate up to 2^30 pairs, cheapest first, for each alive
# structure of the families judged "effectively closed" [E] (tier 1), then for the other alive structures with
# log2 E_Q < 0 or no finite estimate (tier 2).  Restart-safe (both runners resume from their output files).
cd /home/user/Hadamard
CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
[ -s logs/bt/closure_keys.txt ] || python3 logs/bt/closure_keys.py > logs/bt/closure_keys.txt
for tier in 1 2; do
  sed -n "${tier}p" logs/bt/closure_keys.txt > logs/bt/closure_keys_t$tier.txt
  nice -n 10 python3 theory/bt_cost.py logs/bt/cost_v2.jsonl $CEN --top 100000 --min-logE=-1e9 --keys @logs/bt/closure_keys_t$tier.txt
  nice -n 10 python3 theory/bt_gates_run.py logs/bt/cost_v2.jsonl logs/bt/gates_run_closure.jsonl $CEN --cap 30 --keys @logs/bt/closure_keys_t$tier.txt
  echo "closure tier $tier done"
done
