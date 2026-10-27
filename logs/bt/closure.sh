#!/bin/sh
# Exact closure pass [M]: every exact level G/N (N invariant inside K, coarsest first, under a time limit; the pair
# bound of bt_cost.py overestimates badly, so levels are tried under a time limit instead of a bound) for each alive
# structure of the families judged "effectively closed" [E] (tier 1), then for the other alive structures with
# log2 E_Q < 0 or no finite estimate (tier 2).  Each tier runs in three passes:
#   - 30 s per level, the key list split over W runners (--part i/W), each under ulimit -v MEM (98% of the tier-1
#     kills took under 30 s, and a level capped by time costs the whole limit);
#   - 300 s, same split, only on the levels the first pass could not finish (--retry-only);
#   - 300 s, one runner under ulimit -v MEM1, only on the levels refused an allocation under MEM.
# The runners share one memory cgroup (13.4 GiB); on Oct 27 one of three unlimited runners was OOM-killed at 7.0 GB.
# The key lists are regenerated every round (the census keeps adding structures); a round that adds no record ends
# the pass.  Restart-safe (bt_gates_levels.py resumes from its output).
cd /home/user/Hadamard
W=3; MEM=3700000; MEM1=11500000                     # KB
for round in 1 2 3 4 5 6 7 8 9 10; do
  CEN=$(ls logs/bt/census/k4_P4_w*.jsonl logs/bt/census2/*.jsonl | tr '\n' ',' | sed 's/,$//')
  before=$(cat logs/bt/gates_run_closure.jsonl 2>/dev/null | wc -l)
  python3 logs/bt/closure_keys.py > logs/bt/closure_keys.txt
  for tier in 1 2; do
    sed -n "${tier}p" logs/bt/closure_keys.txt > logs/bt/closure_keys_t$tier.txt
    ARGS="logs/bt/gates_run_closure.jsonl $CEN --levels auto --keys @logs/bt/closure_keys_t$tier.txt --seen-all"
    for tmax in 30 300; do
      RETRY=""; [ $tmax -gt 30 ] && RETRY="--retry-only"
      i=0
      while [ $i -lt $W ]; do
        (ulimit -v $MEM; exec nice -n 10 python3 theory/bt_gates_levels.py $ARGS --tmax $tmax $RETRY --part $i/$W) &
        i=$((i + 1))
      done
      wait
      echo "closure round $round tier $tier tmax $tmax done"
    done
    (ulimit -v $MEM1; exec nice -n 10 python3 theory/bt_gates_levels.py $ARGS --tmax 300 --retry-only)
    echo "closure round $round tier $tier memory retry done"
  done
  after=$(cat logs/bt/gates_run_closure.jsonl | wc -l)
  [ "$after" -eq "$before" ] && break
done
