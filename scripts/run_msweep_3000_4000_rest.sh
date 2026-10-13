#!/bin/bash
# restart-safe: run the all-subgroup b = 1 sweep for every order of the 3000-4000 list not yet marked done
cd "$(dirname "$0")/.."
for m in 946 918 750 882 954 970 788 916 932 822 906 838 862; do
  f=logs/msweep34/m$m.jsonl
  grep -q '"done": true' "$f" 2>/dev/null || echo $m
done | xargs -P 3 -I{} sh -c 'python3 theory/p571_msweep.py {} --omin 16 --omax 90 --N 10000 > logs/msweep34/m{}.jsonl 2>&1 && echo "{\"done\": true}" >> logs/msweep34/m{}.jsonl'
