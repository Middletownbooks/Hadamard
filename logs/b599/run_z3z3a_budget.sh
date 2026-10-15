#!/bin/bash
# typical-y stage-2 budgets of the Z3^2-a slice (H' = Z3^2, -1 on one axis), several seeds
cd /home/user/Hadamard
for seed in 1 2 4 5 6 7; do
  timeout 3600 python3 theory/ab_typical_levels.py '{"v":599,"c":6,"b":5,"h":99,"m":594,"Q":[2,3],"H1":[3,3],"q":11,"s":2,"f":5,"mods":[2,3,3,3,11],"nQ":2,"U":[[2,2,2]]}' --seed $seed --N 4000 --memgb 3.0 >> logs/b599/budgets_z3z3a.jsonl 2>> logs/b599/budgets_z3z3a.err
done
