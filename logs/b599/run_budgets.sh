#!/bin/bash
# typical-y stage-2 level budgets (theory/ab_typical_levels.py) for the undecided 599 configurations
cd /home/user/Hadamard
for s in 6 9 12 18; do
  timeout 3600 python3 theory/ab_typical_levels.py "{\"v\":599,\"c\":8,\"b\":7,\"h\":74,\"m\":592,\"Q\":[2,2,2],\"H1\":[2],\"q\":37,\"s\":$s,\"f\":$((36/s)),\"mods\":[2,2,2,2,37],\"nQ\":3,\"U\":[]}" --N 3000 --memgb 3.0 >> logs/b599/budgets_c8.jsonl 2>> logs/b599/budgets_c8.err
done
for U in '[[2,2,8]]' '[]' '[["neg",[1]],[2,2,8]]'; do
  timeout 3600 python3 theory/ab_typical_levels.py "{\"v\":599,\"c\":6,\"b\":5,\"h\":99,\"m\":594,\"Q\":[2,3],\"H1\":[9],\"q\":11,\"s\":5,\"f\":2,\"mods\":[2,3,9,11],\"nQ\":2,\"U\":$U}" --N 3000 --memgb 3.0 >> logs/b599/budgets_c6s5.jsonl 2>> logs/b599/budgets_c6s5.err
done
