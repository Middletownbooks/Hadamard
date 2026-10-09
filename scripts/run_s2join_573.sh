#!/bin/bash
# stage 2 (exact, Z4-class-split join) over recorded 573 stage-1 solutions; one process per npz file
cd "$(dirname "$0")/.."
for f in "$@"; do
  b=$(basename "$f" .npz)
  nohup python3 scripts/paley_s2join.py "$f" 573 11 52 29 --dcomp 4 > "logs/s2join_573_${b}.jsonl" 2>&1 &
done
