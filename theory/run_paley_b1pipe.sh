#!/bin/bash
# exact two-stage b = 1 pipeline on small orders (positive controls), each with a time cap
cd "$(dirname "$0")/.."
for args in "29 7 4" "43 7 6" "57 7 8" "71 7 10" "85 7 12" "99 7 14" "45 11 4" "67 11 6" "89 11 8" "111 11 10" "133 11 12"; do
  echo "== $args"; timeout 3600 python3 theory/paley_b1pipe.py $args 1 2>&1 | tail -1
done
