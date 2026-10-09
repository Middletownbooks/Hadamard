#!/bin/bash
# exact two-stage b = 1 pipeline on small orders (positive controls), each with a time cap
cd "$(dirname "$0")/.."
for args in "43 7 6" "57 7 8" "71 7 10" "67 11 6" "85 7 12" "133 11 12" "99 7 14" "177 11 16"; do
  echo "== $args"; timeout 3600 python3 theory/paley_b1pipe.py $args 1 2>&1 | tail -1
done
