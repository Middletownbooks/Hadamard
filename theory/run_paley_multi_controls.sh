#!/bin/bash
# small doubly-twisted b=1 slices over Z2 x Z_q x Z_p (positive-control search for the exact join)
cd "$(dirname "$0")/.."
for args in "111 1 2 -- 5 11" "131 1 2 -- 5 13" "155 1 2 -- 7 11" "171 1 2 -- 5 17" "183 1 2 -- 7 13" "191 1 2 -- 5 19"; do
  echo "== $args"; timeout 1800 python3 theory/paley_multi.py mitm $args 2>&1 | tail -1
done
