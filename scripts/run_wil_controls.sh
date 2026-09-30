#!/bin/bash
# W2 positive controls (pre-registered), one at a time
cd /home/user/Hadamard
for a in "45 5 9" "51 3 17" "55 5 11" "57 3 19" "63 7 9"; do
  set -- $a
  [ -f logs/wil_n$1.json ] && continue
  PYTHONPATH=scripts timeout 36000 python3 -u scripts/wil_search.py $1 $2 $3 > logs/wil/n$1.log 2>&1
  echo "exit $?" >> logs/wil/n$1.log
done
