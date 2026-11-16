#!/bin/bash
# Large-cap pass, ordered (599 a = 6 slices last); restart-safe: skips keys already present in tcore_bigcap.log as done.
cd /home/user/Hadamard
L=logs/fable_fresh
while IFS= read -r KEY; do
  [ -z "$KEY" ] && continue
  if grep -qF "done $KEY" $L/tcore_bigcap.log 2>/dev/null; then continue; fi
  echo "$(date -u +%T) bigcap $KEY" >> $L/tcore_bigcap.log
  (ulimit -v 11000000; OMP_NUM_THREADS=2 timeout 10800 nice -n 4 python3 theory/fable_tcore.py "$KEY" --stream --maxlift 5e6 --maxpairs 5e7 --cache 4 >> $L/tcore_bigcap.log 2>&1)
  echo "$(date -u +%T) done $KEY" >> $L/tcore_bigcap.log
done < ${TODO:-$L/tcore_bigcap_todo.txt}
echo "$(date -u +%T) bigcap done" >> $L/tcore_bigcap.log
