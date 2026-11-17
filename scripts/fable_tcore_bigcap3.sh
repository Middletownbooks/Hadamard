#!/bin/bash
# Large-cap pass, third ordering: the two fully decidable 719 a = 6 slices first (8 h each), then the partial
# slices, then the classical a = 4 slices that ran into the 3 h limit before (last). Restart-safe like bigcap2:
# skips keys already marked "done" in tcore_bigcap.log. Optionally waits for a running PID ($WAITPID) first.
cd /home/user/Hadamard
L=logs/fable_fresh
if [ -n "$WAITPID" ]; then
  while kill -0 "$WAITPID" 2>/dev/null; do sleep 60; done
  [ -n "$WAITKEY" ] && echo "$(date -u +%T) done $WAITKEY" >> $L/tcore_bigcap.log
fi
n=0
while IFS= read -r KEY; do
  [ -z "$KEY" ] && continue
  n=$((n+1))
  if grep -qF "done $KEY" $L/tcore_bigcap.log 2>/dev/null; then continue; fi
  TL=10800; [ "$n" -le 2 ] && TL=28800
  echo "$(date -u +%T) bigcap $KEY (limit ${TL}s)" >> $L/tcore_bigcap.log
  (ulimit -v 11000000; OMP_NUM_THREADS=2 timeout $TL nice -n 4 python3 theory/fable_tcore.py "$KEY" --stream --maxlift 5e6 --maxpairs 5e7 --cache 4 >> $L/tcore_bigcap.log 2>&1)
  echo "$(date -u +%T) done $KEY" >> $L/tcore_bigcap.log
done < $L/tcore_bigcap_todo3.txt
echo "$(date -u +%T) bigcap done" >> $L/tcore_bigcap.log
