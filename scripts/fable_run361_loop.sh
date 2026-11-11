#!/bin/bash
# Restart-safe multi-seed loop for the skew 361 split-torus slice: fresh Metropolis lists per seed, the top NSIG
# ranked split values per seed (hits on the calibration controls came in the first 1-3 % of the ranked values).
# Usage: fable_run361_loop.sh [threads] [nsig] [first_seed] [last_seed]
TH=${1:-4}; NSIG=${2:-2000}; S0=${3:-2}; S1=${4:-40}
W=/home/user/skew_work; cd /home/user/Hadamard
exec 9>$W/loop.lock; flock -n 9 || { echo "loop already running"; exit 0; }
for s in $(seq $S0 $S1); do
  D=$W/p361_s$s
  [ -f $D/DONE ] && continue
  if ls $W/p361_s*/cert_*.json >/dev/null 2>&1; then echo "$(date -u +%T) hit present, stopping" >> logs/fable_fresh/p361_loop.log; exit 0; fi
  mkdir -p $D
  echo "$(date -u +%F' '%T) seed $s start" >> logs/fable_fresh/p361_loop.log
  SKEW_WORK=$W python3 theory/skew_kt.py theory/specs/p361_torus_a.json --L 3e7,3e7,3e7 --enummax 0 --gensecs 5400 --T 0.6 --dedup --nsig $NSIG --K 4 --seed $s --threads $TH --work $D >> logs/fable_fresh/p361_s$s.log 2>&1
  rc=$?
  echo "$(date -u +%F' '%T) seed $s exit $rc: $(grep -c 'kt4 split values' logs/fable_fresh/p361_s$s.log) chunks, $(grep -c '^..:..:.. SOL' logs/fable_fresh/p361_s$s.log) SOL lines" >> logs/fable_fresh/p361_loop.log
  if [ $rc -eq 0 ]; then
    touch $D/DONE
    ls $D/cert_*.json >/dev/null 2>&1 || rm -f $D/list?.npz $D/kt4_in.bin $D/gen?.bin $D/gen?.out
  else
    sleep 30
  fi
done
