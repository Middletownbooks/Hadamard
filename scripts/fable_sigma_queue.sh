#!/bin/bash
# sigma-invariant GENERAL border runs for the 2396 slice, one at a time (1 worker, ~9.5 GB with the domain reduction).
# Single-character sigma-invariant borders are impossible [P] (report 7.1), so only general runs remain.
cd /home/user/Hadamard
for job in "general 1 7200 1" "general 4 7200 1" "general 7 7200 1"; do
  set -- $job
  if ls logs/fable_fresh/q599_sigma_border_general_*_case0.json > /dev/null 2>&1; then echo "border found, queue stopped" >> logs/fable_fresh/sigma_queue.log; exit 0; fi
  echo "$(date -u +%T) start $job" >> logs/fable_fresh/sigma_queue.log
  nice -n 5 python3 theory/fable_border10.py $1 $2 $3 $4 0 > logs/fable_fresh/sigma_border_$1_f$2.log 2>&1
  echo "$(date -u +%T) end $job: $(tail -1 logs/fable_fresh/sigma_border_$1_f$2.log | cut -c1-150)" >> logs/fable_fresh/sigma_queue.log
done
