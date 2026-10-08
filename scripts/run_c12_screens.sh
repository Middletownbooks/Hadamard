#!/bin/bash
# norm + DP (t=1) screens for 523 c=12 generalized options; one line per option in logs/bge3_general/c12_jobs.txt
cd /home/user/Hadamard
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
i=0
while read -r v Q PSI d r L; do
  PSI=${PSI//\'/}
  ( timeout 3600 python3 scripts/bge3_normscreen_psi.py "$v" "$Q" "$PSI" "$d" "$r" "$L" 1 > logs/bge3_general/c12_norm_$i.log 2>&1
    MAXST=2e7 timeout 7200 python3 scripts/bge3_screen_psi.py "$v" "$Q" "$PSI" "$d" "$r" "$L" 1 > logs/bge3_general/c12_dp_$i.log 2>&1
    echo "{\"done\": $i, \"Q\": \"$Q\", \"psi\": \"$PSI\"}" >> logs/bge3_general/c12_done.log ) &
  i=$((i+1))
  if (( i % 4 == 0 )); then wait; fi
done < logs/bge3_general/c12_jobs.txt
wait
