#!/bin/bash
# restart-safe 223 flat-set probes: 1 worker, nice, two seeds per rho (capture-recapture)
cd /home/user/Hadamard
for rho in $RHOS; do for sd in 1 2; do
  out=/tmp/k223/gen_r${rho}_m${MIN}_s${sd}
  [ -f $out/orbits_23.npy ] && continue
  nice -n 15 python3 scripts/flat_gen_w.py 223 3 $rho $MIN 1 $out --seed $sd --rows 23,11 > logs/k223/gen_r${rho}_m${MIN}_s${sd}.log 2>&1
done; python3 scripts/flat_capture_w.py 223 3 /tmp/k223/gen_r${rho}_m${MIN}_s1 /tmp/k223/gen_r${rho}_m${MIN}_s2 > logs/k223/capture_r${rho}_m${MIN}.txt; done
