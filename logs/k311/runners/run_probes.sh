#!/bin/bash
# restart-safe: each probe skips generation if its pools exist
cd /home/user/Hadamard
SP=${SP:-/tmp/k311}   # scratch dir for pools and caches (large, not in git); set SP to override
for rho in 0.36 0.32 0.40 0.45; do
  for sd in 1 2; do
    out=$SP/probe_r${rho}_s${sd}
    [ -f logs/k311/probe_rho${rho}_s${sd}.log ] && grep -q partial_join logs/k311/probe_rho${rho}_s${sd}.log && continue
    python3 scripts/flat_probe.py 311 5 $rho 1.0 4 $out --profile 21,21,19,1 --P 4096 --kmax 10 --seed $sd > logs/k311/probe_rho${rho}_s${sd}.log 2>&1
  done
done
echo ALLDONE > logs/k311/probes.done
