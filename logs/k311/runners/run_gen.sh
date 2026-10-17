#!/bin/bash
# restart-safe long generation of rho-flat orbit representatives at 311, in 10-minute chunks (one dir per chunk)
# usage: run_gen.sh RHO NCHUNKS
cd /home/user/Hadamard
SP=${SP:-/tmp/k311}   # scratch dir for pools and caches (large, not in git); set SP to override
RHO=$1; N=$2
for c in $(seq 1 $N); do
  out=$SP/gen_r${RHO}_c${c}
  if [ -f $out/orbits_1.npy ] && [ -f $out/orbits_19.npy ] && [ -f $out/orbits_21.npy ]; then continue; fi
  python3 scripts/flat_probe.py 311 5 $RHO 10.0 4 $out --seed $((100 + c)) > logs/k311/gen_r${RHO}_c${c}.log 2>&1
done
echo DONE > logs/k311/gen_r${RHO}.done
