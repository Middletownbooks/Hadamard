#!/bin/bash
# restart-safe flat-set measurements for a ladder rung: per rho, two 1-core generations (seeds 1, 2), capture-recapture,
# then n' (Gaussian + partial joins) on the union.  usage: measure.sh V D PROFILE RHOS MINUTES   (scratch: /tmp/ladder)
cd /home/user/Hadamard
V=$1; D=$2; PROF=$3; RHOS=$4; MIN=$5; ROWS=$(echo $PROF | tr ',' '\n' | sort -un | paste -sd,)
mkdir -p logs/k$V /tmp/ladder
for rho in $RHOS; do
  for sd in 1 2; do
    out=/tmp/ladder/v${V}_r${rho}_m${MIN}_s${sd}
    [ -f $out/orbits_$(echo $ROWS | cut -d, -f1).npy ] && continue
    nice -n 15 python3 scripts/flat_gen_w.py $V $D $rho $MIN 1 $out --seed $sd --rows $ROWS > logs/k$V/gen_r${rho}_m${MIN}_s${sd}.log 2>&1
  done
  [ -s logs/k$V/capture_r${rho}_m${MIN}.txt ] || python3 scripts/flat_capture_w.py $V $D /tmp/ladder/v${V}_r${rho}_m${MIN}_s1 /tmp/ladder/v${V}_r${rho}_m${MIN}_s2 > logs/k$V/capture_r${rho}_m${MIN}.txt
  grep -q '"nprime"' logs/k$V/nprime_r${rho}_m${MIN}.log 2>/dev/null || nice -n 15 python3 scripts/flat_nprime_w.py $V $D $PROF /tmp/ladder/v${V}_r${rho}_m${MIN}_s1 /tmp/ladder/v${V}_r${rho}_m${MIN}_s2 --kmax 12 --big 30000 > logs/k$V/nprime_r${rho}_m${MIN}.log 2>&1
done
