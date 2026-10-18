#!/bin/bash
# restart-safe converging capture-recapture at 283, rho = 0.30: two independent 60-minute 4-worker generations, then
# the Lincoln-Petersen estimate (scripts/flat_capture_w.py) and n' on the union
cd /home/user/Hadamard
for sd in 11 12; do
  out=/tmp/ladder/v283_r0.30_m60w4_s$sd
  [ -f $out/orbits_19.npy ] || python3 scripts/flat_gen_w.py 283 3 0.30 60 4 $out --seed $sd --rows 7,19 > logs/k283/gen_r0.30_m60w4_s$sd.log 2>&1
done
python3 scripts/flat_capture_w.py 283 3 /tmp/ladder/v283_r0.30_m60w4_s11 /tmp/ladder/v283_r0.30_m60w4_s12 > logs/k283/capture_r0.30_m60w4.txt
echo DONE > logs/k283/capture_r0.30_m60w4.done
