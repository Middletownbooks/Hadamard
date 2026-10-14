#!/bin/bash
# control chain (runs after the 495 x83 exhaustive run); every job is the preferred OOM victim
cd /home/user/Hadamard
while kill -0 5322 2>/dev/null; do sleep 30; done
for args in "442 55 17 13" "266 83 7 19" "286 87 11 13"; do
  set -- $args
  bash -c "echo 800 > /proc/self/oom_score_adj; exec python3 theory/p647_control.py $1 $2 $3 $4 --dir logs/p647/control --maxjoint 20" > logs/p647/control/exh_$(( $1 + 1 ))_$2.log 2>&1
done
