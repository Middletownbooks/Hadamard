#!/bin/bash
# profile-restricted b>=3 screens; usage: run_screen2_queue.sh jobs.txt [parallel]   (job line: v c Qmods d r L ts...)
cd /home/user/Hadamard
P=${2:-2}
xargs -P "$P" -L 1 bash -c 'MAXST=${MAXST:-1e7} timeout 7200 python3 scripts/bge3_screen2.py "$0" "$@" > "logs/bge3_screen2_$0_c$1_d$3_r$4_t${6}.log" 2>&1; echo "{\"done\": \"$0 c$1 d$3 r$4 t$6\", \"exit\": $?}" >> logs/bge3_screen2_queue.log' < "$1"
