#!/bin/bash
# profile-restricted b>=3 screens for the undecided / cheapest surviving configs (v c Qmods d r L ts...)
cd /home/user/Hadamard
cat > /tmp/claude-0/-home-user-Hadamard/e1244fe7-42e1-5e13-91cc-287618221334/scratchpad/screen2_jobs.txt <<'J'
647 12 12 13 60 30 1
599 20 20 14 60 32 1
719 10 10 14 60 35 1
719 20 2,10 12 66 41 1 5 7
647 18 18 12 74 39 1 5 7
643 14 14 12 72 39 1 3 5
643 14 14 12 74 39 1 3 5
573 14 14 12 76 43 1 2 4 5
659 6 6 9 78 39 1
659 6 6 12 60 39 1
739 10 10 12 66 41 1
739 10 10 12 70 41 1
J
cat /tmp/claude-0/-home-user-Hadamard/e1244fe7-42e1-5e13-91cc-287618221334/scratchpad/screen2_jobs.txt | xargs -P 2 -L 1 bash -c 'MAXST=1e7 timeout 7200 python3 scripts/bge3_screen2.py "$0" "$@" > logs/bge3_screen2_$0_c$1_d$3_r$4.log 2>&1; echo "{\"done\": \"$0 c$1 d$3 r$4\", \"exit\": $?}" >> logs/bge3_screen2_queue.log'
