#!/bin/bash
# queue: v c Q  (rule-picked border option for each open order; all abelian Q of order c)
cd /home/user/Hadamard
cat > /tmp/bge3_jobs.txt <<'J'
523 4 2,2
523 4 4
647 4 2,2
647 4 4
719 6 6
599 6 6
739 4 2,2
739 4 4
659 4 2,2
659 4 4
643 4 2,2
643 4 4
515 4 2,2
515 4 4
669 10 10
573 14 14
571 22 22
721 38 38
J
cat /tmp/bge3_jobs.txt | xargs -P3 -L1 bash -c 'NS=800 STOPAFTER=3 python3 -u scripts/bge3_cost.py $0 $1 $2 30 220 > logs/bge3_cost_open/v$0_c$1_Q${2//,/x}.log 2>&1'
