#!/bin/bash
# wrapper kept for old check-in messages: starts the multi-seed loop (flock makes it a no-op if already running)
nohup /home/user/skew_work/loop361.sh ${1:-4} 2000 2 40 > /dev/null 2>&1 &
echo "loop launched (no-op if already running)"
