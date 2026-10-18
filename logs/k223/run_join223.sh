#!/bin/bash
# restart-safe 223 disk join only (slots already built by run_prod223.sh); resumes from its log
cd /home/user/Hadamard
export NUMPY_MADVISE_HUGEPAGE=0
SP=${SP:-/tmp/k223}; RHO=0.32; MEM=${MEM:-4e9}
python3 scripts/kt_join_disk.py 223 3 23,11,11,11 $SP/slots_r${RHO} $SP/join_r${RHO}_k12 --k 12 --mem $MEM --hours 200 >> logs/k223/join_r${RHO}_k12.out 2>&1
cp $SP/join_r${RHO}_k12/log.jsonl logs/k223/join_r${RHO}_k12_log.jsonl
cp $SP/join_r${RHO}_k12/solutions_v223.txt logs/k223/ 2>/dev/null
