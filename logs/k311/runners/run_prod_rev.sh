#!/bin/bash
# restart-safe 311 production on reversal-canonical closure slots (one block per reversal pair: each core is one
# A<->B pair of tuples per bucket instead of a cluster of 32), slot join k = 9 (resumes from its log)
cd /home/user/Hadamard
export NUMPY_MADVISE_HUGEPAGE=0     # no THP madvise: fresh multi-GB arrays stalled in compaction
SP=${SP:-/tmp/k311}; K=${1:-9}
SL=$SP/slots_r0.32_rev
if [ ! -f $SL/slot_D.npy ]; then
  python3 scripts/build_slot_pools.py 311 5 21,21,19,1 $SL 67108864 0 $SP/gen_r0.32_c1 $SP/gen_r0.32_c2 $SP/gen_r0.32_c3 $SP/gen_r0.32_c4 --closure --revcanon > logs/k311/slots_r0.32_rev.log 2>&1
fi
python3 scripts/kt_join_slots.py 311 5 21,21,19,1 $SL $SP/join_r0.32_rev_k$K --k $K --mem 3e9 --bits 2.75 --hours 200 >> logs/k311/join_r0.32_rev_k$K.out 2>&1
cp $SP/join_r0.32_rev_k$K/log.jsonl logs/k311/join_r0.32_rev_k${K}_log.jsonl
cp $SP/join_r0.32_rev_k$K/solutions_v311.txt logs/k311/solutions_v311_rev.txt 2>/dev/null
echo DONE > logs/k311/prod_rev.done
