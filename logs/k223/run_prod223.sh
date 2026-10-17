#!/bin/bash
# restart-safe 223 production (rung 2): rho = 0.32 generation in 10-minute chunks (4 workers) until the row-sum-11
# orbits reach TARGET (default 3.6M orbits = 2^28 blocks after closure), then reversal-canonical closure slots (A = row sum 23, whole
# closure; B = C = D = the row-sum-11 closure, ~2^27 blocks, shared files), then the disk join (k = 12).
# usage: run_prod223.sh [TARGET_ORBITS] [MAXCHUNKS]      scratch: SP (default /tmp/k223); env GW (workers, 4), NICE (0)
cd /home/user/Hadamard
export NUMPY_MADVISE_HUGEPAGE=0
SP=${SP:-/tmp/k223}; TARGET=${1:-3600000}; MAXC=${2:-400}; RHO=0.32
mkdir -p $SP logs/k223
count() { python3 - "$SP" <<'PY'
import sys, glob, numpy as np
sys.path.insert(0, 'scripts'); import flatw as FW
fs = sorted(glob.glob(sys.argv[1] + '/prod_c*/orbits_11.npy'))
u = FW.unique_w(np.concatenate([np.load(f) for f in fs])) if fs else np.zeros((0, 2), np.uint64)
print(len(u))
PY
}
for c in $(seq 1 $MAXC); do
  out=$SP/prod_c$c
  if [ ! -f $out/orbits_23.npy ]; then
    n=$(count); echo "chunk $c: union rs11 orbits $n" >> logs/k223/gen_prod_r${RHO}.log
    [ "$n" -ge "$TARGET" ] && break
    nice -n ${NICE:-0} python3 scripts/flat_gen_w.py 223 3 $RHO 10 ${GW:-4} $out --seed $((100 + c)) --rows 23,11 > logs/k223/gen_prod_r${RHO}_c$c.log 2>&1
  fi
done
SL=$SP/slots_r${RHO}
if [ ! -f $SL/slot_D.npy ]; then
  python3 scripts/build_slot_pools_w.py 223 3 23,11,11,11 $SL 536870912 0 $SP/prod_c* --closure --revcanon > logs/k223/slots_r${RHO}.log 2>&1   # B = C = D share one pool (symlinks)
fi
python3 scripts/kt_join_disk.py 223 3 23,11,11,11 $SL $SP/join_r${RHO}_k12 --k 12 --mem 6e9 --hours 200 >> logs/k223/join_r${RHO}_k12.out 2>&1
cp $SP/join_r${RHO}_k12/log.jsonl logs/k223/join_r${RHO}_k12_log.jsonl
cp $SP/join_r${RHO}_k12/solutions_v223.txt logs/k223/ 2>/dev/null
echo DONE > logs/k223/prod_r${RHO}.done
