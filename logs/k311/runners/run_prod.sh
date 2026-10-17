#!/bin/bash
# restart-safe 311 production: generation chunks -> slot pools -> slot k-tree join (resumes from its log)
# usage: run_prod.sh RHO NCHUNKS L LA K BITS
cd /home/user/Hadamard
SP=${SP:-/tmp/k311}   # scratch dir for pools and caches (large, not in git); set SP to override
RHO=$1; N=$2; L=$3; LA=$4; K=$5; BITS=$6
$SP/run_gen.sh $RHO $N
SLOTS=$SP/slots_r${RHO}
if [ ! -f $SLOTS/slot_D.npy ]; then
  python3 scripts/build_slot_pools.py 311 5 21,21,19,1 $SLOTS $L $LA $SP/gen_r${RHO}_c* $SP/probe_r${RHO}_s* > logs/k311/slots_r${RHO}.log 2>&1
fi
python3 scripts/kt_join_slots.py 311 5 21,21,19,1 $SLOTS $SP/join_r${RHO} --k $K --mem 3e9 --bits $BITS --hours 200 >> logs/k311/join_r${RHO}.out 2>&1
cp $SP/join_r${RHO}/log.jsonl logs/k311/join_r${RHO}_log.jsonl
cp $SP/join_r${RHO}/solutions_v311.txt logs/k311/ 2>/dev/null
echo DONE > logs/k311/prod_r${RHO}.done
