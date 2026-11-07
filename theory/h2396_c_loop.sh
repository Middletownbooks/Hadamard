#!/bin/bash
# Part C, gate C2 (sampled): stage-1-only instances of the 1291 index-4 plug-in on Z19 (theory/h2396_c_s1cp.py
# --nozq, one per seed), each then tested exactly at the Z_q level (theory/h2396_c_zqtest.py).
# usage: theory/h2396_c_loop.sh seed0 nseeds secs OUTDIR
s0=$1; ns=$2; secs=$3; out=$4; mkdir -p $out
for ((s=s0; s<s0+ns; s++)); do
  python3 theory/h2396_c_s1cp.py 19 17 --secs $secs --workers 1 --seed $s --nozq --out $out/s1only.jsonl > /dev/null 2>&1
  tail -1 $out/s1only.jsonl > $out/last.jsonl
  python3 theory/h2396_c_zqtest.py 19 17 $out/last.jsonl >> $out/zqtest.jsonl 2>&1
done
