#!/bin/bash
# Gate A2/A3 sampler at fixed budget: level 1 = CP-SAT stage-1 instances (theory/h2396_s1cp.py, one per seed, deduped),
# level 2 = the exact stage-2 counter (theory/h2396_stage2.c) on each new instance.  Appends to OUTDIR/inst.jsonl and
# OUTDIR/stage2.log.  usage: theory/h2396_s1cp_loop.sh n ep em seed0 nseeds secs_per_seed stage2_bin OUTDIR [--zq]
n=$1; ep=$2; em=$3; s0=$4; ns=$5; secs=$6; bin=$7; out=$8; zq=$9
mkdir -p $out
for ((s=s0; s<s0+ns; s++)); do
  tmp=$out/seed$s.jsonl
  python3 theory/h2396_s1cp.py $n 13 $ep $em --secs $secs --workers 1 --seed $s $zq --out $tmp > /dev/null 2>&1
  if grep -q FEASIBLE $tmp 2>/dev/null; then
    u=$(python3 -c "import json; d=json.loads(open('$tmp').readline()); print('\n'.join(' '.join(map(str,r)) for r in d['u']))")
    key=$(echo "$u" | md5sum | cut -c1-12)
    if ! grep -q "$key" $out/inst.jsonl 2>/dev/null; then
      res=$(echo "$u" | $bin $n 13 2>/dev/null | grep COUNT)
      echo "{\"seed\": $s, \"key\": \"$key\", \"stage2\": \"$res\", \"wall\": $(python3 -c "import json; print(json.loads(open('$tmp').readline())['wall'])")}" >> $out/inst.jsonl
      echo "$u" > $out/inst_$key.txt
    else
      echo "{\"seed\": $s, \"key\": \"$key\", \"duplicate\": true}" >> $out/inst.jsonl
    fi
  else
    echo "{\"seed\": $s, \"status\": \"$(python3 -c "import json; print(json.loads(open('$tmp').readline())['status'])" 2>/dev/null)\"}" >> $out/inst.jsonl
  fi
done
