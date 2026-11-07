#!/bin/bash
# Gate A2: run the exact stage-2 counter (theory/h2396_stage2.c) on every stage-1 class written by h2396_s1enum.
# usage: theory/h2396_stage2_all.sh n q sols_file stage2_binary   -> lines "weight COUNT ..." ; last line: weighted total
n=$1; q=$2; f=$3; bin=$4
tot=0
while IFS= read -r line; do
  w=${line%% |*}; rest=${line#* | }
  res=$(echo "$rest" | tr '|' '\n' | $bin $n $q 2>/dev/null)
  c=$(echo "$res" | awk '/^COUNT/{print $2}')
  echo "$w $res"
  tot=$(python3 -c "print($tot + $w * ${c:-0})")
done < "$f"
echo "WEIGHTED_TOTAL $tot"
