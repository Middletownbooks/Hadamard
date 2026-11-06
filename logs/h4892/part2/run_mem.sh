#!/bin/bash
# memory-bound B reruns, one at a time, address space 13 GB, 3600 s each (run only with no other heavy job)
cd "$(dirname "$0")/../../.."
for name in 4792_b1 5892_b1 5164_k4 4748_k4 7916_k4; do
  key=$(grep -P "^$name\t" logs/h4892/part2/keys.txt | cut -f2)
  out="logs/h4892/part2/M_$name.jsonl"
  if [ -f "$out" ] && grep -q '"phase": "final"' "$out"; then continue; fi
  start=$(date +%s)
  timeout 3600 python3 theory/open_screen_stage1_exact.py "$key" --methods B --strong --planted 2 --aslimit 13 --out "$out" \
     > "logs/h4892/part2/M_$name.out" 2> "logs/h4892/part2/M_$name.err"
  echo "{\"name\": \"M_$name\", \"rc\": $?, \"secs\": $(( $(date +%s) - start ))}" >> logs/h4892/part2/done.jsonl
done
