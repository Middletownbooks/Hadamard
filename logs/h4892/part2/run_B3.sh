#!/bin/bash
# B-only reruns (method A first ate the hour): open_screen_stage1_exact.py KEY --methods B --strong --planted 2, 3600 s
cd "$(dirname "$0")/../../.."
while IFS=$'\t' read -r name key; do
  out="logs/h4892/part2/B_$name.jsonl"
  if [ -f "$out" ] && grep -q '"phase": "final"' "$out"; then continue; fi
  start=$(date +%s)
  timeout 3600 python3 theory/open_screen_stage1_exact.py "$key" --methods B --strong --planted 2 --out "$out" \
     > "logs/h4892/part2/B_$name.out" 2> "logs/h4892/part2/B_$name.err"
  echo "{\"name\": \"B_$name\", \"rc\": $?, \"secs\": $(( $(date +%s) - start ))}" >> logs/h4892/part2/done.jsonl
done < logs/h4892/part2/keys_B3.txt
