#!/bin/bash
# h0 + qnorm screens (scripts/paley_bge3.py screen1) on the c = 4 options that pass the w-screen, smallest n' first
cd "$(dirname "$0")/.."
for vp in "751 17" "839 19" "991 19" "955 17" "751 11" "839 11" "991 13" "883 11" "955 7"; do
  set -- $vp
  grep -q "\"v\": $1, \"c\": 4, .*\"p\": $2," logs/c4_qnorm_3000_4000.jsonl 2>/dev/null && continue
  timeout 7200 python3 scripts/paley_bge3.py screen1 $1 4 $2 >> logs/c4_qnorm_3000_4000.jsonl 2>&1 || echo "{\"v\": $1, \"c\": 4, \"p\": $2, \"timeout\": true}" >> logs/c4_qnorm_3000_4000.jsonl
done
