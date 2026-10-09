#!/bin/bash
# exact h0 + qnorm screens for every b>=3 Paley-coset option (logs/paley_bge3_jobs.txt), 3 in parallel
cd "$(dirname "$0")/.."
xargs -P 3 -L 1 bash -c 'timeout 3600 python3 scripts/paley_bge3.py screen1 "$0" "$1" "$2" >> logs/paley_bge3_screen.jsonl || echo "{\"v\": $0, \"c\": $1, \"p\": $2, \"timeout\": true}" >> logs/paley_bge3_screen.jsonl' < logs/paley_bge3_jobs.txt
