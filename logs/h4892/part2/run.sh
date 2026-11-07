#!/bin/bash
# part-2 queue: theory/open_screen_stage1_exact.py KEY --methods A,B --strong --planted 2, 3600 s each, two at a time.
# Resumable: a key whose .jsonl already holds a phase "final" record is skipped.
cd "$(dirname "$0")/../../.."
one() {
  name="$1"; key="$2"; out="logs/h4892/part2/$name.jsonl"
  if [ -f "$out" ] && grep -q '"phase": "final"' "$out"; then echo "skip $name"; return; fi
  start=$(date +%s)
  timeout 3600 python3 theory/open_screen_stage1_exact.py "$key" --methods A,B --strong --planted 2 --out "$out" \
     > "logs/h4892/part2/$name.out" 2> "logs/h4892/part2/$name.err"
  rc=$?
  echo "{\"name\": \"$name\", \"rc\": $rc, \"secs\": $(( $(date +%s) - start ))}" >> logs/h4892/part2/done.jsonl
}
export -f one
while IFS=$'\t' read -r name key; do printf '%s\0%s\0' "$name" "$key"; done < logs/h4892/part2/keys.txt |
  xargs -0 -n 2 -P 2 bash -c 'one "$0" "$1"'
