#!/bin/bash
# second method (theory/h4892_stage1_any.py) over the part-2 keys, one at a time, 3600 s each, resumable
cd "$(dirname "$0")/../../.."
tail -n +7 logs/h4892/part2/keys.txt | while IFS=$'\t' read -r name key; do
  out="logs/h4892/part2/any_$name.jsonl"
  [ -s "$out" ] && { echo "skip $name"; continue; }
  python3 -c "
import sys, json; sys.path.insert(0, 'theory')
import h4892_stage1_any as A, open_screen_stage1_exact as X
k = json.loads(sys.argv[1]); a = A.setup(k); b = X.b_setup(k)
print(json.dumps(dict(name=sys.argv[2], T1_agrees=a['T1'] == b['T1'],
      orbits_agree=sorted(map(sorted, a['orbs'])) == sorted(map(sorted, b['orbits'])))))" "$key" "$name" >> logs/h4892/part2/any_setup_check.jsonl
  timeout 3600 python3 theory/h4892_stage1_any.py "$key" --planted 2 --cap 5e6 --out "$out" > "logs/h4892/part2/any_$name.out" 2>&1
  echo "{\"name\": \"$name\", \"rc\": $?}" >> logs/h4892/part2/any_done.jsonl
done
