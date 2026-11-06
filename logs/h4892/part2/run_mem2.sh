#!/bin/bash
# second memory pass (decision printed to stdout first): keys whose M_ run raised MemoryError; after run_mem.sh ends
cd "$(dirname "$0")/../../.."
while pgrep -f run_mem.sh > /dev/null; do sleep 20; done
for name in $(python3 -c "
import json
d = {}
for l in open('logs/h4892/part2/done.jsonl'):
    r = json.loads(l); d[r['name']] = r['rc']
print(' '.join(n[2:] for n, rc in d.items() if n.startswith('M_') and rc == 1))"); do
  key=$(grep -P "^$name\t" logs/h4892/part2/keys.txt | cut -f2)
  start=$(date +%s)
  timeout 3600 python3 theory/open_screen_stage1_exact.py "$key" --methods B --strong --planted 2 --aslimit 12 --out "logs/h4892/part2/M2_$name.jsonl" \
     > "logs/h4892/part2/M2_$name.out" 2> "logs/h4892/part2/M2_$name.err"
  echo "{\"name\": \"M2_$name\", \"rc\": $?, \"secs\": $(( $(date +%s) - start ))}" >> logs/h4892/part2/done.jsonl
done
