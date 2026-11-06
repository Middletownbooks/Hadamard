#!/bin/bash
# step 3: deficit estimates for every screen survivor (restart-safe; resumes from the output shards)
cd "$(dirname "$0")/../.."
N=${1:-2}
IN=$(ls logs/open_screen/screens_s*.jsonl | paste -sd,)
for i in $(seq 0 $((N-1))); do
  nohup python3 theory/open_screen_deficit.py --batch logs/open_screen/deficit_s$i.jsonl $IN --shard $i/$N --secs 900 \
    > logs/open_screen/deficit_s$i.err 2>&1 &
  echo "deficit $i $!" >> logs/open_screen/pids.txt
done
