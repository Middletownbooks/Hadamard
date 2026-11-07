#!/bin/bash
# step 1-2: exact screens A, B over all enumerated slices (restart-safe; resumes from the shard files)
cd "$(dirname "$0")/../.."
for i in 0 1 2 3; do
  nohup python3 theory/open_screen_enum.py --run logs/open_screen/screens_s$i.jsonl --shard $i/4 --secs 300 \
    > logs/open_screen/screens_s$i.err 2>&1 &
  echo $! >> logs/open_screen/pids.txt
done
