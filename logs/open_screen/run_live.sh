#!/bin/bash
# live runners write to a work directory W (default: the session scratchpad); logs/open_screen/sync_live.sh copies
# them into logs/open_screen/ at each batch commit.  Restart-safe: the deficit batch resumes from W/deficit_s*.jsonl.
cd "$(dirname "$0")/../.."
W=${W:-/tmp/claude-0/-home-user-Hadamard/23d0cbad-e05e-56bb-b7b0-b45b94a4a068/scratchpad/live}
IN=$(ls logs/open_screen/screens_s*.jsonl | paste -sd,)
for i in 0 1 2 3; do
  OS_MEMGB=2.8 nohup python3 theory/open_screen_deficit.py --batch $W/deficit_s$i.jsonl $IN --shard $i/4 --secs 900 \
    > $W/deficit_s$i.err 2>&1 &
done
[ -s $W/prime_rest.txt ] && nohup python3 scripts/b0_open_screen.py $(cat $W/prime_rest.txt) 20000 >> $W/b0_prime_direct.jsonl 2>> $W/b0_prime_direct.err &
