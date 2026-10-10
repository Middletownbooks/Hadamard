#!/bin/bash
# all-subgroup b = 1 slice sweep (theory/p571_msweep.py) for the remaining open orders below 3000, 4 at a time
cd "$(dirname "$0")/.."
printf '%s\n' 522 536 598 642 646 658 668 718 738 | xargs -P 4 -I{} sh -c 'python3 theory/p571_msweep.py {} --omin 16 --omax 90 --N 10000 > logs/msweep/m{}.jsonl 2>&1'
