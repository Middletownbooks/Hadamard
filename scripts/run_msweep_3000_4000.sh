#!/bin/bash
# all-subgroup b = 1 slice sweep for the 3000-4000 open list (user's list, 2026-10-10), 4 at a time, most 573-like first
cd "$(dirname "$0")/.."
printf '%s\n' 858 990 946 918 750 882 954 970 788 916 932 822 906 838 862 | xargs -P 4 -I{} sh -c 'python3 theory/p571_msweep.py {} --omin 16 --omax 90 --N 10000 > logs/msweep34/m{}.jsonl 2>&1'
