#!/bin/bash
# resume after a container restart: the orders not finished (750 rerun in full)
cd "$(dirname "$0")/.."
printf '%s\n' 946 918 750 882 954 970 788 916 932 822 906 838 862 | xargs -P 4 -I{} sh -c 'python3 theory/p571_msweep.py {} --omin 16 --omax 90 --N 10000 > logs/msweep34/m{}.jsonl 2>&1'
