#!/bin/bash
cd /home/user/Hadamard
for v in 61 73 79 97 103 109; do
  N=4096; [ $v -le 73 ] && N=2048
  python3 scripts/w1_null.py $v 3 $N 240 300 > logs/w1_v${v}_d3.log 2>&1
done
