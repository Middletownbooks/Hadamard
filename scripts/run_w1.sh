#!/bin/bash
cd /home/user/Hadamard
( for v in 43 61 73 79 97 103 109; do python3 scripts/w1_null.py $v 3 8192 150 300 > logs/w1_v${v}_d3.log 2>&1; done ) &
( python3 scripts/w1_null.py 131 5 8192 1500 300 > logs/w1_v131_d5.log 2>&1 ) &
( python3 scripts/w1_null.py 151 5 8192 1500 300 > logs/w1_v151_d5.log 2>&1 ) &
wait
