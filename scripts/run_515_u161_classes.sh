#!/bin/bash
# full yc-class enumeration for u=161 sols 3,4,5 (sequential, 4 threads)
cd /home/user/Hadamard
for s in 3 4 5; do python3 -u scripts/bge3_twosided_classes.py 85 5 161 5 $s > logs/twosided/u161_sol${s}_classes.log 2>&1; done
touch logs/twosided/done_u161_classes
