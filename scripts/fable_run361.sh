#!/bin/bash
# restart-safe launcher for the 361 split-torus skew slice (fresh lists, seed 2)
cd /home/user/Hadamard
if pgrep -f "skew_kt.py theory/specs/p361_torus_a.json" > /dev/null; then echo "already running"; exit 0; fi
SKEW_WORK=/home/user/skew_work nohup python3 theory/skew_kt.py theory/specs/p361_torus_a.json --L 3e7,3e7,3e7 --enummax 0 --gensecs 5400 --T 0.6 --dedup --nsig 200000 --K 4 --seed 2 --threads ${1:-4} --work /home/user/skew_work/p361_s2 >> logs/fable_fresh/p361_s2.log 2>&1 &
echo started
