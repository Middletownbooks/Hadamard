#!/bin/sh
# Exact re-evaluation of the prefiltered 839 k8 structures (theory/bt_close.py), repeated while that census adds more.
cd /home/user/Hadamard
while true; do
  sh -c "ulimit -v 6000000; nice -n 15 python3 theory/bt_close.py logs/bt/census2/close.jsonl logs/bt/census2/k8c_839.jsonl" >> logs/bt/census2/close_k8c.log 2>&1
  pgrep -f "bt_census2.py logs/bt/census2/k8c_839" > /dev/null || break
  sleep 600
done
