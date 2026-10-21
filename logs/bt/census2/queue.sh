#!/bin/sh
# census queue for worker $1 of $2 (runs families in order; each command is restart-safe)
cd /home/user/Hadamard
W=$1; NW=$2
python3 theory/bt_census2.py logs/bt/census2/k8_w$W.jsonl --fam k8 --ns 751,599,719,863,839 --rmin 12 --rmax 110 --part $W/$NW > logs/bt/census2/k8_w$W.log 2>&1
python3 theory/bt_census2.py logs/bt/census2/k6_w$W.jsonl --fam k6 --ns 599,719,839,863,947,971 --rmin 12 --rmax 110 --part $W/$NW > logs/bt/census2/k6_w$W.log 2>&1
python3 theory/bt_census971.py logs/bt/census2/k4_971_w$W.jsonl --fam k4 --rmin 12 --rmax 110 --part $W/$NW > logs/bt/census2/k4_971_w$W.log 2>&1
python3 theory/bt_census2.py logs/bt/census2/k4_947_w$W.jsonl --fam k4 --ns 947 --rmin 12 --rmax 110 --part $W/$NW > logs/bt/census2/k4_947_w$W.log 2>&1
python3 theory/bt_census2.py logs/bt/census2/k4one_w$W.jsonl --fam k4one --ns 599,719,751,839,859,863,907,955,947 --rmin 12 --rmax 90 --part $W/$NW > logs/bt/census2/k4one_w$W.log 2>&1
