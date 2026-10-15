#!/bin/bash
# 599 lead: stage-1 solutions with every block symmetric (order-3 or dim-15 involution), bordered classes first
cd /home/user/Hadamard
P=$1
python3 theory/b599_symblock.py --n 9 --q 11 --types 3:10,3:18,2:15 --classes 0,7 --part $P/3 --out logs/b599/sym/sym15_599_p$P.jsonl >> logs/b599/sym/sym15_599_p$P.log 2>&1
python3 theory/b599_symblock.py --n 9 --q 11 --types 3:10,3:18,2:15 --part $P/3 --out logs/b599/sym/sym15_599_p$P.jsonl >> logs/b599/sym/sym15_599_p$P.log 2>&1
