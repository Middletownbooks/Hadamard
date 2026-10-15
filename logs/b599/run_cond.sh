#!/bin/bash
# conditional stage-1 estimates (theory/b599_cond.py) for 599 c = 6 variants
cd /home/user/Hadamard
timeout 7200 python3 theory/b599_cond.py '{"v":599,"c":6,"b":5,"h":99,"m":594,"Q":[2,3],"H1":[3,3],"q":11,"s":2,"f":5,"mods":[2,3,3,3,11],"nQ":2,"U":[["neg",[2,3]]]}' --N 2000 --out logs/b599/cond_599_z3z3_neg23.jsonl > logs/b599/cond_599_z3z3_neg23.log 2>&1
timeout 7200 python3 theory/b599_cond.py '{"v":599,"c":6,"b":5,"h":99,"m":594,"Q":[2,3],"H1":[3,3],"q":11,"s":2,"f":5,"mods":[2,3,3,3,11],"nQ":2,"U":[[2,2,2]]}' --N 2000 --out logs/b599/cond_599_z3z3_x2a2.jsonl > logs/b599/cond_599_z3z3_x2a2.log 2>&1
