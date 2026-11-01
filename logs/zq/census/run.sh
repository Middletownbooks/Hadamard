#!/bin/bash
# census runner, in order of 4n; restart-safe (each run resumes from its output file)
cd /home/user/Hadamard
PART=${1:-0/1}
for x in "z3all 669" "z5all 1115" "z7all 1169" "z3all 1257" "z3all 1437" "z3all 1473" "z3all 1713" "z3all 1719" "z9all 1719" "z9x3all 1719" "z5all 1795" "z3all 1893" "z3all 1929" "z11all 1969" "z7all 1981" "z5all 2095" "z13all 2171" "z5all 2215" "z5all 2315" "z13all 2327" "z7all 2429"; do
  set -- $x
  python3 theory/zq_census.py logs/zq/census/$1_$2.jsonl --fam $1 --ns $2 --part $PART >> logs/zq/census/$1_$2.p${PART%/*}.log 2>&1
done
echo ALLDONE >> logs/zq/census/run.p${PART%/*}.done
