#!/bin/bash
# exact finer levels (theory/bt_gates_levels.py, every invariant N <= K, coarsest first) on the structures alive
# after the census, per order once its census is complete; 10 s per level, then 60 s for the levels capped at 10 s.
# Restart-safe.  usage: logs/zq/gates/run.sh
cd /home/user/Hadamard
ulimit -v 3600000
for x in "z7all 2429" "z13all 2327" "z5all 2315" "z5all 2215" "z13all 2171" "z5all 2095" "z7all 1981" "z11all 1969" "z3all 1929" "z3all 1893" "z5all 1795" "z9x3all 1719" "z9all 1719" "z3all 1719"; do
  set -- $x; tag=$1_$2
  until python3 - "$tag" <<'PY'
import sys, json
tag = sys.argv[1]
try:
    tot = json.loads(open(f'logs/zq/census/{tag}.p0.log').readline())['structures']
    done = sum(1 for _ in open(f'logs/zq/census/{tag}.jsonl'))
except Exception:
    sys.exit(1)
sys.exit(0 if done >= tot else 1)
PY
  do sleep 60; done
  for T in 10 60; do
    python3 theory/zq_summary.py --keys $tag logs/zq/gates/$tag.keys > /dev/null
    [ -s logs/zq/gates/$tag.keys ] || break
    while true; do                       # exit code 75: memory held after a refused level, resume fresh
      python3 theory/bt_gates_levels.py logs/zq/gates/$tag.jsonl logs/zq/census/$tag.jsonl --levels auto \
        --keys @logs/zq/gates/$tag.keys --tmax $T >> logs/zq/gates/$tag.log 2>&1
      [ $? -eq 75 ] || break
    done
  done
  grep -q "^$tag done" logs/zq/gates/run.done 2>/dev/null || echo "$tag done" >> logs/zq/gates/run.done
done
