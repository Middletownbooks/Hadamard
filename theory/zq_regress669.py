r"""Regression [M]: the generalised census driver (theory/bt_census2.py after the Z_q generalisation) reproduces the
669 z3all census records of logs/bt/census2/z3all_669.jsonl (key, verdict, r, e, Q-level solutions, log2 E_Q).
usage: python3 theory/zq_regress669.py OLD.jsonl NEW.jsonl"""
import sys, json
old = {json.loads(l)['key']: json.loads(l) for l in open(sys.argv[1])}
new = {json.loads(l)['key']: json.loads(l) for l in open(sys.argv[2])}
fields = ('verdict', 'r', 'e', 'q_solutions', 'log2E_Q', 'q_patterns')
diff = [(k, f, old[k].get(f), new.get(k, {}).get(f)) for k in old for f in fields if old[k].get(f) != new.get(k, {}).get(f)]
print(json.dumps(dict(old=len(old), new=len(new), same_keys=set(old) == set(new), mismatches=len(diff),
                      first=diff[:5])))
