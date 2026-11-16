#!/usr/bin/env python3
"""Triage of the undecided core slices for the large-cap pass: for every slice, the number of y-quadruples whose
largest lift list lies in (1e6, 5e6] (the only ones a 5e6 cap can newly decide), above 5e6, and at most 1e6.
Lift counts are products of per-orbit option counts (no enumeration).  Output: logs/fable_fresh/tcore_triage.jsonl."""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fable_tcore as TC
out = 'logs/fable_fresh/tcore_triage.jsonl'
done = set()
if os.path.exists(out):
    for l in open(out):
        try: done.add(json.loads(l)['key'])
        except Exception: pass
last = {}
for l in open('logs/fable_fresh/tcore.jsonl'):
    r = json.loads(l); last[r['key']] = r
keys = [k for k, r in last.items() if r['verdict'].startswith('UNDECIDED') and k not in done]
fo = open(out, 'a')
for key in keys:
    t0 = time.time()
    p, a, label, T, s, u, uo = json.loads(key)
    try:
        B = TC.build(a, label, T, s, u)
        cnt, sols, tgt = TC.qlevel_solutions(B, p, len(T))
        if isinstance(cnt, str) or not sols:
            fo.write(json.dumps(dict(key=key, qlevel=str(cnt))) + '\n'); fo.flush(); continue
        ys, quads = sols
        usz = [len(o) for o in B['uorbs']]
        memo = {}
        def nl(i):
            if i not in memo:
                tot = 1
                for yA in ys[i]:
                    tot *= sum(1 for e in itertools.product((1, -1), repeat=len(usz)) if sum(x * z for x, z in zip(e, usz)) == yA)
                memo[i] = tot
            return memo[i]
        small = mid = big = 0
        for q in quads:
            m = max(nl(i) for i in q)
            if m <= 1e6: small += 1
            elif m <= 5e6: mid += 1
            else: big += 1
        fo.write(json.dumps(dict(key=key, quads=len(quads), le1e6=small, mid=mid, gt5e6=big, secs=round(time.time() - t0, 1))) + '\n'); fo.flush()
    except Exception as e:
        fo.write(json.dumps(dict(key=key, error=repr(e)[:200])) + '\n'); fo.flush()
