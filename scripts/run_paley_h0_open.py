"""Fast exact h0/row-sum screen of the Paley-coset ansatz over all open-order options (resumable: skips done lines)."""
import sys, json, os
sys.path.insert(0, 'scripts')
from paley_ansatz import exact_h0_screen_fast
OUT = 'logs/paley_h0_open.jsonl'
def fac(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        r = json.loads(l); done.add((r['v'], r['b'], r['p']))
jobs = []
for v in [523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739]:
    for b in (0, 1):
        m = v - b
        for p, e in fac(m).items():
            if p >= 3 and e == 1: jobs.append((m // p, v, b, p))
for n, v, b, p in sorted(jobs):
    if (v, b, p) in done: continue
    r = exact_h0_screen_fast(v, b, n, p)
    with open(OUT, 'a') as f: f.write(json.dumps(r) + '\n')
    print(json.dumps(r), flush=True)
