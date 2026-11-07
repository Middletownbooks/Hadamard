r"""Brute-force stage-1 count for the plug-in slice on Z_n x F_q (gate A2 control; shares no code with
theory/h2396_s1enum.c): every u in L^n per block (L = {+-1, +-(q-2), +-q}), dictionary 2+2 join on the PAF vector
(lags 0..n/2); counts ordered quadruples, split by the allocation (#letters +-q, #letters +-(q-2)).
usage: python3 theory/h2396_s1brute.py n q"""
import sys, itertools
from collections import Counter
n, q = int(sys.argv[1]), int(sys.argv[2]); m = q * n
L = [1, -1, q - 2, -(q - 2), q, -q]
T = tuple([4 * m + 4 - 4 * q] + [-4 * q] * (n // 2))
def paf(u): return tuple(sum(u[k] * u[(k + d) % n] for k in range(n)) for d in range(n // 2 + 1))
blocks = Counter()
for u in itertools.product(L, repeat=n):
    if sum(u) ** 2 > 4: continue                           # rows (+-2, 0, 0, 0) / (+-1)^4
    if paf(u)[0] > T[0]: continue
    blocks[(paf(u), sum(1 for x in u if abs(x) == q), sum(1 for x in u if abs(x) == q - 2), sum(u) ** 2)] += 1
keys = list(blocks.items())
AB = Counter()
for (ka, ca) in keys:
    for (kb, cb) in keys:
        AB[(tuple(a + b for a, b in zip(ka[0], kb[0])), ka[1] + kb[1], ka[2] + kb[2], ka[3] + kb[3])] += ca * cb
# probe properly
tot = Counter()
ABd = {}
for (k, a13, a11, rr), w in AB.items(): ABd.setdefault(k, []).append((a13, a11, rr, w))
for (kc, cc) in keys:
    for (kd, cd) in keys:
        need = tuple(t - c - d for t, c, d in zip(T, kc[0], kd[0]))
        for (a13, a11, rr, w) in ABd.get(need, []):
            if rr + kc[3] + kd[3] != (4 if n % 2 == 0 else 4): continue
            tot[(a13 + kc[1] + kd[1], a11 + kc[2] + kd[2])] += w * cc * cd
print({f"e_plus={a},e_minus={b}": v for (a, b), v in sorted(tot.items())})
