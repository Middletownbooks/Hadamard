r"""End-to-end control for gate A1/A2 (shares no code with the K-level tools): the number of Legendre-quadruple cores in
the Paley plug-in slice on Z_n x F_q, counted at the G level by brute force.  Every block (x0, x1, x2) in
({+-1}^n)^3 is expanded to G = Z_n x F_q (x0 on Z_n x 0, x1 on the squares, x2 on the non-squares), its PAF over G is
computed by the definition, and a dictionary 2+2 join counts ordered quadruples with sum PAF = 4m at 0 and -4 elsewhere.
usage: python3 theory/h2396_gbrute.py n q"""
import sys, itertools
from collections import Counter
n, q = int(sys.argv[1]), int(sys.argv[2]); m = n * q
sq = {(z * z) % q for z in range(1, q)}
def expand(x0, x1, x2):
    return [x0[k] if z == 0 else (x1[k] if z in sq else x2[k]) for k in range(n) for z in range(q)]
def paf(f):
    return tuple(sum(f[k * q + z] * f[((k + h) % n) * q + (z + w) % q] for k in range(n) for z in range(q))
                 for h in range(n) for w in range(q))
blocks = Counter()
for bits in itertools.product((1, -1), repeat=3 * n):
    f = expand(bits[:n], bits[n:2 * n], bits[2 * n:])
    if sum(f) ** 2 > 4: continue
    blocks[paf(f)] += 1
T = tuple(4 * m if i == 0 else -4 for i in range(m))
AB = Counter()
for ka, ca in blocks.items():
    for kb, cb in blocks.items():
        AB[tuple(a + b for a, b in zip(ka, kb))] += ca * cb
tot = sum(c * AB.get(tuple(t - x for t, x in zip(T, s)), 0) for s, c in AB.items())
print(f"n={n} q={q}: blocks {sum(blocks.values())} (distinct PAF {len(blocks)}); ordered cores {tot}")
