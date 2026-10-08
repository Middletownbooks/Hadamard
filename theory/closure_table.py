"""Composition-closure arithmetic for the open v (H(4v) unknown, 4v in 2000..3000).

Every hypothesis below is a number-theoretic condition that is checked exactly here; the literature
statement of each theorem is in report/theory_arm_2000_3000.md (tags [L] there).  Output: JSON lines
and a markdown table.
"""
import json
from math import isqrt

from sympy import factorint, isprime, primerange

V = [523, 537, 571, 573, 599, 643, 647, 659, 669, 719, 739]


def fac(n):
    return ' * '.join(f'{p}^{e}' if e > 1 else f'{p}' for p, e in sorted(factorint(n).items()))


def is_pp(n):
    f = factorint(n)
    return n > 1 and len(f) == 1


def sum2sq(n):
    """n is a sum of two squares iff every prime = 3 mod 4 divides n to an even power"""
    return all(e % 2 == 0 for p, e in factorint(n).items() if p % 4 == 3)


def prime_powers_upto(N):
    out = set()
    for p in primerange(2, N + 1):
        q = p
        while q <= N:
            out.add(q)
            q *= p
    return out


PP = prime_powers_upto(10000)


def twin_pp(n):
    """n = q(q+2) with q, q+2 prime powers"""
    q = isqrt(n + 1) - 1
    return q >= 2 and q * (q + 2) == n and q in PP and q + 2 in PP


def golay(n):
    """n = 2^a 10^b 26^c"""
    for c in range(0, 4):
        for b in range(0, 6):
            m = n
            d = 26 ** c * 10 ** b
            if m % d:
                continue
            m //= d
            if m & (m - 1) == 0:
                return True
    return False


def abelian_groups(n):
    """isomorphism types of abelian groups of order n, as lists of cyclic factors (primary)"""
    def partitions(k, maxpart=None):
        if maxpart is None:
            maxpart = k
        if k == 0:
            yield []
            return
        for first in range(min(k, maxpart), 0, -1):
            for rest in partitions(k - first, first):
                yield [first] + rest
    import itertools
    per = []
    for p, e in sorted(factorint(n).items()):
        per.append([[p ** a for a in part] for part in partitions(e)])
    out = []
    for combo in itertools.product(*per):
        out.append([x for part in combo for x in part])
    return out


rows = []
for v in V:
    f = factorint(v)
    r = dict(v=v, order=4 * v, fac=fac(v), v_mod8=v % 8)
    r['paleyI_4v-1'] = (fac(4 * v - 1), is_pp(4 * v - 1))
    r['paleyII_2v-1'] = (fac(2 * v - 1), is_pp(2 * v - 1) and (2 * v - 1) % 4 == 1)
    r['conference_2v_two_squares'] = sum2sq(2 * v - 1)
    r['miyamoto_v_pp_1mod4'] = is_pp(v) and v % 4 == 1
    # Whiteman / Mukhopadhyay type orders q^r (q+1)/2 with q a prime power, r >= 1
    wm = [(q, rr) for q in PP for rr in range(1, 4) if q ** rr * (q + 1) == 2 * v]
    r['q^r(q+1)/2'] = wm
    r['q^2+q+1'] = [q for q in PP if q * q + q + 1 == v]
    r['(q+1)/2'] = [q for q in PP if q + 1 == 2 * v]
    r['twin_pp_4v-1'] = twin_pp(4 * v - 1)
    r['twin_pp_2v-1'] = twin_pp(2 * v - 1)
    # Legendre-pair families at l = 2v-1: prime power, 2^k-1, twin prime powers
    l = 2 * v - 1
    r['LP_family'] = dict(prime_power=is_pp(l), mersenne=(l + 1) & l == 0, twin=twin_pp(l))
    r['LP_groups_l'] = abelian_groups(l)
    # Cooper-Wallis / Baumert-Hall splits v = t*w, t,w > 1
    splits = [(t, v // t) for t in range(3, v) if v % t == 0 and t < v]
    r['tw_splits'] = splits
    # T-sequence sources for t in {v} U {t | v}: BS(n+1,n) needs n=(t-1)/2; Golay needs t-1 = 2^a10^b26^c;
    # Turyn-type TT(n) -> T-sequences of length 3n-1 (n even): n=(t+1)/3
    tseq = {}
    for t in sorted({v} | {t for t, w in splits} | {w for t, w in splits}):
        tseq[t] = dict(BS_n=(t - 1) // 2, golay_t_minus_1=golay(t - 1),
                       TT_n=(t + 1) // 3 if (t + 1) % 3 == 0 and ((t + 1) // 3) % 2 == 0 else None)
    r['tseq'] = tseq
    # bordered b=1 core: abelian groups of order v-1; inflated b>=3: even c>=4 with c | v+1
    r['b1_groups_v-1'] = abelian_groups(v - 1)
    r['bge3_c'] = [c for c in range(4, v + 2, 2) if (v + 1) % c == 0 and (v + 1) // c - 1 >= 2]
    # Legendre(v) + circulant D-optimal pair of length v (prime v = 3 mod 4): needs 2v-1 = sum of two squares
    r['LegendreDopt'] = bool(isprime(v) and v % 4 == 3 and sum2sq(2 * v - 1))
    r['Ito_Q8t_paley'] = is_pp(2 * v - 1) or is_pp(4 * v - 1)
    rows.append(r)
    print(json.dumps(r, default=str))

print()
print('| v | 4v | v | 4v-1 | 2v-1 | conf(2v): 2v-1=a^2+b^2 | t*w splits | groups of order v-1 | groups of order 2v-1 | b>=3 c | Leg+Dopt |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    g1 = ', '.join('x'.join(f'Z{c}' for c in G) for G in r['b1_groups_v-1'])
    gl = ', '.join('x'.join(f'Z{c}' for c in G) for G in r['LP_groups_l'])
    print(f"| {r['v']} | {r['order']} | {r['fac']} | {r['paleyI_4v-1'][0]} | {r['paleyII_2v-1'][0]} | "
          f"{'yes' if r['conference_2v_two_squares'] else 'no'} | {r['tw_splits'] or '-'} | {g1} | {gl} | "
          f"{r['bge3_c'] or '-'} | {'yes' if r['LegendreDopt'] else 'no'} |")
