#!/usr/bin/env python3
"""Partial character borders over Z_2^a x Z_s: Parseval-level census at the open orders 4p (p prime).

For p = 2^a s + t (t = p mod 2^a, s > 1) and a cyclic multiplier group U = <g> on Z_s (acting trivially on Z_2^a):
 - r     = orbits per block = 2^a * (number of U-orbits on Z_s);
 - Y     = the possible cell sums (sums of a U-invariant +-1 function over a coset of Z_s);
 - E     = sum of y^2 over the 4*2^a cells, forced by Parseval: ((2^a - t) n + 4 t^2) / 2^a;
 - the number of ways to write E as a sum of 4*2^a squares of elements of Y is counted by a DP (0 = dead at the
   Parseval level [P]; otherwise the least and greatest possible number of cells with |y| > 1 are reported).
Only slices with 24 <= r <= 130 are listed.  This is a necessary condition only.  [M]
Usage: fable_tcensus.py [amax]
"""
import sys, json, math, itertools
def orbits(s, g):
    seen = [False] * s; sizes = []
    for x in range(s):
        if seen[x]: continue
        y = x; c = 0
        while not seen[y]:
            seen[y] = True; y = (y * g) % s; c += 1
        sizes.append(c)
    return sizes
def cyc_subgroups(s):
    units = [x for x in range(1, s) if math.gcd(x, s) == 1]
    seen = set(); out = []
    for g in units:
        H = []; y = 1
        while True:
            H.append(y); y = (y * g) % s
            if y == 1: break
        key = frozenset(H)
        if key not in seen:
            seen.add(key); out.append((g, len(H)))
    return out
def main():
    amax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    op = json.load(open('data/open_orders_lt10000.json'))['open']
    isprime = lambda m: m > 1 and all(m % d for d in range(2, int(m ** .5) + 1))
    primes = [n // 4 for n in op if n % 4 == 0 and isprime(n // 4)]
    rows = []
    for p in primes:
        n = 4 * p
        for a in range(4, amax + 1):
            k = 2 ** a; t = p % k; s = (p - t) // k
            if s < 3: continue
            num = (k - t) * n + 4 * t * t
            assert num % k == 0
            E = num // k; ncell = 4 * k
            for g, u in cyc_subgroups(s):
                sizes = orbits(s, g)
                r = k * len(sizes)
                if not (24 <= r <= 130): continue
                # achievable cell sums
                Y = {0}
                for z in sizes: Y = {y + z for y in Y} | {y - z for y in Y}
                sq = sorted({y * y for y in Y})
                if sq[0] != 1 and s % 2 == 1: pass
                # DP: reach[e] = (min big, max big) number of cells with y^2 > smallest, using exactly ncell cells
                base = sq[0]
                rem = E - ncell * base
                if rem < 0: rows.append((p, a, t, s, u, r, 'dead: E too small')); continue
                incs = [q - base for q in sq[1:]]
                INF = 10 ** 9
                lo = [INF] * (rem + 1); hi = [-1] * (rem + 1); lo[0] = 0; hi[0] = 0
                for e in range(1, rem + 1):
                    for d in incs:
                        if d <= e and lo[e - d] < INF:
                            lo[e] = min(lo[e], lo[e - d] + 1); hi[e] = max(hi[e], hi[e - d] + 1)
                if lo[rem] > ncell: rows.append((p, a, t, s, u, r, 'dead at Parseval'))
                else: rows.append((p, a, t, s, u, r, f'cells above the minimum: {lo[rem]}..{min(hi[rem], ncell)} of {ncell}; |y| in {sorted({abs(y) for y in Y})[:6]}'))
    alive = [x for x in rows if not x[6].startswith('dead')]
    print(f'{len(rows)} slices with 24 <= r <= 130; {len(rows) - len(alive)} dead at the Parseval level; {len(alive)} pass')
    for x in sorted(alive, key=lambda x: (x[5], x[0]))[:60]:
        print('p=%d a=%d |T|=%d s=%d |U|=%d r=%d: %s' % x)
    json.dump(rows, open('logs/fable_fresh/tcensus.json', 'w'))
if __name__ == '__main__':
    main()
