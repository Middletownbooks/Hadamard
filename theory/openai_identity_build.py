#!/usr/bin/env python3
"""Reference *rejected* tensor candidates; no construction of new cores claimed.

Class T: each Q-cell is a product of balanced punctured factors off the coordinate
hyperplanes, with arbitrary binary filling on those hyperplanes. The report proves
that every member fails. The CLI emits examples for the independent exact checker.
All arrays use mixed-radix order, Q first, last K coordinate fastest.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def prime_factors(m):
    out = []
    p = 3
    while p * p <= m:
        while m % p == 0:
            out.append(p)
            m //= p
        p += 2
    if m > 1:
        out.append(m)
    return out


def balanced_factor(d):
    """Quadratic character for primes; a balanced punctured vector otherwise."""
    assert d >= 3 and d % 2
    if all(d % j for j in range(2, math.isqrt(d) + 1)):
        return [0] + [1 if pow(x, (d - 1) // 2, d) == 1 else -1 for x in range(1, d)]
    return [0] + [1] * ((d - 1) // 2) + [-1] * ((d - 1) // 2)


def tensor_candidate(factors, qtype, filling):
    """filling=constant attains the sharp lower bound; varied uses arbitrary signs."""
    coords = list(itertools.product(*(range(d) for d in factors)))
    hs = [balanced_factor(d) for d in factors]
    rows = []
    for i in range(4):
        row = []
        for q in range(4):
            for co in coords:
                if all(co):
                    z = math.prod(h[x] for h, x in zip(hs, co))
                    row.append(z * (-1 if (i + q) % 2 else 1))
                else:
                    bit = i + q
                    if filling == 'varied':
                        bit += sum((j + 1 + i) * x + (q + 1) * x * x for j, x in enumerate(co))
                    row.append(1 if bit % 2 == 0 else -1)
        rows.append(row)
    return dict(name=f'tensor_{factors}_{qtype}_{filling}', q=qtype,
                factors=factors, m=math.prod(factors), filling=filling,
                class_name='balanced-punctured-tensor', sequences=rows)


def suite():
    out = []
    for m in range(3, 100, 2):
        decompositions = {tuple(prime_factors(m)), (m,)}
        for factors in sorted(decompositions):
            for qt, fill in itertools.product(('v4', 'z4'), ('constant', 'varied')):
                out.append(tensor_candidate(list(factors), qt, fill))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--suite', action='store_true')
    ap.add_argument('--factors', default='3,5')
    ap.add_argument('--q', choices=['v4', 'z4'], default='v4')
    ap.add_argument('--filling', choices=['constant', 'varied'], default='constant')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    data = suite() if a.suite else [tensor_candidate([int(v) for v in a.factors.split(',')], a.q, a.filling)]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(data) + '\n')
    print(json.dumps(dict(candidates=len(data), output=a.out)))


if __name__ == '__main__':
    main()
