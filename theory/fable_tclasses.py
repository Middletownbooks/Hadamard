#!/usr/bin/env python3
"""Partial character borders (report/fable_fresh_answer.md section 7): shapes and where they occur.
1. Affine equivalence classes of t-subsets T of F_2^a (translation = twisting by a character, GL(a,2) = automorphisms
   of Q); the classical shapes are a point, a flat and a flat minus a point.
2. For every open order 4p with p prime: the pairs (a, |T|) with p = 2^a s + |T|, s > 1, and whether a non-classical
   class exists.   [M]"""
import itertools, json, sys
def affine_group(a):
    n = 2 ** a; out = []
    for cols in itertools.product(range(1, n), repeat=a):
        img = [0] * n
        for x in range(n):
            y = 0
            for i in range(a):
                if (x >> i) & 1: y ^= cols[i]
            img[x] = y
        if len(set(img)) == n:
            for b in range(n): out.append(tuple(img[x] ^ b for x in range(n)))
    return out
def is_flat(S):
    S = set(S); b = min(S); L = {x ^ b for x in S}
    return all((x ^ y) in L for x in L for y in L)
def kind(S, n):
    if len(S) == 1: return 'point'
    if is_flat(S): return 'flat'
    for z in range(n):
        if z not in S and is_flat(set(S) | {z}): return 'flat minus point'
    return 'NEW'
def classes(a, t):
    n = 2 ** a; G = affine_group(a); seen = set(); reps = []
    for S in itertools.combinations(range(n), t):
        fs = frozenset(S)
        if fs in seen: continue
        orb = set(frozenset(g[x] for x in S) for g in G)
        seen |= orb; reps.append((list(S), len(orb), kind(S, n)))
    return reps
if __name__ == '__main__':
    for a, t in ((3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (4, 3), (4, 7), (4, 11), (4, 15)):
        print(f'a={a} |T|={t}:', classes(a, t))
    op = json.load(open('data/open_orders_lt10000.json'))['open']
    isprime = lambda m: m > 1 and all(m % d for d in range(2, int(m ** .5) + 1))
    primes = [n // 4 for n in op if n % 4 == 0 and isprime(n // 4)]
    for r, mod, t in ((7, 16, 7), (11, 16, 11), (15, 32, 15)):
        print(f'open primes = {r} mod {mod} (|T| = {t} over Z_2^{mod.bit_length()-1}):', [p for p in primes if p % mod == r])
