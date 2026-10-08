"""Interpret Traces generators of Aut(graph(H)) as signed-permutation pairs (P, Q) and verify P H Q^T = H.

Graph convention (hadaut_traces.c): rows r+ = r, r- = n+r; cols c+ = 2n+c, c- = 3n+c.
tau = the antipodal map (-I, -I), which is always an automorphism.
A Williamson-array form needs a subgroup Q8 <= Aut whose central involution is tau and which acts
fixed-point-freely on rows; in particular an element g with g^2 = tau.
usage: check_autgens.py n generator_file
"""
import gzip
import sys

import numpy as np

n = int(sys.argv[1])
f = sys.argv[2]
H = np.loadtxt(gzip.open(f'/home/user/Hadamard/data/hadamard_{n}.csv.gz'), delimiter=',', dtype=np.int64)
gens = [np.array(list(map(int, l.split()[1:]))) for l in open(f) if l.startswith('GEN')]
order_line = [l for l in open(f) if l.startswith('ORDER')]


def to_signed(perm):
    P = np.zeros((n, n), dtype=np.int64)
    Q = np.zeros((n, n), dtype=np.int64)
    for r in range(n):
        t = perm[r]
        assert t < 2 * n
        P[t % n, r] = 1 if t < n else -1
    for c in range(n):
        t = perm[2 * n + c] - 2 * n
        assert 0 <= t < 2 * n
        Q[t % n, c] = 1 if t < n else -1
    return P, Q


def perm_order(g):
    x = np.arange(len(g))
    k = 0
    while True:
        x = g[x]
        k += 1
        if np.array_equal(x, np.arange(len(g))):
            return k


tau = np.concatenate([np.arange(n, 2 * n), np.arange(n), np.arange(3 * n, 4 * n), np.arange(2 * n, 3 * n)])
print('group order line:', order_line[0].strip() if order_line else '?')
# close the generated group (it is small) to list all elements
elems = {tuple(np.arange(4 * n))}
frontier = [np.arange(4 * n)]
while frontier:
    nf = []
    for x in frontier:
        for g in gens:
            y = g[x]
            t = tuple(y)
            if t not in elems:
                elems.add(t)
                nf.append(y)
    frontier = nf
print('closure size', len(elems))
for e in elems:
    g = np.array(e)
    P, Q = to_signed(g)
    ok = np.array_equal(P @ H @ Q.T, H)
    print(f'  elem: automorphism {ok}; order {perm_order(g)}; is_tau {np.array_equal(g, tau)}; '
          f'square_is_tau {np.array_equal(g[g], tau)}; rows fixed up to sign {sum(1 for r in range(n) if g[r] % n == r)}')
