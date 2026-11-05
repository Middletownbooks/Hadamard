"""Exact quotient and lag-zero gates for block-exchanging symmetry. No core search.
Run from repo root: python3 theory/openai_deficit_symmetry.py --m 5
"""
import argparse
import itertools as it
import json
import math
from collections import Counter, defaultdict
from functools import lru_cache

SW = (0, 2, 1, 3)
W = tuple(tuple((-1) ** ((a & b).bit_count() % 2) for b in range(4)) for a in range(4))


def q_gate(m):
    # Each nontrivial Q PSD is 12. Thus each individual Walsh coefficient
    # is in {-2,0,2}; K has odd size, and each Q cell sum is odd.
    v = 4 * m * 43 + 3
    blocks = []
    for a0 in range(0, math.isqrt(4 * v) + 1, 2):
        for tail in it.product((-2, 0, 2), repeat=3):
            a = (a0,) + tail
            z4 = tuple(sum(W[j][k] * a[j] for j in range(4)) for k in range(4))
            if any(t % 4 or (t // 4) % 2 == 0 for t in z4):
                continue
            blocks.append((a, tuple(t // 4 for t in z4)))
    fixed = [b for b in blocks if b[0][1] == b[0][2]]
    pairs = defaultdict(list)
    for ia, a in enumerate(fixed):
        for b in fixed[ia:]:
            key = tuple(a[0][j] ** 2 + b[0][j] ** 2 for j in range(4))
            pairs[key].append((a[1], b[1]))
    target = (4 * v, 12, 12, 12)
    sol = []
    for a, z in blocks:
        key = tuple(target[j] - a[j] ** 2 - a[SW[j]] ** 2 for j in range(4))
        for za, zb in pairs.get(key, []):
            sol.append((za, zb, z, tuple(z[SW[j]] for j in range(4))))
    return sol


def lag0_gate(m):
    # x0,x1 fixed by tau; x3=tau*x2. Independent H-cell orbits have
    # weights |orbit| on blocks 0,1 and weight 2 on every block-2 cell.
    tau = [SW[q] * m + (-u) % m for q in range(4) for u in range(m)]
    orbits = sorted({tuple(sorted({k, tau[k]})) for k in range(4 * m)})
    weights = [len(o) for _ in range(2) for o in orbits] + [2] * (4 * m)
    alphabet = sorted({e + 7 * s for e in (-1, 1) for s in range(-6, 7, 2)})
    target = 16 * m * 43 - 12 * 42 - 16 * m
    # Exact bounded DP, with a witness; each independent cell uses one letter.
    states = {0: []}
    for w in weights:
        new = {}
        for en, witness in states.items():
            for a in alphabet:
                en2 = en + w * (a * a - 1)
                if en2 <= target and en2 not in new:
                    new[en2] = witness + [a]
        states = new
    return dict(weights=dict(Counter(weights)), excess_energy=target,
                reachable=target in states, witness=states.get(target),
                independent_cells=len(weights), weights_list=weights)


def joint_gate(m, solutions):
    """Exact energy reachability conditional on each Q tuple. No off-zero H equations."""
    cap = 16 * m * 42 - 504
    alphabet = sorted({e + 7*s for e in (-1, 1) for s in range(-6, 7, 2)})

    @lru_cache(None)
    def energies(total, weights):
        states = {(0, 0)}
        for w in weights:
            states = {(s+w*a, en+w*(a*a-1)) for s, en in states for a in alphabet
                      if en+w*(a*a-1) <= cap}
        return {en for s, en in states if s == total}

    def add(a, b):
        return {x+y for x in a for y in b if x+y <= cap}

    @lru_cache(None)
    def block(z, fixed):
        if fixed:
            ans = energies(z[0], (1,) + (2,)*((m-1)//2))
            ans = add(ans, energies(z[3], (1,) + (2,)*((m-1)//2)))
            ans = add(ans, {2*e for e in energies(z[1], (1,)*m) if 2*e <= cap})
        else:
            ans = {0}
            for t in z:
                ans = add(ans, energies(t, (1,)*m))
            ans = {2*e for e in ans if 2*e <= cap}
        return ans
    survivors = []
    for zs in solutions:
        en = add(add(block(zs[0], True), block(zs[1], True)), block(zs[2], False))
        if cap in en:
            survivors.append(zs)
    return survivors


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--m', type=int, default=5)
    args = p.parse_args()
    sol = q_gate(args.m)
    lag = lag0_gate(args.m)
    joint = joint_gate(args.m, sol)
    print(json.dumps(dict(m=args.m, v=4*args.m*43+3, seed=None,
                         q_solutions=len(sol), q_example=sol[:1], lag0=lag,
                         joint_survivors=len(joint), joint_examples=joint[:3])))


if __name__ == '__main__':
    main()
