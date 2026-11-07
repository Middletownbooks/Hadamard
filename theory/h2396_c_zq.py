r"""Part C of report/prompt_599_plugin.md: the index-4 Szolloesi plug-in for n = 1291 (core of order N = 323 on
G = Z19 x F17), gate C2 -- the Z_q level (compression along K = Z19 to F17) [P, M].

Credit: the index-r plug-in class, its balance identity and the 10^8.8 estimate are K. Kumar's (vatskishan03/hadamard,
search/NOTES.md section 82; Apache-2.0); the core conventions are those of his verify/verify.py
(_szollosi_zero_diagonal_core): A, B symmetric, A(0) = B(0) = 0, +-1 elsewhere; C, D +-1; sum PAF = 0 off 0, so
sum PAF = (4N - 2) delta = 1290 delta.  This file is new.

Plug-in (index r = 4): H = {1} x S4, S4 = {1, 4, 13, 16} (the fourth powers of F17^*, -1 in S4); a block is x0 on
K x {0} and x_j on K x C_j, C_j = g^j S4 (g = 3).  Compression along K: w(0) = a = sum_k x0(k), w = W_j on C_j,
W_j = sum_k x_j(k).  With the cyclotomic numbers (i, j) = #{y in C_i : y + 1 in C_j} [P]:
    PAF_w(0) = a^2 + 4 sum_j W_j^2,
    PAF_w(z) = 2 a W_l + sum_{j,j'} W_j W_j' (j - l, j' - l)        (z in C_l),
and the Z17-level system is  sum_i PAF_{w_i}(0) = 1290,  sum_i PAF_{w_i}(z) = 0 for z in C_0..C_3 (5 equations).
Ranges [P]: C, D: a, W_j odd in [-19, 19].  A, B (symmetric, vanishing at (0, 0)): a = 2 * (sum of 9 signs) (so
a = 2 mod 4, |a| <= 18) and W_j = x_j(0) + 2 * (sum of 9 signs) (odd, |W_j| <= 19).
Each solution is checked against a direct PAF on F17 of the expanded w.
usage: python3 theory/h2396_c_zq.py [OUT.json]"""
import sys, json, itertools
from collections import defaultdict
import numpy as np

q, nK = 17, 19
g = 3
S4 = sorted({pow(x, 4, q) for x in range(1, q)})
assert S4 == [1, 4, 13, 16]
cls = {}
for j in range(4):
    for s in S4: cls[(pow(g, j, q) * s) % q] = j
cyc = [[sum(1 for y in range(1, q) if cls[y] == i and (y + 1) % q and cls[(y + 1) % q] == j) for j in range(4)] for i in range(4)]

def paf_vec(a, W):
    """(PAF(0), PAF(C_0..C_3)) of w = a delta + sum W_j 1_{C_j}, by the cyclotomic formula"""
    out = [a * a + 4 * sum(x * x for x in W)]
    for l in range(4):
        out.append(2 * a * W[l] + sum(W[j] * W[jj] * cyc[(j - l) % 4][(jj - l) % 4] for j in range(4) for jj in range(4)))
    return tuple(out)

def paf_direct(a, W):
    w = [a] + [W[cls[y]] for y in range(1, q)]
    P = [sum(w[y] * w[(y + z) % q] for y in range(q)) for z in range(q)]
    rep = [0] + [next(y for y in range(1, q) if cls[y] == l) for l in range(4)]
    assert all(P[z] == P[rep[cls[z] + 1]] for z in range(1, q)), "PAF not constant on classes"
    return tuple(P[r] for r in rep)

def main(outp=None):
    odd19 = list(range(-19, 20, 2))
    aAB = [a for a in range(-18, 19) if a % 4 == 2 or a % 4 == -2 % 4]
    blocksCD = defaultdict(int); blocksAB = defaultdict(int)
    for a in odd19:
        for W in itertools.product(odd19, repeat=4):
            P = paf_vec(a, W)
            if P[0] <= 1290: blocksCD[P] += 1
    for a in aAB:
        for W in itertools.product(odd19, repeat=4):
            P = paf_vec(a, W)
            if P[0] <= 1290: blocksAB[P] += 1
    # spot-check the cyclotomic formula against the direct PAF
    rng = np.random.default_rng(1)
    for _ in range(200):
        a = int(rng.choice(odd19)); W = [int(x) for x in rng.choice(odd19, 4)]
        assert paf_vec(a, W) == paf_direct(a, W)
    tgt = (1290, 0, 0, 0, 0)
    AB = defaultdict(int)
    for k1, c1 in blocksAB.items():
        for k2, c2 in blocksAB.items():
            s = tuple(x + y for x, y in zip(k1, k2))
            if s[0] <= 1290: AB[s] += c1 * c2
    tot = 0; ex = []
    for k3, c3 in blocksCD.items():
        for k4, c4 in blocksCD.items():
            s = tuple(x + y for x, y in zip(k3, k4))
            need = tuple(t - x for t, x in zip(tgt, s))
            if need in AB:
                tot += AB[need] * c3 * c4
                if len(ex) < 5: ex.append(dict(C=k3, D=k4, AB=need))
    rec = dict(cyclotomic_numbers_order4_q17=cyc, blocks_AB=len(blocksAB), blocks_CD=len(blocksCD),
               ordered_Zq_level_solutions=tot, examples=ex)
    print(json.dumps(rec))
    if outp: json.dump(rec, open(outp, 'w'), indent=1)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
