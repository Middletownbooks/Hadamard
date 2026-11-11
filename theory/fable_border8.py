#!/usr/bin/env python3
"""Transplant of Kumar's k = 8 border (certificate of his H(3356), Apache-2.0, vatskishan03/hadamard) to a Q-level
solution of the 599 slice (fable_q599.py).

His border rows are, block by block, signed nontrivial characters of Z_2^3; a single-character border depends only on
the number of characters and on the 4 x 4 block Gamma = GS(d), so the 7 labels can be mapped to any 7 characters.
Here: E := his Rhat with the labels mapped to our T (characters of Z_2^4); for every signed permutation of our four
blocks and K in {+-X, +-X^T}:  Z := K^T E M^T / b; accept if Z is a +-1 matrix.  Then check that K is Hadamard and
that the compressed matrix [[K, sqrt(s) E], [sqrt(s) F, M]], F = -Z^T, has Gram n I.  [M]
Usage: fable_border8.py kumar_certificate.json solutions.json index s
"""
import sys, json, itertools
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2

def main():
    cert = json.load(open(sys.argv[1])); d = json.load(open(sys.argv[2])); idx = int(sys.argv[3]); s_ = int(sys.argv[4])
    pm = lambda row: np.array([1 if c == '+' else -1 for c in row])
    X = np.array([pm(r) for r in cert['border']['X']]); Rh = np.array([pm(r) for r in cert['border']['Rhat']])
    b = 28; K_ = 16; t = 7; n = 4 * (K_ * s_ + t)
    # characters of Z_2^3 in every ordering convention: identify each 8-block of Rhat as +- a nontrivial character
    ch8 = np.array([[(-1) ** bin(m & q).count('1') for q in range(8)] for m in range(8)])
    lab = np.zeros((b, 4), dtype=int); sg = np.zeros((b, 4), dtype=int)
    for r in range(b):
        for a in range(4):
            v = Rh[r, a * 8:(a + 1) * 8]; hit = [(m, e) for m in range(1, 8) for e in (1, -1) if np.array_equal(v, e * ch8[m])]
            assert len(hit) == 1, 'block of Rhat is not a signed nontrivial character in this ordering'
            lab[r, a], sg[r, a] = hit[0]
    print('Rhat: every block is a signed nontrivial character of Z_2^3', flush=True)
    # our data
    G = B2.Grp2(4, 1); ch = G.chars()
    sol = d['solutions'][idx]; cellorb = d['cell_orbits']
    pos = [G.Qel.index(tuple((q >> (3 - i)) & 1 for i in range(4))) for q in range(K_)]
    ys = []
    for c in sol['cells']:
        y = np.zeros(K_, dtype=np.int64)
        for j, ob in enumerate(cellorb):
            for q in ob: y[pos[q]] = c[j]
        ys.append(y)
    low = [m for m in range(K_) if sum(int(ch[m] @ y) ** 2 for y in ys) == 4 * t]
    assert len(low) == 7
    E = np.zeros((b, 4 * K_), dtype=np.int64)
    for r in range(b):
        for a in range(4):
            E[r, a * K_:(a + 1) * K_] = sg[r, a] * ch[low[lab[r, a] - 1]]
    assert np.array_equal(E @ E.T, 4 * K_ * np.eye(b, dtype=int)), 'E not orthogonal'
    found = 0
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((1, -1), repeat=4):
            yv = np.array([signs[i] * ys[perm[i]] for i in range(4)])
            M = FB.gs_array(G, yv)
            N_ = E @ M.T
            for nm, Kc in (('X', X), ('X^T', X.T), ('-X', -X), ('-X^T', -X.T)):
                Zn = Kc.T @ N_
                if np.all(np.abs(Zn) == b):
                    Z = Zn // b; F = -Z.T
                    top = Kc @ Kc.T + s_ * (E @ E.T); cross = Kc @ F.T + E @ M.T; bot = s_ * (F @ F.T) + M @ M.T
                    good = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * K_, dtype=int))
                    found += 1
                    if found == 1:
                        print(f'border found: blocks {perm} signs {signs}, K = {nm}; compressed matrix of order {b + 4*K_} has Gram {n} I: {good}', flush=True)
                        if good:
                            json.dump(dict(T=d['T'], s=s_, cells=sol['cells'], block_perm=perm, block_signs=signs,
                                           y=yv.tolist(), K=Kc.tolist(), E=E.tolist(), Z=Z.tolist()),
                                      open(f'logs/fable_fresh/q599_border_{idx}.json', 'w'))
    print('variants with a +-1 Z:', found)
if __name__ == '__main__':
    main()
