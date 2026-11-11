#!/usr/bin/env python3
"""Single-character border for a Q-level solution (no core needed): fable_border4's CP-SAT model applied to the
compressed array M of a solution written by fable_q599.py.  A single-character border depends only on |T| and on the
4 x 4 blocks Gamma_rho, not on the shape of T (the characters enter as labels), so for |T| = 7 and
Gamma = GS(4,2,2,2) it is the same problem as Kumar's k = 8 border.  Check: K Hadamard and the compressed matrix
[[K, sqrt(s) E], [sqrt(s) F, M]] has Gram n I.  [M]
Usage: fable_border7.py solutions.json index s time_limit workers
"""
import sys, json, time
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2, fable_border4 as B4, fable_border6 as B6

def main():
    path = sys.argv[1]; idx = int(sys.argv[2]); s_ = int(sys.argv[3]); tl = float(sys.argv[4]); w = int(sys.argv[5])
    d, sol, M, F, t = B6.load(path, idx)
    K_ = 16; b = 4 * t; n = 4 * (K_ * s_ + t)
    G = B2.Grp2(4, 1); ch = G.chars()
    ys = [M[0, j * K_:(j + 1) * K_] for j in range(4)]
    low = [m for m in range(K_) if sum(int(ch[m] @ y) ** 2 for y in ys) == 4 * t]
    print('T (table of fable_border2):', low, ' Walsh quadruples on T:', [[int(ch[m] @ y) for y in ys] for m in low], flush=True)
    t0 = time.time()
    st, res = B4.solve(G, low, M, tl, workers=w, seed=idx)
    if res is None:
        print(f'single-character border: {st} ({time.time()-t0:.0f}s)', flush=True); return
    Km, Em, Zm = res
    okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=int))
    Fm = -Zm.T
    top = Km @ Km.T + s_ * (Em @ Em.T); cross = Km @ Fm.T + Em @ M.T; bot = s_ * (Fm @ Fm.T) + M @ M.T
    good = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * K_, dtype=int))
    print(f'single-character border: {st} ({time.time()-t0:.0f}s); K Hadamard {okK}; compressed matrix of order {b + 4*K_}: Gram = {n} I is {good}', flush=True)
    json.dump(dict(T=d['T'], s=s_, cells=sol['cells'], M=M.tolist(), K=Km.tolist(), E=Em.tolist(), Z=Zm.tolist()),
              open(f'logs/fable_fresh/q599_border_{idx}.json', 'w'))
if __name__ == '__main__':
    main()
