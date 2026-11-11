#!/usr/bin/env python3
"""Single-character borders with a cyclic symmetry of order t on the labels (t = |T| prime), exact search.

Setting (report section 7): Q = Z_2^a, character set T with |T| = t, and the compressed core acts on the border space
as Gamma (x) I_t (every block's Walsh value constant on T; true for 2592 of the 2688 Q-level solutions of the 599
slice).  A single-character row is (label l_c in Z_t, sign s_c) for each block c = 0..3.  Rows come in orbits of the
shift l -> l + 1.  Conditions, for representatives e (rows of E) and z (rows of Z):
  orthogonality of E:  for every delta in Z_t,  sum_{c : l_c - l'_c = delta} s_c s'_c = 0   (orbits of two reps;
                       inside one orbit automatic),  same for Z;
  K entries:           for every delta in Z_t,  sum_{(c, c') : l_c - l'_{c'} = delta} s_c s'_{c'} Gamma[c][c'] = +-4 * (norm)
Search: all reps (t^4 2^4 / (2t)), compatibility graphs, then 4 pairwise orthogonal E-reps with a common set of
compatible Z-reps containing 4 pairwise orthogonal ones.  Every border found is mapped back to +-1 matrices over
Z_2^4, and K Hadamard and the Gram of the compressed matrix are checked.  [M]
Usage: fable_border9.py solutions.json index s [max_borders]
"""
import sys, json, itertools, time
import numpy as np
sys.path.insert(0, 'theory')
import fable_border as FB, fable_border2 as B2, fable_border6 as B6

def main():
    path = sys.argv[1]; idx = int(sys.argv[2]); s_ = int(sys.argv[3]); maxb = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    d, sol, M, F, t = B6.load(path, idx)
    K_ = 16; b = 4 * t; n = 4 * (K_ * s_ + t)
    G = B2.Grp2(4, 1); ch = G.chars()
    ys = [M[0, j * K_:(j + 1) * K_] for j in range(4)]
    low = [m for m in range(K_) if sum(int(ch[m] @ y) ** 2 for y in ys) == 4 * t]
    # Gamma per character of T: Gamma_m[c][c'] = (1/k) chi_m(block c) M^T chi_m(block c')
    def vec(c, m):
        v = np.zeros(4 * K_, dtype=np.int64); v[c * K_:(c + 1) * K_] = ch[m]; return v
    Gam = [np.array([[int(vec(c, m) @ M.T @ vec(c2, m)) // K_ for c2 in range(4)] for c in range(4)]) for m in low]
    assert all(np.array_equal(Gam[0], g) for g in Gam), 'Gamma not constant on T for this solution'
    Gm = Gam[0]
    # cross terms between different characters vanish
    assert int(vec(0, low[0]) @ M.T @ vec(1, low[1])) == 0
    print('Gamma =', Gm.tolist(), ' Gamma Gamma^T =', (Gm @ Gm.T).tolist(), flush=True)
    reps = []
    seen = set()
    for lab in itertools.product(range(t), repeat=4):
        for sg in itertools.product((1, -1), repeat=4):
            key = None
            orb = []
            for sh in range(t):
                for gs in (1, -1):
                    orb.append((tuple((l + sh) % t for l in lab), tuple(gs * x for x in sg)))
            k0 = min(orb)
            if k0 in seen: continue
            seen.add(k0); reps.append(k0)
    R = len(reps); print(R, 'orbit representatives per side', flush=True)
    L = np.array([r[0] for r in reps]); S = np.array([r[1] for r in reps])
    # orthogonality between reps
    dl = (L[:, None, :] - L[None, :, :]) % t                    # (R, R, 4) same block c
    ss = S[:, None, :] * S[None, :, :]
    orth = np.ones((R, R), dtype=bool)
    for delta in range(t):
        orth &= ((ss * (dl == delta)).sum(2) == 0)
    np.fill_diagonal(orth, False)
    # compatibility e-rep x z-rep
    comp = np.ones((R, R), dtype=bool)
    for delta in range(t):
        acc = np.zeros((R, R), dtype=np.int64)
        for c in range(4):
            for c2 in range(4):
                if Gm[c, c2] == 0: continue
                acc += Gm[c, c2] * (S[:, None, c] * S[None, :, c2]) * (((L[:, None, c] - L[None, :, c2]) % t) == delta)
        comp &= (np.abs(acc) == 4 * t) if False else (np.abs(acc) == 4)
    print('orthogonal rep pairs:', int(orth.sum()) // 2, ' compatible (e, z) pairs:', int(comp.sum()), flush=True)
    found = 0; t0 = time.time()
    def cliques4(mask_nodes, adj):
        idx = np.nonzero(mask_nodes)[0]
        for i in idx:
            m1 = mask_nodes & adj[i]; m1[:i + 1] = False
            for j in np.nonzero(m1)[0]:
                m2 = m1 & adj[j]; m2[:j + 1] = False
                for k in np.nonzero(m2)[0]:
                    m3 = m2 & adj[k]; m3[:k + 1] = False
                    for l in np.nonzero(m3)[0]:
                        yield (i, j, k, l)
    tried = 0
    for Eq in cliques4(np.ones(R, dtype=bool), orth):
        tried += 1
        zmask = comp[Eq[0]] & comp[Eq[1]] & comp[Eq[2]] & comp[Eq[3]]
        if zmask.sum() < 4: continue
        for Zq in cliques4(zmask, orth):
            # build +-1 matrices over Z_2^4
            def rows(q):
                out = []
                for r in q:
                    for sh in range(t):
                        v = np.zeros(4 * K_, dtype=np.int64)
                        for c in range(4):
                            v[c * K_:(c + 1) * K_] = S[r, c] * ch[low[(L[r, c] + sh) % t]]
                        out.append(v)
                return np.array(out)
            E = rows(Eq); Z = rows(Zq)
            Km = (E @ M.T @ Z.T) // (4 * K_)
            okK = np.array_equal(np.abs(Km), np.ones_like(Km)) and np.array_equal(Km @ Km.T, b * np.eye(b, dtype=int))
            Fm = -Z.T
            top = Km @ Km.T + s_ * (E @ E.T); cross = Km @ Fm.T + E @ M.T; bot = s_ * (Fm @ Fm.T) + M @ M.T
            good = np.array_equal(top, n * np.eye(b, dtype=int)) and not cross.any() and np.array_equal(bot, n * np.eye(4 * K_, dtype=int))
            found += 1
            print(f'border {found}: E reps {[reps[i] for i in Eq]}, Z reps {[reps[i] for i in Zq]}; K Hadamard {okK}; '
                  f'compressed matrix of order {b + 4*K_} has Gram {n} I: {good} ({time.time()-t0:.0f}s)', flush=True)
            if good and found == 1:
                json.dump(dict(T=d['T'], s=s_, cells=sol['cells'], M=M.tolist(), K=Km.tolist(), E=E.tolist(), Z=Z.tolist()),
                          open(f'logs/fable_fresh/q599_border_{idx}.json', 'w'))
            if found >= maxb: return
            break
    print(f'done: {found} borders; {tried} E-quadruples tried ({time.time()-t0:.0f}s)')
if __name__ == '__main__':
    main()
