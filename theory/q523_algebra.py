r"""Checks of the exact stage-2 algebra (report/paley_523_answer.md, section 1) [M].

1. eta-basis arithmetic of O_K (K = degree-s subfield of Q(zeta_q)) against complex numerics.
2. P(h) = sum_g conj(Z(g)) Z(g + h) equals the eta-coordinates (R(h, l) - PAF_x(h, 0))_l of the 2-D PAF (FFT on Z_F x Z_q),
   on random slice configurations with the 523 stage-1 y (F = 18, q = 29, s = 4).
3. On the verified H(164) core (v = 41, Paley(5) on Z_8, s = 2): sum_k P_k = 4v delta exactly, and stage 1 holds.
4. Coset alphabets: alpha is a function of y (s = 4, q = 29 and q = 13); quadratic/quartic coordinates
   e = c0 - c1 + c2 - c3, e1 = (c0 - c2) + i (c1 - c3) per letter; the per-coset Parseval identity
   y^2 + f sum_t |Z_t|^2 = q^2.
usage: python3 theory/q523_algebra.py"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_lib import Cyclo, zcoords, pz, expand_x, paf2, target_P
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def check_arith(cy):
    e = cy.eta
    assert abs(e.sum() + 1) < 1e-9
    for i in range(cy.s):
        for j in range(cy.s):
            assert abs(cy.mult[i, j] @ e - e[i] * e[j]) < 1e-8
            assert abs(cy.T[i, j] @ e - e[i] * np.conj(e[j])) < 1e-8
    for j in range(cy.s): assert abs(np.conj(e[j]) - e[(j + cy.cneg) % cy.s]) < 1e-9
    return True

def random_config(cy, Y, rng):
    AC = np.zeros((4, len(Y[0]), 1 + cy.s), np.int64)
    for k in range(4):
        for g, y in enumerate(Y[k]):
            al = cy.alphabet(y); AC[k, g] = al[rng.integers(len(al))]
    return AC

def check_P_vs_fft(cy, Y, n=20, seed=1):
    rng = np.random.default_rng(seed); F = len(Y[0])
    for _ in range(n):
        AC = random_config(cy, Y, rng)
        Z = zcoords(AC); P = pz(Z, cy.T).sum(0)
        X = expand_x(AC, cy); Q = paf2(X)
        R = np.stack([Q[:, cy.reps[l]] for l in range(cy.s)], 1)
        for l in range(cy.s): assert (np.stack([Q[:, z] for z in cy.C[l]], 1) == R[:, [l]]).all()   # H-invariance
        assert (P == R - Q[:, [0]]).all()
        # stage 1: coset sums
        assert (Q.sum(1) == sum(np.array([[sum(int(y[g]) * int(y[(g + h) % F]) for g in range(F)) for h in range(F)] for y in Y]))).all()
    return True

def check_v41():
    L = open(os.path.join(ROOT, 'data', 'core_v41_b1_paley_p5_symX_control.txt')).read().split('\n')
    X = np.array([[1 if ch == '+' else -1 for ch in l.split('=', 1)[1].strip()] for l in L if l[:2] in ('A=', 'B=', 'C=', 'D=')])
    F, q, s, v = 8, 5, 2, 41
    X = X.reshape(4, F, q); cy = Cyclo(q, s)
    AC = np.zeros((4, F, 1 + s), np.int64); AC[:, :, 0] = X[:, :, 0]
    for j in range(s):
        cols = X[:, :, cy.C[j]]
        assert (cols == cols[:, :, :1]).all()                       # slice: constant on classes
        AC[:, :, 1 + j] = cols[:, :, 0]
    P = pz(zcoords(AC), cy.T).sum(0)
    Y = X.sum(2)
    S1 = sum(np.array([[int((y * np.roll(y, -h)).sum()) for h in range(F)] for y in Y]))
    T1 = np.full(F, -4 * q); T1[0] = 4 * F * q - 4 * (q - 1)
    return dict(P_ok=bool((P == target_P(F, v, s)).all()), stage1_ok=bool((S1 == T1).all()), y=Y.tolist())

def alphabet_table(q, s):
    cy = Cyclo(q, s); rows = []
    for y in cy.letters():
        al = cy.alphabet(y)
        alphas = sorted(set(al[:, 0].tolist()))
        ent = []
        for r in al:
            c = r[1:]
            e = int(c[0] - c[1] + c[2] - c[3]) if s == 4 else None
            e1 = (int(c[0] - c[2]), int(c[1] - c[3])) if s == 4 else None
            Zt = [cy.num(r[1:] - r[0], t) for t in range(s)]
            assert abs(y * y + cy.f * sum(abs(z) ** 2 for z in Zt) - q * q) < 1e-6        # Parseval per coset
            ent.append(dict(c=c.tolist(), e=e, e1=e1, Z=(r[1:] - r[0]).tolist()))
        rows.append(dict(y=y, n=len(al), alpha=alphas, configs=ent))
    return rows

if __name__ == "__main__":
    out = {}
    for q in (5, 13, 17, 29, 41):
        for s in (2, 4):
            if (q - 1) % s: continue
            cy = Cyclo(q, s); check_arith(cy)
            out[f"q{q}_s{s}"] = dict(g0=cy.g0, class_of_minus1=cy.cneg, eta=[complex(round(z.real, 6), round(z.imag, 6)) for z in cy.eta].__repr__())
    print(json.dumps(dict(arith_ok=True, fields=out)))
    Y = json.load(open(os.path.join(ROOT, 'logs', 'cyclo', 'stage1_cpsat_523.json')))['y']
    print(json.dumps(dict(P_vs_fft_523=check_P_vs_fft(Cyclo(29, 4), Y))))
    print(json.dumps(dict(v41=check_v41())))
    for q in (29, 13):
        for r in alphabet_table(q, 4):
            print(json.dumps(dict(q=q, y=r['y'], n=r['n'], alpha=r['alpha'],
                                  e_e1=sorted({(c['e'], c['e1']) for c in r['configs']}))))
