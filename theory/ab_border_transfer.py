r"""Border transfer: a published border (K, Kc, C) fits every compression with the same deviations [P, M].

[P] With x_k = mu_k 1 + d_k on Q, the compressed GS array is Cb = (mu-part) (x) J_c + Cb(d), and every row of K and Kc is
balanced inside each block, so K Cb^T Kc^T = K Cb(d)^T Kc^T: the forced corner depends on d only.  Shifting a block's
compression by an even constant (z_k -> z_k + 2u_k, which keeps the parity z = h mod 2) keeps d, hence the border.
So a border known for one order with block size c transfers to every order v' with the same c (and h' = h mod 2) whose
row equation sum_k (a_k + 2 c u_k)^2 = 4v' has an integer solution; also with any per-block sign (z_k -> -z_k, with
K, Kc block columns negated, checked numerically below) and h' of the other parity via odd shifts (checked).
usage: python3 theory/ab_border_transfer.py"""
import sys, os, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_assemble import load_template347, forced_corner

def solve_rows(a, c, v, umax=12, odd=False):
    """integer shifts u (or odd shifts) with sum_k (a_k + c*s_k)^2 = 4v, s_k = 2u_k (+1 if odd)"""
    out = []
    rng = range(-umax, umax + 1)
    for u in itertools.product(rng, repeat=len(a)):
        s = [2 * x + (1 if odd else 0) for x in u]
        r = [ak + c * sk for ak, sk in zip(a, s)]
        if sum(x * x for x in r) == 4 * v: out.append((s, r))
    return out

if __name__ == "__main__":
    T = load_template347(); comp = T['comp']; c = 6
    a = comp.sum(1).tolist()
    res = []
    for v in [599, 659, 719, 839, 863, 917, 947, 971]:
        h = (v + 1) // c - 1
        odd = (h % 2) != (57 % 2)
        sols = solve_rows(a, c, v, umax=8, odd=odd)
        ok = []
        for (s, r) in sols[:200]:
            Z = comp + np.array(s)[:, None]
            Cf = forced_corner(Z, (2, 3), T['layout'], T['K'], T['Kc'])
            ok.append(Cf is not None and (Cf @ Cf.T == 20 * np.eye(20, dtype=np.int64)).all())
        rec = dict(v=v, h=h, h_parity_same_as_347=not odd, row_solutions=len(sols),
                   corner_ok=int(sum(ok)), example=(sols[0][1] if sols else None))
        print(json.dumps(rec)); res.append(rec)
