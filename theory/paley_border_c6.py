r"""Which Q-compressions does the PUBLISHED 347 border (c = 6, b = 5) accept, for the open v = 599 (c = 6, b = 5)?

The border rows K, columns Kc and layout come from data/template347.npz (scripts/bge3_template347.py).  For a compression
xbar (4 blocks x 6 cosets of Q = Z2 x Z3), the corner is forced: C = -K Cb^T Kc^T / (4c) (scripts/bge3_border_check.py,
verified exactly on 347).  A compression is admissible iff every entry of K Cb^T Kc^T is +-4c (then C is +-1 and
sqrt(4b)-orthogonal automatically).  We enumerate compressions "constant + deviation": xbar_k = c_k + dev_k with
dev_k in {0, 2}^6 or {0, -2}^6 patterns of small weight, sum_k |dev^_k(xi)|^2 = 4b on Q^ \ 1, sum_k r_k^2 = 4v,
all values odd (|H| = 99 odd), and report which pass the corner test.
usage: python3 theory/paley_border_c6.py [v]"""
import sys, os, json, itertools as it
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from bge3_border_check import compressed_core, T as TPL

def corner_ok(xbar, c=6):
    K, Kc = TPL['K'], TPL['Kc']
    Cb = compressed_core(np.asarray(xbar), (2, 3))
    M = K @ Cb.T @ Kc.T
    return bool((np.abs(M) == 4 * c).all())

def main(v=599):
    b, c = 5, 6
    comp347 = TPL['comp']
    out = {"v": v, "control_347_ok": corner_ok(comp347)}
    # characters of Q = Z2 x Z3 (index = 3a + b)
    co = np.array([[a, bb] for a in range(2) for bb in range(3)])
    w3 = np.exp(2j * np.pi / 3)
    chars = np.array([[(-1) ** (x * a) * w3 ** (y * bb) for (a, bb) in co] for (x, y) in co])
    devs = [np.array(d) for d in it.product((0, 2), repeat=6) if sum(d) <= 6] + [-np.array(d) for d in it.product((0, 2), repeat=6) if 0 < sum(d) <= 6]
    devs = list({tuple(d) for d in devs}); devs = [np.array(d) for d in devs]
    spec = {tuple(d): np.round(np.abs(chars @ d) ** 2, 6) for d in devs}
    found = []
    # deviation 4-tuples with sum_k |dev^(xi)|^2 = 20 on xi != 1 (use the spectrum of the deviation itself)
    D = [(tuple(d), spec[tuple(d)]) for d in devs]
    for (d1, s1), (d2, s2), (d3, s3), (d4, s4) in it.combinations_with_replacement(D, 4):
        tot = s1 + s2 + s3 + s4
        if not np.allclose(tot[1:], 4 * b): continue
        for perm in set(it.permutations([d1, d2, d3, d4])):
            dv = [np.array(t) for t in perm]
            for cs in it.product(range(-15, 16, 2), repeat=4):
                xbar = [cs[k] + dv[k] for k in range(4)]
                if any((x % 2 == 0).any() for x in xbar): continue
                if sum(int(x.sum()) ** 2 for x in xbar) != 4 * v: continue
                if corner_ok(xbar): found.append([x.tolist() for x in xbar])
            if len(found) > 200: break
        if len(found) > 200: break
    out.update(admissible_found=len(found), examples=found[:6])
    print(json.dumps(out))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 599)
