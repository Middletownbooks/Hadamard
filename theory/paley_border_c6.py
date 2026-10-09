r"""Which Q-compressions does the PUBLISHED 347 border (c = 6, b = 5) accept, for the open v = 599 (c = 6, b = 5)?

Border rows K / columns Kc and the core layout come from data/template347.npz (scripts/bge3_template347.py).  For a
compression xbar (4 blocks x 6 cosets of Q = Z2 x Z3) the corner is forced, C = -K Cb^T Kc^T / (4c)
(scripts/bge3_border_check.py, exact on 347).  Admissible iff every entry of K Cb^T Kc^T is +-4c.
Fact used [P]: in 347's layout every rho, gam is the trivial character, and K, Kc are balanced inside each block, so a
block-constant compression is annihilated: admissibility depends on xbar only modulo block constants.  Hence we
(1) enumerate deviation patterns dev_k in {0, 2}^6 or -{0,2}^6 (up to block constants) with
    sum_k |dev^_k(xi)|^2 = 4b on Q^ \ 1, (2) test the corner once per pattern, (3) list the block constants that give
    sum_k r_k^2 = 4v with odd entries (|H| = 99 odd) -- every such compression is admissible.
usage: python3 theory/paley_border_c6.py [v]"""
import sys, os, json, math, itertools as it
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from bge3_border_check import compressed_core, T as TPL

KB, KCB = np.array(TPL['K']), np.array(TPL['Kc'])               # cached (NpzFile re-reads on every access)

def corner_ok(xbar, c=6):
    Cb = compressed_core(np.asarray(xbar), (2, 3))
    M = KB @ Cb.T @ KCB.T
    return bool((np.abs(M) == 4 * c).all())

def main(v=599):
    b, c = 5, 6
    comp347 = np.array(TPL['comp'])
    # constants are annihilated: check on 347 itself and on random constant shifts
    rng = np.random.default_rng(0)
    shift_ok = all(corner_ok(comp347 + 2 * rng.integers(-5, 6, (4, 1))) for _ in range(20))
    co = [(a, bb) for a in range(2) for bb in range(3)]
    w3 = np.exp(2j * np.pi / 3)
    chars = np.array([[(-1) ** (x * a) * w3 ** (y * bb) for (a, bb) in co] for (x, y) in co])
    devs = sorted({tuple(int(s * x) for x in d) for d in it.product((0, 2), repeat=6) for s in (1, -1)})
    spec = {d: np.round(np.abs(chars @ np.array(d)) ** 2, 6)[1:] for d in devs}
    # pair-sum hashing on the spectrum (xi != 1) to find 4-tuples with total 4b
    byspec = defaultdict(list)
    for d in devs: byspec[tuple(spec[d])].append(d)
    keys = list(byspec)
    pair = defaultdict(list)
    for i, k1 in enumerate(keys):
        for k2 in keys[i:]:
            s = tuple(np.round(np.array(k1) + np.array(k2), 6))
            if all(x <= 4 * b + 1e-9 for x in s): pair[s].append((k1, k2))
    target = tuple([float(4 * b)] * 5)
    quads = set()
    for s, lst in pair.items():
        need = tuple(np.round(np.array(target) - np.array(s), 6))
        if need in pair:
            for (k1, k2) in lst:
                for (k3, k4) in pair[need]:
                    quads.add(tuple(sorted((k1, k2, k3, k4))))
    admissible = []; tested = 0
    for qd in quads:
        for spec_perm in set(it.permutations(qd)):
            for dv in it.product(*[byspec[k] for k in spec_perm]):
                tested += 1
                if corner_ok([np.array(x) for x in dv]): admissible.append(dv)
    # constants giving sum r^2 = 4v (odd entries): xbar_k = const_k + dev_k
    sols = []
    for dv in admissible:
        sd = [sum(d) for d in dv]
        for cs in it.product(range(-17, 18, 2), repeat=4):
            if sum((6 * cs[k] + sd[k]) ** 2 for k in range(4)) == 4 * v:
                sols.append(dict(dev=[[int(x) for x in d] for d in dv], const=[int(x) for x in cs], rows=[int(6 * cs[k] + sd[k]) for k in range(4)]))
                break
    # alpha bits at fixed compression, assuming a balanced gadget (y_Q = alpha_Q, as in 359): sum log2 C(|H'|, #minus)
    hs = 9
    for sl in sols:
        ys = [[sl['const'][k] + x for x in sl['dev'][k]] for k in range(4)]
        sl['alpha_bits'] = round(sum(math.log2(math.comb(hs, (hs - y) // 2)) for row in ys for y in row if abs(y) <= hs), 1)
    ab = sorted(sl['alpha_bits'] for sl in sols)
    out = dict(v=v, control_347_ok=corner_ok(comp347), constant_shifts_annihilated=shift_ok,
               deviation_patterns_tested=tested, admissible_patterns=len(admissible),
               admissible_with_rowsums_4v=len(sols),
               alpha_bits_min=ab[0] if ab else None, alpha_bits_median=ab[len(ab) // 2] if ab else None,
               alpha_bits_max=ab[-1] if ab else None, examples=sols[:4])
    print(json.dumps(out))
    with open(os.path.join(ROOT, 'theory', 'logs', f'paley_border_c6_{v}_admissible.json'), 'w') as fh:
        json.dump(sols, fh)

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 599)
