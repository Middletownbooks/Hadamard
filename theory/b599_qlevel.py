r"""All Q-level solutions of the c = 6 family (Q = Z_2 x Z_3) up to symmetry [M].

Q-level (report/alpoge_design_answer.md 1.2): sum_k PAF_Q(z_k) = T_Q, i.e. sum_k |z_k^(psi)|^2 = 20 (psi != 1) and
sum_k r_k^2 = 4v.  Every z_k is odd-valued, so z_k = c_k 1 + 2 u_k with c_k odd and u_k integral (min u_k = 0):
sum_k |u_k^(psi)|^2 = 5 for every psi != 1.  The three classes of psi != 1 (psi ~ conj psi):
  P1 = (eps = -1, trivial on Z_3),  P2 = (eps = 1, omega),  P3 = (eps = -1, omega).
Symmetries used for classes: per-block translation by Q and negation (z -> -z: u -> max u - u, c -> -c - 2 max u),
block permutations, and the automorphism t -> -t of Z_3 (all are symmetries of stage 1; report/b599_answer.md 2).
usage: python3 theory/b599_qlevel.py [--v 599]"""
import sys, os, json, itertools, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

QE = [(e, t) for e in range(2) for t in range(3)]          # q index = 3 e + t


def prof(u):
    u = np.asarray(u, float)
    w = np.exp(2j * np.pi / 3)
    def ch(eps, a):
        return sum(u[3 * e + t] * (eps ** e) * (w ** (a * t)) for e, t in QE)
    return (round(abs(ch(-1, 0)) ** 2), round(abs(ch(1, 1)) ** 2), round(abs(ch(-1, 1)) ** 2))


def tr(u, s):
    """translate u by s in Q: u'(q) = u(q + s)"""
    se, st = divmod(s, 3)
    return tuple(u[3 * ((e + se) % 2) + (t + st) % 3] for e, t in QE)


def negt(u):
    return tuple(u[3 * e + (-t) % 3] for e, t in QE)


def canon_u(u):
    """canonical representative of u under translation, negation (max - u) and t -> -t"""
    best = None
    for v in (tuple(u), tuple(max(u) - x for x in u)):
        for w in (v, negt(v)):
            for s in range(6):
                x = tr(w, s); m = min(x); x = tuple(a - m for a in x)
                if best is None or x < best: best = x
    return best


def block_shapes(maxval=3):
    """all u (min 0) with p1 + 2 p2 + 2 p3 <= 25, grouped by canonical class"""
    shapes = {}
    for u in itertools.product(range(maxval + 1), repeat=6):
        if min(u) != 0: continue
        p = prof(u)
        if p[0] + 2 * p[1] + 2 * p[2] > 25: continue
        if max(p) > 5: continue
        shapes.setdefault(canon_u(u), p)
    return shapes


def qlevel_classes():
    """multisets of 4 canonical shapes whose profiles sum to (5, 5, 5)"""
    S = block_shapes(); keys = sorted(S)
    out = []
    for comb in itertools.combinations_with_replacement(range(len(keys)), 4):
        tot = tuple(sum(S[keys[i]][j] for i in comb) for j in range(3))
        if tot == (5, 5, 5): out.append([keys[i] for i in comb])
    return S, out


def realisations(cls_shapes):
    """for a class (4 canonical shapes): all translates/negations of each block that keep the profile sum exact
    (profiles are translation and negation invariant, so every combination works); return one representative"""
    return [list(u) for u in cls_shapes]


def rows_for(us, v):
    """odd c_k with sum_k (6 c_k + 2 sum u_k)^2 = 4 v; up to the per-block negation (u -> max u - u) the sign of the
    row is free, so list all"""
    U = [sum(u) for u in us]
    lim = int(math.isqrt(4 * v) // 6) + 3
    sols = []
    for c in itertools.product(range(-lim - (1 - lim % 2), lim + 2, 2), repeat=4):
        if all(x % 2 for x in c) and sum((6 * c[k] + 2 * U[k]) ** 2 for k in range(4)) == 4 * v:
            sols.append(c)
    return sols


if __name__ == "__main__":
    a = sys.argv[1:]
    v = int(a[a.index('--v') + 1]) if '--v' in a else 599
    S, C = qlevel_classes()
    print(json.dumps(dict(shapes=len(S), classes=len(C))))
    for cl in C:
        rs = rows_for(cl, v)
        print(json.dumps(dict(shapes=cl, profiles=[S[u] for u in cl], sizes=[sum(u) for u in cl], row_solutions=len(rs))))
