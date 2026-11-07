r"""H(4892) slice: independent check of the stage-2 obstruction at the Z13 compression [M].

Shares no code with theory/h4892_stage2*.py, h4892_slevel.py, h4892_stage1.py (reads only data/h4892_stage1.json).
If x_0..x_3 on G = Z47 x Z2 x Z13 satisfy sum_k PAF_G(x_k) = 4888 delta - 4, then the compressions
w_k(z) = sum_{g' in Z47 x Z2} x_k(g', z) satisfy sum_k PAF_{Z13}(w_k)(t) = -4*94 = -376 for every t != 0 [P: PAF of
a compression is the sum of PAF over the cosets].  For a slice block with U-orbit values y (sizes 1, 1, 23, 23, 23,
23 on the six orbits of Z47 x Z2 -- recomputed here), w(0) = sum_O |O| alpha(O) and w(z) = sum_O |O| c_j(O) for
z in C_j = 2^j {1, 3, 9}, where (c_0..c_3)(O) runs over the sign 4-tuples with sum sigma(O), y = alpha + 3 sigma.
For each block, the set of reachable (w(C_0), .., w(C_3)) is built exactly (product over the orbits, as a set);
each is expanded to a length-13 integer sequence and its PAF computed by the definition (all 13 lags); the four
blocks are joined 2+2 on the full 12-vector of nonzero-lag PAF values (ordered, dictionary).  A stage-1 class with
no joined quadruple has no stage-2 lift [M].  Planted: a random reachable quadruple's own PAF sum as the target.
usage: python3 theory/h4892_z13check.py [--planted N]"""
import sys, json, itertools, random
from collections import defaultdict

q, p = 13, 47
sq = {(a * a) % p for a in range(1, p)}
# orbits of Z47 x Z2 under the squares, by index 2a + b, as sets
orbs = []
for b in (0, 1):
    orbs.append({b}); orbs.append({2 * a + b for a in sq}); orbs.append({2 * a + b for a in range(1, p) if a not in sq})
cos = []
for j in range(4):
    g = pow(2, j, q); cos.append(sorted({g * m % q for m in (1, 3, 9)}))
pos_of = {}
for j, c in enumerate(cos):
    for z in c: pos_of[z] = j


def alpha_sigma(v):
    for al in (1, -1):
        if (v - al) % 6 == 0 or ((v - al) % 3 == 0 and ((v - al) // 3) % 2 == 0): return al, (v - al) // 3
    raise ValueError(v)


def reach(y):
    """y: the block as 94 values; returns (w0, set of reachable 4-tuples)"""
    acc = {(0, 0, 0, 0)}; w0 = 0
    for O in orbs:
        vals = {y[g] for g in O}; assert len(vals) == 1
        al, sg = alpha_sigma(vals.pop()); m = len(O); w0 += m * al
        pats = [t for t in itertools.product((1, -1), repeat=4) if sum(t) == sg]
        acc = {tuple(a + m * c for a, c in zip(x, t)) for x in acc for t in pats}
    return w0, acc


def seq(w0, w):
    return [w0] + [w[pos_of[z]] for z in range(1, q)]


def paf(s):
    return tuple(sum(s[i] * s[(i + t) % q] for i in range(q)) for t in range(1, q))


def join(lists, target):
    """lists: four lists of PAF 12-vectors (with multiplicity irrelevant); returns number of (distinct-vector) matches"""
    D = [defaultdict(int) for _ in range(4)]
    for k in range(4):
        for P in lists[k]: D[k][P] += 1
    left = defaultdict(int)
    for a, na in D[0].items():
        for b, nb in D[1].items(): left[tuple(x + y for x, y in zip(a, b))] += na * nb
    cnt = 0
    for c, nc in D[2].items():
        for d, nd in D[3].items():
            need = tuple(t - x - y for t, x, y in zip(target, c, d))
            cnt += left.get(need, 0) * nc * nd
    return cnt


def main():
    a = sys.argv[1:]; npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    D = json.load(open('data/h4892_stage1.json'))
    Ys = D['canonical']; T = tuple([-4 * 94] * 12)
    cache = {}
    def lst(y):
        k = tuple(y)
        if k not in cache:
            w0, R = reach(y); cache[k] = (w0, sorted(R), [paf(seq(w0, w)) for w in sorted(R)])
        return cache[k]
    rnd = random.Random(13)
    pl = []
    for t in range(npl):
        Y = Ys[rnd.randrange(len(Ys))]; L = [lst(y) for y in Y]
        pick = [rnd.randrange(len(l[2])) for l in L]
        Tp = tuple(sum(L[k][2][pick[k]][i] for k in range(4)) for i in range(12))
        pl.append(join([l[2] for l in L], Tp) > 0)
    print('planted recovered', pl, flush=True)
    res = []
    for i, Y in enumerate(Ys):
        L = [lst(y) for y in Y]
        # lag-0 consistency: sum_k sum_z w_k(z) = sum_k r_k and the energy of w is determined at lag 0 by the rest
        res.append(join([l[2] for l in L], T))
    print('classes', len(res), 'with a Z13-level quadruple', sum(1 for r in res if r), 'reachable sizes',
          sorted({len(l[1]) for l in cache.values()}))
    json.dump(dict(script='theory/h4892_z13check.py', classes=len(res), z13_matches=res, planted=pl,
                   reachable_sizes={str(k[:2]): len(v[1]) for k, v in cache.items()}),
              open('logs/h4892/z13check.json', 'w'))


if __name__ == "__main__":
    main()
