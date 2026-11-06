r"""Independent necessary check of stage 2 at the Z_q compression, for any b1 slice key and a stage-1 list [M].

Generalises theory/h4892_z13check.py; shares no code with the stage-2 files or the stage-1 tool.  If a b1 core
x_0..x_3 on G = G' x Z_q exists, the compressions w_k(z) = sum_{g'} x_k(g', z) satisfy
sum_k PAF_{Z_q}(w_k)(t) = -4 |G'| for t != 0 [P].  For a slice block (S = U x M, M of index s in Z_q^*, f = (q-1)/s)
with U-orbit values y(O) = alpha + f sigma: w(0) = sum_O |O| alpha(O), w(z) = sum_O |O| c_j(O) for z in the j-th
coset C_j = g^j M (g a generator of Z_q^*), (c_0..c_{s-1})(O) a sign s-tuple with sum sigma(O).
Per block the exact reachable set of (w(C_0), .., w(C_{s-1})) is built (as a set, orbit by orbit), each element
expanded to a length-q sequence, its nonzero-lag PAF computed on one lag per coset (PAF of an M-invariant sequence is
M-invariant), and the four blocks joined 2+2 on that s-vector (dictionary with multiplicities).
Input list: JSON with 'orbits' (U-orbits of G' as index lists) and 'solutions' (canonical classes as orbit-value
vectors), as written by theory/h4892_stage1_any.py --full.
usage: python3 theory/h4892_zqcheck.py 'KEY' LIST [--planted N] [--out FILE] [--maxreach R]"""
import sys, json, itertools, random
from collections import defaultdict


def main():
    a = sys.argv[1:]; key = json.loads(a[0]); L = json.load(open(a[1]))
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    outp = a[a.index('--out') + 1] if '--out' in a else None
    maxr = int(float(a[a.index('--maxreach') + 1])) if '--maxreach' in a else 3 * 10 ** 6
    fam, v, Q, H1, q, s, U = key; assert fam == 'b1'
    f = (q - 1) // s
    assert f >= 2, 'alpha is not determined by y when f = 1; not supported'
    g = next(x for x in range(2, q) if all(pow(x, (q - 1) // p_, q) != 1 for p_ in range(2, q) if (q - 1) % p_ == 0 and all(p_ % d for d in range(2, p_))))
    Mg = sorted({pow(g, s * e, q) for e in range(f)})
    cos = [sorted({pow(g, j, q) * m % q for m in Mg}) for j in range(s)]
    assert sorted(sum(cos, [])) == list(range(1, q))
    pos = {z: j for j, c in enumerate(cos) for z in c}
    sizes = [len(o) for o in L['orbits']]; nG = sum(sizes)
    pats = {sg: [t for t in itertools.product((1, -1), repeat=s) if sum(t) == sg] for sg in range(-s, s + 1, 2)}
    lagrep = [c[0] for c in cos]
    def decode(y):
        for al in (1, -1):
            if (y - al) % f == 0 and ((y - al) // f + s) % 2 == 0 and abs((y - al) // f) <= s: return al, (y - al) // f
        raise ValueError(y)
    cache = {}
    def reach(y):
        k = tuple(y)
        if k in cache: return cache[k]
        acc = {tuple([0] * s)}; w0 = 0
        for val, m in zip(y, sizes):
            al, sg = decode(val); w0 += m * al
            acc = {tuple(x + m * c for x, c in zip(A, t)) for A in acc for t in pats[sg]}
            if len(acc) > maxr: raise OverflowError(f'reachable set > {maxr}')
        P = defaultdict(int)
        for w in acc:
            sq = [w0] + [w[pos[z]] for z in range(1, q)]
            P[tuple(sum(sq[i] * sq[(i + t) % q] for i in range(q)) for t in lagrep)] += 1
        cache[k] = (w0, len(acc), P); return cache[k]
    def join(Ps, T):
        left = defaultdict(int)
        for x, nx in Ps[0].items():
            for y_, ny in Ps[1].items(): left[tuple(p + r for p, r in zip(x, y_))] += nx * ny
        c = 0
        for x, nx in Ps[2].items():
            for y_, ny in Ps[3].items():
                c += left.get(tuple(t - p - r for t, p, r in zip(T, x, y_)), 0) * nx * ny
        return c
    T = tuple([-4 * nG] * s)
    sols = L['solutions']; rnd = random.Random(q * 1000 + s); pl = []
    out = dict(key=key, script='theory/h4892_zqcheck.py', classes=len(sols))
    try:
        for t in range(npl):
            Y = sols[rnd.randrange(len(sols))]; R = [reach(y) for y in Y]
            pick = [rnd.choice(list(r[2])) for r in R]
            Tp = tuple(sum(p[i] for p in pick) for i in range(s))
            pl.append(join([r[2] for r in R], Tp) > 0)
        res = []
        for Y in sols:
            res.append(join([reach(y)[2] for y in Y], T))
        out.update(planted=pl, classes_with_zq_match=sum(1 for r in res if r), matches=res,
                   reach_sizes=sorted({c[1] for c in cache.values()}), distinct_paf=sorted({len(c[2]) for c in cache.values()}))
    except OverflowError as e:
        out.update(undecided=repr(e), planted=pl)
    print(json.dumps({k: v for k, v in out.items() if k != 'matches'}))
    if outp: json.dump(out, open(outp, 'w'))


if __name__ == "__main__":
    main()
