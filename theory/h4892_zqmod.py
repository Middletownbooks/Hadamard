r"""Necessary check of stage 2 at the Z_q compression, reduced modulo m [M] -- for slices whose U-orbits have sizes 1
and multiples of m.  Shares no code with the other h4892 files or the stage-1 tool.

With w_k the Z_q compression (theory/h4892_zqcheck.py docstring), write w = s1 + m s2, where s1 collects the orbits
whose size is not divisible by m (here: the singletons) and s2 the rest.  Then PAF(w)(t) = PAF(s1)(t) (mod m) [P],
so a core needs sum_k PAF(s1_k)(t) = -4|G'| (mod m) for every t != 0.  Per block, the reachable s1 (exactly, from the
sign s-tuples of the small orbits with sum sigma(O), and alpha) give a set of PAF vectors mod m on one lag per coset
of M; the four blocks are joined 2+2 on these vectors mod m.  A class with no joined quadruple has no lift.
Input list as for h4892_zqcheck.py.  Planted: random reachable s1 of the four blocks of a random class, target their
PAF sum mod m.
Mode --all: w itself is reduced mod m over ALL orbits (PAF(w) mod m depends only on w mod m [P]); the reachable
set of w mod m is built exactly, orbit by orbit (at most m^s states).
usage: python3 theory/h4892_zqmod.py 'KEY' LIST M [--all] [--planted N] [--out FILE]"""
import sys, json, itertools, random
from collections import defaultdict


def main():
    a = sys.argv[1:]; key = json.loads(a[0]); L = json.load(open(a[1])); m = int(a[2])
    npl = int(a[a.index('--planted') + 1]) if '--planted' in a else 4
    outp = a[a.index('--out') + 1] if '--out' in a else None
    fam, v, Q, H1, q, s, U = key; assert fam == 'b1'
    f = (q - 1) // s
    assert f >= 2, 'alpha is not determined by y when f = 1; not supported'
    primes = [p_ for p_ in range(2, q) if (q - 1) % p_ == 0 and all(p_ % d for d in range(2, p_))]
    g = next(x for x in range(2, q) if all(pow(x, (q - 1) // p_, q) != 1 for p_ in primes))
    Mg = sorted({pow(g, s * e, q) for e in range(f)})
    cos = [sorted({pow(g, j, q) * u % q for u in Mg}) for j in range(s)]
    pos = {z: j for j, c in enumerate(cos) for z in c}
    sizes = [len(o) for o in L['orbits']]; nG = sum(sizes)
    small = list(range(len(sizes))) if '--all' in a else [i for i, z in enumerate(sizes) if z % m]
    pats = {sg: [t for t in itertools.product((1, -1), repeat=s) if sum(t) == sg] for sg in range(-s, s + 1, 2)}
    lagrep = [c[0] for c in cos]
    def decode(y):
        for al in (1, -1):
            if (y - al) % f == 0 and ((y - al) // f + s) % 2 == 0 and abs((y - al) // f) <= s: return al, (y - al) // f
        raise ValueError(y)
    cache = {}
    def pafset(y):
        k = tuple(y)
        if k in cache: return cache[k]
        acc = {tuple([0] * s)}; w0 = 0
        for i in small:
            al, sg = decode(y[i]); w0 += sizes[i] * al
            acc = {tuple((x + sizes[i] * c) % m for x, c in zip(A, t)) for A in acc for t in pats[sg]}
        for i in range(len(sizes)):
            if i not in small: w0 += sizes[i] * decode(y[i])[0]
        P = set()
        for w in acc:
            sq = [w0 % m] + [w[pos[z]] for z in range(1, q)]
            P.add(tuple(sum(sq[i] * sq[(i + t) % q] for i in range(q)) % m for t in lagrep))
        cache[k] = P; return P
    def join(Ps, T):
        left = {tuple((x + y_) % m for x, y_ in zip(A, B)) for A in Ps[0] for B in Ps[1]}
        for C in Ps[2]:
            for D in Ps[3]:
                if tuple((t - x - y_) % m for t, x, y_ in zip(T, C, D)) in left: return True
        return False
    T = tuple([(-4 * nG) % m] * s)
    sols = L['solutions']; rnd = random.Random(7 * q + m); pl = []
    for t in range(npl):
        Y = sols[rnd.randrange(len(sols))]; Ps = [pafset(y) for y in Y]
        pick = [rnd.choice(sorted(P)) for P in Ps]
        pl.append(join(Ps, tuple(sum(p[i] for p in pick) % m for i in range(s))))
    res = [join([pafset(y) for y in Y], T) for Y in sols]
    out = dict(key=key, script='theory/h4892_zqmod.py', modulus=m, mode='all' if '--all' in a else 'small', small_orbits=small, classes=len(sols),
               classes_passing=sum(res), planted=pl, passing=[i for i, r in enumerate(res) if r][:200],
               pafset_sizes=sorted({len(P) for P in cache.values()}))
    print(json.dumps({k: v for k, v in out.items() if k != 'passing'}))
    if outp: json.dump(out, open(outp, 'w'))


if __name__ == "__main__":
    main()
