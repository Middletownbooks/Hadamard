r"""Gate A2 diagnostic: the stage-2 quotient ladder for one stage-1 solution y of the plug-in slice on Z_n x F_13 [M].
For a divisor m of n (quotient Z_m = K/L, |L| = n/m), the compression t -> T(j) = sum_{k = j mod m} t(k) is a ring map,
so a completion t gives T with
    sum_i PAF_{T_i}(d) = T2_m(d) := sum_{h = d mod m} T2(h)     (T2(0) = #clean, T2(h != 0) = -(1/13) sum PAF_s(h)),
    sum_i [CCF(S_i, T_i) + CCF(T_i, S_i)](d) = 0,             S_i(j) = sum_{k = j mod m} s_i(k),
with T_i(j) in {-c, -c + 2, ..., c}, c = #clean positions of block i in the coset j.  Each level is decided exactly by
a dictionary 2+2 join over the per-block value vectors (m small).  m = 1 is the Z_13 level (the trivial character of
K, report/h4892_answer.md section 5); m = n/2 is the S-level of theory/h2396_stage2.c.  Prints, per m, the number of
compressed solutions (as distinct value quadruples).
usage: python3 theory/h2396_s2ladder.py n FILE [m1,m2,...]"""
import sys, itertools
from collections import defaultdict


def decode(u, q=13):
    s, cl = [], []
    for row in u:
        s.append([v if abs(v) == 1 else (0 if abs(v) == q else -2 * (v // (q - 2))) for v in row])
        cl.append([1 if abs(v) == 1 else 0 for v in row])
    return s, cl


def level(n, u, m, q=13):
    s, cl = decode(u, q)
    T2 = [sum(cl[i][x] for i in range(4) for x in range(n))]
    for h in range(1, n):
        p = sum(s[i][x] * s[i][(x + h) % n] for i in range(4) for x in range(n))
        assert p % q == 0; T2.append(-p // q)
    Tm = [sum(T2[h] for h in range(n) if h % m == d) for d in range(m)]
    hl = m // 2
    target = tuple(Tm[:hl + 1]) + (0,) * (hl + 1)
    lists = []
    for i in range(4):
        S = [sum(s[i][k] for k in range(n) if k % m == j) for j in range(m)]
        c = [sum(cl[i][k] for k in range(n) if k % m == j) for j in range(m)]
        keys = set()
        for T in itertools.product(*[range(-cj, cj + 1, 2) for cj in c]):
            P = tuple(sum(T[j] * T[(j + d) % m] for j in range(m)) for d in range(hl + 1))
            Q = tuple(sum(S[j] * T[(j + d) % m] + T[j] * S[(j + d) % m] for j in range(m)) for d in range(hl + 1))
            keys.add(P + Q)
        lists.append(keys)
    AB = defaultdict(int)
    for a in lists[0]:
        for b in lists[1]: AB[tuple(x + y for x, y in zip(a, b))] += 1
    tot = 0
    for c_ in lists[2]:
        for d_ in lists[3]:
            tot += AB.get(tuple(t - x - y for t, x, y in zip(target, c_, d_)), 0)
    return tot, [len(L) for L in lists]


if __name__ == "__main__":
    n = int(sys.argv[1]); u = [list(map(int, l.split())) for l in open(sys.argv[2]) if l.strip()]
    ms = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [m for m in range(1, n) if n % m == 0 and m <= 6]
    for m in ms:
        tot, sizes = level(n, u, m)
        print(f"n={n} quotient Z_{m}: compressed stage-2 solutions {tot}  (distinct per-block keys {sizes})", flush=True)
