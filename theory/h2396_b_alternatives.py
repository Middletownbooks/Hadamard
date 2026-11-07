r"""Part B, gate B1: which multiplier groups H <= GL(2, Z4) x Z37^* give 70 orbits on G = Z4^2 x Z37, and which
7-sets T = {1} u (6 order-4 characters, closed under conjugation) does each preserve [M].  Cyclic H = <(A, m)> and
two-generator H = <(A, m), (B, m')> with B, m' ranging over generators of the factors (A, B of order dividing 12).
usage: python3 theory/h2396_b_alternatives.py"""
import itertools, json
E = [(a, b) for a in range(4) for b in range(4)]
def matmul(A, B): return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) % 4 for j in range(2)) for i in range(2))
def act(A, x): return tuple(sum(A[i][k] * x[k] for k in range(2)) % 4 for i in range(2))
GL = [((a, b), (c, d)) for a, b, c, d in itertools.product(range(4), repeat=4) if ((a * d - b * c) % 2) == 1]
I2 = ((1, 0), (0, 1))
def order(A):
    o, X = 1, A
    while X != I2: X = matmul(X, A); o += 1
    return o
p = 37
def ordm(m):
    o, x = 1, m
    while x != 1: x = x * m % p; o += 1
    return o
units37 = list(range(1, p))
def orbits(gens):
    seen = set(); n = 0
    for x in E:
        for k in range(p):
            if (x, k) in seen: continue
            n += 1; st = [(x, k)]; seen.add((x, k))
            while st:
                y, l = st.pop()
                for (A, m) in gens:
                    z = (act(A, y), l * m % p)
                    if z not in seen: seen.add(z); st.append(z)
    return n
# characters chi_a(x) = i^{a.x} (standard pairing on Z4^2); order-4 iff a not in 2Z4^2
ord4 = [a for a in E if a[0] % 2 or a[1] % 2]
def dual_act(A, a):            # chi_a o A = chi_{A^T a}
    At = ((A[0][0], A[1][0]), (A[0][1], A[1][1])); return act(At, a)
res = []
reps = {}
for A in GL:                    # conjugacy is not reduced; small enough
    oA = order(A)
    for m in units37:
        if 36 % ordm(m): continue
        n = orbits([(A, m)])
        if n == 70: res.append((A, m, oA, ordm(m)))
print("cyclic H with 70 orbits:", len(res), "orders (A, m):", sorted({(r[2], r[3]) for r in res}))
# T-sets preserved, for one representative per (order A, order m)
out = []
for (oA, om) in sorted({(r[2], r[3]) for r in res}):
    A, m = next((r[0], r[1]) for r in res if r[2] == oA and r[3] == om)
    # 6-sets of order-4 characters closed under negation and A-invariant
    pairs = sorted({tuple(sorted([a, ((-a[0]) % 4, (-a[1]) % 4)])) for a in ord4})
    good = []
    for ch in itertools.combinations(pairs, 3):
        S = {a for pr in ch for a in pr}
        if all(dual_act(A, a) in S for a in S): good.append(sorted(S))
    out.append(dict(A=A, m=m, order_A=oA, order_m=om, invariant_T6=len(good), example=good[:2]))
print(json.dumps(out, indent=1))
