"""Exact Z_5-quotient test for the 955 b = 0 'BC' slice u = (-1 on Z_5, x5 of order 19 on Z_191) [P + M].

Compress every block along Z_191 onto Z_5: y(x) = sum_b X(x, b).  A, D are u-invariant, so y_A, y_D are symmetric
(y(x) = y(-x)); C(x, b) = B(-x, mu^-1 b) gives y_C(x) = y_B(-x), so PAF_{y_C} = PAF_{y_B}.  Every compressed value is
(value at b = 0) + 19 * (sum of 10 orbit signs), i.e. lies in +-1 + 38Z.  The quotient condition is
    PAF_{y_A}(j) + PAF_{y_D}(j) + 2 PAF_{y_B}(j) = 3820 [j = 0]   (j in Z_5),
with row sums (s_A, s_B, s_D) in {(+-37, +-35, +-1), (+-1, +-35, +-37)} (the slice's admissible patterns).
Enumerates every admissible value tuple; prints the solutions (if any)."""
import itertools, json
L = [v for v in range(-191, 192) if (v - 1) % 38 == 0 or (v + 1) % 38 == 0]
L = [v for v in L if 2 * v * v <= 3820]                      # every value appears at least twice in the energy, except
L0 = [v for v in range(-191, 192) if ((v - 1) % 38 == 0 or (v + 1) % 38 == 0) and v * v <= 3820]
def paf5(y):
    return [sum(y[x] * y[(x + j) % 5] for x in range(5)) for j in range(5)]
sym = []                                                       # symmetric vectors (a0, a1, a2, a2, a1)
for a0 in L0:
    for a1 in L:
        for a2 in L:
            y = (a0, a1, a2, a2, a1); sym.append((y, paf5(y), sum(y)))
gen = []
for y in itertools.product(L, repeat=5):
    if 2 * sum(t * t for t in y) > 3820: continue
    gen.append((y, paf5(y), sum(y)))
print('values', L0, 'symmetric vectors', len(sym), 'B vectors', len(gen))
pats = [(sa, sb, sd) for sa in (37, -37, 1, -1) for sb in (35, -35) for sd in (37, -37, 1, -1) if sa * sa + 2 * sb * sb + sd * sd == 3820]
B = {}
for y, p, s in gen: B.setdefault((s, tuple(p)), []).append(y)
sols = []
for ya, pa, sa in sym:
    for yd, pd, sd in sym:
        for sb in (35, -35):
            if (sa, sb, sd) not in pats: continue
            need = []
            for j in range(5):
                t = (3820 if j == 0 else 0) - pa[j] - pd[j]
                if t % 2: break
                need.append(t // 2)
            else:
                for yb in B.get((sb, tuple(need)), []): sols.append((ya, yd, yb))
print(json.dumps(dict(patterns=pats, solutions=len(sols), examples=sols[:5])))
