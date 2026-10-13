r"""Exact screen for generalised-Paley slices (b = 1): lag 0 + row sums + Z_2 compression of stage 1 [P, M].

Slice: G = Z_F x Z_q, Z_F fixed pointwise, index-s cyclotomic classes on Z_q, f = (q - 1)/s.  Stage-1 coset sums
y_k(g) take letters a + f * (sum of s signs).  Necessary conditions, all exact consequences of
sum_k PAF_{Z_F}(y_k) = T_J (T_J(0) = 4m - 4(q - 1), T_J(h != 0) = -4q):
  (L0) lag 0: sum_{k,g} y_k(g)^2 = 4m - 4(q - 1)          (fixes the multiset of |y|: a "composition")
  (R)  trivial character: sum_k r_k^2 = 4, r_k = sum_g y_k(g)  =>  rows (+-2, 0, 0, 0) up to order
  (A)  order-2 character: sum_k A_k^2 = 4v, A_k = sum_g (-1)^g y_k(g)
       (sum_h chi(h) T_J(h) = T_J(0) + 4q = 4(m + 1) for every nontrivial chi of Z_F)
For each composition the large letters are distributed over the 4 blocks and their two parity halves with signs; the
|y| = 1 entries of a half of size F/2 with j large letters contribute any sum in {-(F/2 - j), ..., F/2 - j} (step 2).
The slice is exactly empty if no composition admits (R) and (A) together.  Reports, per composition, whether (R) alone
and (R) + (A) are satisfiable.
usage: python3 theory/gp_z2screen.py F q s [F q s ...]"""
import sys, os, json, itertools as it
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_stage1 import compositions

def block_options(F, bigs):
    """set of (|r|, A^2) reachable by one block of length F holding the large letters `bigs` (magnitudes, any signs)"""
    n = F // 2; st = {(0, 0, 0, 0)}                              # (j0, B0, j1, B1) over the large letters
    for a in bigs:
        st = {t for (j0, B0, j1, B1) in st for t in ((j0 + 1, B0 + a, j1, B1), (j0 + 1, B0 - a, j1, B1),
                                                     (j0, B0, j1 + 1, B1 + a), (j0, B0, j1 + 1, B1 - a))
              if t[0] <= n and t[2] <= n}
    R = set(); RA = set()
    for (j0, B0, j1, B1) in st:
        u0, u1 = n - j0, n - j1
        for E in range(B0 - u0, B0 + u0 + 1, 2):
            for O in range(max(B1 - u1, -2 - E), min(B1 + u1, 2 - E) + 1):
                if (O - B1 + u1) % 2: continue
                R.add(abs(E + O)); RA.add((abs(E + O), (E - O) ** 2))
    return R, RA

def distributions(c):
    """multiset distributions of the large letters of composition c over 4 unlabelled blocks"""
    mags = [a for a in sorted(c) if a > 1 and c[a]]
    def splits(n):                                                # labelled 4-tuples summing to n
        return [t for t in it.product(range(n + 1), repeat=4) if sum(t) == n]
    out = set()
    for combo in it.product(*[splits(c[a]) for a in mags]):
        per = tuple(sorted(tuple(a for a, cnt in zip(mags, [combo[i][k] for i in range(len(mags))]) for _ in range(cnt))
                           for k in range(4)))
        out.add(per)
    return out

def screen(F, q, s):
    v = F * q + 1
    comps, lw = compositions(F, q, s)
    res = []
    for c in comps:
        big = sorted(a for a in c if a > 1 for _ in range(c[a]))
        okR = okRA = False; witness = None; cache = {}
        for per in distributions(c):
            opts = []
            for bl in per:
                if bl not in cache: cache[bl] = block_options(F, bl)
                opts.append(cache[bl])
            # (R): exactly one block with |r| = 2
            for two in range(4):
                if all((2 if k == two else 0) in opts[k][0] for k in range(4)): okR = True
                # (R) + (A): sum of A^2 over blocks = 4v
                sets = [{a2 for (r, a2) in opts[k][1] if r == (2 if k == two else 0)} for k in range(4)]
                acc = {0}
                for S in sets: acc = {x + y for x in acc for y in S if x + y <= 4 * v}
                if 4 * v in acc:
                    okRA = True; witness = dict(blocks=[list(b) for b in per], two=two)
            if okRA: break
        res.append(dict(comp={str(a): n for a, n in c.items() if n}, log2_weight=round(lw(c), 1), rows_ok=okR,
                        rows_and_Z2_ok=okRA, witness=witness))
    return dict(F=F, q=q, s=s, v=v, compositions=len(comps), survivors=sum(r["rows_and_Z2_ok"] for r in res),
                empty=not any(r["rows_and_Z2_ok"] for r in res), detail=res)

if __name__ == "__main__":
    a = list(map(int, sys.argv[1:]))
    for i in range(0, len(a), 3):
        print(json.dumps(screen(*a[i:i + 3])), flush=True)
