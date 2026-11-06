r"""A closed-form energy obstruction at the Z_q compression, for b1 slices [P + M].

If a b1 core exists on G = G' x Z_q, its Z_q compressions w_k satisfy (lag 0 of sum_k PAF_{Z_q}(w_k), using
sum_t sum_k PAF(w_k)(t) = sum_k r_k^2 = sum_h T1(h)):
    sum_k E(w_k) = sum_h T1(h) + (q - 1) * 4 |G'|,        E(w) = w(0)^2 + f sum_j W_j^2,
where, for a slice block with U-orbit values y(O) = alpha(O) + f sigma(O):
    w(0) = sum_O |O| alpha(O)  (fixed by y),   W_j = w(C_j),   sum_j W_j = sum_O |O| sigma(O) =: S  (fixed by y),
    W_j = sum_O |O| c_j(O) = sum_O |O|  (mod 2) for every j.
Relaxation [P]: E(w) >= LB := w(0)^2 + f * min{ sum_j W_j^2 : W in Z^s, sum W_j = S, W_j = par (mod 2) }, the min
attained by the most balanced split (computed exactly below by moving 2 between coordinates).  So a stage-1 class
with sum_k LB_k > the right-hand side has no lift.  This uses only y, the orbit sizes and (q, s): no enumeration.
usage: python3 theory/h4892_energybound.py 'KEY' LIST [--out FILE]   (LIST as written by h4892_stage1_any --full)"""
import sys, json, itertools


def minsq(S, s, par):
    # s integers of parity par with sum S, minimising the sum of squares: as equal as possible in steps of 2
    base = [par] * s; rem = S - par * s
    assert rem % 2 == 0, (S, s, par)
    k, r = divmod(rem // 2, s)                     # add 2k to all, and 2 more to r of them (floor division handles sign)
    W = [b + 2 * k + (2 if i < r else 0) for i, b in enumerate(base)]
    assert sum(W) == S
    return sum(x * x for x in W)


def main():
    a = sys.argv[1:]; key = json.loads(a[0]); L = json.load(open(a[1]))
    fam, v, Q, H1, q, s, U = key; assert fam == 'b1'
    f = (q - 1) // s; sizes = [len(o) for o in L['orbits']]; nG = sum(sizes)
    par = sum(sizes) % 2
    rhs = 4 + (q - 1) * 4 * nG                     # b1: sum_h T1(h) = sum_k r_k^2 = 4
    def lb(y):
        # every valid decoding y = alpha + f sigma per orbit (unique for f >= 2; for f = 1 both alpha are valid),
        # minimised over all combinations
        opts = [[(x, (val - x) // f) for x in (1, -1) if (val - x) % f == 0 and abs((val - x) // f) <= s
                 and ((val - x) // f + s) % 2 == 0] for val in y]
        best = None
        for comb in itertools.product(*opts):
            w0 = sum(m * al for (al, sg), m in zip(comb, sizes)); S = sum(m * sg for (al, sg), m in zip(comb, sizes))
            b = w0 * w0 + f * minsq(S, s, par); best = b if best is None else min(best, b)
        return best
    tot = [sum(lb(y) for y in Y) for Y in L['solutions']]
    out = dict(key=key, script='theory/h4892_energybound.py', rhs=rhs, classes=len(tot),
               classes_killed=sum(1 for t in tot if t > rhs), min_lower_bound=min(tot), max_lower_bound=max(tot))
    print(json.dumps(out))
    if outp := (a[a.index('--out') + 1] if '--out' in a else None): json.dump(dict(out, per_class=tot), open(outp, 'w'))


if __name__ == "__main__":
    main()
