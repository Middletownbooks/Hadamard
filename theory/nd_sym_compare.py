r"""Cross-check theory/nd_sym_stage1.py (alpha-congruence + S/D split) against theory/nd_sym_bruteforce.py [M].

nd_sym_stage1 returns stage-1 solutions with canonical alpha's (blocks 1..3 sorted); every solution is mapped here to
the brute-force convention (block 0 the lexicographically least of its shift-by-F/2 pair; blocks 1..3 least under
shift and negation; blocks 1..3 in every order) and the distinct quadruples are counted.  The two counts must agree.
usage: python3 theory/nd_sym_compare.py F q s"""
import sys, os, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nd_sym_stage1 as S

def canon(y, neg):
    F = len(y); n = F // 2; vs = [tuple(y), tuple(np.roll(y, -n))]
    if neg: vs += [tuple(-y), tuple(-np.roll(y, -n))]
    return min(vs)

def main(F, q, s):
    out, sols = S.solve(F, q, s, log=lambda x: None)
    keys = set()
    for x in sols:
        Y = [np.array(y) for y in x['y']]
        c0 = canon(Y[0], False); rest = [canon(y, True) for y in Y[1:]]
        for p in itertools.permutations(rest): keys.add((c0,) + p)
    res = dict(F=F, q=q, s=s, solver_solutions=len(sols), all_verified=out['all_verified'],
               canonical_ordered=len(keys), alpha_quads=out['alpha_quads'], S_solutions=out['S_solutions'])
    print(json.dumps(res), flush=True)

if __name__ == "__main__":
    main(*map(int, sys.argv[1:4]))
