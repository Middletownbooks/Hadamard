"""H(4892) step 1: compare the two complete stage-1 enumerations as sets of explicit quadruples [M].
(a) logs/h4892/stage1_toolB_full.json (method B of open_screen_stage1_exact.py, up to block order and sign);
(b) data/h4892_stage1.json (theory/h4892_stage1.py, every ordered quadruple).
Checks: canonical sets equal; the ordered set of (b) equals the full expansion of the canonical set by the 24 block
orders and 16 block signs; each of the 384 classes has 4 distinct blocks (so 384 * 384 ordered quadruples)."""
import json, itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h4892_stage1 as S
A = json.load(open('data/h4892_stage1.json')); B = json.load(open('logs/h4892/stage1_toolB_full.json'))
can = lambda Y: tuple(sorted(tuple(y) if y[0] > 0 else tuple(-x for x in y) for y in Y))
sa = {can(Y) for Y in A['canonical']}; sb = {can(Y) for Y in B['list']}
tab = [tuple(S.expand(v)) for v in A['block_table']]
ordset = {tuple(tab[i] for i in Q) for Q in A['ordered']}
exp = {tuple(tuple(sg[k] * x for x in Y[p[k]]) for k in range(4))
       for Y in sa for p in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)}
out = dict(canonical_A=len(sa), canonical_B=len(sb), canonical_equal=sa == sb, ordered=len(ordset),
           expansion=len(exp), ordered_equals_expansion=ordset == exp,
           all_blocks_distinct=all(len(set(Y)) == 4 for Y in sa),
           distinct_blocks_up_to_sign=len({y for Y in sa for y in Y}),
           rows=sorted({tuple(sorted(abs(sum(y)) for y in Y)) for Y in sa}))
print(json.dumps(out)); json.dump(out, open('logs/h4892/stage1_compare.json', 'w'))
