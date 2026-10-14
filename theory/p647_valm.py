r'''Validation of join4m on real data: one class-merged unit (logs/p647/prod/units.json, unit 0) on lists restricted to sigma with
multiplicity >= 4 and the two central energy classes; the merged join must equal the union of the separate join4p runs of
its (left, right) joint solutions. Output: logs/p647/validation_join4m_unit0.json.'''
import sys, json, time; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np, p647_common as PC, p647_cross as P, p647_join as J
from nd_cross import orbit_maps
X, ys = PC.load_setup(); rec, sols = PC.joint_solutions(X, ys, 0)
reps, oi, maps = orbit_maps(X); NH = P.lag_tensor(X, oi)
unit = json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))), 'logs', 'p647', 'prod', 'units.json')))['units'][0]
JJ = unit['joint']; nl, nr = len(JJ), len(JJ[0])
cache = {}
def blocks_of(e):
    bl = P.build(X, ys[0], X.J.ya_full(sols[e]['layers']), maps, NH, oi)
    return bl
B = {}
for i in range(nl):
    for j in range(nr):
        B[(i, j)] = blocks_of(JJ[i][j])
keep = lambda Bk: Bk.Ssig[Bk.mult >= 4]
lefts = [(keep(B[(i, 0)][0]), keep(B[(i, 0)][1])) for i in range(nl)]
rights = [(keep(B[(0, j)][2]), keep(B[(0, j)][3])) for j in range(nr)]
# consistency: left lists of (i, j) do not depend on j, right lists not on i
for i in range(nl):
    for j in range(nr):
        assert np.array_equal(keep(B[(i, j)][0]), lefts[i][0]) and np.array_equal(keep(B[(i, j)][1]), lefts[i][1])
        assert np.array_equal(keep(B[(i, j)][2]), rights[j][0]) and np.array_equal(keep(B[(i, j)][3]), rights[j][1])
els = [319, 343]
t = time.time(); M, st = J.join4m(lefts, rights, X.TS, els=els); tm = time.time() - t
sep = set(); tsep = 0
for i in range(nl):
    for j in range(nr):
        t = time.time(); s_, _ = J.join4p([lefts[i][0], lefts[i][1], rights[j][0], rights[j][1]], X.TS, order=[0, 1, 2, 3], els=els); tsep += time.time() - t
        sep |= {(i, j) + tuple(x) for x in s_}
print(json.dumps(dict(unit_shape=[nl, nr], merged=len(M), separate=len(sep), identical=set(M) == sep, merged_pairs=st['left'] + st['right'],
                      merged_secs=round(tm, 1), separate_secs=round(tsep, 1), per_pair=[sorted({(a, b) for a, b, *_ in M})])))
