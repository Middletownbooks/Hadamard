r"""Part C, gate C2: the exact Z_q-level test of a stage-1 solution y of the index-4 Szolloesi plug-in on Z_nK x F_q
(as report/h4892_answer.md section 2.4: a lift exists only if the compression along K to F_q can be met) [M].
Given y (4 x nK; y = x0 + L sigma, x0 = 0 at k = 0 for A, B), every block's class sums W_j = sum_k x_j(k) range over the
sums of per-position class vectors (x_0..x_3)(k) in {+-1}^4 with sum sigma(k) (A, B: one vector per pair {k, -k},
counted twice).  The Z_q level is the five equations of theory/h2396_c_zq.py with a_i = sum_k x0_i(k) fixed by y.
Per block the reachable W sets are built exactly; then a 2+2 join on the five PAF-class values.
Input: JSON lines with "y" (4 lists), or a theory/h2396_c_s1cp.py output file (its "solution" field).
usage: python3 theory/h2396_c_zqtest.py nK q FILE"""
import sys, json, itertools
from collections import defaultdict
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from h2396_c_s1cp import classes


def zq_test(nK, q, y):
    S, cls, cyc = classes(q); L = (q - 1) // 4; T0 = 4 * nK * q - 2
    vecs = {s: [v for v in itertools.product((1, -1), repeat=4) if sum(v) == s] for s in (-4, -2, 0, 2, 4)}
    keys = []
    for i in range(4):
        x0 = []; sig = []
        for k in range(nK):
            v = y[i][k]
            a = 0 if (i < 2 and k == 0) else (1 if (v - 1) % (2 * L) == 0 else -1)
            s = (v - a) // L; assert a + L * s == v; x0.append(a); sig.append(s)
        a_i = sum(x0)
        reach = {(0, 0, 0, 0)}
        pos = [k for k in range(nK) if not (i < 2 and k > nK // 2)]
        for k in pos:
            mult = 2 if (i < 2 and k != 0) else 1
            new = set()
            for W in reach:
                for v in vecs[sig[k]]:
                    new.add(tuple(W[j] + mult * v[j] for j in range(4)))
            reach = new
        ks = set()
        for W in reach:
            P = [a_i * a_i + L * sum(w * w for w in W)]
            for l in range(4):
                P.append(2 * a_i * W[l] + sum(W[j] * W[jj] * cyc[(j - l) % 4][(jj - l) % 4] for j in range(4) for jj in range(4)))
            if P[0] <= T0: ks.add(tuple(P))
        keys.append(ks)
    AB = defaultdict(int)
    for a in keys[0]:
        for b in keys[1]: AB[tuple(x + z for x, z in zip(a, b))] += 1
    tgt = (T0, 0, 0, 0, 0); hits = 0
    for c in keys[2]:
        for d in keys[3]:
            hits += AB.get(tuple(t - x - z for t, x, z in zip(tgt, c, d)), 0)
    return hits, [len(k) for k in keys]


if __name__ == "__main__":
    nK, q = int(sys.argv[1]), int(sys.argv[2])
    S, cls, cyc = classes(q); L = (q - 1) // 4
    for line in open(sys.argv[3]):
        d = json.loads(line)
        if 'y' in d: y = d['y']
        elif d.get('solution'):
            sol = d['solution']; y = [[sol['x0'][i][k] + L * sum(sol['xj'][i][k]) for k in range(nK)] for i in range(4)]
        else: continue
        hits, sizes = zq_test(nK, q, y)
        print(json.dumps(dict(nK=nK, q=q, Zq_level_quadruples=hits, reachable_W_keys=sizes, status=d.get('status'))), flush=True)
