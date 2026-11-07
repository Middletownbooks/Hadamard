r"""Reference stage-2 counter for the plug-in slice (gate A2 control; shares no code with theory/h2396_stage2.c).
Given a stage-1 solution u (4 x n), decode s and the clean sets, enumerate every t_i in {+-1}^clean_i, and count the
quadruples with sum_i PAF_{t_i}(h) = -(1/q) sum_i PAF_{s_i}(h) (h != 0) and sum_i [CCF(s_i,t_i) + eps CCF(t_i,s_i)](h) = 0
by a dictionary 2+2 join on the full key (all lags).  Also builds planted instances for the C counter.
usage: python3 theory/h2396_stage2_ref.py count n q 'u0;u1;u2;u3'        (u_i comma-separated)
       python3 theory/h2396_stage2_ref.py plant n q E seed              (prints a planted instance for h2396_stage2)"""
import sys, itertools, random
from collections import Counter


def decode(u, q):
    s, cl = [], []
    for row in u:
        sr, cr = [], []
        for v in row:
            if abs(v) == 1: sr.append(v); cr.append(1)
            elif abs(v) == q: sr.append(0); cr.append(0)
            elif abs(v) == q - 2: sr.append(-2 * (v // (q - 2))); cr.append(0)
            else: raise ValueError(v)
        s.append(sr); cl.append(cr)
    return s, cl


def key(t, s, n, eps):
    P = tuple(sum(t[x] * t[(x + h) % n] for x in range(n)) for h in range(1, n))
    Q = tuple(sum(s[x] * t[(x + h) % n] + eps * t[x] * s[(x + h) % n] for x in range(n)) for h in range(n))
    return P + Q


def count(n, q, u, T2=None, QT=None):
    eps = -1 if ((q - 1) // 2) % 2 else 1
    s, cl = decode(u, q)
    if T2 is None:
        T2 = []
        for h in range(1, n):
            p = sum(s[i][x] * s[i][(x + h) % n] for i in range(4) for x in range(n))
            if p % q: return 0
            T2.append(-p // q)
        QT = [0] * n
    target = tuple(T2) + tuple(QT)
    lists = []
    for i in range(4):
        pos = [x for x in range(n) if cl[i][x]]
        C = Counter()
        for signs in itertools.product((1, -1), repeat=len(pos)):
            t = [0] * n
            for x, sg in zip(pos, signs): t[x] = sg
            C[key(t, s[i], n, eps)] += 1
        lists.append(C)
    AB = Counter()
    for ka, ca in lists[0].items():
        for kb, cb in lists[1].items():
            AB[tuple(a + b for a, b in zip(ka, kb))] += ca * cb
    tot = 0
    for kc, cc in lists[2].items():
        for kd, cd in lists[3].items():
            need = tuple(t - c - d for t, c, d in zip(target, kc, kd))
            tot += AB.get(need, 0) * cc * cd
    return tot


def plant(n, q, E, seed):
    """random s / clean sets with E defects (letters +-q or +-(q-2)), random t; targets = the planted sums"""
    rng = random.Random(seed); eps = -1 if ((q - 1) // 2) % 2 else 1
    u = [[rng.choice((1, -1)) for _ in range(n)] for _ in range(4)]
    cells = rng.sample([(i, x) for i in range(4) for x in range(n)], E)
    for (i, x) in cells: u[i][x] = rng.choice((q, -q, q - 2, -(q - 2)))
    s, cl = decode(u, q)
    t = [[rng.choice((1, -1)) if cl[i][x] else 0 for x in range(n)] for i in range(4)]
    T2 = [sum(t[i][x] * t[i][(x + h) % n] for i in range(4) for x in range(n)) for h in range(1, n)]
    QT = [sum(s[i][x] * t[i][(x + h) % n] + eps * t[i][x] * s[i][(x + h) % n] for i in range(4) for x in range(n)) for h in range(n)]
    return u, t, T2, QT


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == 'count':
        n, q = int(sys.argv[2]), int(sys.argv[3])
        u = [list(map(int, r.split(','))) for r in sys.argv[4].split(';')]
        print(count(n, q, u))
    elif cmd == 'plant':
        n, q, E, seed = map(int, sys.argv[2:6])
        u, t, T2, QT = plant(n, q, E, seed)
        for r in u: print(' '.join(map(str, r)))
        print('PLANT'); print(' '.join(map(str, T2))); print(' '.join(map(str, QT)))
        print('#T ' + ' | '.join(' '.join(map(str, r)) for r in t), file=sys.stderr)
    elif cmd == 'countplant':
        n, q, E, seed = map(int, sys.argv[2:6])
        u, t, T2, QT = plant(n, q, E, seed)
        print(count(n, q, u, T2, QT))
