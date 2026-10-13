r"""Stage split of a slice on G = Z_F x F_q (q = p^2): compressed stage 1 (along F_q) and its Gaussian [E, counts exact].

For a multiplier group M = <(u_i, gamma^{e_i})> (theory/nd_fq_sweep.py conventions) the F_q-compression
y(g) = sum_z x(g, z) is M-invariant on Z_F; on the M-orbit of g it ranges over the signed sums of the sizes of the
M_g-orbits partitioning the fibre {g} x F_q (exact), and
    sum_k PAF_{Z_F}(y_k)(0) = 4m - 4(q - 1),  sum_k PAF_{Z_F}(y_k)(h) = -4q (h != 0),  rows (2, 0, 0, 0).
Stage-1 estimate: exact count of quadruples with the lag-0 energy and the rows (DP over the Z_F-orbits), times the
lattice Gaussian of the remaining lags (one per <M_F, -1>-orbit) under uniform sampling from that set.
Stage-2 per y = whole slice (theory/nd_fq_sweep.py) - stage 1 (an average; per-y values need real y's).
Also writes the compressed problem (orbits, letters, targets) for the exact solvers.
usage: python3 theory/nd_fq_split.py F p "u,e[;u,e]" [--N 20000] [--dump FILE]"""
import sys, os, json, math
import numpy as np
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits
from nd_fq_sweep import Group, closure, orbits_of

def compressed(F, p, gens):
    Gr = Group(F, p); n = Gr.n; q = Gr.q
    M = closure([Gr.act(u, e) for (u, e) in gens], n)
    lab, orbs = orbits_of(M, n)
    # Z_F orbits of the projected action
    MF = sorted({int(pm[g * q] // q) for pm in M for g in [1 % F]})   # image of 1 under M on Z_F
    UF = sorted({u % F for (u, e) in gens})
    Fo = -np.ones(F, np.int64); freps = []; fsz = []
    for g in range(F):
        if Fo[g] >= 0: continue
        o = sorted({int(pm[g * q]) // q for pm in M}); Fo[o] = len(freps); freps.append(g); fsz.append(len(o))
    letters = []
    for g in freps:
        parts = defaultdict(int)
        for z in range(q): parts[int(lab[g * q + z])] += 1
        L = {0}
        for s_ in parts.values(): L = {a + e * s_ for a in L for e in (1, -1)}
        letters.append(sorted(L))
    # lag orbits on Z_F under <M_F, -1>
    Lo = -np.ones(F, np.int64); lreps = []
    for h in range(F):
        if Lo[h] >= 0: continue
        o = set()
        for pm in M:
            x = int(pm[h * q]) // q; o.add(x); o.add((-x) % F)
        Lo[sorted(o)] = len(lreps); lreps.append(h)
    T = [4 * Gr.n - 4 * (q - 1)] + [-4 * q] * (len(lreps) - 1)
    return dict(F=F, p=p, q=q, m=Gr.n, gens=gens, Fo=Fo.tolist(), freps=freps, fsz=fsz, letters=letters, lreps=lreps, T=T)

def s1_estimate(C, N=20000, seed=1):
    rng = np.random.default_rng(seed); F = C['F']; T = np.array(C['T'], np.int64)
    L = C['letters']; w = C['fsz']; allL = sorted({x for l in L for x in l})
    g = 0
    for y in allL: g = math.gcd(g, y * y - 1)
    unit = lambda y: (y * y - 1) // g
    E = (int(T[0]) - 4 * F) // g
    nr = len(w)
    tabs = [None] * (nr + 1); tabs[nr] = {(0, 0): 1}
    for i in range(nr - 1, -1, -1):
        nd = defaultdict(int)
        for (e, r), c in tabs[i + 1].items():
            for y in L[i]:
                e2 = e + w[i] * unit(y)
                if e2 <= E: nd[(e2, r + w[i] * y)] += c
        tabs[i] = dict(nd)
    c0 = defaultdict(int); c2 = defaultdict(int)
    for (e, r), c in tabs[0].items():
        if r == 0: c0[e] += c
        if r == 2: c2[e] += c
    parts = []; Q = 0
    for e0 in range(E + 1):
        for e1 in range(E + 1 - e0):
            for e2 in range(E + 1 - e0 - e1):
                e3 = E - e0 - e1 - e2
                ww = c2.get(e0, 0) * c0.get(e1, 0) * c0.get(e2, 0) * c0.get(e3, 0)
                if ww: parts.append((e0, e1, e2, e3, ww)); Q += ww
    out = dict(energy_units=E, unit=g, quadruples_log2=round(math.log2(Q), 2) if Q else None, energy_splits=len(parts))
    if not Q: out['stage1_log2E'] = -math.inf; out['note'] = 'no quadruple with lag-0 energy and rows'; return out
    W = np.array([p_[4] for p_ in parts], float); W /= W.sum()
    Fo = np.array(C['Fo']); lr = C['lreps'][1:]
    def sample(e, r):
        vals = []
        for i in range(nr):
            opts, wts = [], []
            for y in L[i]:
                c = tabs[i + 1].get((e - w[i] * unit(y), r - w[i] * y), 0)
                if c: opts.append(y); wts.append(c)
            wts = np.array(wts, float); y = opts[rng.choice(len(opts), p=wts / wts.sum())]
            vals.append(y); e -= w[i] * unit(y); r -= w[i] * y
        return np.array(vals)[Fo]
    V = np.zeros((N, len(lr)), np.int64)
    for t in range(N):
        e0, e1, e2, e3, _ = parts[rng.choice(len(parts), p=W)]
        acc = np.zeros(len(lr), np.int64)
        for e, r in ((e0, 2), (e1, 0), (e2, 0), (e3, 0)):
            y = sample(e, r)
            acc += np.array([int((y * np.roll(y, -h)).sum()) for h in lr])
        V[t] = acc
    gb = gauss_bits(V, T[1:]) if len(lr) else dict(bits=0.0, rank=0, in_lattice=True)
    out.update(lag_bits=gb, stage1_log2E=round(math.log2(Q) + 3 - gb['bits'], 2) if gb.get('in_lattice') else -math.inf)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "20000", "--dump": ""}; pos = []; i = 0
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(a[i]); i += 1
    F, p = int(pos[0]), int(pos[1]); gens = [tuple(map(int, t.split(','))) for t in pos[2].split(';')]
    C = compressed(F, p, gens)
    r = s1_estimate(C, int(opt["--N"]))
    print(json.dumps(dict(F=F, p=p, gens=gens, Z_F_orbits=len(C['freps']), orbit_sizes=C['fsz'],
                          letter_sets=[l if len(l) <= 12 else (l[:3] + ['...'] + l[-3:]) for l in C['letters']],
                          lag_classes=len(C['lreps']) - 1, **r)), flush=True)
    if opt["--dump"]: json.dump(C, open(opt["--dump"], "w"))
