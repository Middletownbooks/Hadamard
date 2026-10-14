r"""Stage-1 Gaussian for a symmetric generalised-Paley slice (739 x215 and its controls) [E, counts exact].

Slice: G = Z_F x Z_q, M = <-1 on Z_F> x <index-s subgroup on Z_q> (for 739 x215: F = 18, q = 41, s = 8; the multiplier
215 is -1 mod 18 and has order 5 mod 41).  Stage 1 = the Z_q-coset sums y_k on Z_F: symmetric (y(-g) = y(g)), letters
alpha + f*(sum of s signs) (s even: alpha + 2f*t, |t| <= s/2), with
    sum_k PAF(y_k)(0) = 4m - 4(q - 1),  sum_k PAF(y_k)(h) = -4q (h != 0),  rows (+-2, 0, 0, 0).
Here:
  * the number Q of symmetric quadruples with the lag-0 energy and the rows is EXACT (DP over the orbit positions);
  * the probability that the remaining lags (h = 1 .. F/2) hit the target is the lattice-coordinate Gaussian of
    theory/q523_gauss.py fitted to quadruples drawn uniformly from that set (energy split drawn with its exact weight,
    then each block uniform in its (energy, row) class by backward sampling of the DP).
usage: python3 theory/nd_sym_s1est.py F q s [--N 200000] [--seed 1] [--nosym]   (q may be a prime power)"""
import sys, os, math, json
import numpy as np
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q523_gauss import gauss_bits

def letters(q, s):
    f = (q - 1) // s
    L = sorted({a + f * c for a in (1, -1) for c in range(-s, s + 1, 2)})
    return L, f

SYM = True

def orbit_positions(F):
    """negation orbits of Z_F (or single points with --nosym): representative and size"""
    if not SYM: return [(g, 1) for g in range(F)]
    seen = set(); reps = []
    for g in range(F):
        if g in seen: continue
        o = {g, (-g) % F}; seen |= o; reps.append((g, len(o)))
    return reps

def block_dp(F, L, cap, unit):
    """DP tables over the orbit positions: tab[i][(e, r)] = number of completions of positions i.. with energy e and
    row r (weighted by orbit size)"""
    reps = orbit_positions(F)
    tabs = [None] * (len(reps) + 1)
    tabs[len(reps)] = {(0, 0): 1}
    for i in range(len(reps) - 1, -1, -1):
        w = reps[i][1]; nd = defaultdict(int)
        for (e, r), c in tabs[i + 1].items():
            for y in L:
                e2 = e + w * unit[y]
                if e2 <= cap: nd[(e2, r + w * y)] += c
        tabs[i] = dict(nd)
    return reps, tabs

def sample_block(reps, tabs, L, unit, e, r, rng):
    """uniform configuration (orbit values) with energy e and row r"""
    vals = []
    for i, (g, w) in enumerate(reps):
        opts = []; wts = []
        for y in L:
            e2 = e - w * unit[y]; r2 = r - w * y
            c = tabs[i + 1].get((e2, r2), 0)
            if c: opts.append(y); wts.append(c)
        wts = np.array(wts, float); y = opts[rng.choice(len(opts), p=wts / wts.sum())]
        vals.append(y); e -= w * unit[y]; r -= w * y
    assert e == 0 and r == 0
    return vals

def expand(F, reps, vals):
    y = np.zeros(F, np.int64)
    for (g, w), v in zip(reps, vals):
        y[g] = v
        if SYM: y[(-g) % F] = v
    return y

def main(F, q, s, N=200000, seed=1):
    rng = np.random.default_rng(seed)
    L, f = letters(q, s)
    m = F * q; T0 = 4 * m - 4 * (q - 1); Th = -4 * q
    # energy units: (y^2 - 1) is divisible by g = gcd over letters
    g = 0
    for y in L: g = math.gcd(g, y * y - 1)
    unit = {y: (y * y - 1) // g for y in L}
    E = (T0 - 4 * F) // g                                         # total energy in units
    assert (T0 - 4 * F) % g == 0
    reps, tabs = block_dp(F, L, E, unit)
    base = tabs[0]
    c0 = defaultdict(int); c2 = defaultdict(int)
    for (e, r), c in base.items():
        if r == 0: c0[e] += c
        if r == 2: c2[e] += c
    # exact count of quadruples: block 0 row +2, blocks 1..3 row 0, energies summing to E
    parts = []; Q = 0
    for e0 in range(E + 1):
        for e1 in range(E + 1 - e0):
            for e2 in range(E + 1 - e0 - e1):
                e3 = E - e0 - e1 - e2
                w = c2.get(e0, 0) * c0.get(e1, 0) * c0.get(e2, 0) * c0.get(e3, 0)
                if w: parts.append((e0, e1, e2, e3, w)); Q += w
    Qbits = math.log2(Q) if Q else -math.inf
    out = dict(F=F, q=q, s=s, f=f, letters=L, energy_units=E, unit=g, quadruples_log2=round(Qbits, 2),
               energy_splits=len(parts))
    if not Q:
        out['note'] = 'no quadruple with the lag-0 energy and the rows'; print(json.dumps(out)); return out
    W = np.array([p[4] for p in parts], float); W /= W.sum()
    H = F // 2
    V = np.zeros((N, H), np.int64)
    for n in range(N):
        e0, e1, e2, e3, _ = parts[rng.choice(len(parts), p=W)]
        acc = np.zeros(H, np.int64)
        for e, r in ((e0, 2), (e1, 0), (e2, 0), (e3, 0)):
            y = expand(F, reps, sample_block(reps, tabs, L, unit, e, r, rng))
            for h in range(1, H + 1): acc[h - 1] += int((y * np.roll(y, -h)).sum())
        V[n] = acc
    T = np.full(H, Th, np.int64)
    gb = gauss_bits(V, T)
    # x 8: which block carries the row (4) and its sign (2)
    est = Qbits + 3 - gb['bits'] if gb['in_lattice'] else -math.inf
    out.update(lag_bits=gb, stage1_log2E=round(est, 2), N=N)
    print(json.dumps(out), flush=True)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; opt = {"--N": "200000", "--seed": "1"}; pos = []; i = 0
    if "--nosym" in a: SYM = False; a.remove("--nosym")
    while i < len(a):
        if a[i] in opt: opt[a[i]] = a[i + 1]; i += 2
        else: pos.append(int(a[i])); i += 1
    main(*pos, N=int(opt["--N"]), seed=int(opt["--seed"]))
