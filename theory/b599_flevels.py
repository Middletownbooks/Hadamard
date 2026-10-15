r"""Stage-1 level budget for the c = 6 family (theory/b599_fam.py) given a Q-compression quadruple z [E, counts exact].

As theory/b599_levels.py, for any (n, q): y-quadruples with Q-compressions z (4 x 6) and sum_k E_k = T1(0) are counted
exactly (DP over cosets) and sampled uniformly; the lattice Gaussian (theory/q523_gauss.gauss_bits) of the remaining
stage-1 equations gives log2E(stage 1 | z) = H1 - C_full; the S-level (compression along Z_2) and the D-freedom h_D
(exact per sample) give the S/D budget.  Also: per-class totals over all Q-level solutions (theory/b599_qlevel.py).
usage: python3 theory/b599_flevels.py --n 9 --q 11 [--classes] [--N 3000] [--out FILE]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b599_fam import Fam
from q523_gauss import gauss_bits
import b599_qlevel as QL


class BlockSampler:
    """uniform sampling of block vectors (family F) with Q-compression z and energy E"""
    def __init__(self, F, z, Emax):
        self.F = F; self.z = np.asarray(z, np.int64); J = F.J
        TUP = np.array(list(itertools.product(F.LET.tolist(), repeat=J)), np.int64)
        TSUM = TUP @ F.WJ; TEN = (TUP * TUP) @ F.WJ
        self.TUP = TUP
        self.cos = []
        for q in range(6):
            d = defaultdict(list)
            for i in np.nonzero(TSUM == self.z[q])[0]: d[int(TEN[i])].append(i)
            self.cos.append({e: np.array(ix) for e, ix in d.items()})
        self.cnt = [None] * 7; self.cnt[6] = {0: 1}
        for q in range(5, -1, -1):
            nd = defaultdict(int)
            for E, c in self.cnt[q + 1].items():
                for e, ix in self.cos[q].items():
                    if E + e <= Emax: nd[E + e] += c * len(ix)
            self.cnt[q] = dict(nd)
        self.N = self.cnt[0]

    def draw(self, E, rng):
        J = self.F.J; out = np.zeros(6 * J, np.int64); R = E
        for q in range(6):
            es = [e for e in self.cos[q] if (R - e) in self.cnt[q + 1]]
            w = np.array([len(self.cos[q][e]) * self.cnt[q + 1][R - e] for e in es], float)
            e = es[rng.choice(len(es), p=w / w.sum())]
            ix = self.cos[q][e]; out[J * q:J * q + J] = self.TUP[ix[rng.integers(len(ix))]]; R -= e
        return out


def energy_law(samplers, T0):
    N = [s.N for s in samplers]
    P01 = defaultdict(list)
    for e0, c0 in N[0].items():
        for e1, c1 in N[1].items():
            if e0 + e1 <= T0: P01[e0 + e1].append((e0, e1, c0 * c1))
    splits = []; wts = []
    for e2, c2 in N[2].items():
        for e3, c3 in N[3].items():
            for (e0, e1, c01) in P01.get(T0 - e2 - e3, []):
                splits.append((e0, e1, e2, e3)); wts.append(c01 * c2 * c3)
    tot = sum(wts)
    if not tot: return [], None, None
    return splits, np.array([w / tot for w in wts], float), math.log2(tot)


def levels_of(F, P):
    """P: (n, nl) summed class signatures -> dict of level vectors"""
    J = F.J; n = P.shape[0]
    S = P[:, :2 * J] + P[:, 2 * J:]
    return dict(full=P, S=S)


def hD_block(F, v, z, E):
    """log2 #block vectors with the same sigma, Q-compression z and energy E"""
    J = F.J; sg, _ = F.sd_of(v); LET = F.LET.tolist()
    total = {0: 1}
    for t in range(3):
        st = {(0, 0): 1}
        for j in range(J):
            c = J * t + j; w = 1 if j == 0 else 2; s = int(sg[0, c])
            opts = [(a, b) for a in LET for b in LET if a + b == 2 * s]
            nd = defaultdict(int)
            for (S0, Ee), m in st.items():
                for (a, b) in opts: nd[(S0 + w * a, Ee + w * (a * a + b * b))] += m
            st = nd
        per = defaultdict(int)
        for (S0, Ee), m in st.items():
            if S0 == int(z[t]): per[Ee] += m
        nt = defaultdict(int)
        for E1, m1 in total.items():
            for E2, m2 in per.items(): nt[E1 + E2] += m1 * m2
        total = nt
    return math.log2(total.get(int(E), 0))


def budget(F, z, N=3000, seed=1, nh=150):
    rng = np.random.default_rng(seed); t0 = time.time()
    T0 = int(F.T1SIG[0])
    samp = [BlockSampler(F, z[k], T0) for k in range(4)]
    splits, w, H1 = energy_law(samp, T0)
    if H1 is None: return dict(H1=None, note="energy infeasible")
    pick = rng.choice(len(splits), size=N, p=w)
    Y = np.zeros((N, 4, F.no), np.int64)
    for i, si in enumerate(pick):
        for k in range(4): Y[i, k] = samp[k].draw(splits[si][k], rng)
    P = F.sig(Y.reshape(4 * N, F.no)).reshape(N, 4, F.nl).sum(1)
    V = levels_of(F, P); Tt = levels_of(F, F.T1SIG[None, :])
    out = dict(H1=round(H1, 2))
    for name in ('full', 'S'):
        gb = gauss_bits(V[name], Tt[name][0])
        out[name] = dict(bits=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'))
    hs = [sum(hD_block(F, Y[i, k], z[k], splits[pick[i]][k]) for k in range(4)) for i in range(min(nh, N))]
    out['h_D'] = round(float(np.mean(hs)), 2)
    if out['full']['in_lattice'] and out['S']['in_lattice']:
        Cf, CS = out['full']['bits'], out['S']['bits']
        out['log2E'] = round(H1 - Cf, 2); out['D_def'] = round(Cf - CS - out['h_D'], 2)
        out['log2_S_sols'] = round(H1 - out['h_D'] - CS, 2)
    else:
        out['log2E'] = -math.inf
    out['secs'] = round(time.time() - t0, 1)
    return out


def class_quads(F, maxc=None):
    """for every Q-level class: every row solution -> z-quadruple (shapes in class order, c_k odd)"""
    S, C = QL.qlevel_classes()
    out = []
    for ci, cl in enumerate(C):
        for c in QL.rows_for(cl, F.v):
            z = np.array([[c[k] + 2 * cl[k][q] for q in range(6)] for k in range(4)], np.int64)
            out.append((ci, cl, c, z))
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    n, q = int(opt('--n', '9')), int(opt('--q', '11'))
    N = int(opt('--N', '3000')); outp = opt('--out', '')
    F = Fam(n, q)
    done = set()
    if outp and os.path.exists(outp):
        for l in open(outp):
            try: d = json.loads(l); done.add((d['class'], tuple(d['c'])))
            except Exception: pass
    fo = open(outp, 'a') if outp else None
    for ci, cl, c, z in class_quads(F):
        if (ci, tuple(c)) in done: continue
        r = budget(F, z, N=N)
        r.update(n=n, q=q, v=F.v, **{'class': ci}, c=list(c), rows=[int(x) for x in z.sum(1)], shapes=[list(u) for u in cl])
        s = json.dumps(r); print(s, flush=True)
        if fo: fo.write(s + "\n"); fo.flush()
