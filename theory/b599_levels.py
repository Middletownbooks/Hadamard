r"""Stage-1 level budget of the v = 599 lead slice for a bordered Q-compression quadruple [E, counts exact].

Given z = (z_0, ..., z_3) (Q-compressions, 4 x 6), the y-quadruples with these compressions and lag-0 energy
sum_k E_k = 2176 are counted exactly (DP over the six cosets of each block: N_k(E)); a uniform sample of them
(energy split drawn with weight prod_k N_k(E_k), each block uniform given (z_k, E_k)) gives the lattice Gaussian
(theory/q523_gauss.gauss_bits) of
  full : sum_k PAF_{G'}(y_k) at the 20 lag classes (stage 1);
  S    : the S-level (compression of the PAF along the Z_2 axis: 10 classes of G'' = Z_3 x Z_9);
  L3   : the Q x Z_3 level (characters of order dividing 3 on Z_9: compression along 3 Z_9);
  S3   : S-level of L3 (characters trivial on Z_2 and of order dividing 3 on Z_9).
h-values (exact, per sample, averaged):  H1 = log2 #quadruples;  H_S = log2 #distinct sigma-quadruples,
h_D = E[log2 #y with the sampled sigma] (the D-freedom once S is fixed, z and energies included).
Then  log2E(stage 1 | z) = H1 - C_full,  S-solutions per z = H_S - C_S,  D_def = C(D|S) - h_D.
usage: python3 theory/b599_levels.py [--zi I] [--N 6000] [--seed 1]"""
import sys, os, json, math, itertools, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b599_lib as B
from q523_gauss import gauss_bits

# all 6^5 coset tuples (y_0..y_4), their coset sums and energies
TUP = np.array(list(itertools.product(B.LET.tolist(), repeat=5)), np.int64)
WJ = np.array([1, 2, 2, 2, 2], np.int64)
TSUM = TUP @ WJ
TEN = (TUP * TUP) @ WJ


class BlockSampler:
    """uniform sampling of block vectors with Q-compression z (6 values) and energy E"""
    def __init__(self, z, Emax=2176):
        self.z = np.asarray(z, np.int64)
        self.cos = []                                  # per coset: dict energy -> array of tuples
        for q in range(6):
            m = TSUM == self.z[q]
            d = defaultdict(list)
            for i in np.nonzero(m)[0]: d[int(TEN[i])].append(i)
            self.cos.append({e: np.array(ix) for e, ix in d.items()})
        # backward DP: cnt[q][E] = number of ways for cosets q..5 with total energy E
        self.cnt = [None] * 7; self.cnt[6] = {0: 1}
        for q in range(5, -1, -1):
            nd = defaultdict(int)
            for E, c in self.cnt[q + 1].items():
                for e, ix in self.cos[q].items():
                    if E + e <= Emax: nd[E + e] += c * len(ix)
            self.cnt[q] = dict(nd)
        self.N = self.cnt[0]

    def draw(self, E, rng, n=1):
        out = np.zeros((n, 30), np.int64)
        for i in range(n):
            R = E
            for q in range(6):
                es = [e for e in self.cos[q] if (R - e) in self.cnt[q + 1]]
                w = np.array([len(self.cos[q][e]) * self.cnt[q + 1][R - e] for e in es], float)
                e = es[rng.choice(len(es), p=w / w.sum())]
                ix = self.cos[q][e]; t = ix[rng.integers(len(ix))]
                out[i, 5 * q:5 * q + 5] = TUP[t]; R -= e
            assert R == 0
        return out


def quad_energy_law(samplers, T0=2176):
    """weights of energy splits (E_0..E_3) with sum T0: prod_k N_k(E_k); returns (list of splits, weights, log2 total)"""
    N = [s.N for s in samplers]
    P01 = defaultdict(list)
    for e0, c0 in N[0].items():
        for e1, c1 in N[1].items():
            if e0 + e1 <= T0: P01[e0 + e1].append((e0, e1, c0 * c1))
    splits = []; wts = []
    for e2, c2 in N[2].items():
        for e3, c3 in N[3].items():
            rest = T0 - e2 - e3
            for (e0, e1, c01) in P01.get(rest, []):
                splits.append((e0, e1, e2, e3)); wts.append(c01 * c2 * c3)
    tot = sum(wts)
    return splits, np.array([w / tot for w in wts], float), (math.log2(tot) if tot else None)


def level_vectors(Y):
    """Y: (n, 4, 30) -> dict of level vectors (sum over blocks)"""
    n = Y.shape[0]
    P = B.sig(Y.reshape(n * 4, 30)).reshape(n, 4, 20).sum(1)            # 20 classes
    full = P
    # S-level: PAF compressed along Z_2: classes of G'' (t in {0,1}, j in 0..4): P(0,t,j) + P(1,t,j)
    S = P[:, :10] + P[:, 10:]
    # L3 (Q x Z_3): sum over the coset j + 3 Z_9 of the full PAF: needs the full PAF; use class sizes
    # full PAF at lag (e, t, j9) equals class value; compress j9 -> j9 mod 3
    L3 = np.zeros((n, 2, 2, 3), np.int64)          # (e, t>0, j mod 3) sum over j9 with that residue
    for e in range(2):
        for t in range(2):
            for j9 in range(9):
                L = 10 * e + 5 * t + min(j9, 9 - j9)
                L3[:, e, t, j9 % 3] += P[:, L]
    L3 = L3.reshape(n, 12)
    S3 = (L3.reshape(n, 2, 6)).sum(1)
    return dict(full=full, S=S, L3=L3, S3=S3)


def targets():
    P = B.T1SIG[None, :]
    out = level_vectors_from_sig(P)
    return {k: v[0] for k, v in out.items()}


def level_vectors_from_sig(P):
    n = P.shape[0]
    S = P[:, :10] + P[:, 10:]
    L3 = np.zeros((n, 2, 2, 3), np.int64)
    for e in range(2):
        for t in range(2):
            for j9 in range(9):
                L = 10 * e + 5 * t + min(j9, 9 - j9)
                L3[:, e, t, j9 % 3] += P[:, L]
    L3 = L3.reshape(n, 12)
    return dict(full=P, S=S, L3=L3, S3=L3.reshape(n, 2, 6).sum(1))


# sigma multiplicities: pairs (y0, y1) -> sigma = (y0 + y1)/2
PAIRS = [(a, b) for a in B.LET.tolist() for b in B.LET.tolist()]


def h_values(Y, samplers, Es):
    """exact H_S-type quantities for one quadruple: log2 #y with the same sigma (z and E_k fixed per block)"""
    hD = 0.0
    for k in range(4):
        v = Y[k]; sg, dl = B.sd_of(v)
        z = samplers[k].z
        # count block vectors with the same sigma, same Q-compression z_k and same energy E_k: per cell pairs (y0, y1)
        # with (y0 + y1)/2 = sigma[c]; constraint: coset sums of e = 0 and e = 1 layers and energy -> DP over cells
        # state: (sum of layer-0 coset t, energy) per t handled by DP over cells ordered by t
        st = {(0, 0): 1}        # per t: (layer-0 coset sum, energy) ; layer-1 sum = 2 sigma-sum - layer-0
        total = {0: 1}
        for t in range(3):
            st = {(0, 0): 1}
            for j in range(5):
                c = 5 * t + j; w = int(B.W15[c]); s = int(sg[0, c])
                opts = [(a, b) for (a, b) in PAIRS if a + b == 2 * s]
                nd = defaultdict(int)
                for (S0, E), m in st.items():
                    for (a, b) in opts:
                        nd[(S0 + w * a, E + w * (a * a + b * b))] += m
                st = nd
            # layer-0 coset sum must equal z[(0, t)] = z[t]; layer-1 = z[3 + t]
            per = defaultdict(int)
            for (S0, E), m in st.items():
                if S0 == int(z[t]): per[E] += m
            nt = defaultdict(int)
            for E1, m1 in total.items():
                for E2, m2 in per.items(): nt[E1 + E2] += m1 * m2
            total = nt
        cntk = total.get(int(Es[k]), 0)
        hD += math.log2(cntk)
    return hD


def run(z, N=6000, seed=1, log=print):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    samplers = [BlockSampler(z[k]) for k in range(4)]
    splits, w, H1 = quad_energy_law(samplers)
    if H1 is None: return dict(H1=None)
    pick = rng.choice(len(splits), size=N, p=w)
    Y = np.zeros((N, 4, 30), np.int64)
    for i, si in enumerate(pick):
        for k in range(4): Y[i, k] = samplers[k].draw(splits[si][k], rng)[0]
    V = level_vectors(Y); Tt = level_vectors_from_sig(B.T1SIG[None, :])
    out = dict(H1=round(H1, 2), secs_sample=round(time.time() - t0, 1))
    for name in ('full', 'S', 'L3', 'S3'):
        gb = gauss_bits(V[name], Tt[name][0])
        out[name] = dict(bits=gb.get('bits'), rank=gb.get('rank'), in_lattice=gb.get('in_lattice'), mahal=gb.get('mahal'))
    # h_D over a subsample
    hs = [h_values(Y[i], samplers, splits[pick[i]]) for i in range(min(200, N))]
    out['h_D'] = round(float(np.mean(hs)), 2)
    if out['full']['in_lattice'] and out['S']['in_lattice']:
        Cf = out['full']['bits']; CS = out['S']['bits']
        out['log2E_stage1_given_z'] = round(H1 - Cf, 2)
        out['C_D_given_S'] = round(Cf - CS, 2)
        out['D_def'] = round(Cf - CS - out['h_D'], 2)
        out['log2_S_solutions_given_z'] = round(H1 - out['h_D'] - CS, 2)
    out['secs'] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    zi = int(opt('--zi', '0')); N = int(opt('--N', '6000')); seed = int(opt('--seed', '1'))
    import ab_border_transfer2 as BT
    T = BT.load(1388); comp = T['comp']
    a_ = [int(comp[k].sum()) for k in range(4)]
    sols = []
    for s in itertools.product(range(-40, 41, 2), repeat=4):
        r = [a_[k] + 6 * s[k] for k in range(4)]
        if sum(x * x for x in r) == 2396: sols.append(s)
    s = sols[zi]
    z = np.array([comp[k] + s[k] for k in range(4)])
    r = run(z, N, seed)
    r.update(zi=zi, s=list(s), rows=[int(x) for x in z.sum(1)], z=z.tolist())
    print(json.dumps(r))
