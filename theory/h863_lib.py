r"""Shared helpers for the v = 863 (order 3452) lift study: S7 host 2e93be9718324148 of the k4 bordered family
(EKLO Theorem 7 in Kumar's form, t = 3), G = Z2^2 x Z5 x Z43, K = Z5 x Z43, Q = Z2^2, multiplier x4 on Z43 (order 7).

The construction set-up, the certificate data/h3452_stage1.json (the Z43-compressed solution y) and the border are by
OpenAI Codex, directed by the repository owner (branch codex/order-3452-compression).  The bordered construction is
Kumar's / EKLO's (theory/bt_engine.py, theory/bt_border.py).  This module is new and imports neither bt_engine nor the
Codex scripts: it re-derives the orbit structure from scratch.

Conventions (checked against the certificate by check_y):
  H = Z2^2 x Z5, cell index k = q*5 + u, q = 2*q1 + q2 in 0..3, u in Z5.
  Z43^* = <4> x {coset reps}; CLS[t] = class index 0..5 of t != 0 (order of first appearance: reps 1, 2, 3, 6, 7, 9).
  Block i on cell h: e_i(h) at (h, 0), c_ij(h) on h x C_j.  y_i(h) = e_i(h) + 7 sum_j c_ij(h), S_i(h) = sum_j c_ij(h).
  G/Z5 level (Z2^2 x Z43): z_i(q, 0) = E_i(q) = sum_u e_i(q, u), z_i(q, t) = w_i(q, CLS[t]), w_i(q, j) = sum_u c_ij(q,u).
  Targets: full level F = 3440 at 0, -12 on K \ 0, +4 off K; H level 2936 / -516 / 172; G/Z5 level 3392 / -60 / 20.
"""
import json, math, itertools as it
import numpy as np

P = 43
ORD = 7
HSUB = [pow(4, k, P) for k in range(ORD)]
CLS = [-1] * P
REPS = []
for _t in range(1, P):
    if CLS[_t] < 0:
        for _h in HSUB:
            CLS[_t * _h % P] = len(REPS)
        REPS.append(_t)
assert REPS == [1, 2, 3, 6, 7, 9], REPS
NEG = [CLS[(P - r) % P] for r in REPS]            # class of -C_j
assert NEG == [1, 0, 3, 2, 5, 4], NEG              # pairs (1,2), (3,6), (7,9)
PAIRS = [(0, 1), (2, 3), (4, 5)]                    # class-index pairs {C, -C}
PAIR_REP = [REPS[a] for a, b in PAIRS]              # one lag per pair: 1, 3, 7


def load_cert(path='data/h3452_stage1.json'):
    d = json.load(open(path))
    Y = np.array(d['y'], dtype=np.int64)
    E = np.where(np.mod(Y, 14) == 1, 1, -1).astype(np.int64)
    assert ((Y - E) % 14 == 0).all()
    S = (Y - E) // 7
    assert (np.abs(S) <= 6).all() and (S % 2 == 0).all()
    return Y, E, S


def hpaf(y, a, b, s):
    """PAF of a length-20 sequence on H = Z2^2 x Z5 at lag (a, b, s)"""
    return sum(y[(q1 * 2 + q2) * 5 + u] * y[((q1 ^ a) * 2 + (q2 ^ b)) * 5 + (u + s) % 5]
               for q1 in (0, 1) for q2 in (0, 1) for u in range(5))


def check_y(Y):
    for a, b, s in it.product((0, 1), (0, 1), range(5)):
        v = sum(hpaf(y, a, b, s) for y in Y)
        tgt = 2936 if (a, b, s) == (0, 0, 0) else (-516 if (a, b) == (0, 0) else 172)
        assert v == tgt, (a, b, s, v, tgt)
    return True


# ------------------------------------------------------------------ 5 x 6 sign tables with fixed row sums
def table_options(rowsums, ncol=6):
    """all column-sum vectors w (odd entries) of +-1 tables with the given row sums (len(rowsums) rows, ncol columns),
    with exact multiplicities (number of tables).  Returns (W: (m, ncol) int64, mult: (m,) python ints as object)."""
    nplus = [(s + ncol) // 2 for s in rowsums]
    subsets = {n: [c for c in it.combinations(range(ncol), n)] for n in set(nplus)}
    states = {tuple([0] * ncol): 1}
    for n in nplus:
        new = {}
        for st, m in states.items():
            for sub in subsets[n]:
                a = list(st)
                for j in sub: a[j] += 1
                a = tuple(a)
                new[a] = new.get(a, 0) + m
        states = new
    R = len(rowsums)
    keys = sorted(states)
    W = np.array([[2 * a - R for a in k] for k in keys], dtype=np.int64)
    mult = np.array([states[k] for k in keys], dtype=object)
    return W, mult


def cell_lists(S):
    """per (i, q): (W, mult) for the 5 x 6 table of block i, Z2^2-cell q (rows u = 0..4 with sums S_i(q, u))"""
    out = {}
    for i in range(4):
        for q in range(4):
            out[i, q] = table_options([int(S[i, q * 5 + u]) for u in range(5)])
    return out


# ------------------------------------------------------------------ G/Z5 level (Z2^2 x Z43)
def g5_seq(Ei, wi):
    """explicit length-172 sequence z(q, t), index q*43 + t, from E_i(q) (4,) and w_i (4, 6)"""
    z = np.empty(4 * P, dtype=np.int64)
    for q in range(4):
        z[q * P] = Ei[q]
        for t in range(1, P):
            z[q * P + t] = wi[q][CLS[t]]
    return z


def g5_paf(z, qa, t0):
    s = 0
    for q in range(4):
        for t in range(P):
            s += int(z[q * P + t]) * int(z[(q ^ qa) * P + (t + t0) % P])
    return s


def g5_target(qa, t0):
    return 3392 if (qa, t0) == (0, 0) else (-60 if qa == 0 else 20)


LAGS16 = [(qa, t0) for qa in range(4) for t0 in [0] + PAIR_REP]


def g5_check(E4, w4):
    """direct check of the G/Z5 level on explicit sequences: E4 (4,4) singleton sums, w4 (4,4,6); all 172 lags"""
    zs = [g5_seq(E4[i], w4[i]) for i in range(4)]
    bad = []
    for qa in range(4):
        for t0 in range(P):
            v = sum(g5_paf(z, qa, t0) for z in zs)
            if v != g5_target(qa, t0): bad.append((qa, t0, v))
    return bad


# cyclotomic numbers for the class-level PAF: for t0 in class k, #{t in C_j : t + t0 in C_l}
def _cyc_counts():
    M = np.zeros((6, 6, 6), dtype=np.int64)       # M[k, j, l] with t0 = REPS[k]
    for k in range(6):
        t0 = REPS[k]
        for t in range(1, P):
            s = (t + t0) % P
            if s: M[k, CLS[t], CLS[s]] += 1
    return M
CYC = _cyc_counts()


def g5_sig_fast(E4, w4):
    """the 16 G/Z5 lag values (LAGS16 order) of the summed PAF, from class data (vectorised over leading dims).
    E4: (..., 4, 4) [block, q]; w4: (..., 4, 4, 6) [block, q, class].  Returns (..., 16) int64."""
    E4 = np.asarray(E4, np.int64); w4 = np.asarray(w4, np.int64)
    out = []
    for qa in range(4):
        qq = [q ^ qa for q in range(4)]
        Eb = E4[..., qq]; wb = w4[..., qq, :]
        # t0 = 0: E E' + 7 sum_j w_j w'_j
        v0 = (E4 * Eb).sum((-1, -2)) + 7 * (w4 * wb).sum((-1, -2, -3))
        out.append(v0)
        for k in (0, 2, 4):
            t0 = REPS[k]
            # t = 0 term: E(q) * w'(q^qa, cls(t0)); t = -t0 term: w(q, cls(-t0)) * E'(q^qa)
            v = (E4 * wb[..., k]).sum((-1, -2)) + (w4[..., NEG[k]] * Eb).sum((-1, -2))
            # t, t + t0 both nonzero: sum_{j,l} CYC[k, j, l] w(q, j) w'(q^qa, l)
            v = v + np.einsum('jl,...bqj,...bql->...', CYC[k], w4, wb)
            out.append(v)
    return np.stack(out, -1)


TARGET16 = np.array([g5_target(qa, t0) for qa, t0 in LAGS16], dtype=np.int64)


# ------------------------------------------------------------------ Z43-level (class-constant) PAF, vectorised
def z43_paf3(z0, zc):
    """PAF at the three class-pair lags (reps 1, 3, 7) of Z43 sequences constant on the six classes.
    z0: (...,) value at 0; zc: (..., 6) class values.  Returns (..., 3) int64."""
    z0 = np.asarray(z0, np.int64); zc = np.asarray(zc, np.int64)
    out = []
    for k in (0, 2, 4):
        v = z0 * (zc[..., k] + zc[..., NEG[k]]) + np.einsum('jl,...j,...l->...', CYC[k], zc, zc)
        out.append(v)
    return np.stack(out, -1)


def z43_paf0(z0, zc):
    z0 = np.asarray(z0, np.int64); zc = np.asarray(zc, np.int64)
    return z0 * z0 + 7 * (zc * zc).sum(-1)


# Z2^2 characters psi = (a1, a2): psi(q) = (-1)^(a1 q1 + a2 q2), q = 2 q1 + q2.  Index psi = 2 a1 + a2.
PSI = np.array([[(-1) ** (a1 * (q >> 1) + a2 * (q & 1)) for q in range(4)] for a1 in (0, 1) for a2 in (0, 1)],
               dtype=np.int64)
# per-psi Z43 targets (lag form, divided by nothing): psi0: 3452 delta; psi != 0: 3372 at 0, -80 off 0
PSI_TGT3 = np.array([[0, 0, 0]] + [[-80, -80, -80]] * 3, dtype=np.int64)
PSI_TGT0 = np.array([3452, 3372, 3372, 3372], dtype=np.int64)


def psi_sig(E4, w4):
    """(..., 4 psi, 3) Z43 PAF values of z_{i,psi} summed over blocks, from E4 (...,4,4), w4 (...,4,4,6)"""
    z0 = np.einsum('pq,...iq->...ip', PSI, E4)               # (..., 4 blocks, 4 psi)
    zc = np.einsum('pq,...iqj->...ipj', PSI, w4)             # (..., 4, 4, 6)
    return z43_paf3(z0, zc).sum(-3)


# ------------------------------------------------------------------ full level G = Z2^2 x Z5 x Z43
# index g = ((q * 5 + u) * 43 + t), q = 2 q1 + q2.  Lag classes under S (t -> 4t) and negation (q, u, t) -> (q, -u, -t).
def _lag_classes():
    rep = {}
    reps = []
    for q in range(4):
        for u in range(5):
            for t in range(P):
                g = (q, u, t)
                if g in rep: continue
                orb = set()
                for s in HSUB:
                    orb.add((q, u, t * s % P)); orb.add((q, (-u) % 5, (-t * s) % P))
                for h in orb: rep[h] = len(reps)
                reps.append(g)
    return rep, reps
LAGCLS, LAGREPS = _lag_classes()
assert len(LAGREPS) == 72
LAGREP_IDX = np.array([(q * 5 + u) * P + t for q, u, t in LAGREPS], dtype=np.int64)
LAGCLS_SIZE = np.bincount([LAGCLS[(q, u, t)] for q in range(4) for u in range(5) for t in range(P)], minlength=72)


def full_target():
    T = np.empty(72, np.int64)
    for c, (q, u, t) in enumerate(LAGREPS):
        T[c] = 3440 if (q, u, t) == (0, 0, 0) else (-12 if q == 0 else 4)
    return T
FULL_TARGET = full_target()

# character classes (psi, theta, phi-class) for the "mixed" bookkeeping: a lag class is H-level if its t = 0 ... (the
# equations are read in lag form; which combinations are fixed by y / w is found by the lattice reduction)


def assemble(e, c):
    """explicit block sequences. e: (..., 20) singleton signs; c: (..., 20, 6) class signs -> (..., 860) int8"""
    e = np.asarray(e); c = np.asarray(c)
    lead = e.shape[:-1]
    x = np.empty(lead + (20, P), np.int8)
    x[..., 0] = e
    x[..., 1:] = c[..., np.array(CLS[1:])]
    return x.reshape(lead + (20 * P,))


def class_paf(x, chunk=4096):
    """class PAF of block sequences x (..., 860): returns (..., 72) int64 via 4-D FFT (chunked to bound memory)"""
    x = np.asarray(x)
    lead = x.shape[:-1]
    xf = x.reshape(-1, 20 * P)
    out = np.empty((len(xf), 72), np.int64)
    for a in range(0, len(xf), chunk):
        xc = xf[a:a + chunk].astype(np.float64).reshape(-1, 2, 2, 5, P)
        F = np.fft.rfftn(xc, axes=(-4, -3, -2, -1))
        pa = np.fft.irfftn(F.real ** 2 + F.imag ** 2, s=(2, 2, 5, P), axes=(-4, -3, -2, -1))
        out[a:a + chunk] = np.rint(pa.reshape(len(xc), 20 * P)[:, LAGREP_IDX]).astype(np.int64)
    return out.reshape(lead + (72,))


# ------------------------------------------------------------------ samplers
_SUBS = {n: np.array([[1 if j in s else 0 for j in range(6)] for s in it.combinations(range(6), n)], np.int8)
         for n in range(7)}


def sample_y_fibre(S, N, rng):
    """uniform class signs c (N, 4, 20, 6) with row sums S (4, 20)"""
    c = np.empty((N, 4, 20, 6), np.int8)
    for i in range(4):
        for k in range(20):
            n = (int(S[i, k]) + 6) // 2
            sub = _SUBS[n]
            c[:, i, k] = 2 * sub[rng.integers(len(sub), size=N)] - 1
    return c


class TableSampler:
    """uniform +-1 R x C tables (R = 5, C = 6) with given row sums and column sums (exact DP over rows)"""
    _cache = {}

    def __init__(self, rowsums, colsums):
        self.R, self.C = len(rowsums), len(colsums)
        self.rp = [(s + self.C) // 2 for s in rowsums]
        cp = [(s + self.R) // 2 for s in colsums]
        base = self.R + 1
        self.base = base
        self.pw = np.array([base ** j for j in range(self.C)], np.int64)
        start = int(sum(a * b for a, b in zip(cp, self.pw)))
        self.start = start
        nst = base ** self.C
        # count[r][state] = ways to fill rows r..R-1 exactly consuming 'state' (remaining column demands)
        digits = (np.arange(nst)[:, None] // self.pw[None, :]) % base          # (nst, C)
        self.digits = digits
        cnt = [None] * (self.R + 1)
        z = np.zeros(nst, object); z[0] = 1
        cnt[self.R] = z
        self.nxt = []
        for r in range(self.R - 1, -1, -1):
            sub = _SUBS[self.rp[r]]                                            # (m, C)
            ok = (digits[:, None, :] >= sub[None, :, :]).all(-1)               # (nst, m)
            nx = np.arange(nst)[:, None] - (sub.astype(np.int64) @ self.pw)[None, :]
            nx = np.where(ok, nx, 0)
            c = np.zeros(nst, object)
            for s in range(len(sub)):
                vals = cnt[r + 1][nx[:, s]]
                c = c + np.where(ok[:, s], vals, 0)
            cnt[r] = c
            self.nxt.insert(0, (sub, ok, nx))
        self.cnt = cnt
        self.total = int(cnt[0][start])
        self.cntf = [np.array([float(v) for v in cc]) for cc in cnt]

    def sample(self, N, rng):
        T = np.empty((N, self.R, self.C), np.int8)
        st = np.full(N, self.start, np.int64)
        for r in range(self.R):
            sub, ok, nx = self.nxt[r]
            w = np.where(ok[st], self.cntf[r + 1][nx[st]], 0.0)               # (N, m)
            cs = np.cumsum(w, 1); u = rng.random(N) * cs[:, -1]
            k = (cs < u[:, None]).sum(1)
            T[:, r] = 2 * sub[k] - 1
            st = nx[st, k]
        assert (st == 0).all()
        return T
