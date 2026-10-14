r"""Exact search for b >= 3 borders (templates) with a Singer-type symmetry [M/P].

Setting (report section 3.1).  Monomial compressions z_k = a_k delta (+ constants), sum a_k^2 = 4b, b = c - 1.  When
one a_k = 0 (or Q is an elementary abelian 2-group) every permutation of Q acts on the GS array Cb (applied as sigma or
rho sigma rho blockwise), so a group Gamma of order b acting regularly on Q \ {0} fixes Cb [P].  In "Singer
coordinates" (core block J = {inf} u Gamma, border orbits = 4 regular Gamma-orbits) Cb = Q (x) I_c with a 4 x 4
integer matrix Q, Q Q^T = 4b I.  Gamma-invariant templates are then exactly:
  K  = [e  (x) 1 | circ(f)],   Kc = [e' (x) 1 | circ(g)]       (e, e' in H(4); f, g: 4 x 4 arrays of +-1 functions on
                                                               Gamma with sum(f_IJ) = -e_IJ, sum(g_IJ) = -e'_IJ)
  frames:   every row-block (f_I0..f_I3) has sum_J PAF(f_IJ) = -4 off 0, and two row-blocks are orthogonal iff their
            cross-correlation sum_J [e_IJ e_I'J + PCF(f_IJ, f_I'J)(d)] vanishes for every d (shift-invariant);
  corner:   C[(I,i),(I',l)] = -(1/4c) sum_{J,J',x} K_(I,i)(J,x) Q[J',J] Kc_(I',l)(J',x) must be +-1, and block (I, I')
            depends only on row-block I of K and row-block I' of Kc.
Search: (A) first Kc row-block up to Gamma-affine maps and column signed permutations; (B) its partners by a hashed
join on (PAF sums, cross-correlation), then 3-cliques -> Kc frames; (C) for each frame and each Q (signed row
permutations of the class representatives), all K row-blocks with a flat corner against the frame (hashed join, keys
mod 8c and PAF sums), then 4-cliques of mutually orthogonal ones -> templates, each verified exactly.
The published c = 8 template (359) is of this form with Gamma = Z_7 (theory/ab_singer.py --c 8 is the control).
usage: python3 theory/ab_singer.py --c 10 --group cyc|z3z3 [--a 4,4,2,0] [--max-first N] [--part i/W] [--out DIR]"""
import sys, os, json, time, itertools
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------- the group Gamma (order b) and its affine maps
def group(kind, b):
    if kind == 'cyc':
        add = np.array([[(x + y) % b for y in range(b)] for x in range(b)])
        auts = [np.array([(m * x) % b for x in range(b)]) for m in range(1, b) if np.gcd(m, b) == 1]
    elif kind == 'z3z3':
        assert b == 9
        el = [(x // 3, x % 3) for x in range(9)]
        idx = lambda u, v: 3 * (u % 3) + (v % 3)
        add = np.array([[idx(el[x][0] + el[y][0], el[x][1] + el[y][1]) for y in range(9)] for x in range(9)])
        auts = []
        for (p, q, r, s) in itertools.product(range(3), repeat=4):
            if (p * s - q * r) % 3 == 0: continue
            auts.append(np.array([idx(p * el[x][0] + q * el[x][1], r * el[x][0] + s * el[x][1]) for x in range(9)]))
    elif kind == 'z2z2z2':
        assert b == 7
        raise ValueError("use cyc for b = 7")
    else:
        raise ValueError(kind)
    neg = np.array([int(np.nonzero(add[x] == 0)[0][0]) for x in range(b)])
    return add, neg, auts


# ---------------------------------------------------------------- sequences, PAF, cross-correlation
def sequences(b):
    return np.array([s for s in itertools.product((1, -1), repeat=b) if abs(sum(s)) == 1], np.int64)


def paf(S, add):
    b = add.shape[0]
    return np.stack([(S * S[:, add[:, d]]).sum(1) for d in range(1, b)], 1)


def pcf_matrix(h, add):
    """M[d, y] with (M g)[d] = sum_x h(x) g(x + d)"""
    b = len(h); M = np.zeros((b, b), np.int64)
    for d in range(b):
        for x in range(b): M[d, add[x, d]] += h[x]
    return M


def hkey(A, w):
    return (A.astype(np.int64) * w[None, :A.shape[1]]).sum(1)


def join(LK, RK):
    """pairs (i, j) with LK[i] == RK[j] (int64 hash keys)"""
    o = np.argsort(RK, kind='stable'); Rs = RK[o]
    lo = np.searchsorted(Rs, LK, 'left'); hi = np.searchsorted(Rs, LK, 'right')
    I = np.repeat(np.arange(len(LK)), hi - lo)
    J = np.concatenate([o[a:b_] for a, b_ in zip(lo, hi) if b_ > a]) if (hi > lo).any() else np.zeros(0, np.int64)
    return I, J


# ---------------------------------------------------------------- (A) first row-blocks up to symmetry
def first_actions(S, add, auts):
    """rep[i]: index of the sum(+1) representative of +-S[i]; acts: for every affine map x -> mu(x) + a of Gamma, the
    induced map on sequence indices followed by rep (column negations absorbed)"""
    b = add.shape[0]
    idx = {tuple(s): i for i, s in enumerate(S.tolist())}
    rep = np.array([i if sum(s) == 1 else idx[tuple(-x for x in s)] for i, s in enumerate(S.tolist())])
    acts = []
    for mu in auts:
        inv = np.argsort(mu)
        for a in range(b):
            na = int(np.nonzero(add[a] == 0)[0][0])
            pinv = np.array([inv[add[y, na]] for y in range(b)])
            acts.append(np.array([rep[idx[tuple(S[i][pinv].tolist())]] for i in range(len(S))]))
    return rep, np.array(acts)


def canon_first(rep, acts, g):
    """canonical form of a row-block (sequence indices) under affine maps of Gamma and column signed permutations"""
    r = rep[np.asarray(g)]
    return min(tuple(sorted(row.tolist())) for row in acts[:, r])


def first_blocks(S, P, add, auts, chunk=1 << 15):
    """multisets {i0 <= i1 <= i2 <= i3} of sum(+1) sequences with sum PAF = -4, up to affine maps of Gamma and column
    signed permutations (vectorised: every match is mapped by every affine map, sorted, encoded, minimised)"""
    rep, acts = first_actions(S, add, auts)
    pos = np.nonzero(S.sum(1) == 1)[0]
    Pp = P[pos]; n = len(pos)
    ii, jj = np.triu_indices(n)
    PS = Pp[ii] + Pp[jj]
    w = np.random.default_rng(1).integers(1, 2 ** 61, size=PS.shape[1], dtype=np.int64)
    LK = hkey(PS, w); RK = hkey(-4 - PS, w)
    I, J = join(LK, RK)
    Qd = np.stack([pos[ii[I]], pos[jj[I]], pos[ii[J]], pos[jj[J]]], 1)
    Qd = np.sort(Qd, 1)
    Qd = np.unique(Qd, axis=0)
    base = len(S)
    out = set()
    for c0 in range(0, len(Qd), chunk):
        q = Qd[c0:c0 + chunk]
        img = np.sort(acts[:, q], axis=2)                       # (n_acts, m, 4)
        code = ((img[..., 0] * base + img[..., 1]) * base + img[..., 2]) * base + img[..., 3]
        mn = code.min(0)
        out.update(np.unique(mn).tolist())
    res = []
    for code in sorted(out):
        t = []
        for _ in range(4): t.append(code % base); code //= base
        res.append(tuple(int(x) for x in reversed(t)))
    return sorted(res)


# ---------------------------------------------------------------- (B) partners and Kc frames
class Ctx:
    def __init__(self, kind, b):
        self.kind, self.b = kind, b
        self.add, self.neg, self.auts = group(kind, b)
        self.S = sequences(b); self.P = paf(self.S, self.add); self.sums = self.S.sum(1)
        n = len(self.S)
        self.pi = np.repeat(np.arange(n), n); self.pj = np.tile(np.arange(n), n)
        self.w = np.random.default_rng(7).integers(1, 2 ** 61, size=64, dtype=np.int64)
        self.PP = self.P[self.pi] + self.P[self.pj]
        self.hPL = hkey(self.PP, self.w[:32]); self.hPR = hkey(-4 - self.PP, self.w[:32])
        self.rep, self.acts = first_actions(self.S, self.add, self.auts)


def cross(ctx, h, gidx):
    """X[d] = sum_J [e'_hJ e'_gJ + PCF(h_J, g_J)(d)] for a block h (4 x b) and blocks given by sequence indices"""
    b = ctx.b; X = np.zeros((len(gidx), b), np.int64)
    for J in range(4):
        M = pcf_matrix(h[J], ctx.add); g = ctx.S[gidx[:, J]]
        X += g @ M.T + (-h[J].sum()) * (-ctx.sums[gidx[:, J]])[:, None]
    return X


def partners(ctx, h):
    """all row-blocks (index quadruples) g with sum PAF = -4 and cross(h, g) = 0"""
    b = ctx.b
    M = [pcf_matrix(h[J], ctx.add) for J in range(4)]; eh = -h.sum(1)
    S, sm = ctx.S, ctx.sums
    XL = S[ctx.pi] @ M[0].T + S[ctx.pj] @ M[1].T + (eh[0] * -sm[ctx.pi] + eh[1] * -sm[ctx.pj])[:, None]
    XR = S[ctx.pi] @ M[2].T + S[ctx.pj] @ M[3].T + (eh[2] * -sm[ctx.pi] + eh[3] * -sm[ctx.pj])[:, None]
    LK = hkey(np.concatenate([ctx.PP, XL], 1), ctx.w); RK = hkey(np.concatenate([-4 - ctx.PP, -XR], 1), ctx.w)
    I, J = join(LK, RK)
    G = np.stack([ctx.pi[I], ctx.pj[I], ctx.pi[J], ctx.pj[J]], 1)
    ok = ((ctx.P[G].sum(1) == -4).all(1)) & (cross(ctx, h, G) == 0).all(1)
    return G[ok]


def canon_block(ctx, g):
    """canonical form of a row-block (index quadruple) under shifts of Gamma and global negation"""
    b = ctx.b; best = None
    seqs = ctx.S[g]
    for a in range(b):
        sh = seqs[:, ctx.add[:, a]]
        for sg in (1, -1):
            t = tuple((sg * sh).ravel().tolist())
            if best is None or t > best: best = t
    return best


def frames_from(ctx, h, G):
    """3-cliques of mutually orthogonal partner classes -> list of frames (4 x 4 x b arrays)"""
    cls = {}
    for g in G:
        k = canon_block(ctx, g)
        if k not in cls: cls[k] = g
    reps = list(cls.values())
    R = np.array(reps) if reps else np.zeros((0, 4), np.int64)
    blocks = [ctx.S[g] for g in R]
    n = len(R); adj = np.zeros((n, n), bool)
    for i in range(n):
        if i + 1 < n:
            X = cross(ctx, blocks[i], R[i + 1:])
            adj[i, i + 1:] = (X == 0).all(1)
    adj |= adj.T
    out = []
    for i in range(n):
        for j in np.nonzero(adj[i, i + 1:])[0] + i + 1:
            for k in np.nonzero(adj[i, j + 1:] & adj[j, j + 1:])[0] + j + 1:
                out.append((np.stack([h, blocks[i], blocks[j], blocks[k]]), (R[i], R[j], R[k])))
    return out, n


# ---------------------------------------------------------------- (C) admissible K row-blocks and K frames
def frame_tables(ctx, Kc):
    """Q-independent part of the corner rows: U[J'][f, I' b + l] = sum_t f(t) g_I'J'(t - l) + e(f) e'_I'J'"""
    b = ctx.b; S, sm = ctx.S, ctx.sums
    ep = -Kc.sum(2)
    U = np.zeros((4, len(S), 4 * b), np.int64)
    for Jp in range(4):
        for Ip in range(4):
            for l in range(b):
                v = Kc[Ip, Jp][ctx.add[:, ctx.neg[l]]]          # v[t] = g_I'J'(t - l)
                U[Jp, :, Ip * b + l] = S @ v + (-sm) * ep[Ip, Jp]
    return U


def corner_parity_ok(Kc, Q):
    """necessary condition at theta = 1: C(1) = -e Q^T e'^T / 4 must be odd for some e in H(4); with M = Q^T e'^T the
    entries of e M are = 4 mod 8 iff, in every column of M, the number of entries = 6 mod 8 has the parity of the
    row of e (H(4) rows all have one parity), so all columns of M need the same such parity"""
    ep = -Kc.sum(2)                                   # e'[I', J']
    M = Q.T @ ep.T                                    # M[J, I'] = sum_J' Q[J', J] e'[I', J']
    par = ((np.mod(M, 8) == 6).sum(0)) % 2
    return bool((par == par[0]).all())


def admissible(ctx, Kc, Q, c, U=None, nkey=12):
    """all K row-blocks f (4 x b, with e = -sums) such that the corner rows against the Kc frame are +-1:
    with T_J = sum_J' Q[J', J] U_J' the corner entry is -(1/4c) sum_J T_J(f_J); hashed join on the PAF sums and
    nkey of the residues mod 8c, then every entry checked exactly"""
    if U is None: U = frame_tables(ctx, Kc)
    m = 8 * c
    T = np.einsum('aj,afk->jfk', Q, U)                  # T[J] = sum_J' Q[J', J] U[J']
    sel = np.arange(nkey)
    L = (T[0][:, sel][ctx.pi] + T[1][:, sel][ctx.pj])
    R = (T[2][:, sel][ctx.pi] + T[3][:, sel][ctx.pj])
    LK = ctx.hPL + hkey(np.mod(L, m), ctx.w[32:])
    RK = ctx.hPR + hkey(np.mod(4 * c - R, m), ctx.w[32:])
    I, J = join(LK, RK)
    if len(I) == 0: return np.zeros((0, 4), np.int64)
    tot = T[0][ctx.pi[I]] + T[1][ctx.pj[I]] + T[2][ctx.pi[J]] + T[3][ctx.pj[J]]
    ok = (np.abs(tot) == 4 * c).all(1) & ((ctx.PP[I] + ctx.PP[J]) == -4).all(1)
    return np.stack([ctx.pi[I[ok]], ctx.pj[I[ok]], ctx.pi[J[ok]], ctx.pj[J[ok]]], 1)


try:
    import numba as _nb

    @_nb.njit(cache=True)
    def _pair_keys(TA, TB, hP, w, m, shift, neg):
        """key[i*n + j] = hP[i*n + j] + sum_k r_k w_k, r_k = (+-(TA[i,k] + TB[j,k]) + shift) mod m"""
        n = TA.shape[0]; K = TA.shape[1]
        out = np.empty(n * n, np.uint64)
        for i in range(n):
            for j in range(n):
                acc = np.uint64(hP[i * n + j])
                for k in range(K):
                    v = TA[i, k] + TB[j, k]
                    if neg: v = -v
                    v = (v + shift) % m
                    if v < 0: v += m
                    acc += np.uint64(v) * np.uint64(w[k])
                out[i * n + j] = acc
        return out

    @_nb.njit(cache=True)
    def _hash_join(LK, RK, bits):
        size = 1 << bits; mask = np.uint64(size - 1)
        head = -np.ones(size, np.int64); nxt = -np.ones(RK.shape[0], np.int64)
        for j in range(RK.shape[0]):
            hsh = (RK[j] * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(64 - bits)
            nxt[j] = head[hsh]; head[hsh] = j
        cnt = 0
        for i in range(LK.shape[0]):
            hsh = (LK[i] * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(64 - bits)
            j = head[hsh]
            while j >= 0:
                if RK[j] == LK[i]: cnt += 1
                j = nxt[j]
        I = np.empty(cnt, np.int64); J = np.empty(cnt, np.int64); t = 0
        for i in range(LK.shape[0]):
            hsh = (LK[i] * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(64 - bits)
            j = head[hsh]
            while j >= 0:
                if RK[j] == LK[i]:
                    I[t] = i; J[t] = j; t += 1
                j = nxt[j]
        return I, J
    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False


def admissible_fast(ctx, Kc, Q, c, U=None, nkey=12):
    """same result as admissible(), with the key computation and the join compiled (numba)"""
    if not HAVE_NUMBA: return admissible(ctx, Kc, Q, c, U, nkey)
    if U is None: U = frame_tables(ctx, Kc)
    m = 8 * c
    T = np.einsum('aj,afk->jfk', Q, U)
    sel = np.arange(nkey)
    if not hasattr(ctx, 'wu'):
        ctx.wu = ctx.w[32:32 + 32].astype(np.uint64)
        ctx.hPLu = ctx.hPL.astype(np.uint64); ctx.hPRu = ctx.hPR.astype(np.uint64)
    T0, T1, T2, T3 = (np.ascontiguousarray(T[J][:, sel]) for J in range(4))
    LK = _pair_keys(T0, T1, ctx.hPLu, ctx.wu, m, 0, False)
    RK = _pair_keys(T2, T3, ctx.hPRu, ctx.wu, m, 4 * c, True)
    bits = max(10, int(np.ceil(np.log2(len(RK)))) + 1)
    I, J = _hash_join(LK, RK, bits)
    if len(I) == 0: return np.zeros((0, 4), np.int64)
    tot = T[0][ctx.pi[I]] + T[1][ctx.pj[I]] + T[2][ctx.pi[J]] + T[3][ctx.pj[J]]
    ok = (np.abs(tot) == 4 * c).all(1) & ((ctx.PP[I] + ctx.PP[J]) == -4).all(1)
    return np.stack([ctx.pi[I[ok]], ctx.pj[I[ok]], ctx.pi[J[ok]], ctx.pj[J[ok]]], 1)


def k_frames(ctx, F, limit=10):
    """4 mutually orthogonal admissible row-blocks (up to shift/negation classes)"""
    cls = {}
    for g in F:
        k = canon_block(ctx, g)
        if k not in cls: cls[k] = g
    R = np.array(list(cls.values())) if cls else np.zeros((0, 4), np.int64)
    n = len(R); adj = np.zeros((n, n), bool)
    blocks = [ctx.S[g] for g in R]
    for i in range(n - 1):
        X = cross(ctx, blocks[i], R[i + 1:]); adj[i, i + 1:] = (X == 0).all(1)
    adj |= adj.T
    out = []
    for i in range(n):
        for j in np.nonzero(adj[i, i + 1:])[0] + i + 1:
            cand = np.nonzero(adj[i] & adj[j])[0]; cand = cand[cand > j]
            for k in cand:
                for l in cand[cand > k]:
                    if adj[k, l]:
                        out.append(np.stack([blocks[i], blocks[j], blocks[k], blocks[l]]))
                        if len(out) >= limit: return out, n
    return out, n


# ---------------------------------------------------------------- assembly and exact check (Singer coordinates)
def assemble(ctx, Kblocks, Kcblocks, Q, c):
    """K, Kc (4b x 4c, columns (J, x) with x = inf first), Cb = Q (x) I_c, forced corner C; returns (ok, N)"""
    b = ctx.b
    def frame(blocks):
        M = np.zeros((4 * b, 4 * c), np.int64)
        for I in range(4):
            for i in range(b):
                for J in range(4):
                    f = blocks[I][J]
                    M[I * b + i, J * c] = -f.sum()
                    for t in range(b): M[I * b + i, J * c + 1 + t] = f[ctx.add[t, ctx.neg[i]]]   # f(t - i)
        return M
    K = frame(Kblocks); Kc = frame(Kcblocks)
    Cb = np.kron(Q, np.eye(c, dtype=np.int64))
    num = -(K @ Cb.T @ Kc.T)
    ok = bool((K @ K.T == 4 * c * np.eye(4 * b, dtype=np.int64)).all() and (Kc @ Kc.T == 4 * c * np.eye(4 * b, dtype=np.int64)).all())
    ok &= all((K[:, J * c:(J + 1) * c].sum(1) == 0).all() and (Kc[:, J * c:(J + 1) * c].sum(1) == 0).all() for J in range(4))
    ok &= bool((num % (4 * c) == 0).all())
    C = num // (4 * c)
    ok &= bool((np.abs(C) == 1).all() and (C @ C.T == 4 * b * np.eye(4 * b, dtype=np.int64)).all())
    N = np.block([[C, K], [Kc.T, Cb]])
    return ok, N


# ---------------------------------------------------------------- Q matrices
def q_matrices(a):
    """all 4 x 4 integer matrices with rows signed permutations of a and Q Q^T = (sum a^2) I, up to right signed
    permutations (the K side is searched in full), each given by a representative"""
    n2 = sum(x * x for x in a)
    rows = sorted({tuple(s * v for s, v in zip(sg, p)) for p in itertools.permutations(a) for sg in itertools.product((1, -1), repeat=4)})
    rows = [np.array(r) for r in rows]
    out = []
    def rec(cur):
        if len(cur) == 4:
            out.append(np.array(cur)); return
        for r in rows:
            if all(int(r @ x) == 0 for x in cur): rec(cur + [r])
    rec([])
    # reduce by right signed permutations: canonical form = min over column signed perms
    colops = [(p, sg) for p in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)]
    seen = {}; reps = []
    for Q in out:
        key = min(tuple((Q[:, list(p)] * np.array(sg)[None, :]).ravel().tolist()) for p, sg in colops)
        if key not in seen:
            seen[key] = Q; reps.append(Q)
    return reps


def main():
    a_ = sys.argv[1:]
    opt = lambda k, d: a_[a_.index(k) + 1] if k in a_ else d
    c = int(opt('--c', '10')); kind = opt('--group', 'cyc'); b = c - 1
    a = [int(x) for x in opt('--a', '4,4,2,0' if c == 10 else '2,2,2,4').split(',')]
    maxfirst = int(opt('--max-first', '1000000000'))
    part, nparts = (int(x) for x in opt('--part', '0/1').split('/'))
    outdir = os.path.join(ROOT, opt('--out', 'logs/ab/singer')); os.makedirs(outdir, exist_ok=True)
    tag = f"c{c}_{kind}_a{''.join(str(abs(x)) for x in a)}"
    logf = os.path.join(outdir, f'{tag}.jsonl' if nparts == 1 else f'{tag}_p{part}of{nparts}.jsonl')
    done = set()
    for fn in os.listdir(outdir):
        if not (fn.startswith(tag) and fn.endswith('.jsonl')): continue
        for l in open(os.path.join(outdir, fn)):
            try:
                r = json.loads(l)
                if r.get('stage') == 'first': done.add(tuple(r['first']))
            except Exception: pass
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    t0 = time.time()
    ctx = Ctx(kind, b)
    Qs = q_matrices(a)
    firsts = first_blocks(ctx.S, ctx.P, ctx.add, ctx.auts)
    log(dict(stage='setup', c=c, b=b, group=kind, a=a, sequences=len(ctx.S), Q_reps=len(Qs), first_blocks=len(firsts),
             secs=round(time.time() - t0, 1)))
    # Q's paired with Kc frames: Kc columns were normalized (sorted, sum +1), so Q runs over left signed permutations
    rowops = [(p, sg) for p in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)]
    Qall = []; seenQ = set()
    colops = [(p, sg) for p in itertools.permutations(range(4)) for sg in itertools.product((1, -1), repeat=4)]
    for Q0 in Qs:
        for p, sg in rowops:
            Q = (np.array(sg)[:, None] * Q0[list(p), :])
            key = min(tuple((Q[:, list(pp)] * np.array(ss)[None, :]).ravel().tolist()) for pp, ss in colops)
            if key in seenQ: continue
            seenQ.add(key); Qall.append(Q)
    log(dict(stage='Q', left_classes=len(Qall)))
    nfr = 0; ntpl = 0
    for fi, first in enumerate(firsts[:maxfirst]):
        if fi % nparts != part or tuple(first) in done: continue
        t1 = time.time()
        h = ctx.S[list(first)]
        G = partners(ctx, h)
        frs, npc = frames_from(ctx, h, G)
        rec = dict(stage='first', first=list(first), partners=int(len(G)), partner_classes=int(npc), frames=len(frs))
        found = []; nskip = 0
        for fr, gidx in frs:
            if any(canon_first(ctx.rep, ctx.acts, g) < tuple(first) for g in gidx):
                nskip += 1; continue                    # processed from its minimal row-block
            nfr += 1
            U = frame_tables(ctx, fr)
            for qi, Q in enumerate(Qall):
                if not corner_parity_ok(fr, Q): continue
                F = admissible_fast(ctx, fr, Q, c, U)
                if len(F) < 4: continue
                kf, nk = k_frames(ctx, F, limit=2)
                for kb in kf:
                    ok, N = assemble(ctx, list(kb), list(fr), Q, c)
                    if ok:
                        ntpl += 1
                        f = os.path.join(outdir, f'template_{tag}_f{fi}_{ntpl}.npz')
                        np.savez_compressed(f, N=N.astype(np.int8), Q=Q, Kc_blocks=fr, K_blocks=kb)
                        found.append(dict(Q=qi, admissible=int(len(F)), saved=os.path.relpath(f, ROOT)))
        rec.update(frames_processed=len(frs) - nskip, templates=found, secs=round(time.time() - t1, 1))
        log(rec)
    log(dict(stage='end', frames=nfr, templates=ntpl, secs=round(time.time() - t0, 1)))


if __name__ == "__main__":
    main()
