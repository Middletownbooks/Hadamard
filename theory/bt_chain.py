r"""Exact compression chain for bordered-GS structures (theory/bt_engine.BStruct) [M].

The method of Kumar's order-3964 search (search/NOTES.md, METHOD.md "Bordered structures"; Apache-2.0), written anew:
- level lists: for a quotient Gbar = G/N, every S_i-invariant compressed vector of block i (reduced coordinates =
  S_i-orbits of Gbar) with its multiplicity (number of orbit sign vectors), filtered by the row sums allowed and by
  the per-character bound |y^(psi)|^2 <= sum_u F_N(u) psi(u) (every block's energy at psi is at most the total);
- refinement: the vectors of a finer level inside the fibre of fixed coarser vectors, by a dynamic program over the
  cells of a base level (orbits grouped by their base-level coordinate), with range pruning on the fixed coordinates;
- joins: exact 4-list joins of the compressed PAF vectors (64-bit linear hash of the pair sums, sort + binary
  search, every key match re-checked exactly);
- expansion: all orbit sign vectors in a fibre (DP-guided backtracking), for the final join on every class.
"""
import math, time, itertools
import numpy as np
from numba import njit
import bt_engine as E


# ------------------------------------------------------------------------------------------------ characters
def char_matrix(qt):
    """Complex character table of Gbar = G/N: characters of G trivial on N, evaluated on the coset reps.
    Returns (Ch (nchar x nb), targets (nchar,)) with targets[psi] = sum_u F_N(u) psi(u) (real)."""
    G = qt.st.G
    m = G.marr.astype(float)
    # xi in G with <xi, n> integral for all n in N
    Nco = G.co[qt.N]
    phase = (G.co[:, None, :] * Nco[None, :, :] / m).sum(-1)          # (n, |N|)
    ok = np.all(np.abs(phase - np.round(phase)) < 1e-9, axis=1)
    xis = np.flatnonzero(ok)
    assert len(xis) == qt.nb, (len(xis), qt.nb)
    ph = (G.co[xis][:, None, :] * G.co[qt.rep][None, :, :] / m).sum(-1)
    Ch = np.exp(2j * np.pi * ph)
    tg = (Ch @ qt.FN).real
    return Ch, tg


# ------------------------------------------------------------------------------------------------ level lists
def coord_values(qt, i):
    """Per reduced coordinate U of block i: dict value -> number of orbit sign assignments (subset sums of +-w_O
    over the G-orbits O mapping to U)."""
    st = qt.st
    nbo = len(qt.bsize[i])
    out = []
    for U in range(nbo):
        ws = qt.w[i][qt.omap[i] == U]
        d = {0: 1}
        for w in ws:
            nd = {}
            for s, c in d.items():
                nd[s + w] = nd.get(s + w, 0) + c
                nd[s - w] = nd.get(s - w, 0) + c
            d = nd
        out.append(d)
    return out


def level_list(qt, i, rows=None, chars=None, maxlen=5 * 10 ** 7, cvals=None):
    """All reduced compressed vectors of block i at quotient qt with row sum in `rows` (None: any) and every
    per-character energy at most the total target.  Returns (Y int64 (L x nbo), cnt float64 (L,)).
    Coordinates are combined largest-range first with row-sum range pruning."""
    cv = coord_values(qt, i) if cvals is None else cvals
    sz = qt.bsize[i]
    nbo = len(sz)
    vals = [np.array(sorted(d), np.int64) for d in cv]
    cnts = [np.array([d[v] for v in sorted(d)], np.float64) for d in cv]
    order = np.argsort([-(sz[U] * (vals[U].max() - vals[U].min())) for U in range(nbo)], kind='stable')
    # remaining row-sum range after position k
    lo = np.array([sz[U] * vals[U].min() for U in order]); hi = np.array([sz[U] * vals[U].max() for U in order])
    rlo = np.concatenate([np.cumsum(lo[::-1])[::-1], [0]]); rhi = np.concatenate([np.cumsum(hi[::-1])[::-1], [0]])
    rows_arr = None if rows is None else np.array(sorted(set(int(r) for r in rows)), np.int64)
    Y = np.zeros((1, 0), np.int64); C = np.ones(1); R = np.zeros(1, np.int64)
    for k, U in enumerate(order):
        v, c = vals[U], cnts[U]
        Rn = (R[:, None] + sz[U] * v[None, :]).ravel()
        keep = np.ones(len(Rn), bool)
        if rows_arr is not None:
            # some allowed row in [Rn + rlo[k+1], Rn + rhi[k+1]]
            a = np.searchsorted(rows_arr, Rn + rlo[k + 1], 'left')
            b = np.searchsorted(rows_arr, Rn + rhi[k + 1], 'right')
            keep = b > a
        kept = np.flatnonzero(keep)
        if len(kept) > maxlen:
            raise MemoryError(f"level list > {maxlen}")
        ri, vi = kept // len(v), kept % len(v)
        Y = np.concatenate([Y[ri], v[vi][:, None]], 1)
        C = C[ri] * c[vi]; R = Rn[kept]
        if len(Y) > maxlen:
            raise MemoryError(f"level list > {maxlen}")
    inv = np.argsort(order)
    Y = Y[:, inv]
    if rows_arr is not None:
        ok = np.isin(R, rows_arr); Y, C = Y[ok], C[ok]
    if chars is not None:
        Ch, tg = chars
        keep = np.ones(len(Y), bool)
        for s in range(0, len(Y), 200000):
            yf = Y[s:s + 200000][:, qt.borb[i]].astype(float)
            en = np.abs(yf @ Ch.T) ** 2
            keep[s:s + 200000] = np.all(en <= tg[None, :] + 1e-6, axis=1)
        Y, C = Y[keep], C[keep]
    return Y, C


def nested_map(qc, qf, i):
    """For N_f <= N_c: (cell of every fine reduced coordinate of block i, weight |U_f| / |U_c|)."""
    # fine coordinate U_f -> coarse coordinate of the image of its representative point
    st = qf.st
    rep_f = qf.breps[i]                                  # Gbar_f points representing the fine coordinates
    g = qf.rep[rep_f]                                    # G elements
    cell = qc.borb[i][qc.pi[g]]
    wt = qf.bsize[i] // qc.bsize[i][cell]
    assert (wt * qc.bsize[i][cell] == qf.bsize[i]).all()
    return cell, wt


def refine_list(qc, yc, qf, i, chars=None, cvals=None, maxlen=5 * 10 ** 7):
    """All reduced vectors of block i at the finer level qf (N_f <= N_c) compressing to yc at qc, with
    multiplicities, filtered by the per-character bounds of qf.  Product over the coarse cells of the fine
    combinations meeting each cell value."""
    cell, wt = nested_map(qc, qf, i)
    cv = coord_values(qf, i) if cvals is None else cvals
    nbo = len(qf.bsize[i])
    parts = []                                            # per coarse cell: (fine coord indices, Y, cnt)
    for c in range(len(qc.bsize[i])):
        idx = np.flatnonzero(cell == c)
        Y = np.zeros((1, 0), np.int64); C = np.ones(1); S = np.zeros(1, np.int64)
        # remaining range
        lo = np.array([wt[u] * min(cv[u]) for u in idx]); hi = np.array([wt[u] * max(cv[u]) for u in idx])
        rlo = np.concatenate([np.cumsum(lo[::-1])[::-1], [0]]); rhi = np.concatenate([np.cumsum(hi[::-1])[::-1], [0]])
        tgt = int(yc[c])
        for k, u in enumerate(idx):
            v = np.array(sorted(cv[u]), np.int64); cc = np.array([cv[u][x] for x in v], float)
            Sn = (S[:, None] + wt[u] * v[None, :]).ravel()
            keep = (Sn + rlo[k + 1] <= tgt) & (Sn + rhi[k + 1] >= tgt)
            kept = np.flatnonzero(keep)
            if len(kept) > maxlen:
                raise MemoryError(f"refined cell list > {maxlen}")
            ri, vi = kept // len(v), kept % len(v)
            Y = np.concatenate([Y[ri], v[vi][:, None]], 1)
            C = C[ri] * cc[vi]; S = Sn[kept]
        ok = S == tgt
        parts.append((idx, Y[ok], C[ok]))
        if ok.sum() == 0:
            return np.zeros((0, nbo), np.int64), np.zeros(0)
    total = np.prod([float(len(p[1])) for p in parts])
    if total > maxlen:
        raise MemoryError(f"refined list {total:.3g} > {maxlen}")
    Y = np.zeros((1, nbo), np.int64); C = np.ones(1)
    for idx, Yp, Cp in parts:
        Yn = np.repeat(Y, len(Yp), 0)
        Yn[:, idx] = np.tile(Yp, (len(Y), 1))
        C = (C[:, None] * Cp[None, :]).ravel(); Y = Yn
    if chars is not None:
        Ch, tg = chars
        keep = np.ones(len(Y), bool)
        for s in range(0, len(Y), 200000):
            yf = Y[s:s + 200000][:, qf.borb[i]].astype(float)
            en = np.abs(yf @ Ch.T) ** 2
            keep[s:s + 200000] = np.all(en <= tg[None, :] + 1e-6, axis=1)
        Y, C = Y[keep], C[keep]
    return Y, C


def paf_reduced(qt, i, Y):
    """class PAF vectors (L x ncls) of reduced vectors Y of block i"""
    out = np.zeros((len(Y), qt.ncls), np.int64)
    for s in range(0, len(Y), 100000):
        yf = Y[s:s + 100000][:, qt.borb[i]]
        out[s:s + 100000] = qt.paf_classes(yf)
    return out


# ------------------------------------------------------------------------------------------------ 4-list join
@njit(cache=True)
def _hash(P, R):
    """linear 64-bit hash (int64 wrap-around arithmetic): h(P + P') = h(P) + h(P') mod 2^64"""
    n, d = P.shape
    h = np.zeros(n, np.int64)
    for a in range(n):
        s = np.int64(0)
        for j in range(d):
            s += P[a, j] * R[j]
        h[a] = s
    return h


@njit(cache=True)
def _join_core(h1, h2, h3, h4, ht, P1, P2, P3, P4, T, maxout, chunk):
    """all (a,b,c,d) with P1[a]+P2[b]+P3[c]+P4[d] == T (exact), via sorted pair-sum keys of (1,2) in chunks"""
    n1, n2, n3, n4 = len(h1), len(h2), len(h3), len(h4)
    d = P1.shape[1]
    out = np.empty((maxout, 4), np.int64); no = 0
    rows_per = max(1, chunk // max(n2, 1))
    a0 = 0
    while a0 < n1:
        a1 = min(n1, a0 + rows_per)
        m = (a1 - a0) * n2
        keys = np.empty(m, np.int64); ia = np.empty(m, np.int64); ib = np.empty(m, np.int64)
        t = 0
        for a in range(a0, a1):
            for b in range(n2):
                keys[t] = h1[a] + h2[b]; ia[t] = a; ib[t] = b; t += 1
        o = np.argsort(keys)
        keys = keys[o]; ia = ia[o]; ib = ib[o]
        for c in range(n3):
            for e in range(n4):
                k = ht - h3[c] - h4[e]
                lo = np.searchsorted(keys, k)
                while lo < m and keys[lo] == k:
                    a = ia[lo]; b = ib[lo]
                    ok = True
                    for j in range(d):
                        if P1[a, j] + P2[b, j] + P3[c, j] + P4[e, j] != T[j]:
                            ok = False; break
                    if ok:
                        if no < maxout:
                            out[no, 0] = a; out[no, 1] = b; out[no, 2] = c; out[no, 3] = e
                        no += 1
                    lo += 1
        a0 = a1
    return out[:min(no, maxout)], no


_RNG = np.random.default_rng(12345)
_R = _RNG.integers(1, 2 ** 62, size=4096, dtype=np.int64) * 2 + 1


def _wrap64(x):
    x %= 2 ** 64
    return np.int64(x - 2 ** 64 if x >= 2 ** 63 else x)


def join4(P, T, maxout=10 ** 7, chunk=4 * 10 ** 7):
    """All index quadruples with P[0][a] + P[1][b] + P[2][c] + P[3][d] == T.  The two pair sides are chosen to
    balance the work.  Returns (array (k x 4) in the ORIGINAL list order, total count)."""
    P = [np.ascontiguousarray(p, dtype=np.int64) for p in P]
    T = np.asarray(T, np.int64)
    if min(len(p) for p in P) == 0:
        return np.zeros((0, 4), np.int64), 0
    d = len(T)
    R = _R[:d]
    sizes = [len(p) for p in P]
    best = None
    for perm in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
        a, b, c, e = perm
        s1, s2 = sizes[a] * sizes[b], sizes[c] * sizes[e]
        if s1 > s2:
            perm = (c, e, a, b); s1, s2 = s2, s1
        cost = s1 * math.log2(s1 + 2) + s2 * math.log2(s1 + 2)
        if best is None or cost < best[0]:
            best = (cost, perm)
    perm = best[1]
    hs = [_hash(p, R) for p in P]
    ht = _wrap64(sum(int(T[j]) * int(R[j]) for j in range(d)))
    out, no = _join_core(hs[perm[0]], hs[perm[1]], hs[perm[2]], hs[perm[3]], ht,
                         P[perm[0]], P[perm[1]], P[perm[2]], P[perm[3]], T, maxout, chunk)
    res = np.empty_like(out)
    for k in range(4):
        res[:, perm[k]] = out[:, k]
    return res, no


# ------------------------------------------------------------------------------------------------ fibres
class Cells:
    """Orbits of block i grouped by their coordinate at a base level (qb): cell U holds the orbits O with
    omap[O] = U.  valid(U, value) lists the sign assignments (as +-1 arrays) meeting the base value."""

    def __init__(self, qb, i):
        self.qb, self.i = qb, i
        self.members = [np.flatnonzero(qb.omap[i] == U) for U in range(len(qb.bsize[i]))]
        self.w = [qb.w[i][m] for m in self.members]
        self._cache = {}

    def assignments(self, U, value):
        key = (U, int(value))
        if key not in self._cache:
            m = self.members[U]; w = self.w[U]; k = len(m)
            if k > 26:
                raise MemoryError("cell too large")
            # meet in the middle on the cell's own sum
            bits = ((np.arange(2 ** k)[:, None] >> np.arange(k)[None, :]) & 1) if k <= 20 else None
            if bits is None:
                h = k // 2
                A = ((np.arange(2 ** h)[:, None] >> np.arange(h)[None, :]) & 1) * 2 - 1
                B = ((np.arange(2 ** (k - h))[:, None] >> np.arange(k - h)[None, :]) & 1) * 2 - 1
                sa = A @ w[:h]; sb = B @ w[h:]
                o = np.argsort(sb); sb_s = sb[o]
                res = []
                for ai in range(len(A)):
                    need = value - sa[ai]
                    lo = np.searchsorted(sb_s, need, 'left'); hi = np.searchsorted(sb_s, need, 'right')
                    for bj in o[lo:hi]:
                        res.append(np.concatenate([A[ai], B[bj]]))
                eps = np.array(res, np.int64).reshape(-1, k)
            else:
                eps = bits * 2 - 1
                eps = eps[eps @ w == value]
            self._cache[key] = eps
        return self._cache[key]


def fibre_dp(st, i, qb, yb, qs_fixed, ys_fixed, qnew, maxstates=5 * 10 ** 6):
    """Distinct reduced vectors of block i at level qnew inside the fibre {x : compression at qb = yb, and at each
    qs_fixed[j] = ys_fixed[j]}, with multiplicities.  DP over the cells of qb.
    Returns (Ynew (L x nbo_new) int64, counts float64)."""
    cells = Cells(qb, i)
    nbo_new = len(qnew.bsize[i])
    fixed_dims = [len(q.bsize[i]) for q in qs_fixed]
    D = sum(fixed_dims) + nbo_new
    # contribution matrix of every orbit to (fixed coords..., new coords)
    Cmat = np.concatenate([q.contrib(i) for q in qs_fixed] + [qnew.contrib(i)], axis=1)
    target_fixed = np.concatenate([np.asarray(y, np.int64) for y in ys_fixed]) if qs_fixed else np.zeros(0, np.int64)
    nf = len(target_fixed)
    # remaining range per fixed coordinate after each cell
    ncell = len(cells.members)
    absC = np.abs(Cmat[:, :nf])
    rem = np.zeros((ncell + 1, nf), np.int64)
    for U in range(ncell - 1, -1, -1):
        rem[U] = rem[U + 1] + absC[cells.members[U]].sum(0)
    S = np.zeros((1, D), np.int64); cnt = np.ones(1)
    vt = np.dtype((np.void, 8 * D))
    for U in range(ncell):
        eps = cells.assignments(U, yb[U])
        if len(eps) == 0:
            return np.zeros((0, nbo_new), np.int64), np.zeros(0)
        contrib = eps @ Cmat[cells.members[U]]                       # (a x D)
        cu, inv = np.unique(contrib, axis=0, return_inverse=True)
        cc = np.bincount(inv.ravel(), minlength=len(cu)).astype(float)
        Sn = (S[:, None, :] + cu[None, :, :]).reshape(-1, D)
        Cn = (cnt[:, None] * cc[None, :]).ravel()
        if nf:
            ok = np.all(np.abs(target_fixed[None, :] - Sn[:, :nf]) <= rem[U + 1][None, :], axis=1)
            Sn, Cn = Sn[ok], Cn[ok]
        key = np.ascontiguousarray(Sn).view(vt).ravel()
        u, idx, inv = np.unique(key, return_index=True, return_inverse=True)
        S = Sn[idx]; cnt = np.bincount(inv.ravel(), weights=Cn, minlength=len(u))
        if len(S) > maxstates:
            raise MemoryError(f"fibre DP states > {maxstates}")
    if nf:
        ok = np.all(S[:, :nf] == target_fixed[None, :], axis=1)
        S, cnt = S[ok], cnt[ok]
    return S[:, nf:], cnt


def _cartesian(per, cells_idx):
    """enumerate the product of the cells' assignment lists: (index matrix (n x len(cells_idx)), partial sums)"""
    nf = per[0][1].shape[1]
    I = np.zeros((1, 0), np.int64); S = np.zeros((1, nf), np.int64)
    for U in cells_idx:
        a = len(per[U][0])
        I = np.concatenate([np.repeat(I, a, 0), np.tile(np.arange(a), len(I))[:, None]], 1)
        S = (S[:, None, :] + per[U][1][None, :, :]).reshape(-1, nf)
    return I, S


def expand_fibre(st, i, qb, yb, qs_fixed, ys_fixed, maxn=2 * 10 ** 7, maxhalf=3 * 10 ** 7):
    """All orbit sign vectors (L x r_i, int8) of block i with compression yb at qb and ys_fixed[j] at qs_fixed[j]:
    per base cell the assignments meeting yb, then a meet in the middle over the cells on the fixed coordinates."""
    cells = Cells(qb, i)
    r = st.r[i]
    Cmat = np.concatenate([q.contrib(i) for q in qs_fixed], axis=1) if qs_fixed else np.zeros((r, 0), np.int64)
    target = np.concatenate([np.asarray(y, np.int64) for y in ys_fixed]) if qs_fixed else np.zeros(0, np.int64)
    nf = len(target)
    ncell = len(cells.members)
    per = []
    for U in range(ncell):
        eps = cells.assignments(U, yb[U])
        if len(eps) == 0:
            return np.zeros((0, r), np.int8)
        per.append((eps, eps @ Cmat[cells.members[U]] if nf else np.zeros((len(eps), 0), np.int64)))
    # balanced split of the cells by log size
    logs = np.array([math.log(len(p[0])) for p in per])
    order = np.argsort(-logs)
    L, R = [], []; sl = sr = 0.0
    for U in order:
        if sl <= sr:
            L.append(int(U)); sl += logs[U]
        else:
            R.append(int(U)); sr += logs[U]
    if math.exp(max(sl, sr)) > maxhalf:
        raise MemoryError("fibre half too large")
    IL, SL = _cartesian(per, L)
    IR, SR = _cartesian(per, R)
    Rh = _R[:max(nf, 1)]
    hR = SR @ Rh[:nf] if nf else np.zeros(len(SR), np.int64)
    hL = (target @ Rh[:nf] if nf else 0) - (SL @ Rh[:nf] if nf else np.zeros(len(SL), np.int64))
    o = np.argsort(hR, kind='stable'); hRs = hR[o]
    lo = np.searchsorted(hRs, hL, 'left'); hi = np.searchsorted(hRs, hL, 'right')
    cnt = hi - lo
    tot = int(cnt.sum())
    if tot > maxn:
        raise MemoryError(f"fibre > {maxn}")
    li = np.repeat(np.arange(len(SL)), cnt)
    ri = o[np.concatenate([np.arange(a, b) for a, b in zip(lo[cnt > 0], hi[cnt > 0])])] if tot else np.zeros(0, np.int64)
    if nf and tot:
        ok = np.all(SL[li] + SR[ri] == target[None, :], axis=1)
        li, ri = li[ok], ri[ok]
    out = np.empty((len(li), r), np.int8)
    for k, U in enumerate(L):
        out[:, cells.members[U]] = per[U][0][IL[li, k]]
    for k, U in enumerate(R):
        out[:, cells.members[U]] = per[U][0][IR[ri, k]]
    return out


def class_paf_orbits(st, i, EPS):
    """class PAF vectors (L x ncls) of orbit sign vectors EPS (L x r_i) of block i: PAF on G by FFT over the
    factor axes (exact after rounding), read at the class representatives."""
    G = st.G
    out = np.zeros((len(EPS), st.ncls), np.int64)
    reps = st.cls_rep
    B = max(1, 2 ** 22 // max(G.n, 1))
    for s in range(0, len(EPS), B):
        X = EPS[s:s + B].astype(np.float64)[:, st.orb[i]].reshape((-1,) + G.mods)
        axes = tuple(range(1, G.k + 1))
        F = np.fft.fftn(X, axes=axes)
        Pf = np.fft.ifftn(np.abs(F) ** 2, axes=axes).real.reshape(len(X), G.n)
        out[s:s + B] = np.rint(Pf[:, reps]).astype(np.int64)
    return out


def final_join(st, EPS, qC=None, maxout=10 ** 6, log=None, joiner=None):
    """Exact 4-list join of orbit sign vectors on every class of the structure, optionally split first by the
    compressions at qC (distinct level-C vectors joined on the level-C equations, then the full join inside each
    compatible quadruple).  Returns list of index quadruples into EPS."""
    joiner = joiner or join4b
    P = [class_paf_orbits(st, i, EPS[i]) for i in range(4)]
    if qC is None:
        res, no = joiner(P, st.cls_target, maxout=maxout)
        return [tuple(x) for x in res]
    groups, keysP = [], []
    for i in range(4):
        yC = EPS[i].astype(np.int64) @ qC.contrib(i)
        u, inv = np.unique(yC, axis=0, return_inverse=True)
        inv = inv.ravel()
        order = np.argsort(inv, kind='stable')
        bounds = np.searchsorted(inv[order], np.arange(len(u) + 1))
        groups.append((order, bounds))
        keysP.append(paf_reduced(qC, i, u))
    kq, nk = join4(keysP, qC.cls_target, maxout=10 ** 7)
    out = []
    for q in kq:
        idx = [groups[i][0][groups[i][1][q[i]]:groups[i][1][q[i] + 1]] for i in range(4)]
        res, no = joiner([P[i][idx[i]] for i in range(4)], st.cls_target, maxout=maxout)
        for x in res:
            out.append(tuple(int(idx[i][x[i]]) for i in range(4)))
    if log is not None:
        log['level_C_quadruples'] = int(nk)
    return out


# ------------------------------------------------------------------------------------------------ bucketed join
@njit(cache=True)
def _groups(vals):
    """sort order of vals and the boundaries of equal-value groups: (order, uniq values, starts (len+1))"""
    o = np.argsort(vals, kind='mergesort')
    sv = vals[o]
    n = len(sv)
    nu = 0
    for k in range(n):
        if k == 0 or sv[k] != sv[k - 1]:
            nu += 1
    u = np.empty(nu, np.int64); st = np.empty(nu + 1, np.int64)
    t = 0
    for k in range(n):
        if k == 0 or sv[k] != sv[k - 1]:
            u[t] = sv[k]; st[t] = k; t += 1
    st[nu] = n
    return o, u, st


@njit(cache=True)
def _bucket_join(P1, P2, P3, P4, h1, h2, h3, h4, ht, T, j, maxout, maxbuf):
    o1, u1, s1 = _groups(P1[:, j]); o2, u2, s2 = _groups(P2[:, j])
    o3, u3, s3 = _groups(P3[:, j]); o4, u4, s4 = _groups(P4[:, j])
    d = P1.shape[1]
    # value -> group index lookups for lists 2 and 4
    lo2, hi2 = u2[0], u2[-1]; lo4, hi4 = u4[0], u4[-1]
    map2 = -np.ones(hi2 - lo2 + 1, np.int64); map4 = -np.ones(hi4 - lo4 + 1, np.int64)
    for k in range(len(u2)):
        map2[u2[k] - lo2] = k
    for k in range(len(u4)):
        map4[u4[k] - lo4] = k
    out = np.empty((maxout, 4), np.int64); no = 0
    smin = u1[0] + u2[0]; smax = u1[-1] + u2[-1]
    keys = np.empty(maxbuf, np.int64); ia = np.empty(maxbuf, np.int64); ib = np.empty(maxbuf, np.int64)
    overflow = 0
    for s in range(smin, smax + 1):
        sd = T[j] - s
        # any (3,4) pair with P3 + P4 = sd ?
        n34 = 0
        for g3 in range(len(u3)):
            v4 = sd - u3[g3]
            if v4 < lo4 or v4 > hi4: continue
            g4 = map4[v4 - lo4]
            if g4 < 0: continue
            n34 += (s3[g3 + 1] - s3[g3]) * (s4[g4 + 1] - s4[g4])
        if n34 == 0: continue
        m = 0
        for g1 in range(len(u1)):
            v2 = s - u1[g1]
            if v2 < lo2 or v2 > hi2: continue
            g2 = map2[v2 - lo2]
            if g2 < 0: continue
            for x in range(s1[g1], s1[g1 + 1]):
                a = o1[x]
                for y in range(s2[g2], s2[g2 + 1]):
                    b = o2[y]
                    if m < maxbuf:
                        keys[m] = h1[a] + h2[b]; ia[m] = a; ib[m] = b
                    m += 1
        if m == 0: continue
        if m > maxbuf:
            overflow += 1
            m = maxbuf
        oo = np.argsort(keys[:m])
        ks = keys[:m][oo]; ias = ia[:m][oo]; ibs = ib[:m][oo]
        for g3 in range(len(u3)):
            v4 = sd - u3[g3]
            if v4 < lo4 or v4 > hi4: continue
            g4 = map4[v4 - lo4]
            if g4 < 0: continue
            for x in range(s3[g3], s3[g3 + 1]):
                c = o3[x]
                for y in range(s4[g4], s4[g4 + 1]):
                    e = o4[y]
                    k = ht - h3[c] - h4[e]
                    lo = np.searchsorted(ks, k)
                    while lo < m and ks[lo] == k:
                        a = ias[lo]; b = ibs[lo]
                        ok = True
                        for jj in range(d):
                            if P1[a, jj] + P2[b, jj] + P3[c, jj] + P4[e, jj] != T[jj]:
                                ok = False; break
                        if ok:
                            if no < maxout:
                                out[no, 0] = a; out[no, 1] = b; out[no, 2] = c; out[no, 3] = e
                            no += 1
                        lo += 1
    return out[:min(no, maxout)], no, overflow


def join4b(P, T, maxout=10 ** 7, maxbuf=6 * 10 ** 7):
    """Exact 4-list join, pairs bucketed by one class value (the most spread) and hashed on the rest.
    Raises if a bucket overflows maxbuf (then nothing is lost silently)."""
    P = [np.ascontiguousarray(p, dtype=np.int64) for p in P]
    T = np.asarray(T, np.int64)
    if min(len(p) for p in P) == 0:
        return np.zeros((0, 4), np.int64), 0
    d = len(T)
    sizes = [len(p) for p in P]
    best = None
    for perm in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
        a, b, c, e = perm
        s1, s2 = sizes[a] * sizes[b], sizes[c] * sizes[e]
        if s1 > s2:
            perm = (c, e, a, b)
        cost = max(s1, s2) + 0.3 * min(s1, s2)
        if best is None or cost < best[0]:
            best = (cost, perm)
    perm = best[1]
    Q = [P[k] for k in perm]
    # bucket coordinate: largest spread of the stored pair sums
    spread = [np.ptp(Q[0][:, j]) + np.ptp(Q[1][:, j]) for j in range(d)]
    j = int(np.argmax(spread))
    R = _R[:d]
    hs = [_hash(p, R) for p in Q]
    ht = _wrap64(sum(int(T[k]) * int(R[k]) for k in range(d)))
    out, no, ovf = _bucket_join(Q[0], Q[1], Q[2], Q[3], hs[0], hs[1], hs[2], hs[3], ht, T, j, maxout, maxbuf)
    if ovf:
        raise MemoryError(f"join bucket overflow ({ovf} buckets > {maxbuf})")
    res = np.empty_like(out)
    for k in range(4):
        res[:, perm[k]] = out[:, k]
    return res, no


def char_orbit_matrix(st, i, maxchars=None):
    """Phi (r_i x nc) with Phi[O, psi] = sum_{g in O} psi(g) for one character psi per dual orbit of <S_i, conj>
    (|x^(psi)| is constant on these orbits for S_i-invariant x), and the per-character targets (4t or 4n)."""
    G = st.G
    m = G.marr.astype(float)
    # dual action: psi_xi o s = psi_{s^T xi}; orbits of xi under S_i^T and -1
    perms = []
    for M in st.bgens[i]:
        M = np.array(M, np.int64)
        # psi_xi o M = psi_xi' with xi'_k = sum_j xi_j M_jk m_k / m_j  (integral since M is well defined)
        W = (M * G.marr[None, :]) // G.marr[:, None]
        perms.append(G.idx(G.co @ W))
    perms.append(G.neg())
    o, reps, _ = E.perm_closure_orbits(G.n, perms)
    ph = (G.co[reps][:, None, :] * G.co[None, :, :] / m).sum(-1)
    Ch = np.exp(2j * np.pi * ph)                                       # (nc x n)
    Phi = np.zeros((st.r[i], len(reps)), complex)
    for O in range(st.r[i]):
        Phi[O] = Ch[:, st.orb[i] == O].sum(1)
    tg = (Ch @ st.F).real
    return Phi, tg


def char_filter(st, i, EPS, Phi_tg):
    Phi, tg = Phi_tg
    keep = np.ones(len(EPS), bool)
    for s in range(0, len(EPS), 100000):
        en = np.abs(EPS[s:s + 100000].astype(float) @ Phi) ** 2
        keep[s:s + 100000] = np.all(en <= tg[None, :] + 1e-6, axis=1)
    return keep


def q_level_list(qQ, i, row, chars, cvals=None):
    """Q-level list of block i for row sum `row`: z(u) - row/q is bounded pointwise by
    D = sqrt(sum_{chi != 1} bound(chi) / q) (Parseval on Q), so only coordinate values within D of row/q are
    combined; then the exact row and character filters.  Exact (no vector meeting the bounds is lost)."""
    Ch, tg = chars
    q = qQ.nb
    # bound(chi) for chi != 1: the per-block bound tg; trivial character is the row sum itself
    triv = int(np.argmax(np.all(np.abs(Ch - 1) < 1e-9, axis=1)))
    D = math.sqrt(max(0.0, sum(tg[k] for k in range(len(tg)) if k != triv)) / q)
    cv = coord_values(qQ, i) if cvals is None else cvals
    sz = qQ.bsize[i]
    mean = row / q
    vals, cnts = [], []
    for U in range(len(sz)):
        v = np.array([x for x in sorted(cv[U]) if abs(x - mean) <= D + 1e-9], np.int64)
        if len(v) == 0:
            return np.zeros((0, len(sz)), np.int64), np.zeros(0)
        vals.append(v); cnts.append(np.array([cv[U][x] for x in v], np.float64))
    Y = np.zeros((1, 0), np.int64); C = np.ones(1)
    for U in range(len(sz)):
        Y = np.concatenate([np.repeat(Y, len(vals[U]), 0), np.tile(vals[U], len(Y))[:, None]], 1)
        C = (C[:, None] * cnts[U][None, :]).ravel()
    ok = (Y * sz[None, :]).sum(1) == row
    Y, C = Y[ok], C[ok]
    yf = Y[:, qQ.borb[i]].astype(float)
    en = np.abs(yf @ Ch.T) ** 2
    keep = np.all(en <= tg[None, :] + 1e-6, axis=1)
    return Y[keep], C[keep]


def q_level_list_two_step(st, qQ, i, row, chars, cvals=None, maxlen=2 * 10 ** 6):
    """Q-level list when some nontrivial characters of Q have the loose bound 4n (product families): the loose
    characters are those trivial on Q0 = the intersection of their kernels; z = (Q/Q0-compression w) + deviation,
    and the deviation energy sum_u (z(u) - w(c(u))/|Q0|)^2 <= (1/q) sum_{tight chi} bound(chi) (Parseval) [exact].
    Enumerates w (coarse level, row and loose-character filters), then per Q0-coset the fine combinations within
    the deviation budget, then the exact character filter."""
    Ch, tg = chars
    q = qQ.nb
    n4 = 4 * st.nH
    triv = int(np.argmax(np.all(np.abs(Ch - 1) < 1e-9, axis=1)))
    loose = [k for k in range(len(tg)) if k != triv and tg[k] > 4 * st.t + 1e-6]
    tight = [k for k in range(len(tg)) if k != triv and k not in loose]
    if not loose:
        return q_level_list(qQ, i, row, chars, cvals)
    # Q0 = common kernel of the loose characters, as a set of cosets (indices into Q = cosets of K)
    ker = np.all(np.abs(Ch[loose] - 1) < 1e-9, axis=0)
    Q0 = np.flatnonzero(ker)
    # coarse cosets of Q0 in Q
    addQ = qQ.add
    cc = -np.ones(q, np.int64); nc = 0
    for u in range(q):
        if cc[u] < 0:
            for v in Q0:
                cc[addQ[u, v]] = nc
            nc += 1
    Dsq = sum(tg[k] for k in tight) / q
    cv = coord_values(qQ, i) if cvals is None else cvals
    sz = qQ.bsize[i]
    # group the fine coordinates (S-orbits U on Q) by the S-orbit of their coarse cosets ("coarse orbits"); every U
    # meets each coset of its coarse orbit in |U| / |O_c| points (S-invariance), so the coarse value per coset is
    # w = sum_U (|U| / |O_c|) z_U and the deviation energy of the group is sum_U |U| (z_U - w / m0)^2
    pts = [np.flatnonzero(qQ.borb[i] == U) for U in range(len(sz))]
    cos_of = [sorted(set(cc[pts[U]].tolist())) for U in range(len(sz))]
    # union-find coarse cosets that share a fine orbit
    parent = list(range(nc))
    def fd(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for cs in cos_of:
        for c2 in cs[1:]:
            ra, rb = fd(cs[0]), fd(c2)
            if ra != rb: parent[rb] = ra
    groups = {}
    for c in range(nc):
        groups.setdefault(fd(c), []).append(c)
    gl = list(groups.values())
    m0 = len(Q0)
    out_Y, out_C = [], []
    per = []                                                  # per coarse orbit: (fine coords, {row contribution: [...]})
    for cos_list in gl:
        oc = len(cos_list)
        idx = [U for U in range(len(sz)) if cos_of[U][0] in cos_list]
        vals_sorted = [np.array(sorted(cv[U])) for U in idx]
        res = {}
        def rec2(j, vals, cnt):
            if j == len(idx):
                tot = sum(sz[idx[k]] * vals[k] for k in range(len(idx)))   # row contribution of the group
                w = tot / oc
                dev = sum(sz[idx[k]] * (vals[k] - w / m0) ** 2 for k in range(len(idx)))
                if dev <= Dsq + 1e-9:
                    res.setdefault(int(tot), []).append((tuple(vals), cnt, dev))
                return
            lo_ = min(vals) if vals else None; hi_ = max(vals) if vals else None
            for v in vals_sorted[j]:
                if vals and (max(hi_, v) - min(lo_, v)) > 2 * Dsq ** 0.5 + 1e-9:
                    continue
                rec2(j + 1, vals + [int(v)], cnt * cv[idx[j]][int(v)])
        rec2(0, [], 1.0)
        per.append((idx, res))
    nc = len(per)
    # combine coarse cosets: row sum and total deviation <= Dsq (pruned merge, lists sorted by deviation)
    keys = [sorted(p[1].keys()) for p in per]
    for c in range(nc):
        for k in per[c][1]:
            per[c][1][k].sort(key=lambda t: t[2])
    mins = [{k: per[c][1][k][0][2] for k in per[c][1]} for c in range(nc)]
    def rec3(c, wsum, partial, cnt, dev):
        if c == nc:
            if wsum != row:
                return
            y = np.zeros(len(sz), np.int64)
            for cc_, combo in enumerate(partial):
                for k, U in enumerate(per[cc_][0]):
                    y[U] = combo[k]
            out_Y.append(y); out_C.append(cnt)
            if len(out_Y) > maxlen:
                raise MemoryError("two-step Q-level list too large")
            return
        rest_min = 0.0
        for k in keys[c]:
            if c == nc - 1 and wsum + k != row:
                continue
            base = mins[c][k]
            if dev + base > Dsq + 1e-9:
                continue
            for vals, ct, dv in per[c][1][k]:
                if dev + dv > Dsq + 1e-9:
                    break
                rec3(c + 1, wsum + k, partial + [vals], cnt * ct, dev + dv)
    rec3(0, 0, [], 1.0, 0.0)
    if not out_Y:
        return np.zeros((0, len(sz)), np.int64), np.zeros(0)
    Y = np.array(out_Y); Cn = np.array(out_C)
    en = np.abs(Y[:, qQ.borb[i]].astype(float) @ Ch.T) ** 2
    keep = np.all(en <= tg[None, :] + 1e-6, axis=1)
    return Y[keep], Cn[keep]
