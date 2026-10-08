"""Exhaustive solver for symmetric 'slices' of periodic complementary problems over a finite abelian group G.

Problem types (target T for sum_k PAF_k(g), g != 0):
  'lp' : 2 sequences, sum PAF = -2   (Legendre pair; H(2|G|+2) via two-circulant core)  row sums +-1
  'b1' : 4 sequences, sum PAF = -4   (bordered Goethals-Seidel; H(4(|G|+1)))            row sums (+-2,0,0,0) or (+-1,..)
  'b0' : 4 sequences, sum PAF =  0   (Goethals-Seidel; H(4|G|))                          sum a_k^2 = 4|G|
All sequences are invariant under the symmetry group S = <auts> <= Aut(G) (one sign per S-orbit).
Method: Gray-code enumeration of sign vectors per orbit (numba), filter by row sum and the PSD bound
|hat X(xi)|^2 <= (total PSD), exact integer PAF on lag-class representatives from the PSD, then an exact
hash join (2-list or 4-list).  Complete for the slice: if it reports 0, no S-invariant solution exists.
"""
import time

import numba as nb
import numpy as np

from abgroup import slice_structure


def prepare(G, auts):
    orb, lag, frq = slice_structure(G, auts)
    r = int(orb.max()) + 1
    sizes = np.bincount(orb, minlength=r).astype(np.int64)
    assert sizes[orb[0]] == 1 and orb[0] == 0
    # frequency classes (exclude the class of the trivial character, index 0)
    fl = frq.copy()
    f0 = fl[0]
    fclasses = [c for c in np.unique(fl) if c != f0]
    frep = np.array([np.flatnonzero(fl == c)[0] for c in fclasses])
    fsize = np.array([np.sum(fl == c) for c in fclasses])
    # lag classes (exclude {0})
    l0 = lag[0]
    lclasses = [c for c in np.unique(lag) if c != l0]
    lrep = np.array([np.flatnonzero(lag == c)[0] for c in lclasses])
    lsize = np.array([np.sum(lag == c) for c in lclasses])
    mods = np.array(G.mods, dtype=float)
    # eta[i, c] = sum_{x in O_i} e(<x, xi_c>)
    ph = (G.coords[:, None, :] * G.coords[frep][None, :, :] / mods).sum(-1)  # n x Lf
    E = np.exp(2j * np.pi * ph)
    eta = np.zeros((r, len(frep)), dtype=np.complex128)
    for c in range(len(frep)):
        eta[:, c] = np.bincount(orb, weights=E[:, c].real, minlength=r) + 1j * np.bincount(orb, weights=E[:, c].imag, minlength=r)
    # W[c, j] = sum_{xi in freq class c} e(-<g_j, xi>) for lag reps g_j  (real by -1 symmetry)
    W = np.zeros((len(frep), len(lrep)))
    for c, cl in enumerate(fclasses):
        xs = np.flatnonzero(fl == cl)
        php = (G.coords[lrep][:, None, :] * G.coords[xs][None, :, :] / mods).sum(-1)
        W[c] = np.cos(2 * np.pi * php).sum(1)
    return dict(G=G, orb=orb, r=r, sizes=sizes, frep=frep, fsize=fsize, lrep=lrep, lsize=lsize, eta=eta, W=W,
                n=G.n)


@nb.njit(cache=True)
def _enum(eta, sizes, bound, allowed, maxout):
    r, L = eta.shape
    F = np.zeros(L, dtype=np.complex128)
    for c in range(L):
        for i in range(r):
            F[c] += eta[i, c]
    a = np.ones(r, dtype=np.int64)
    rs = 0
    for i in range(r):
        rs += sizes[i]
    masks = np.zeros(maxout, dtype=np.int64)
    psd = np.zeros((maxout, L))
    rsums = np.zeros(maxout, dtype=np.int64)
    nout = 0
    mask = 0
    total = 1 << (r - 1)
    order = np.arange(L)
    for k in range(total):
        if k > 0:
            j = 0
            kk = k
            while (kk & 1) == 0:
                kk >>= 1
                j += 1
            i = j + 1
            a[i] = -a[i]
            mask ^= (1 << i)
            s = 2 * a[i]
            for c in range(L):
                F[c] += s * eta[i, c]
            rs += s * sizes[i]
        okr = False
        for t in range(allowed.shape[0]):
            if rs == allowed[t]:
                okr = True
                break
        if not okr:
            continue
        ok = True
        for cc in range(L):
            c = order[cc]
            v = F[c].real * F[c].real + F[c].imag * F[c].imag
            if v > bound + 1e-6:
                ok = False
                break
        if not ok:
            continue
        if nout < maxout:
            masks[nout] = mask
            rsums[nout] = rs
            for c in range(L):
                psd[nout, c] = F[c].real * F[c].real + F[c].imag * F[c].imag
        nout += 1
    return masks, rsums, psd, nout


def enumerate_blocks(P, bound, allowed, maxout=4_000_000):
    masks, rsums, psd, nout = _enum(P['eta'], P['sizes'], float(bound), np.array(sorted(allowed), dtype=np.int64),
                                    maxout)
    if nout > maxout:
        raise RuntimeError(f'too many survivors {nout} > {maxout}')
    return masks[:nout], rsums[:nout], psd[:nout]


def exact_paf(P, rsums, psd):
    """integer PAF at lag-class reps from PSD on frequency classes"""
    n = P['n']
    val = (rsums.astype(float)[:, None] ** 2 + psd @ P['W']) / n
    R = np.rint(val)
    err = np.abs(val - R).max() if len(val) else 0.0
    assert err < 1e-4, err
    return R.astype(np.int64)


def mask_to_seq(P, mask):
    r = P['r']
    a = np.array([-1 if (int(mask) >> i) & 1 else 1 for i in range(r)])
    return a[P['orb']]


def _hash(paf, R):
    return (paf.astype(np.uint64) * R[None, :]).sum(1)  # wraps mod 2^64: linear hash


def solve(P, kind, max_solutions=10, verbose=True):
    """returns list of solutions; each solution = list of +-1 sequences on G (length n)"""
    n = P['n']
    t0 = time.time()
    rng = np.random.default_rng(12345)
    Lg = len(P['lrep'])
    Rh = rng.integers(1, 2 ** 62, size=Lg, dtype=np.int64).astype(np.uint64)
    sols = []
    if kind == 'lp':
        bound = 2 * n + 2
        masks, rs, psd = enumerate_blocks(P, bound, [1, -1])
        paf = exact_paf(P, rs, psd)
        target = np.full(Lg, -2, dtype=np.int64)
        h = _hash(paf, Rh)
        ht = _hash(target[None, :], Rh)[0]
        order = np.argsort(h)
        hs = h[order]
        need = ht - h
        pos = np.searchsorted(hs, need)
        for i in range(len(h)):
            j = pos[i]
            while j < len(hs) and hs[j] == need[i]:
                b = order[j]
                if b >= i and np.array_equal(paf[i] + paf[b], target):
                    sols.append([mask_to_seq(P, masks[i]), mask_to_seq(P, masks[b])])
                    if len(sols) >= max_solutions:
                        break
                j += 1
            if len(sols) >= max_solutions:
                break
        info = dict(survivors=len(masks))
    elif kind in ('b1', 'b0'):
        if kind == 'b1':
            bound = 4 * (n + 1)
            T = -4
            allowed = [2, -2, 0] if n % 2 == 0 else [1, -1]
        else:
            bound = 4 * n
            T = 0
            amax = int(np.floor(np.sqrt(4 * n)))
            allowed = [a for a in range(-amax, amax + 1) if (a - n) % 2 == 0]
        masks, rs, psd = enumerate_blocks(P, bound, allowed)
        paf = exact_paf(P, rs, psd)
        h = _hash(paf, Rh)
        target = np.full(Lg, T, dtype=np.int64)
        ht = _hash(target[None, :], Rh)[0]
        info = dict(survivors=len(masks))
        try:
            sols = _join4(P, kind, masks, rs, psd, paf, h, ht, target, bound, max_solutions)
        except JoinTooLarge as e:
            sols = []
            info['undecided'] = str(e)
    else:
        raise ValueError(kind)
    info['secs'] = round(time.time() - t0, 2)
    info['nsol'] = len(sols)
    if verbose:
        print(f"  [{kind}] r={P['r']} L={Lg} survivors={info['survivors']} sols={len(sols)} {info['secs']}s", flush=True)
    return sols, info


MAX_CHECK = 3_000_000_000   # pair checks (before the PSD pair filter) allowed per side
MAX_STORE = 150_000_000     # filtered (C,D) pairs kept in memory


class JoinTooLarge(Exception):
    pass


def _pairs(I, J, psd, h, bound, L, same):
    """PSD-filtered pairs (i in I, j in J) [i<=j if same], returned as (keys, pairs) sorted by hash sum"""
    if len(I) * len(J) > MAX_CHECK:
        raise JoinTooLarge(f'{len(I)}x{len(J)} pair checks')
    chunk = max(1, int(4e7 // max(1, len(J) * max(L, 1))))
    keys, pairs = [], []
    pj = psd[J]
    stored = 0
    for s in range(0, len(I), chunk):
        A = I[s:s + chunk]
        ok = (psd[A][:, None, :] + pj[None, :, :] <= bound + 1e-6).all(-1)
        ii, jj = np.nonzero(ok)
        if same:
            keep = jj >= ii + s
            ii, jj = ii[keep], jj[keep]
        stored += len(ii)
        if stored > MAX_STORE:
            raise JoinTooLarge(f'>{MAX_STORE} filtered pairs')
        keys.append(h[A[ii]] + h[J[jj]])
        pairs.append(np.stack([A[ii], J[jj]], 1))
    if not keys:
        return np.zeros(0, dtype=np.uint64), np.zeros((0, 2), dtype=np.int64)
    keys = np.concatenate(keys)
    pairs = np.concatenate(pairs)
    o = np.argsort(keys)
    return keys[o], pairs[o]


def _probe(IA, IB, keys, pairs, psd, h, ht, paf, target, bound, L, masks, P, sols, max_solutions):
    if len(keys) == 0:
        return False
    if len(IA) * len(IB) > MAX_CHECK:
        raise JoinTooLarge(f'{len(IA)}x{len(IB)} probe checks')
    chunk = max(1, int(4e7 // max(1, len(IB) * max(L, 1))))
    pb = psd[IB]
    for s in range(0, len(IA), chunk):
        A = IA[s:s + chunk]
        ok = (psd[A][:, None, :] + pb[None, :, :] <= bound + 1e-6).all(-1)
        ii, jj = np.nonzero(ok)
        need = ht - h[A[ii]] - h[IB[jj]]
        pos = np.searchsorted(keys, need)
        hit = (pos < len(keys)) & (keys[np.minimum(pos, len(keys) - 1)] == need)
        for q in np.flatnonzero(hit):
            a, b = A[ii[q]], IB[jj[q]]
            j = pos[q]
            while j < len(keys) and keys[j] == need[q]:
                c, d = pairs[j]
                if np.array_equal(paf[a] + paf[b] + paf[c] + paf[d], target):
                    sols.append([mask_to_seq(P, masks[x]) for x in (a, b, c, d)])
                    if len(sols) >= max_solutions:
                        return True
                j += 1
    return False


def _join4(P, kind, masks, rs, psd, paf, h, ht, target, bound, max_solutions):
    """complete 4-list join: all (a,b,c,d) with sum PAF = target and the row-sum condition
    (b1, n even: one block |a|=2, three blocks a=0; b1, n odd: all |a|=1; b0: sum a^2 = 4n).
    Raises JoinTooLarge when the pair lists exceed the caps (slice reported as undecided)."""
    n = P['n']
    L = psd.shape[1] if psd.ndim == 2 else 0
    sols = []
    idx = np.arange(len(masks))
    if kind == 'b1':
        if n % 2 == 0:
            I2 = idx[np.abs(rs) == 2]
            I0 = idx[rs == 0]
            plans = [(I2, I0, I0, I0, True)]
        else:
            plans = [(idx, idx, idx, idx, True)]
    else:
        vals = sorted(set(int(x) for x in np.abs(rs)))
        byabs = {a: idx[np.abs(rs) == a] for a in vals}
        plans = []
        for a in vals:
            for b in vals:
                if b < a:
                    continue
                for c in vals:
                    if c < b:
                        continue
                    d2 = 4 * n - a * a - b * b - c * c
                    if d2 < c * c:
                        continue
                    d = int(round(np.sqrt(d2)))
                    if d * d == d2 and d in vals:
                        plans.append((byabs[a], byabs[b], byabs[c], byabs[d], c == d))
    for IA, IB, IC, ID, same in plans:
        if min(len(IA), len(IB), len(IC), len(ID)) == 0:
            continue
        keys, pairs = _pairs(IC, ID, psd, h, bound, L, same and IC is ID)
        if _probe(IA, IB, keys, pairs, psd, h, ht, paf, target, bound, L, masks, P, sols, max_solutions):
            return sols
    return sols
