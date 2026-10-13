r"""Independent exhaustive check of symmetric GP stage 1 (no alpha-congruence, no S/D split) [M].

Same problem as theory/nd_sym_stage1.py: y_k symmetric on Z_F with letters alpha + 2f t (|t| <= s/2), rows
(2, 0, 0, 0), sum_k PAF(y_k)(0) = 4Fq - 4(q - 1), sum_k PAF(y_k)(h) = -4q (h = 1..F/2).
Method: per block, every configuration (orbit values) with energy (y^2 - 1)/gcd <= E_total and the block's row is
listed with its full PAF signature (lags 0..F/2); blocks 1..3 are listed only as canonical representatives under shift
by F/2 and negation (they have row 0), block 0 under shift by F/2; then an exhaustive 2+2 join by left energy class
with a sorted-key match on the whole signature.  The count is the number of solutions with block 0 of row +2 and
blocks 1..3 canonical (blocks 1..3 in every order).
usage: python3 theory/nd_sym_bruteforce.py F q s [--T "t0,t1,..."]   (default: the slice target)"""
import sys, os, json, time, math, itertools
import numpy as np
from numba import njit

def letters(q, s):
    f = (q - 1) // s
    return sorted({a + 2 * f * t for a in (1, -1) for t in range(-(s // 2), s // 2 + 1)}), f

def configs(F, L, E, unit, row):
    """all orbit-value vectors (reps 0..F/2) with energy <= E and row sum = row"""
    reps = list(range(F // 2 + 1)); w = [1 if (g == 0 or 2 * g == F) else 2 for g in reps]
    out = []
    def rec(i, e, r, cur):
        if i == len(reps):
            if r == row: out.append(list(cur))
            return
        # remaining max |row| contribution bound
        rest = sum(w[i:]) * max(abs(x) for x in L)
        for y in L:
            e2 = e + w[i] * unit[y]
            if e2 > E: continue
            r2 = r + w[i] * y
            cur.append(y); rec(i + 1, e2, r2, cur); cur.pop()
    rec(0, 0, 0, [])
    return np.array(out, np.int64).reshape(-1, len(reps))

def expand(F, V):
    Y = np.zeros((V.shape[0], F), np.int64)
    for i, g in enumerate(range(F // 2 + 1)):
        Y[:, g] = V[:, i]; Y[:, (-g) % F] = V[:, i]
    return Y

def canon_mask(F, Y, neg):
    """True for the canonical representative (lexicographically least) under shift by F/2 (and negation)"""
    n = F // 2; keep = np.ones(len(Y), bool)
    variants = [np.roll(Y, -n, axis=1)]
    if neg: variants += [-Y, -np.roll(Y, -n, axis=1)]
    for Z in variants:
        # Y is canonical iff Y <= Z lexicographically for every variant
        diff = Z - Y; first = np.argmax(diff != 0, axis=1); has = (diff != 0).any(1)
        keep &= ~has | (diff[np.arange(len(Y)), first] > 0)
    return keep

def sigs(Y):
    F = Y.shape[1]
    return np.stack([(Y * np.roll(Y, -h, axis=1)).sum(1) for h in range(F // 2 + 1)], 1)

@njit(cache=True)
def count_join(lk, rk):
    """number of (i, j) with lk[i] == rk[j]; both sorted"""
    i = 0; j = 0; c = 0
    while i < lk.shape[0] and j < rk.shape[0]:
        if lk[i] < rk[j]: i += 1
        elif lk[i] > rk[j]: j += 1
        else:
            a = i
            while i < lk.shape[0] and lk[i] == lk[a]: i += 1
            b = j
            while j < rk.shape[0] and rk[j] == rk[b]: j += 1
            c += (i - a) * (j - b)
    return c

def keyof(S, base, off):
    return ((S + off) * base[None, :]).sum(1)

def run(F, q, s, T=None):
    L, f = letters(q, s); g = 0
    for y in L: g = math.gcd(g, y * y - 1)
    unit = {y: (y * y - 1) // g for y in L}
    if T is None: T = np.array([4 * F * q - 4 * (q - 1)] + [-4 * q] * (F // 2), np.int64)
    E = (int(T[0]) - 4 * F) // g
    t0 = time.time()
    B0 = expand(F, configs(F, L, E, unit, 2)); B0 = B0[canon_mask(F, B0, False)]
    B1 = expand(F, configs(F, L, E, unit, 0)); B1 = B1[canon_mask(F, B1, True)]
    S0 = sigs(B0); S1 = sigs(B1)
    e0 = (S0[:, 0] - F) // g; e1 = (S1[:, 0] - F) // g
    total = 0; npairs = 0
    for EL in range(E + 1):
        # left: block 0 (energy a) + block 1 (energy EL - a); right: blocks 2, 3 (energies summing to E - EL);
        # inside a class the lag-0 component is fixed, so rows are matched exactly on lags 1..F/2 (np.unique ids)
        lk = []
        for a in range(EL + 1):
            I = np.nonzero(e0 == a)[0]; J = np.nonzero(e1 == EL - a)[0]
            if len(I) and len(J):
                lk.append((S0[I][:, None, 1:] + S1[J][None, :, 1:]).reshape(-1, S0.shape[1] - 1))
        if not lk: continue
        rk = []
        ER = E - EL
        for a in range(ER + 1):
            I = np.nonzero(e1 == a)[0]; J = np.nonzero(e1 == ER - a)[0]
            if len(I) and len(J):
                rk.append((T[None, None, 1:] - S1[I][:, None, 1:] - S1[J][None, :, 1:]).reshape(-1, S0.shape[1] - 1))
        if not rk: continue
        lk = np.concatenate(lk); rk = np.concatenate(rk)
        npairs += len(lk) + len(rk)
        _, inv = np.unique(np.concatenate([lk, rk]), axis=0, return_inverse=True)
        inv = inv.ravel(); nl = len(lk)
        cl = np.bincount(inv[:nl], minlength=inv.max() + 1); cr = np.bincount(inv[nl:], minlength=inv.max() + 1)
        total += int((cl * cr).sum())
    out = dict(F=F, q=q, s=s, energy_units=E, block0_configs=int(len(B0)), row0_canonical=int(len(B1)),
               pairs=int(npairs), solutions_ordered_canonical=int(total), secs=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    return out

if __name__ == "__main__":
    a = sys.argv[1:]; T = None
    if "--T" in a:
        i = a.index("--T"); T = np.array(list(map(int, a[i + 1].split(","))), np.int64); a = a[:i] + a[i + 2:]
    run(*map(int, a), T=T)
