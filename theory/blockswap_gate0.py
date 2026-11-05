r"""Gate 0 of report/prompt_blockswap_863.md: independent re-derivation of the block-exchanging ansatz screens.

Credit.  The ansatz (x0, x1, x2, x3) = (A, B, C, tau C) and its cost formula are by the OpenAI model
(report/openai_deficit_answer.md, branch codex/order-3452-compression, commit 5b26dae).  The k4 set-up at v = 863
(G = Z2^2 x Z5 x Z43, S = x4 on Z43) and its certificate are by OpenAI Codex.  The bordered construction is
Kumar's form of Theorem 7 of Egan-Kotsireas-London-O Cathain [EKLO].  This file is new and shares no code with
theory/openai_deficit_*.py (nor with bt_engine): it works from the definitions only.

What it does, for m in {5, 3} (v = 172 m + 3, G = Z2^2 x Z_m x Z43, K = Z_m x Z43):
  1. [M] checks R_{tau C}(g) = R_C(tau g) on random +-1 sequences over G by explicit integer PAFs (all |G| lags), and
     that tau is an automorphism of G that commutes with S and maps K onto K, and that the target F is tau-invariant;
  2. [M] enumerates the Q-level tuples directly in the cell sums z(q) (no Walsh transform): z odd, A and B with
     z(01) = z(10), the fourth block's z equal to C's with 01 and 10 exchanged, and the integer Q-level PAF target
     sum_i PAF_{z_i}(a) = (|K| - 1)(-12) + 4(|G| - |K|)... i.e. 872 at a = 0 and 860 at a != 0 for m = 5;
     normalised by block sign (block total > 0, or lexicographically positive if the total is 0) and by the
     interchange of A and B (unordered pair).  C <-> tau C is NOT identified (as in the codex count);
  3. [M] the joint Q-row-sum / energy gate: per tuple, the set of reachable excess energies sum w (y^2 - 1) of every
     block with its Q row sums, by exhaustive enumeration of rows (numpy, not a DP), and an exact join of the three
     independent blocks' sets against the excess 16 m 43 - 504 - 16 m ... = 2936 - 80 = 2856 at m = 5.
  4. [M] the plain lag-0 subset-sum (weights four 1s and 2(m-1)+... twos) for the record.
Output: logs/blockswap/gate0_m{m}.json
"""
import sys, os, json, time, itertools as it
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
ALPH = np.array(sorted(s * v for v in (1, 13, 15, 27, 29, 41, 43) for s in (1, -1)), np.int64)


# ---------------------------------------------------------------- 1. the group, tau and the PAF identity
def group(m):
    els = [(q1, q2, u, t) for q1 in range(2) for q2 in range(2) for u in range(m) for t in range(43)]
    idx = {g: i for i, g in enumerate(els)}
    return els, idx


def add(g, h, m):
    return ((g[0] + h[0]) % 2, (g[1] + h[1]) % 2, (g[2] + h[2]) % m, (g[3] + h[3]) % 43)


def tau(g, m):
    return (g[1], g[0], (-g[2]) % m, g[3])


def smul(g):
    return (g[0], g[1], g[2], (4 * g[3]) % 43)


def paf_full(x, els, idx, m):
    """R_x(g) = sum_h x(h) x(h + g), integer, for every g (explicit index tables)."""
    n = len(els)
    R = np.zeros(n, np.int64)
    for gi, g in enumerate(els):
        perm = np.array([idx[add(h, g, m)] for h in els])
        R[gi] = int((x * x[perm]).sum())
    return R


def target_full(m, els):
    out = []
    for g in els:
        if g == (0, 0, 0, 0):
            out.append(4 * (len(els)))
        elif g[0] == 0 and g[1] == 0:
            out.append(-12)
        else:
            out.append(4)
    return np.array(out, np.int64)


def check_identity(m, trials=3, seed=863):
    els, idx = group(m)
    rng = np.random.default_rng(seed)
    tperm = np.array([idx[tau(g, m)] for g in els])
    res = {}
    # automorphism: tau(g + h) = tau g + tau h on all pairs (small group, exhaustive on a sample of g)
    auto = all(tau(add(g, h, m), m) == add(tau(g, m), tau(h, m), m) for g in els[::7] for h in els)
    comm = all(tau(smul(g), m) == smul(tau(g, m)) for g in els)
    kk = all((tau(g, m)[0], tau(g, m)[1]) == (0, 0) for g in els if (g[0], g[1]) == (0, 0))
    inv = all(tau(tau(g, m), m) == g for g in els)
    F = target_full(m, els)
    res.update(tau_automorphism=auto, tau_commutes_with_S=comm, tau_maps_K_to_K=kk, tau_involution=inv,
               target_tau_invariant=bool((F == F[tperm]).all()))
    ok = True
    for _ in range(trials):
        C = rng.choice([-1, 1], size=len(els)).astype(np.int64)
        tC = C[tperm]                        # (tau C)(g) = C(tau g)
        RC, RtC = paf_full(C, els, idx, m), paf_full(tC, els, idx, m)
        ok &= bool((RtC == RC[tperm]).all())
        # and with an S-invariant C (constant on <4>-orbits in Z43, singleton t = 0)
    res['identity_R_tauC_eq_RC_tau_random'] = ok
    res['identity_trials'] = trials
    return res


# ---------------------------------------------------------------- 2. Q-level tuples from the definitions
def qpaf(z):
    """PAF on Q = Z2^2 of z = (z00, z01, z10, z11), cell index q = 2 q1 + q2; lag a: sum_q z(q) z(q ^ a)."""
    return tuple(sum(z[q] * z[q ^ a] for q in range(4)) for a in range(4))


def sign_normal(z):
    s = sum(z)
    if s > 0:
        return True
    if s < 0:
        return False
    return next(v for v in z if v != 0) > 0


def q_tuples(m):
    K = m * 43
    n = 4 * K
    tgt = (4 * n + (K - 1) * (-12), 4 * K, 4 * K, 4 * K)
    lim = int(np.sqrt(tgt[0])) + 1
    odd = [v for v in range(-lim, lim + 1) if v % 2]
    fixed = []          # tau-invariant blocks: z01 = z10
    for a, b, c in it.product(odd, repeat=3):
        z = (a, b, b, c)
        if a * a + 2 * b * b + c * c <= tgt[0] and sign_normal(z):
            fixed.append(z)
    free = []
    for z in it.product(odd, repeat=4):
        if 2 * sum(v * v for v in z) <= tgt[0] and sign_normal(z):
            free.append(z)
    pafs = {z: qpaf(z) for z in fixed}
    # contribution of C and tau C
    cont = {}
    for z in free:
        tz = (z[0], z[2], z[1], z[3])
        cont.setdefault(tuple(x + y for x, y in zip(qpaf(z), qpaf(tz))), []).append(z)
    sols = []
    for i, za in enumerate(fixed):
        for zb in fixed[i:]:
            key = tuple(t - x - y for t, x, y in zip(tgt, pafs[za], pafs[zb]))
            for zc in cont.get(key, []):
                sols.append((za, zb, zc, (zc[0], zc[2], zc[1], zc[3])))
    return tgt, sols


# ---------------------------------------------------------------- 3. joint row-sum / energy gate (exhaustive rows)
def row_table(sw, ew):
    """all rows with cell sum-weights sw (a cell of sum-weight w stands for w equal entries of the row) and energy
    weights ew (how often the cell's value occurs in the lag-0 energy of the whole quadruple); returns
    {row_sum: set(excess)} with row_sum = sum sw*y and excess = sum ew*(y^2 - 1)."""
    k = len(sw)
    grid = np.array(np.meshgrid(*([ALPH] * k), indexing='ij')).reshape(k, -1)
    w = np.array(sw, np.int64)[:, None]
    we = np.array(ew, np.int64)[:, None]
    s = (w * grid).sum(0)
    e = (we * (grid * grid - 1)).sum(0)
    out = {}
    for sv, ev in set(zip(s.tolist(), e.tolist())):
        out.setdefault(sv, set()).add(ev)
    return out


def sumset(a, b, cap):
    return {x + y for x in a for y in b if x + y <= cap}


def block_energy_sets(m, cap):
    """returns two functions: fixed(z) and free(z) -> set of excess energies (already weighted, C counted twice)."""
    half = (m - 1) // 2
    # q = 00 or 11 row of a tau-invariant block: cell u=0 and the pairs {u, -u} (both in the row and in the energy)
    t_fix = row_table([1] + [2] * half, [1] + [2] * half)
    # q = 01 row (the q = 10 row is its mirror image): m cells, row sum unweighted, energy counted twice
    t_mir = row_table([1] * m, [2] * m)
    # a row of C: row sum unweighted, energy counts twice (C and tau C)
    t_c = row_table([1] * m, [2] * m)

    def fixed(z):
        r = {0}
        for part, tab in ((z[0], t_fix), (z[1], t_mir), (z[3], t_fix)):
            r = sumset(r, tab.get(part, set()), cap)
        return r

    def free(z):
        r = {0}
        for part in z:
            r = sumset(r, t_c.get(part, set()), cap)
        return r
    return fixed, free


def joint_gate(m, sols):
    cap = 16 * m * 43 - 12 * 42 - 16 * m       # sum of the H-level lag-0 excess: |G|*4 + 42*(-12) - 4*4m
    fixed, free = block_energy_sets(m, cap)
    surv = []
    detail = []
    for za, zb, zc, ztc in sols:
        ea, eb, ec = fixed(za), fixed(zb), free(zc)
        ok = any((cap - x - y) in ec for x in ea for y in eb)
        detail.append(dict(z=[za, zb, zc, ztc], minA=min(ea) if ea else None, minB=min(eb) if eb else None,
                           minC=min(ec) if ec else None, survives=ok))
        if ok:
            surv.append((za, zb, zc, ztc))
    return cap, surv, detail


def lag0_only(m):
    cap = 16 * m * 43 - 12 * 42 - 16 * m
    half = (m - 1) // 2
    weights = ([1, 1] + [2] * (2 * half) + [2] * m) * 2 + [2] * (4 * m)
    reach = {0}
    for w in weights:
        reach = {r + w * (int(a) * int(a) - 1) for r in reach for a in ALPH[ALPH > 0] if r + w * (int(a) * int(a) - 1) <= cap}
    return dict(cap=cap, n_cells=len(weights), weights={1: weights.count(1), 2: weights.count(2)},
                reachable=cap in reach)


def main(m):
    t0 = time.time()
    out = dict(m=m, v=172 * m + 3, script='theory/blockswap_gate0.py')
    out['identity'] = check_identity(m)
    tgt, sols = q_tuples(m)
    out['q_target_paf'] = tgt
    out['q_target_psd'] = [sum(tgt), tgt[0] - tgt[1]]
    out['q_tuples'] = len(sols)
    out['q_tuples_zero_total_blocks'] = sum(1 for s in sols for z in s[:3] if sum(z) == 0)
    cap, surv, detail = joint_gate(m, sols)
    out['excess_energy'] = cap
    out['joint_survivors'] = len(surv)
    out['lag0_only'] = lag0_only(m)
    out['tuples'] = detail
    # count modulo C <-> tau C as well (for information)
    canon = {(tuple(s[0]), tuple(s[1]), min(tuple(s[2]), tuple(s[3]))) for s in surv}
    out['joint_survivors_mod_C_tauC'] = len(canon)
    out['secs'] = round(time.time() - t0, 2)
    os.makedirs(os.path.join(ROOT, 'logs', 'blockswap'), exist_ok=True)
    with open(os.path.join(ROOT, 'logs', 'blockswap', f'gate0_m{m}.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'tuples'}))
    return out


if __name__ == '__main__':
    for mm in ([int(a) for a in sys.argv[1:]] or [5, 3]):
        main(mm)
