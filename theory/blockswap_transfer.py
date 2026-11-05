r"""Gate 5 of report/prompt_blockswap_863.md: transfer screen of the block-exchanging ansatz (A, B, C, tau C) to the
863-shaped k4 hosts of report/h863_lift_answer.md section 4.

Credit.  Ansatz: the OpenAI model (report/openai_deficit_answer.md, codex/order-3452-compression, 5b26dae).  The
host family (k4, t = 3, G = Z2^2 x Z_m x Z_p, multiplier of order f on Z_p) follows OpenAI Codex's 863 set-up; the
host list is that of theory/h863_transfer_run.py.  Bordered construction: Kumar's form of Theorem 7 of EKLO.
New code (shares nothing with theory/openai_deficit_*.py).

Host (m, p, f): G = Z2^2 x Z_m x Z_p, K = Z_m x Z_p, n = |G| = 4mp, v = n + 3, e = (p - 1)/f classes.  Target
sum_i PAF_i = 4n at 0, -12 on K \ 0, +4 off K.  The involution is tau(q1, q2, u, t) = (q2, q1, -u, t) (an automorphism
commuting with the multiplier and fixing K; the target is tau-invariant [P, as at 863]).  H-level letters
y = eps + f s, eps = +-1, s in {-e, -e+2, ..., e}.
  1. Q-level tuples: odd cell sums, A and B with z01 = z10, the fourth block = C with 01 <-> 10, target
     PAF_Q = 4n - 12(mp - 1) at 0 and 4mp elsewhere; normalised by block sign and A <-> B.
  2. joint row-sum / energy gate: reachable excess energies per block under its Q row sums (capped (sum, excess)
     DP over the cells of a row), joined exactly against the H-level excess 4n - 12(p - 1) - 16m.
  3. the plain H-level lag-0 subset-sum with the ansatz weights (cells counted with their multiplicity), and the
     unrestricted one of theory/h863_lag0_check.py (16m cells of weight 1), for the record.
usage: python3 theory/blockswap_transfer.py [--validate]       -> logs/blockswap/transfer.jsonl
"""
import sys, os, json, time, itertools as it

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

HOSTS = [  # (label, m, p, f)  -- theory/h863_transfer_run.py
    ('analogue 519', 3, 43, 7), ('analogue 583', 5, 29, 7),
    ('863 S7', 5, 43, 7),
    ('1223', 5, 61, 15), ('1223', 5, 61, 5),
    ('1423', 5, 71, 7), ('1583', 5, 79, 13), ('1783', 5, 89, 11), ('2063', 5, 103, 17),
    ('1039', 7, 37, 9), ('1039', 7, 37, 3), ('1319', 7, 47, 23), ('1487', 7, 53, 13), ('1879', 7, 67, 11),
    ('2215', 7, 79, 13),
]


def mult_order_exists(p, f):
    """an element of order f in Z_p^* exists iff f | p - 1 (cyclic group); -1 is in it iff f is even"""
    return (p - 1) % f == 0


def alphabet(p, f):
    e = (p - 1) // f
    return sorted({eps + f * s for eps in (-1, 1) for s in range(-e, e + 1, 2)})


def qpaf(z):
    return tuple(sum(z[q] * z[q ^ a] for q in range(4)) for a in range(4))


def sign_normal(z):
    s = sum(z)
    return s > 0 or (s == 0 and next(v for v in z if v) > 0)


def q_tuples(m, p):
    n = 4 * m * p
    K = m * p
    tgt = (4 * n - 12 * (K - 1), 4 * K, 4 * K, 4 * K)
    lim = int(tgt[0] ** 0.5) + 1
    odd = [v for v in range(-lim, lim + 1) if v % 2]
    fixed = [(a, b, b, c) for a, b, c in it.product(odd, repeat=3)
             if a * a + 2 * b * b + c * c <= tgt[0] and sign_normal((a, b, b, c))]
    cont = {}
    sq = {v: v * v for v in odd}
    for z in it.product(odd, repeat=4):
        if 2 * (sq[z[0]] + sq[z[1]] + sq[z[2]] + sq[z[3]]) > tgt[0] or not sign_normal(z):
            continue
        tz = (z[0], z[2], z[1], z[3])
        key = tuple(x + y for x, y in zip(qpaf(z), qpaf(tz)))
        cont.setdefault(key, []).append(z)
    pf = {z: qpaf(z) for z in fixed}
    sols = []
    for i, za in enumerate(fixed):
        for zb in fixed[i:]:
            key = tuple(t - x - y for t, x, y in zip(tgt, pf[za], pf[zb]))
            for zc in cont.get(key, []):
                sols.append((za, zb, zc))
    return tgt, sols


def row_states(sw, ew, alph, cap):
    """{row_sum: set(excess)} over rows with cell sum-weights sw and energy weights ew; excess capped at cap"""
    st = {(0, 0)}
    for a_w, e_w in zip(sw, ew):
        st = {(s + a_w * y, en + e_w * (y * y - 1)) for s, en in st for y in alph if en + e_w * (y * y - 1) <= cap}
    out = {}
    for s, en in st:
        out.setdefault(s, set()).add(en)
    return out


def sumset(a, b, cap):
    return {x + y for x in a for y in b if x + y <= cap}


def screen(label, m, p, f, validate=False):
    t0 = time.time()
    n = 4 * m * p
    v = n + 3
    e = (p - 1) // f
    alph = alphabet(p, f)
    cap = 4 * n - 12 * (p - 1) - 16 * m
    half = (m - 1) // 2
    out = dict(label=label, m=m, p=p, f=f, e=e, v=v, mult_exists=mult_order_exists(p, f), neg1_in_S=(f % 2 == 0),
               alphabet_pos=[a for a in alph if a > 0], excess=cap)
    # 3. lag-0 subset sums
    w_ans = ([1, 1] + [2] * (2 * half) + [2] * m) * 2 + [2] * (4 * m)
    vals = sorted({y * y - 1 for y in alph})

    def reach(weights):
        r = {0}
        for w in weights:
            r = {x + w * d for x in r for d in vals if x + w * d <= cap}
        return cap in r
    out['lag0_ansatz_weights'] = reach(w_ans)
    out['lag0_unrestricted'] = reach([1] * (16 * m))
    # 1. Q-level
    tgt, sols = q_tuples(m, p)
    out['q_target'] = tgt
    out['q_tuples'] = len(sols)
    # 2. joint gate
    t_fix = row_states([1] + [2] * half, [1] + [2] * half, alph, cap)
    t_mir = row_states([1] * m, [2] * m, alph, cap)

    memo = {}

    def fixed(z):
        k = ('f', z)
        if k not in memo:
            r = {0}
            for part, tab in ((z[0], t_fix), (z[1], t_mir), (z[3], t_fix)):
                r = sumset(r, tab.get(part, set()), cap)
            memo[k] = r
        return memo[k]

    def free(z):
        k = ('c', z)
        if k not in memo:
            r = {0}
            for part in z:
                r = sumset(r, t_mir.get(part, set()), cap)
            memo[k] = r
        return memo[k]
    surv = []
    for za, zb, zc in sols:
        ea, eb, ec = fixed(za), fixed(zb), free(zc)
        if any((cap - x - y) in ec for x in ea for y in eb):
            surv.append(dict(z=[za, zb, zc], minA=min(ea), minB=min(eb), minC=min(ec)))
    out['joint_survivors'] = len(surv)
    canon = {(tuple(s['z'][0]), tuple(s['z'][1]), min(tuple(s['z'][2]), (s['z'][2][0], s['z'][2][2], s['z'][2][1], s['z'][2][3])))
             for s in surv}
    out['joint_survivors_mod_C_tauC'] = len(canon)
    out['survivors'] = surv[:200]
    out['secs'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    validate = '--validate' in sys.argv
    outp = os.path.join(ROOT, 'logs', 'blockswap', 'transfer.jsonl')
    hosts = HOSTS if not validate else [h for h in HOSTS if h[0] in ('analogue 519', '863 S7')]
    for h in hosts:
        r = screen(*h, validate=validate)
        s = json.dumps(r)
        print(json.dumps({k: v for k, v in r.items() if k != 'survivors'}), flush=True)
        if not validate:
            with open(outp, 'a') as fo:
                fo.write(s + '\n')
