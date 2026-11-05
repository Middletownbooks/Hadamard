#!/usr/bin/env python3
r"""[M] Review of Claim 1 (the k4 compression normal form) of report/openai_4p3_answer.md by exhaustive brute force.

Claim 1, its proof and the repair-cost table are by the OpenAI model (report/openai_4p3_answer.md, Proposal 1).  This
script is new (Claude Code review, report/openai_4p3_review.md) and imports no repository code at all, in particular
nothing from theory/openai_4p3_*.py.

Problem.  Q = Z2^2 (addition XOR on 0..3) or Q = Z4; m odd.  Count every ORDERED 4-tuple (z_1, .., z_4) of odd
integer vectors on Q with
        sum_i PAF_{z_i}(d) = 4m + 12 (d = 0),   4m (d != 0),        PAF_z(d) = sum_q z(q) z(q + d).
Method: list every odd 4-vector of energy <= 4m (the other three rows have energy >= 4 each), group the rows by their
PAF vector (the condition depends on the PAF vectors only), and join the classes 2 + 2 on exact integer keys (a
linear balanced-digit encoding of the 4-vector of summed PAFs).  Two exact outputs per (Q, m):

  N(m)   = number of ordered solutions, by the weighted pair convolution sum_k W(k) W(T - k);
  multisets of PAF classes solving the target, each checked for the normal form: every row of every class used is
           constant (a, a, a, a) or a singleton row (three equal entries, the fourth differing by 2), with exactly one
           constant class and three singleton classes per solution.

The normal form predicts N(m) = 4 * 2 * 8^3 * R(m), R(m) = #{(a, k1, k2, k3): a > 0 odd, k_j >= 0,
a^2 + sum_j k_j (k_j + 1) = m} (position of the constant row; +-a; for each singleton row with pronic value k(k+1),
its unique odd t, the position of the exceptional entry and the sign).  N(m) = 4096 R(m), together with the converse
(every normal-form tuple solves the target; checked here directly on every multiset), proves that the solution set is
exactly the normal-form set at that (Q, m) [M].

usage: python3 theory/4p3_review_nf.py --q v4 --mmax 99 --out logs/4p3_review/nf_v4.json [--seconds 1200]
"""
import argparse, json, math, sys, time
import numpy as np

B = 1 << 14                       # digit base of the linear key (|digit| < B/2 always holds for m <= 2000)


def partner(q):
    """partner[d, u] = u + d on Q"""
    if q == 'v4':
        return np.array([[u ^ d for u in range(4)] for d in range(4)], np.int64)
    return np.array([[(u + d) % 4 for u in range(4)] for d in range(4)], np.int64)


def all_rows(E):
    """every odd integer 4-vector with sum of squares <= E"""
    R = int(math.isqrt(max(E - 3, 1)))
    if R % 2 == 0:
        R -= 1
    vals = np.arange(-R, R + 1, 2, dtype=np.int64)
    g = np.stack(np.meshgrid(vals, vals, vals, vals, indexing='ij'), -1).reshape(-1, 4)
    return g[(g * g).sum(1) <= E]


def pafs(Z, part):
    return np.stack([(Z * Z[:, part[d]]).sum(1) for d in range(4)], 1)


def shapes(Z):
    """0 = constant, 1 = singleton (three equal, fourth differs by 2), 2 = other"""
    s = np.sort(Z, 1)
    const = (s[:, 0] == s[:, 3])
    sing = ((s[:, 0] == s[:, 2]) & (s[:, 3] == s[:, 2] + 2)) | ((s[:, 1] == s[:, 3]) & (s[:, 0] == s[:, 1] - 2))
    return np.where(const, 0, np.where(sing, 1, 2)).astype(np.int64)


def enc(P):
    P = np.asarray(P, np.int64)
    return ((P[..., 0] * B + P[..., 1]) * B + P[..., 2]) * B + P[..., 3]


def R_count(m):
    """ordered (a, k1, k2, k3), a > 0 odd, k_j >= 0, a^2 + sum k_j (k_j + 1) = m"""
    pron = []
    k = 0
    while k * (k + 1) <= m:
        pron.append(k * (k + 1)); k += 1
    cnt = 0
    a = 1
    while a * a <= m:
        rest = m - a * a
        for x in pron:
            for y in pron:
                z = rest - x - y
                if z >= 0 and z in pron:
                    cnt += 1
        a += 2
    return cnt


def nf_multisets(m):
    """the normal-form solutions as multisets: (a, sorted (t1, t2, t3)) with t_j the odd t of each pronic value"""
    out = []
    a = 1
    while a * a <= m:
        rest = m - a * a
        ts = []
        t = -1
        cand = [t for t in range(-2 * int(math.isqrt(m)) - 3, 2 * int(math.isqrt(m)) + 4) if t % 2]
        tv = {t: t * (t + 1) for t in cand if 0 <= t * (t + 1) <= rest}
        keys = sorted(tv)
        for i1, t1 in enumerate(keys):
            for i2 in range(i1, len(keys)):
                t2 = keys[i2]
                for i3 in range(i2, len(keys)):
                    t3 = keys[i3]
                    if tv[t1] + tv[t2] + tv[t3] == rest:
                        out.append((a, (t1, t2, t3)))
        a += 2
    return out


def run_m(m, P, cnt, shp_sets, part, cap_matches=10 ** 6):
    E = 4 * m + 12
    T = np.array([E, 4 * m, 4 * m, 4 * m], np.int64)
    keep = np.flatnonzero(P[:, 0] <= 4 * m)
    Pk, ck = P[keep], cnt[keep]
    en = Pk[:, 0]
    order = np.argsort(en, kind='stable')
    Pk, ck, en, keep = Pk[order], ck[order], en[order], keep[order]
    nc = len(Pk)
    # ---- ordered pair convolution (exact int64 weights, Python-int total)
    keys_l, w_l, iu, ju = [], [], [], []
    for i in range(nc):
        lim = E - 8 - en[i]
        jmax = int(np.searchsorted(en, lim, 'right'))
        if jmax == 0:
            continue
        js = np.arange(jmax)
        keys_l.append(enc(Pk[i][None, :] + Pk[js]))
        w_l.append(ck[i] * ck[js])
        iu.append(np.full(jmax, i, np.int32)); ju.append(js.astype(np.int32))
    keys = np.concatenate(keys_l); w = np.concatenate(w_l)
    I = np.concatenate(iu); J = np.concatenate(ju)
    o = np.argsort(keys, kind='stable')
    ks, ws, Is, Js = keys[o], w[o], I[o], J[o]
    starts = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    uk = ks[starts]
    W = np.add.reduceat(ws, starts)
    tk = int(enc(T))
    pos = np.searchsorted(uk, tk - uk)
    pos = np.minimum(pos, len(uk) - 1)
    hit = uk[pos] == tk - uk
    N = sum(int(a) * int(b) for a, b in zip(W[hit], W[pos[hit]]))
    # ---- explicit multisets: unordered pairs (i <= j) matched with (k <= l), j <= k
    ends = np.r_[starts[1:], len(ks)]
    found = set()
    nmatch = 0
    for g in np.flatnonzero(hit):
        if uk[g] > tk - uk[g]:
            continue
        h = pos[g]
        A = [(Is[x], Js[x]) for x in range(starts[g], ends[g]) if Is[x] <= Js[x]]
        Bq = [(Is[x], Js[x]) for x in range(starts[h], ends[h]) if Is[x] <= Js[x]]
        nmatch += len(A) * len(Bq)
        if nmatch > cap_matches:
            raise MemoryError('multiset enumeration cap')
        for a in A:
            for b in Bq:
                found.add(tuple(sorted((int(keep[a[0]]), int(keep[a[1]]), int(keep[b[0]]), int(keep[b[1]])))))
    viol = []
    nfsets = []
    for ms in sorted(found):
        sh = [shp_sets[c] for c in ms]
        Ptot = P[list(ms)].sum(0)
        assert (Ptot == T).all()
        pure = all(len(s) == 1 for s in sh)
        kinds = sorted(next(iter(s)) for s in sh) if pure else None
        if not pure or kinds != [0, 1, 1, 1]:
            viol.append(dict(classes=[P[c].tolist() for c in ms], shapes=[sorted(s) for s in sh]))
        else:
            nfsets.append(ms)
    return dict(m=m, classes_used=int(nc), ordered_pairs=int(len(keys)), N=N, R=R_count(m),
                predicted=4096 * R_count(m), N_equals_prediction=(N == 4096 * R_count(m)),
                multisets=len(found), nf_multisets_predicted=len(nf_multisets(m)), violations=viol)


def naive_count(m, part):
    """second method for small m (3 + 1, no PAF classes, no key encoding): for every ordered triple of rows, the number
    of rows whose PAF is the remainder (dictionary of PAF tuples), and the normal form of every solution row"""
    E = 4 * m + 12
    Z = all_rows(4 * m)
    P = pafs(Z, part)
    sh = shapes(Z)
    d = {}
    for p, s in zip(map(tuple, P.tolist()), sh.tolist()):
        d.setdefault(p, []).append(s)
    T = (E, 4 * m, 4 * m, 4 * m)
    N, bad = 0, 0
    Pl = P.tolist()
    for i in range(len(Z)):
        for j in range(len(Z)):
            e2 = Pl[i][0] + Pl[j][0]
            if e2 > E - 8:
                continue
            for k in range(len(Z)):
                rest = tuple(T[x] - Pl[i][x] - Pl[j][x] - Pl[k][x] for x in range(4))
                if rest in d:
                    N += len(d[rest])
                    for s4 in d[rest]:
                        if sorted((sh[i], sh[j], sh[k], s4)) != [0, 1, 1, 1]:
                            bad += 1
    return N, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', choices=('v4', 'z4'), required=True)
    ap.add_argument('--mmin', type=int, default=3)
    ap.add_argument('--mmax', type=int, default=99)
    ap.add_argument('--naive-mmax', type=int, default=0, help='also run the 3 + 1 second method up to this m')
    ap.add_argument('--seconds', type=float, default=1200)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    t0 = time.time()
    part = partner(a.q)
    naive = []
    for m in range(3, a.naive_mmax + 1, 2):
        t1 = time.time()
        N, bad = naive_count(m, part)
        naive.append(dict(m=m, N=N, predicted=4096 * R_count(m), non_nf_solution_rows=bad, secs=round(time.time() - t1, 1)))
        print(json.dumps(dict(naive=naive[-1])), flush=True)
    Z = all_rows(4 * a.mmax)
    P = pafs(Z, part)
    sh = shapes(Z)
    k = enc(P)
    uk, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    first = np.zeros(len(uk), np.int64)
    first[inv[::-1]] = np.arange(len(inv))[::-1]
    Pc = P[first]
    shp_sets = [set() for _ in range(len(uk))]
    for c, s in set(zip(inv.tolist(), sh.tolist())):
        shp_sets[c].add(s)
    mixed = sum(1 for s in shp_sets if len(s) > 1)
    res = dict(script='theory/4p3_review_nf.py', command=' '.join(sys.argv), seed=None, q=a.q,
               rows_listed=int(len(Z)), paf_classes=int(len(uk)), mixed_shape_classes=int(mixed),
               arithmetic='int64 keys and weights, Python-int totals', naive_3plus1=naive, records=[])
    for m in range(a.mmin | 1, a.mmax + 1, 2):
        if time.time() - t0 > a.seconds:
            res['stopped_at_time_limit_before_m'] = m
            break
        t1 = time.time()
        r = run_m(m, Pc, cnt.astype(np.int64), shp_sets, part)
        r['secs'] = round(time.time() - t1, 2)
        res['records'].append(r)
        print(json.dumps({x: r[x] for x in ('m', 'N', 'predicted', 'N_equals_prediction', 'multisets',
                                             'nf_multisets_predicted', 'secs')} | {'violations': len(r['violations'])}),
              flush=True)
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)
    res['all_ok'] = all(r['N_equals_prediction'] and not r['violations'] and r['multisets'] == r['nf_multisets_predicted']
                        for r in res['records']) and all(x['N'] == x['predicted'] and not x['non_nf_solution_rows']
                                                         for x in naive)
    res['secs'] = round(time.time() - t0, 1)
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(dict(q=a.q, all_ok=res['all_ok'], m_done=[res['records'][0]['m'], res['records'][-1]['m']],
                          secs=res['secs'])))


if __name__ == '__main__':
    main()
