r"""Check the staged b >= 3 (and b = 0) systems on the published cores [M].

For each core: its scalar stabiliser S (computed in theory/ab_lib conventions), the CRT split of the cyclic H-factor so
that a prime q | h is its own axis, and for every such q:
  * stage 1: y = compression along Z_q satisfies sum_k PAF_{G'}(y_k) = T1 exactly (two-valued for b >= 3:
    4m - 4b(q - 1) at 0, -4bq on H' \ 0, +4q off H'), y is S-invariant, and y_k(b) is in the letter set of its fibre;
  * alpha-congruence: with f0 = the largest f with every letter = +-1 (mod 2f), alpha = (+-1 with y = alpha mod 2 f0)
    satisfies sum_k PAF(alpha_k) = T1 (mod 4 f0) at every lag, rows of alpha = rows of y (mod 2 f0);
  * Q-compression: z_k on Q = G'/H' has sum_k |z_k^(phi)|^2 = 4b (phi != 1), 4v (phi = 1);
  * stage 2: sum_k |X_k^(phi, psi)|^2 = 4v for every psi != 1 on Z_q (numerical, max error printed).
Also the b = 0 core of H(892) (Z_223, d = 3) as a one-axis slice.
usage: python3 theory/ab_verify_published.py"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ab_lib import AG, Slice, stage2_check, alpha_congruence_mod, factorint

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_core(path):
    S = {}; hdr = None
    for l in open(path):
        if hdr is None: hdr = l.strip(); continue
        if len(l) > 2 and l[0] in 'ABCD' and l[1] == '=':
            S[l[0]] = np.array([1 if ch == '+' else -1 for ch in l[2:].strip()])
    return hdr, np.array([S[k] for k in 'ABCD'])


def recoordinate(X, mods_old, mods_new, axis_map):
    """X on mods_old (ravel); new coordinates: for each new axis, (old axis, modulus) by CRT reduction"""
    Go, Gn = AG(mods_old), AG(mods_new)
    co_n = np.zeros((Go.n, len(mods_new)), np.int64)
    for j, (oa, mod) in enumerate(axis_map): co_n[:, j] = Go.co[:, oa] % mod
    newidx = Gn.idx(co_n)
    assert len(set(newidx.tolist())) == Go.n
    Y = np.zeros_like(X); Y[:, newidx] = X
    return Y


CASES = [
    # name, file, old mods, new mods, axis map, stabiliser generators (per new axis), q axes to test
    ("479", "data/core_v479_b3.txt", (2, 2, 119), (2, 2, 7, 17), [(0, 2), (1, 2), (2, 7), (2, 17)], [(1, 1, 6, 13)], [2, 3], 479, 3),
    ("347", "data/core_v347_b5.txt", (2, 3, 57), (2, 3, 3, 19), [(0, 2), (1, 3), (2, 3), (2, 19)], [(1, 1, 1, 7), (1, 2, 2, 11)], [3, 2], 347, 5),
    ("359", "data/core_v359_b7.txt", (2, 2, 8, 11), (2, 2, 8, 11), [(0, 2), (1, 2), (2, 8), (3, 11)], [(1, 1, 1, 9)], [3], 359, 7),
]


def check_case(name, path, mo, mn, amap, gens, qaxes, v, b, log):
    hdr, X0 = load_core(os.path.join(ROOT, path))
    X = recoordinate(X0, mo, mn, amap)
    G = AG(mn)
    P = G.paf(X).sum(0)
    Hmask = (P == -4 * b); Hmask[0] = True
    m = G.n; c = b + 1
    assert Hmask.sum() * c == m and set(P[1:].tolist()) == {4, -4 * b}
    out = []
    for qa in qaxes:
        sl = Slice(mn, v, b, gens, qa, Hmask)
        assert sl.is_invariant(X), "core not S-invariant"
        Y = sl.compress(X)
        P1 = sl.Gp.paf(Y).sum(0)
        st1 = bool((P1 == sl.T1).all())
        inv = all((Y[:, sl.Gp.act(u)] == Y).all() for u in sl.Sp)
        lett = all(int(Y[k, bb]) in sl.letters[sl.borb[bb]] for k in range(4) for bb in range(sl.Gp.n))
        f0 = alpha_congruence_mod(sl.letters)
        A = np.where(((Y - 1) % (2 * f0)) == 0, 1, -1)
        assert (((Y - A) % (2 * f0)) == 0).all()
        PA = sl.Gp.paf(A).sum(0)
        acong = bool((((PA - sl.T1) % (4 * f0)) == 0).all())
        rows_ok = bool((((A.sum(1) - Y.sum(1)) % (2 * f0)) == 0).all())
        # Q-compression
        Hp = sl.Hpmask
        cos = -np.ones(sl.Gp.n, np.int64); nq = 0
        for bb in range(sl.Gp.n):
            if cos[bb] < 0:
                for hh in np.nonzero(Hp)[0]: cos[sl.Gp.add(bb, hh)] = nq
                nq += 1
        Z = np.zeros((4, nq), np.int64)
        for k in range(4): np.add.at(Z[k], cos, Y[k])
        # characters of Q via the norm identity and per-Q-character power: use Q-PAF on the coset structure
        PZ = np.zeros(nq, np.int64)          # Q-PAF: sum_k sum_q z(q) z(q + d)
        # build Q addition from G' coset representatives
        rep = [int(np.nonzero(cos == j)[0][0]) for j in range(nq)]
        qadd = np.array([[cos[sl.Gp.add(rep[i], rep[j])] for j in range(nq)] for i in range(nq)])
        for d in range(nq):
            PZ[d] = sum(int(Z[k, i]) * int(Z[k, qadd[i, d]]) for k in range(4) for i in range(nq))
        TQ = np.array([sl.T1[(cos == d)].sum() for d in range(nq)])
        qok = bool((PZ == TQ).all())
        s2 = stage2_check(sl, X)
        rec = dict(core=name, q=sl.q, Gp=sl.Gp.mods, S_order=len(sl.S), orbits_per_block=sl.r, Gp_orbits=len(sl.breps),
                   letter_sets=[list(L) for L in sorted({tuple(L) for L in sl.letters}, key=len)],
                   T1_values=sorted(set(sl.T1.tolist())), stage1_exact=st1, y_invariant=inv, letters_ok=lett,
                   f0=f0, alpha_congruence_ok=acong, alpha_rows_ok=rows_ok, Q_order=nq, Q_compression_ok=qok,
                   Q_norm=int((Z ** 2).sum()), Q_norm_expected=(4 * v + 4 * b * (c - 1)) // c,
                   stage2_max_err=s2, rows_y=Y.sum(1).tolist())
        log(json.dumps(rec)); out.append(rec)
    return out


def check_892(log):
    hdr, X = load_core(os.path.join(ROOT, "data/core_v223_b0.txt"))
    v = 223
    # the record is in the g-builder format: B, C, D stored as g(s) = f(-1 - s); natural coordinates undo it
    X = X.copy(); X[1:] = X[1:, (-1 - np.arange(v)) % v]
    assert (AG((v,)).paf(X).sum(0) == np.r_[4 * v, np.zeros(v - 1, np.int64)]).all()
    # stabiliser of the record (d = 3): the cube roots of unity mod 223 (unique subgroup of order 3)
    g = next(u for u in range(2, v) if pow(u, 3, v) == 1)
    sl = Slice((v,), v, 0, [(g,)], 0, None)
    ok_inv = sl.is_invariant(X)
    Y = sl.compress(X)
    rec = dict(core="892 (b = 0, Z_223, d = 3)", invariant=bool(ok_inv), orbits_per_block=sl.r,
               letters=sl.letters[0][:6] + ['...'], rows=Y[:, 0].tolist(), rows_sq=int((Y ** 2).sum()), four_v=4 * v,
               letters_ok=all(int(Y[k, 0]) in sl.letters[0] for k in range(4)), stage2_max_err=stage2_check(sl, X))
    log(json.dumps(rec)); return rec


if __name__ == "__main__":
    logf = open(os.path.join(ROOT, "logs/ab/verify_published.jsonl"), "w")
    def log(s): print(s, flush=True); logf.write(s + "\n")
    allok = True
    for case in CASES:
        for r in check_case(*case, log=log):
            allok &= r['stage1_exact'] and r['y_invariant'] and r['letters_ok'] and r['alpha_congruence_ok'] and \
                r['alpha_rows_ok'] and r['Q_compression_ok'] and r['Q_norm'] == r['Q_norm_expected'] and r['stage2_max_err'] < 1e-6
    r = check_892(log)
    allok &= r['invariant'] and r['letters_ok'] and r['rows_sq'] == r['four_v'] and r['stage2_max_err'] < 1e-6
    log(json.dumps(dict(all_checks_pass=bool(allok))))
