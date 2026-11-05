#!/usr/bin/env python3
r"""[M] Review of Claim 2 of report/openai_4p3_answer.md (one finite border certificate for every valid k4 core with
odd |K|): independent check of the certificates, and end-to-end assembly of real cores with ONLY that certificate.

Credit.  The universal-border theorem and the certificates (logs/openai_4p3/border_certificates.json) are by the OpenAI
model; the bordered construction is Kumar's form of Theorem 7 of Egan-Kotsireas-London-O Cathain (EKLO); the
frame-search border reduction that found the certificates is the repository's (theory/bt_border.py).  The cores: 1916
(v = 479, data/core_v479_b3.txt, this repository), 3964 (Kumar, https://github.com/vatskishan03/hadamard,
hadamard/3964/certificate.json, Apache-2.0), and elliptic-ansatz cores found by theory/4p3_review_elliptic.py.
This script is new (Claude Code review) and imports no repository code: the group, the Goethals-Seidel array (Kumar's
layout, from its definition), the core check, the normalisation and the assembly are written here from scratch.  The
assembled matrices are verified by theory/4p3_review_hverify.py, a separate process that shares no code with this file.

GS array over an abelian group (Kumar's layout): with A[g, h] = a(h - g), Y2[g, h] = y(g + h), Y2'[g, h] = y(-g - h),
        [ A    B2    C2    D2  ]
        [-B2   A     D2'  -C2' ]
        [-C2  -D2'   A     B2' ]
        [-D2   C2'  -B2'   A   ]
H = [[X, R], [C, M]], R[x, (i, g)] = Rh[x, 4i + q(g)], C[(i, g), x] = Ch[4i + q(g), x], q(g) the coset of g mod K.
Certificate identities (checked here by integer arithmetic): entries +-1; X X^T = 12 I; Rh Rh^T = 16 I;
Ch Ch^T = I_4 (x) Gamma_Q (12 on the diagonal, -4 off it); X Ch^T + Rh Mc(z*)^T = 0 for z* = (0, 2d0, 2d0, 2d0); every
4-entry block of every Rh row sums to 0 (so Rh Mc(z)^T is unchanged by adding constants to the rows of z).

usage: python3 theory/4p3_review_border.py --out logs/4p3_review/border.json [--kumar DIR] [--elliptic FILE]
"""
import argparse, hashlib, json, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


# ----------------------------------------------------------------------------------------------- group (mixed radix)
class Group:
    """G = Z_m1 x ... x Z_mk, element index in mixed radix with the FIRST coordinate slowest; the 2-part is the first
    one or two coordinates ((2, 2) or (4,)), K = the odd coordinates; the coset index q(g) is 2 x1 + x2 or x1."""

    def __init__(self, mods, qtype):
        self.mods = tuple(mods)
        self.n = int(np.prod(mods))
        self.qtype = qtype
        co = np.array(np.unravel_index(np.arange(self.n), self.mods)).T
        self.co = co
        M = np.array(self.mods)

        def idx(c):
            return np.ravel_multi_index(tuple(np.moveaxis(np.mod(c, M), -1, 0)), self.mods)
        self.idx = idx
        self.add = idx(co[:, None, :] + co[None, :, :]).astype(np.int32)
        self.neg = idx(-co)
        if qtype == 'v4':
            assert self.mods[:2] == (2, 2) and all(m % 2 for m in self.mods[2:])
            self.q = 2 * co[:, 0] + co[:, 1]
            self.qadd = np.array([[u ^ v for v in range(4)] for u in range(4)])
        else:
            assert self.mods[0] == 4 and all(m % 2 for m in self.mods[1:])
            self.q = co[:, 0] % 4
            self.qadd = np.array([[(u + v) % 4 for v in range(4)] for u in range(4)])
        self.K = np.flatnonzero(self.q == 0)

    def qrep(self, s):
        """the element of coset s with zero odd part"""
        c = np.zeros(len(self.mods), np.int64)
        if self.qtype == 'v4':
            c[0], c[1] = s // 2, s % 2
        else:
            c[0] = s
        return int(self.idx(c))


def gs_array(seqs, add, neg):
    a, b, c, d = [np.asarray(s, np.int64) for s in seqs]
    n = len(a)
    diff = add[neg[:, None], np.arange(n)[None, :]]          # h - g
    A = a[diff]
    Y = lambda y: y[add]                                       # y(g + h)
    Yr = lambda y: y[neg[add]]                                 # y(-g - h)
    return np.block([[A, Y(b), Y(c), Y(d)],
                     [-Y(b), A, Yr(d), -Yr(c)],
                     [-Y(c), -Yr(d), A, Yr(b)],
                     [-Y(d), Yr(c), -Yr(b), A]])


def q_tables(qtype):
    add = np.array([[(u ^ v) if qtype == 'v4' else (u + v) % 4 for v in range(4)] for u in range(4)])
    neg = np.array([np.flatnonzero(add[u] == 0)[0] for u in range(4)])
    return add, neg


# ----------------------------------------------------------------------------------------------- certificate
def check_certificate(cert):
    q = cert['q']
    X, Rh, Ch = (np.array(cert[k], np.int64) for k in ('X', 'Rh', 'Ch'))
    add, neg = q_tables(q)
    zstar = [np.zeros(4, np.int64)] + [np.array([2, 0, 0, 0], np.int64)] * 3
    Mc = gs_array(zstar, add, neg)
    Gam = np.where(np.eye(4, dtype=bool), 12, -4)
    out = dict(q=q, frame_index_zero_based=cert.get('frame_index_zero_based'),
               shapes=[list(X.shape), list(Rh.shape), list(Ch.shape)],
               entries_pm1=bool(all(np.isin(M, (-1, 1)).all() for M in (X, Rh, Ch))),
               XXt_12I=bool(np.array_equal(X @ X.T, 12 * np.eye(12, dtype=np.int64))),
               RhRht_16I=bool(np.array_equal(Rh @ Rh.T, 16 * np.eye(12, dtype=np.int64))),
               ChCht_I4xGamma=bool(np.array_equal(Ch @ Ch.T, np.kron(np.eye(4, dtype=np.int64), Gam))),
               link_zero=bool(not (X @ Ch.T + Rh @ Mc.T).any()),
               Rh_blocks_balanced=bool(all((Rh[:, 4 * j:4 * j + 4].sum(1) == 0).all() for j in range(4))),
               own_Mc_equals_stored_Mc=bool('Mc' in cert and np.array_equal(Mc, np.array(cert['Mc']))))
    # transfer: adding the constant c_i to row i of z changes Rh Mc(z)^T by c_i Rh Mc(1_i)^T, which must vanish
    tr = []
    for i in range(4):
        z = [np.zeros(4, np.int64) for _ in range(4)]
        z[i] = np.ones(4, np.int64)
        tr.append(bool(not (Rh @ gs_array(z, add, neg).T).any()))
    out['constant_rows_annihilated'] = tr
    out['ok'] = all(out[k] for k in ('entries_pm1', 'XXt_12I', 'RhRht_16I', 'ChCht_I4xGamma', 'link_zero',
                                     'Rh_blocks_balanced')) and all(tr)
    return out, (X, Rh, Ch)


# ----------------------------------------------------------------------------------------------- cores
def paf_sum(xs, G):
    tot = np.zeros(G.n, np.int64)
    for x in xs:
        x = np.asarray(x, np.int64)
        for g in range(G.n):
            tot[g] += int(x @ x[G.add[:, g]])
    return tot


def core_ok(xs, G):
    tot = paf_sum(xs, G)
    want = np.where(G.q == 0, -12, 4).astype(np.int64)
    want[0] = 4 * G.n
    bad = np.flatnonzero(tot != want)
    return dict(ok=bool(len(bad) == 0 and all(np.isin(x, (-1, 1)).all() and len(x) == G.n for x in xs)),
                lags_checked=int(G.n), bad_lags=int(len(bad)))


def compressions(xs, G):
    return [np.bincount(G.q, weights=np.asarray(x, np.int64), minlength=4).astype(np.int64).tolist() for x in xs]


def normalise(xs, G):
    """row permutation (constant row first), row signs and Q-translations putting the exceptional compression entries
    at q = 0 with value t + 2; returns the new sequences and the record of the operations"""
    z = compressions(xs, G)
    const, sing = [], []
    for i, r in enumerate(z):
        if len(set(r)) == 1:
            const.append(i)
            continue
        vals = sorted(set(r))
        cnt = {v: r.count(v) for v in vals}
        if len(vals) != 2 or sorted(cnt.values()) != [1, 3]:
            raise ValueError('row %d compression %s not in normal form' % (i, r))
        e = next(v for v in vals if cnt[v] == 1)
        u = next(v for v in vals if cnt[v] == 3)
        if abs(e - u) != 2:
            raise ValueError('row %d compression %s not in normal form' % (i, r))
        sing.append((i, r.index(e), 1 if e > u else -1))
    if len(const) != 1 or len(sing) != 3:
        raise ValueError('need one constant and three singleton rows: %s' % z)
    ops = [dict(row=const[0], sign=1, q_translation=0)] + [dict(row=i, sign=s, q_translation=sq) for i, sq, s in sing]
    new = []
    for op in ops:
        g0 = G.qrep(op['q_translation'])
        x = np.asarray(xs[op['row']], np.int64)
        new.append(op['sign'] * x[G.add[:, g0]])              # x'(g) = sign * x(g + s)
    zn = compressions(new, G)
    c = [zn[0][0]] + [zn[j][1] for j in (1, 2, 3)]
    assert len(set(zn[0])) == 1 and all(zn[j] == [c[j] + 2, c[j], c[j], c[j]] for j in (1, 2, 3)), zn
    return new, dict(original_compressions=z, operations=ops, normalised_compressions=zn, constants=c)


def assemble(xs, G, X, Rh, Ch):
    M = gs_array(xs, G.add, G.neg)
    n = G.n
    cols = np.concatenate([4 * i + G.q for i in range(4)])     # column (i, g) -> 4i + q(g)
    R = Rh[:, cols]
    C = Ch[cols, :]
    H = np.block([[X, R], [C, M]])
    return H.astype(np.int8)


def load_1916():
    lines = open(os.path.join(ROOT, 'data/core_v479_b3.txt')).read().split('\n')
    assert lines[0].startswith('v=479') and 'group=2x2x119' in lines[0]
    seqs = {}
    for L in lines[1:]:
        if '=' in L:
            k, v = L.split('=', 1)
            seqs[k.strip()] = np.array([1 if ch == '+' else -1 for ch in v.strip()], np.int64)
    return Group((2, 2, 119), 'v4'), [seqs[k] for k in 'ABCD']


def load_3964(kdir):
    c = json.load(open(os.path.join(kdir, 'hadamard/3964/certificate.json')))
    assert c['group'] == [2, 2, 247] and c['subgroup_generators'] == [[0, 0, 1]]
    xs = [np.array([1 if ch == '+' else -1 for ch in s], np.int64) for s in c['sequences']]
    return Group((2, 2, 247), 'v4'), xs, c


def scramble(xs, G, rng):
    """a random equivalent core: x'_i(g) = s_i x_{pi(i)}(phi(g) + g_i) with a random row permutation pi, row signs s_i,
    translations g_i in G and an automorphism phi of G (GL(2, 2) or -1 on the 2-part, a unit on each odd coordinate);
    phi(K) = K, so sum_i PAF is preserved [P]"""
    k = len(G.mods)
    M = np.zeros((k, k), np.int64)
    if G.qtype == 'v4':
        while True:
            A2 = rng.integers(0, 2, (2, 2))
            if (A2[0, 0] * A2[1, 1] - A2[0, 1] * A2[1, 0]) % 2:
                break
        M[:2, :2] = A2
        first_odd = 2
    else:
        M[0, 0] = int(rng.choice((1, 3)))
        first_odd = 1
    for j in range(first_odd, k):
        m = G.mods[j]
        while True:
            u = int(rng.integers(1, m))
            if np.gcd(u, m) == 1:
                break
        M[j, j] = u
    phi = G.idx(G.co @ M.T)
    assert len(set(phi.tolist())) == G.n
    perm = rng.permutation(4)
    signs = rng.choice((-1, 1), 4)
    trans = rng.integers(0, G.n, 4)
    out = []
    for i in range(4):
        x = np.asarray(xs[perm[i]], np.int64)
        out.append(int(signs[i]) * x[G.add[phi, int(trans[i])]])
    return out, dict(phi=M.tolist(), perm=perm.tolist(), signs=signs.tolist(), translations=trans.tolist())


def run_core(name, G, xs, certs, scratch, extra=None):
    t0 = time.time()
    rec = dict(core=name, G=list(G.mods), Q=G.qtype, K_order=int(len(G.K)), order=4 * G.n + 12)
    rec['core_check'] = core_ok(xs, G)
    bad = [np.array(x) for x in xs]
    bad[0][0] *= -1
    rec['core_negative_control_rejected'] = not core_ok(bad, G)['ok']
    if not rec['core_check']['ok']:
        rec['status'] = 'invalid core'
        return rec
    new, nrec = normalise(xs, G)
    rec['normalisation'] = nrec
    rec['normalised_core_check'] = core_ok(new, G)
    X, Rh, Ch = certs[G.qtype]
    H = assemble(new, G, X, Rh, Ch)
    path = os.path.join(scratch, 'H_%s.npy' % name.replace(' ', '_'))
    np.save(path, H)
    rec['H_sha256_int8'] = hashlib.sha256(np.ascontiguousarray(H).tobytes()).hexdigest()
    vout = path + '.verify.json'
    p = subprocess.run([sys.executable, '-I', os.path.join(HERE, '4p3_review_hverify.py'), path, '--out', vout],
                       capture_output=True, text=True)
    rec['verifier_exit'] = p.returncode
    rec['verifier'] = json.load(open(vout))
    # negative control at the matrix level: a corrupted copy must be rejected by the verifier
    Hb = H.copy()
    Hb[5, 17] *= -1
    pb = path.replace('.npy', '_corrupt.npy')
    np.save(pb, Hb)
    q = subprocess.run([sys.executable, '-I', os.path.join(HERE, '4p3_review_hverify.py'), pb], capture_output=True,
                       text=True)
    rec['corrupted_matrix_rejected'] = q.returncode != 0
    os.remove(pb)
    rec['status'] = 'Hadamard' if (p.returncode == 0 and rec['verifier']['ok'] and rec['corrupted_matrix_rejected']
                                   and rec['normalised_core_check']['ok']) else 'FAILED'
    if extra:
        rec.update(extra)
    rec['secs'] = round(time.time() - t0, 1)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--certs', default=os.path.join(ROOT, 'logs/openai_4p3/border_certificates.json'))
    ap.add_argument('--kumar', default='/home/user/vatskishan03/hadamard')
    ap.add_argument('--elliptic', default='', help='JSON list of extra cores {name, p, q, sequences} (q*p+k order)')
    ap.add_argument('--scratch', default=os.environ.get('SCRATCH', '/tmp'))
    ap.add_argument('--only', default='', help='comma list: 1916,3964,elliptic')
    ap.add_argument('--scramble', type=int, default=0, help='also this many random equivalent copies of 1916 and 3964')
    ap.add_argument('--scramble-seed', type=int, default=5)
    a = ap.parse_args()
    t0 = time.time()
    only = set(a.only.split(',')) if a.only else {'1916', '3964', 'elliptic'}
    data = json.load(open(a.certs))
    res = dict(script='theory/4p3_review_border.py', command=' '.join(sys.argv), seed=None, certificates=[], cores=[])
    certs = {}
    for c in data['certificates']:
        chk, mats = check_certificate(c)
        res['certificates'].append(chk)
        certs[c['q']] = mats
        print(json.dumps(chk), flush=True)
    if '1916' in only:
        G, xs = load_1916()
        res['cores'].append(run_core('1916', G, xs, certs, a.scratch))
        print(json.dumps({k: v for k, v in res['cores'][-1].items() if k not in ('verifier',)}), flush=True)
    if '3964' in only:
        G, xs, c = load_3964(a.kumar)
        rec = run_core('3964', G, xs, certs, a.scratch, extra=dict(
            kumar_sha256_field=c.get('sha256'), kumar_reported_coset_sums='a0: -11,-13,-13,-13; a1: -5,-7,-7,-7; '
                                                                         'a2: 1,-1,-1,-1; a3: 7,7,7,7'))
        res['cores'].append(rec)
        print(json.dumps({k: v for k, v in rec.items() if k not in ('verifier',)}), flush=True)
    if a.scramble:
        rng = np.random.default_rng(a.scramble_seed)
        pool = []
        if '1916' in only:
            pool.append(('1916',) + load_1916())
        if '3964' in only:
            G3, x3, _ = load_3964(a.kumar)
            pool.append(('3964', G3, x3))
        for nm, G, xs in pool:
            for r in range(a.scramble):
                ys, how = scramble(xs, G, rng)
                rec = run_core('%s scrambled %d' % (nm, r), G, ys, certs, a.scratch, extra=dict(scramble=how))
                res['cores'].append(rec)
                print(json.dumps({k: v for k, v in rec.items() if k not in ('verifier',)}), flush=True)
    if 'elliptic' in only and a.elliptic:
        for e in json.load(open(a.elliptic)):
            p, qt = e['p'], e['q']
            G = Group((2, 2, p) if qt == 'v4' else (4, p), qt)
            # input in q*p + k order, q = 2 x1 + x2 (v4) or x1 (z4): equals the mixed-radix index here
            xs = [np.array(s, np.int64) for s in e['sequences']]
            res['cores'].append(run_core(e['name'], G, xs, certs, a.scratch,
                                         extra=dict(source=e.get('source'), parameters=e.get('parameters'))))
            print(json.dumps({k: v for k, v in res['cores'][-1].items() if k not in ('verifier', 'parameters')}),
                  flush=True)
    res['all_certificates_ok'] = all(c['ok'] for c in res['certificates'])
    res['all_cores_hadamard'] = all(c.get('status') == 'Hadamard' for c in res['cores'])
    res['secs'] = round(time.time() - t0, 1)
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(dict(all_certificates_ok=res['all_certificates_ok'], cores=[(c['core'], c.get('status'))
                                                                               for c in res['cores']],
                          secs=res['secs'])))


if __name__ == '__main__':
    main()
