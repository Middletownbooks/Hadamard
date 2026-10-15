r"""H(3668) = 4 * 917 by Cooper-Wallis: T-matrices of order 131 (cyclotomic, K of order 5) x Williamson matrices of
order 7, written as a b = 0 Goethals-Seidel core of four +-1 circulants of length 917 [M].

[P] (Cooper-Wallis.)  T_1..T_4: {0,+-1} circulants of order t with disjoint supports covering Z_t and
sum_i PAF(T_i) = t delta; W_1..W_4 (= A, B, C, D): symmetric +-1 circulants of order w with sum_i PAF(W_i) = 4w delta.
On G = Z_t x Z_w put
    x_1 = T1 A + T2 B + T3 C + T4 D,     x_2 = -T1 B + T2 A + T3 D - T4 C,
    x_3 = -T1 C - T2 D + T3 A + T4 B,    x_4 = -T1 D + T2 C - T3 B + T4 A      (products as functions on Z_t x Z_w).
The x_k are +-1 (exactly one T_i is nonzero at each point of Z_t), and sum_k PAF_G(x_k) = 4tw delta: the cross terms
cancel because each W_i is symmetric (PAF-type cross-correlations CC(W_i, W_j) are symmetric in the lag, and they
enter with opposite signs).  With gcd(t, w) = 1, Z_t x Z_w = Z_tw by the CRT, so the x_k are circulants of length tw
and the Goethals-Seidel array gives H(4tw).
T-matrices of order t are equivalently a b = 0 GS core of order t with the product condition a_j b_j c_j d_j = +1
(report/tmat_routes.md 2).  At t = 131 the K-invariant ansatz with K the order-5 subgroup of F_131^* (27 orbits,
row sums (0, +-5, +-5, +-9)) is row-sum feasible and has log2 E = 30.0 (scripts/tmat_estimate.py 131 5); the support
split + k-tree join (scripts/tmat_kt.py 131 5) finds one in seconds.

usage: python3 theory/b917_cw.py [--T logs/b917/T131_k5.npy | --search SEED] [--out FILE]
(default FILE: logs/b917/seq_H3668_Z917.txt, or logs/b917/seq_H3668_Z917_search<SEED>.txt with --search; --search
also saves its T-matrices as logs/b917/T131_k5_search<SEED>.npy)"""
import sys, os, json, itertools, subprocess, hashlib
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T_, W_ = 131, 7


def paf_fft(x):
    """periodic autocorrelation over Z_n1 x ... (exact for integer x)"""
    F = np.fft.fftn(np.asarray(x, float))
    return np.rint(np.fft.ifftn(np.abs(F) ** 2).real).astype(np.int64)


def check_T(T):
    t = T.shape[1]
    assert set(np.unique(T).tolist()) <= {-1, 0, 1} and (np.abs(T).sum(0) == 1).all(), "supports must partition Z_t"
    P = sum(paf_fft(r) for r in T)
    assert P[0] == t and (P[1:] == 0).all(), "sum PAF(T_i) != t delta"
    return True


def williamson(w):
    """all symmetric +-1 sequences of length w (first entry +1); first quadruple (sorted) with sum PAF = 4w delta"""
    h = (w - 1) // 2
    syms = []
    for bits in itertools.product((1, -1), repeat=h):
        syms.append(np.array([1] + list(bits) + list(bits[::-1]), np.int64))
    pafs = [paf_fft(x) for x in syms]
    for q in itertools.combinations_with_replacement(range(len(syms)), 4):
        tot = sum(pafs[i] for i in q)
        if (tot[1:] == 0).all() and tot[0] == 4 * w:
            return np.array([syms[i] for i in q])
    raise ValueError("no Williamson quadruple of order %d" % w)


def cooper_wallis(T, W):
    """the four +-1 functions x_k on Z_t x Z_w (shape (4, t, w))"""
    T1, T2, T3, T4 = T; A, B, C, D = W
    o = np.outer
    X = np.array([o(T1, A) + o(T2, B) + o(T3, C) + o(T4, D),
                  -o(T1, B) + o(T2, A) + o(T3, D) - o(T4, C),
                  -o(T1, C) - o(T2, D) + o(T3, A) + o(T4, B),
                  -o(T1, D) + o(T2, C) - o(T3, B) + o(T4, A)], np.int64)
    assert set(np.unique(X).tolist()) == {-1, 1}
    return X


def to_cyclic(X):
    """CRT Z_tw -> Z_t x Z_w, n -> (n mod t, n mod w)"""
    t, w = X.shape[1:]; n = np.arange(t * w)
    return X[:, n % t, n % w]


def search_T(seed):
    log = os.path.join(ROOT, 'logs/b917', f'tmat_kt_131_k5_s{seed}.log')
    with open(log, 'w') as fh:
        subprocess.run([sys.executable, 'scripts/tmat_kt.py', '131', '5', '--nS', '200', '--seed', str(seed), '--buckets', '200'],
                       cwd=ROOT, stdout=fh, check=True)
    recs = [json.loads(l) for l in open(log) if l.startswith('{"S_index"')]
    hit = [r for r in recs if r.get('FOUND')]
    assert hit, "no T-matrices found"
    return np.array(hit[0]['T'], np.int64)


def main():
    a = sys.argv[1:]
    opt = lambda n, d=None: a[a.index(n) + 1] if n in a else d
    if '--search' in a:
        sd = int(opt('--search'))
        T = search_T(sd)
        np.save(os.path.join(ROOT, f'logs/b917/T131_k5_search{sd}.npy'), T)
        default_out = f'logs/b917/seq_H3668_Z917_search{sd}.txt'
    else:
        T = np.load(os.path.join(ROOT, opt('--T', 'logs/b917/T131_k5.npy'))).astype(np.int64)
        default_out = 'logs/b917/seq_H3668_Z917.txt'
    check_T(T)
    W = williamson(W_)
    X = cooper_wallis(T, W)
    P = sum(paf_fft(x) for x in X)                      # PAF over Z_131 x Z_7
    ok_G = bool(P[0, 0] == 4 * T_ * W_ and (P.ravel()[1:] == 0).all())
    Xc = to_cyclic(X)
    Pc = sum(paf_fft(x) for x in Xc)
    ok_c = bool(Pc[0] == 4 * T_ * W_ and (Pc[1:] == 0).all())
    out = os.path.join(ROOT, opt('--out', default_out))
    rs = Xc.sum(1).tolist()
    with open(out, 'w') as fh:
        fh.write("# Four +-1 circulant sequences of length 917 (cyclic Z_917), b = 0 Goethals-Seidel core of H(3668)\n")
        fh.write(f"# sum_k PAF_k(s) = 0 for s != 0; row sums (A,B,C,D) = ({','.join(map(str, rs))})\n")
        for lab, x in zip('ABCD', Xc): fh.write(f"{lab}=" + ''.join('+' if v > 0 else '-' for v in x) + "\n")
    rec = dict(T_row_sums=T.sum(1).tolist(), T_weights=np.abs(T).sum(1).tolist(), williamson7=W.tolist(),
               paf_ok_Z131xZ7=ok_G, paf_ok_Z917=ok_c, row_sums=rs, sum_sq_row_sums=int(sum(r * r for r in rs)),
               out=os.path.relpath(out, ROOT),
               sha256_T=hashlib.sha256(T.astype(np.int8).tobytes()).hexdigest()[:16])
    print(json.dumps(rec), flush=True)
    return ok_G and ok_c


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
