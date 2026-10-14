"""Write a pub/H<N>/ bundle for a b >= 3 bordered group-developed Goethals-Seidel Hadamard matrix (format of
pub/H2092/, with the border and corner written out).  Shares no code with theory/.
    python3 scripts/make_bundle_bge3.py CORE.npy MODS NQ LAYOUT.json K.npy KC.npy C.npy OUTDIR ["description"]
CORE.npy: 4 x m (+-1) in the mixed-radix order of G = Z_{mods[0]} x ... (last axis fastest); Q = the first NQ axes.
Writes OUTDIR/core_G.txt, OUTDIR/border.json, OUTDIR/H<N>_normalized.txt.gz, OUTDIR/verify.py, runs verify.py."""
import gzip, json, os, sys, subprocess
import numpy as np

VERIFY = r'''"""Standalone verification of H(%(N)d), b = %(b)d, c = %(c)d.  Needs only Python 3 and numpy.

  python3 verify.py     # (1) the core: four +-1 functions on G = %(Gname)s with
                        #     sum_k PAF_k = 4m at 0, -4b on H minus 0, +4 off H (H = the subgroup with Q-coordinates 0)
                        # (2) build H(%(N)d): [[C, K inflated], [Kc^T inflated, GS array of the core]] (border.json)
                        #     and check H H^T = %(N)d I exactly (integer arithmetic)
                        # (3) check the published normalized matrix the same way, and that it is the normalization of (2)
"""
import gzip, json, itertools, numpy as np

def group(mods):
    co = np.array(list(itertools.product(*[range(m) for m in mods])), dtype=np.int64)
    M = np.array(mods, dtype=np.int64)
    w = np.array([int(np.prod(mods[i + 1:])) for i in range(len(mods))], dtype=np.int64)
    idx = lambda c: ((c %% M) * w).sum(-1)
    return co, idx

def load():
    L = [l.strip() for l in open('core_G.txt') if l.strip() and not l.startswith('#')]
    hdr = json.loads(L[0]); rows = [l.split('=', 1)[1] for l in L[1:5]]
    X = np.array([[1 if ch == '+' else -1 for ch in r] for r in rows], dtype=np.int64)
    return hdr, X

def paf_sum(X, co, idx):
    m = X.shape[1]; P = np.zeros(m, dtype=np.int64)
    for t in range(m):
        sh = idx(co + co[t]); P[t] = int((X * X[:, sh]).sum())
    return P

def is_hadamard(H):
    N = H.shape[0]; H = H.astype(np.int64)
    for i in range(0, N, 256):
        G = H[i:i + 256] @ H.T
        E = np.zeros_like(G); E[np.arange(G.shape[0]), np.arange(i, i + G.shape[0])] = N
        if not np.array_equal(G, E): return False
    return True

def build(hdr, X, B):
    mods, nQ, b = hdr['mods'], hdr['nQ'], hdr['b']
    co, idx = group(mods); m = len(co)
    Qm = mods[:nQ]; wq = np.array([int(np.prod(Qm[i + 1:])) for i in range(nQ)], dtype=np.int64)
    qg = (co[:, :nQ] * wq).sum(1) if nQ else np.zeros(m, dtype=np.int64)
    K, Kc, C = (np.array(B[k], dtype=np.int64) for k in ('K', 'Kc', 'C'))
    c = K.shape[1] // 4; off = 4 * b; N = off + 4 * m
    H = np.zeros((N, N), dtype=np.int64)
    H[:off, :off] = C
    for J in range(4):
        H[:off, off + J * m: off + (J + 1) * m] = K[:, J * c + qg]
        H[off + J * m: off + (J + 1) * m, :off] = Kc[:, J * c + qg].T
    sig = {'h-g': lambda g, h: idx(h - g), 'g+h': lambda g, h: idx(g + h), 'g-h': lambda g, h: idx(g - h),
           '-g-h': lambda g, h: idx(-g - h)}
    for key, L in B['layout'].items():
        I, J = int(key[0]), int(key[1])
        S = sig[L['sigma']](co[:, None, :], co[None, :, :])
        H[off + I * m: off + (I + 1) * m, off + J * m: off + (J + 1) * m] = L['sign'] * X[L['k']][S]
    return H

if __name__ == '__main__':
    hdr, X = load(); B = json.load(open('border.json'))
    co, idx = group(hdr['mods']); m = len(co); b = hdr['b']; nQ = hdr['nQ']
    P = paf_sum(X, co, idx)
    inH = (co[:, :nQ] == 0).all(1)
    T = np.where(inH, -4 * b, 4); T[0] = 4 * m
    print('core: |G| =', m, ' sum PAF = target:', bool((P == T).all()), ' row sums:', X.sum(1).tolist())
    assert (P == T).all()
    H = build(hdr, X, B); N = H.shape[0]
    print('built H of order', N, ' entries', sorted(set(np.unique(H).tolist())), ' Hadamard:', is_hadamard(H))
    rows = [l.strip() for l in gzip.open('H%(N)d_normalized.txt.gz', 'rt') if l.strip()]
    Hn = np.array([[1 if ch == '+' else -1 for ch in r] for r in rows], dtype=np.int64)
    Hb = H * H[:, [0]]; Hb = Hb * Hb[[0], :]
    print('published normalized matrix: order', Hn.shape, ' first row/col all +1:', bool((Hn[0] == 1).all() and (Hn[:, 0] == 1).all()),
          ' Hadamard:', is_hadamard(Hn), ' equals normalization of built H:', bool(np.array_equal(Hn, Hb)))
'''


def main(core, mods, nQ, layout, K, Kc, C, outdir, desc=""):
    X = np.load(core).astype(np.int64); mods = [int(x) for x in mods.split(',')]; nQ = int(nQ)
    L = json.load(open(layout)); K = np.load(K).astype(int); Kc = np.load(Kc).astype(int); C = np.load(C).astype(int)
    b = K.shape[0] // 4; c = K.shape[1] // 4; m = X.shape[1]; v = m + b; N = 4 * v
    os.makedirs(outdir, exist_ok=True)
    Gname = ' x '.join(f'Z_{x}' for x in mods)
    with open(os.path.join(outdir, 'core_G.txt'), 'w') as fh:
        fh.write(f"# b = {b} core of H({N}): four +-1 functions on G = {Gname} (|G| = {m}), mixed-radix order, last axis fastest\n")
        fh.write(f"# H = the subgroup with the first {nQ} coordinates 0 (index c = {c}); sum_k PAF_k = 4m at 0, -4b on H\\0, +4 off H\n")
        for line in desc.split('\n'):
            if line: fh.write(f"# {line}\n")
        fh.write(json.dumps(dict(mods=mods, nQ=nQ, b=b, c=c, v=v, N=N)) + "\n")
        for lab, x in zip('ABCD', X): fh.write(f"{lab}=" + ''.join('+' if t > 0 else '-' for t in x) + "\n")
    lay = {k: dict(k=int(v_['k']), sigma=v_['sigma'], sign=int(v_['rho_sign'] * v_['gam_sign'])) for k, v_ in L.items()}
    json.dump(dict(layout=lay, K=K.tolist(), Kc=Kc.tolist(), C=C.tolist(),
                   note="H = [[C, K inflated], [Kc^T inflated, core GS array]]; inflation: column J*m + g of the core part "
                        "uses column J*c + q(g) of K (q(g) = Q-coordinates of g, mixed radix); core block (I, J)[g, h] = "
                        "sign * x_k(sigma(g, h))"), open(os.path.join(outdir, 'border.json'), 'w'))
    open(os.path.join(outdir, 'verify.py'), 'w').write(VERIFY % dict(N=N, b=b, c=c, Gname=Gname))
    # normalized matrix from the standalone builder itself
    sys.path.insert(0, outdir)
    ns = {}
    exec(compile(open(os.path.join(outdir, 'verify.py')).read().split("if __name__ == '__main__':")[0], 'verify', 'exec'), ns)
    cwd = os.getcwd(); os.chdir(outdir)
    hdr, Xl = ns['load'](); H = ns['build'](hdr, Xl, json.load(open('border.json')))
    os.chdir(cwd)
    Hn = H * H[:, [0]]; Hn = Hn * Hn[[0], :]
    with gzip.open(os.path.join(outdir, f'H{N}_normalized.txt.gz'), 'wt') as fh:
        for r in Hn: fh.write(''.join('+' if t > 0 else '-' for t in r) + "\n")
    print(subprocess.run([sys.executable, 'verify.py'], cwd=outdir, capture_output=True, text=True).stdout)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(*a[:8], desc=a[8] if len(a) > 8 else "")
