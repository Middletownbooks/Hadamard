r"""c >= 10 borders: Hadamard completion of the h = 1 template by Douglas-Rachford [M].

A b >= 3 template for c is the same thing as a Hadamard matrix of order 4(2c - 1)
    N = [[C, K], [Kc^T, Cb]],   Cb = the 4c x 4c GS array (347 layout) of four +-1 sequences z_k on Z_c with
                                sum_k PAF_{Z_c}(z_k) = 4 off 0   (the h = 1 Q-level equation)
(report section 3.1: K, Kc are then automatically block-balanced frames and C = -K Cb^T Kc^T/(4c); every template
transfers to every h along the all-ones direction).  We search N with Cb fixed by Douglas-Rachford between
    A = {+-1 matrices with the Cb block fixed}   and   B = {sqrt(n) x orthogonal matrices}
(projection onto B = polar factor), as in Aragon Artacho-Borwein-Tam's DR searches for Hadamard matrices.
Every hit is checked exactly (N N^T = n I) and saved.
usage: python3 theory/ab_border_dr.py --c 10 [--z I] [--iters 20000] [--restarts 50] [--seed S] [--selftest]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS
from ab_border10 import qlevel_h1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAYOUT = AS.load_template347()['layout']


def core_of(z, c):
    return AS.compressed_core(np.asarray(z, np.int64), (c,), LAYOUT)


def polar(M):
    U, s, Vt = np.linalg.svd(M)
    return U @ Vt


def dr(Cb, nb, iters, rng, log_every=0, beta=1.0):
    """Douglas-Rachford for N = [[C, K], [Kc^T, Cb]] (nb = 4b border rows); returns N (+-1, Hadamard) or None"""
    m = Cb.shape[0]; n = nb + m; rt = np.sqrt(n)
    free = np.ones((n, n), bool); free[nb:, nb:] = False
    x = rng.standard_normal((n, n))
    x[nb:, nb:] = Cb
    best = 1e18
    for it in range(iters):
        pa = np.where(x >= 0, 1.0, -1.0); pa[nb:, nb:] = Cb
        r = 2 * pa - x
        pb = rt * polar(r)
        x = x + beta * (pb - pa)
        if it % 10 == 0:
            g = pa @ pa.T
            err = float(np.abs(g - n * np.eye(n)).sum())
            best = min(best, err)
            if err == 0:
                return pa.astype(np.int64), it, best
    return None, iters, best


def main():
    a = sys.argv[1:]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    c = int(opt('--c', '10')); iters = int(opt('--iters', '20000')); restarts = int(opt('--restarts', '50'))
    seed = int(opt('--seed', '0')); zi = opt('--z', None)
    b = c - 1; nb = 4 * b; n = nb + 4 * c
    rng = np.random.default_rng(seed)
    outdir = os.path.join(ROOT, 'logs/ab/border10'); os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, f'dr_c{c}.jsonl')
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    if '--selftest' in a:
        assert c == 6
        T = AS.load_template347(); comp = T['comp']
        z = np.zeros_like(comp)
        for k in range(4):
            lo_ = comp[k].min(); z[k] = np.where(comp[k] == lo_, -1, 1) if len(set(comp[k])) > 1 else -np.ones(6, int)
        Cb = AS.compressed_core(z, (2, 3), T['layout'])
        zs = [('347', Cb)]
    else:
        X, key, sols = qlevel_h1(c)
        zs = []
        for i, (k0, k1, k2, k3) in enumerate(sols):
            if zi is not None and i != int(zi): continue
            z = np.array([X[key[k][0]] for k in (k0, k1, k2, k3)])
            zs.append((i, core_of(z, c)))
    for (name, Cb) in zs:
        P = Cb @ Cb.T
        for r in range(restarts):
            t0 = time.time()
            N, it, best = dr(Cb, nb, iters, rng)
            rec = dict(c=c, z=name, restart=r, found=N is not None, iters=it, best_err=round(best, 1), secs=round(time.time() - t0, 1))
            if N is not None:
                ok = bool((N @ N.T == n * np.eye(n, dtype=np.int64)).all())
                rec['hadamard'] = ok
                np.save(os.path.join(outdir, f'template_c{c}_z{name}_r{r}.npy'), N.astype(np.int8))
            log(rec)
            if N is not None: break


if __name__ == "__main__":
    main()
