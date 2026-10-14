r"""Drivers for the CP-SAT border search (theory/ab_border_sat.py) [M].

Monomial compressions.  z_k = a_k delta_0 (+ constants, which do not matter on balanced frames), sum_k a_k^2 = 4b.
When one a_k is 0 the GS array of z (347 layout) is fixed by sigma acting on the core blocks as
    (sigma, rho sigma rho, rho sigma rho, sigma)  or the analogous split (rho = negation on Q = Z_c),
for EVERY permutation sigma of Q (theory: the type-2 blocks only force pi'_J = rho pi_I rho on the nonzero blocks, and
with one block zero these conditions are consistent for all sigma; checked by the model, which rejects a generator
that does not fix Cb).  The Singer-type ansatz takes sigma = a (c-1)-cycle fixing 0: V_Q is then the regular
representation of Z_(c-1), the frames are 4 regular orbits (border blocks of b = c - 1 rows), the corner is a 4 x 4
array of b x b circulants, as in the published c = 8 template of 359 (Z_7 = Singer cycle of GF(8)).
  --neg: the generator acts as -sigma on the core and negacyclically on the border (order 2(c - 1)).
usage: python3 theory/ab_border_sat_run.py --c 10 --a 4,4,2,0 [--neg] [--time 600] [--seed S] [--workers 4]
       python3 theory/ab_border_sat_run.py --c 8 --a 2,2,2,4 ...      (control: the 359 profile)"""
import sys, os, json, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_assemble as AS, ab_border_sat as BS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAYOUT = AS.load_template347()['layout']


def core(a, c, pos=None):
    z = np.zeros((4, c), np.int64)
    for k in range(4): z[k, 0 if pos is None else pos[k]] = a[k]
    return z, AS.compressed_core(z, (c,), LAYOUT)


def block_perms(a, c, sigma):
    """per-block point permutations (rows and columns) fixing the GS array of a_k delta_0: try all assignments of
    sigma / rho sigma rho to the 4 blocks; return the first that fixes Cb (as a 4c-point permutation)"""
    rho = np.array([(-x) % c for x in range(c)])
    rsr = rho[sigma[rho]]
    _, Cb = core(a, c)
    for pat in itertools.product((0, 1), repeat=4):
        p = np.concatenate([J * c + (sigma if pat[J] == 0 else rsr) for J in range(4)])
        if (Cb[np.ix_(p, p)] == Cb).all(): return p, pat        # Cb[p[r], p[s]] = Cb[r, s]
    return None, None


def singer_gen(a, c, neg=False):
    b = c - 1; nb = 4 * b; n = nb + 4 * c
    sigma = np.array([0] + [(x % b) + 1 for x in range(1, c)])        # 0 fixed, 1 -> 2 -> ... -> b -> 1
    p, pat = block_perms(a, c, sigma)
    if p is None: return None
    pr = np.zeros(n, np.int64); sr = np.ones(n, np.int64)
    for I in range(4):
        for i in range(b):
            pr[I * b + i] = I * b + (i + 1) % b
            if neg and i == b - 1: sr[I * b + i] = -1
    pr[nb:] = nb + p
    if neg: sr[nb:] = -1
    return (pr.tolist(), sr.tolist(), pr.tolist(), sr.tolist()), pat


def main():
    a_ = sys.argv[1:]
    opt = lambda k, d: a_[a_.index(k) + 1] if k in a_ else d
    c = int(opt('--c', '10')); a = [int(x) for x in opt('--a', '4,4,2,0').split(',')]
    neg = '--neg' in a_; tl = float(opt('--time', '600')); seed = int(opt('--seed', '0')); wk = int(opt('--workers', '4'))
    b = c - 1
    assert sum(x * x for x in a) == 4 * b, "monomial profile must have sum a_k^2 = 4b"
    z, Cb = core(a, c)
    g = singer_gen(a, c, neg)
    outdir = os.path.join(ROOT, 'logs/ab/border10'); os.makedirs(outdir, exist_ok=True)
    rec = dict(c=c, a=a, sym=f"Z{b}" + ("neg" if neg else ""), seed=seed, time_limit=tl)
    if g is None:
        rec['status'] = 'no block assignment fixes Cb'
    else:
        gen, pat = g; rec['block_pattern'] = list(pat)
        res, N = BS.solve(Cb, b, c, [gen], time_limit=tl, workers=wk, seed=seed)
        rec.update(res)
        if N is not None and res.get('verified'):
            f = os.path.join(outdir, f"template_c{c}_a{''.join(map(str, a))}_Z{b}{'neg' if neg else ''}_s{seed}.npz")
            np.savez_compressed(f, N=N.astype(np.int8), z=z, a=np.array(a), Cb=Cb)
            rec['saved'] = os.path.relpath(f, ROOT)
    with open(os.path.join(outdir, 'sat_runs.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
