r"""Calibration of the z-conditioned stage-1 count (theory/s1z_lib.py) against exact counts [M vs E].

Per configuration: every class multiset of compression quadruples (the whole Q-level solution set up to symmetry),
its exact lag-0 count N0 and the number of ordered compression quadruples it stands for; per multiset the Gaussian
estimate E = N0 * 2^-bits (remaining lag classes; raw and with the finite-sample log-det correction) and, where the
block lists fit, the exact number of stage-1 solutions with those compressions.  The exact solutions are also
classified by symmetry: a block is "symmetric" if an affine map g -> phi(g) + t of G' (phi a per-axis unit scalar,
t in G') that moves some U-orbit fixes it.
Output JSONL, one line per multiset, then a summary line: totals over all ordered compression quadruples.
usage: python3 theory/s1z_calib.py CFG(.json|json) OUT.jsonl [--N 2000] [--exact] [--max-ms 400] [--seed 1]
                                   [--max-pairs 6e8]"""
import sys, os, json, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import ZOrb


def affine_maps(Z):
    """affine maps of G' (per-axis unit scalars, any translation) that move some U-orbit, as index permutations"""
    Gp = Z.Gp; mods = Gp.mods
    opts = [[u for u in range(1, m) if math.gcd(u, m) == 1] if m > 1 else [0] for m in mods]
    maps = []
    borb = Z.borb
    for us in itertools.product(*opts):
        A = Gp.act(us)
        for t in range(Gp.n):
            P = Gp.idx(Gp.co[A] + Gp.co[t])           # g -> phi(g) + t
            if (borb[P] == borb).all(): continue       # acts trivially on U-invariant functions
            # the map must send U-orbits to U-orbits (otherwise it never fixes a U-invariant y in a useful way)
            maps.append(P)
    return maps


def sym_share(Z, Ys, maps):
    nall = 0; nblk = 0
    for Y in Ys:
        flags = []
        for k in range(4):
            y = Y[k][Z.borb]
            flags.append(any((y[P] == y).all() for P in maps))
        nblk += sum(flags); nall += all(flags)
    return nall, nblk


def main():
    a = sys.argv[1:]
    src = a[0]; outp = a[1]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    N = int(opt('--N', '2000')); do_exact = '--exact' in a; maxms = int(opt('--max-ms', '400'))
    seed = int(opt('--seed', '1')); max_pairs = int(float(opt('--max-pairs', '6e8')))
    cfg = SD.cfg_from_json(open(src).read() if src.endswith('.json') else src)
    t0 = time.time()
    Z = ZOrb(cfg)
    rows = Z.rows(); C = Z.compressions(rows); cls, idx = Z.classes(C); qd = Z.class_quads(cls)
    cache = {}
    for q in qd: q['n0'] = Z.n0(cls, q['m'], cache)
    qd = [q for q in qd if q['n0'] > 0]
    head = dict(stage='setup', cfg={k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()}, L1=Z.L1,
                rest_lags=len(Z.rest), compressions=len(C), classes=len(cls), multisets=len(qd),
                ordered_zquads=sum(q['ordered'] for q in qd),
                log2_N0_total=round(math.log2(sum(q['ordered'] * q['n0'] for q in qd)), 3) if qd else None,
                secs=round(time.time() - t0, 1))
    print(json.dumps(head), flush=True)
    fo = open(outp, 'w'); fo.write(json.dumps(head) + "\n"); fo.flush()
    if len(qd) > maxms:
        rng = np.random.default_rng(seed)
        w = np.array([q['ordered'] * q['n0'] for q in qd], float); w /= w.sum()
        pick = sorted(set(rng.choice(len(qd), size=maxms, replace=True, p=w).tolist()))
    else:
        pick = list(range(len(qd)))
    maps = affine_maps(Z) if do_exact else None
    tot_e = 0.0; tot_e_raw = 0.0; tot_x = 0; n_ex = 0; tot_n0 = 0; nsym = 0; nsol = 0
    for t in pick:
        q = qd[t]
        zq = [cls[i]['rep'] for i in q['m']]
        r = Z.estimate(zq, N=N, seed=seed + t)
        rec = dict(stage='ms', i=t, m=list(q['m']), z=[list(z) for z in zq], ordered=q['ordered'],
                   log2N0=round(math.log2(q['n0']), 3), **{k: r.get(k) for k in ('bits', 'bits_raw', 'rank', 'mahal', 'log2E', 'in_lattice')})
        if r.get('log2E') is not None and np.isfinite(r['log2E']):
            tot_e += q['ordered'] * 2.0 ** r['log2E']
            tot_e_raw += q['ordered'] * 2.0 ** (math.log2(q['n0']) - r['bits_raw'])
        tot_n0 += q['ordered'] * q['n0']
        if do_exact:
            try:
                ex, sizes, Ys = Z.exact_count(zq, max_pairs=max_pairs)
                ns, nb = sym_share(Z, Ys, maps)
                rec.update(exact=ex, log2_exact=round(math.log2(ex), 3) if ex else None, block_lists=sizes,
                           sym_all=ns, sym_blocks=nb)
                tot_x += q['ordered'] * ex; n_ex += 1; nsym += q['ordered'] * ns; nsol += q['ordered'] * ex
            except (RuntimeError, MemoryError) as e:
                rec['exact_skipped'] = repr(e)[:120]
        rec['secs'] = round(time.time() - t0, 1)
        print(json.dumps(rec), flush=True); fo.write(json.dumps(rec) + "\n"); fo.flush()
    summ = dict(stage='summary', multisets_evaluated=len(pick), multisets=len(qd), sampled=len(qd) > maxms,
                log2_total_est=round(math.log2(tot_e), 3) if tot_e > 0 else None,
                log2_total_est_raw=round(math.log2(tot_e_raw), 3) if tot_e_raw > 0 else None,
                log2_N0_evaluated=round(math.log2(tot_n0), 3) if tot_n0 else None)
    if do_exact:
        summ.update(exact_multisets=n_ex, total_exact=tot_x, log2_total_exact=round(math.log2(tot_x), 3) if tot_x else None,
                    sym_share=round(nsym / nsol, 4) if nsol else None)
    summ['secs'] = round(time.time() - t0, 1)
    print(json.dumps(summ), flush=True); fo.write(json.dumps(summ) + "\n"); fo.close()


if __name__ == '__main__':
    main()
