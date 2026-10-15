r"""Exact decision of small Q-automorphism GP slices (theory/ab_qaut.py), end to end [M].

Per configuration: exact stage 1 on G' (per-block letter vectors under the lag-0 energy, the Q-level bound
|z^(chi)|^2 <= 4b and the row set of screen A's row equation; 2+2 hashed join on the full PAF over G', every match
re-checked), the monomial border of each solution (theory/ab_mono_border.border_for; Singer templates for c = 4, 6, 8),
exact stage 2 per bordered y (theory/ab_small_exact.solve_y, block fibres <= MAXFIB), assembly and H H^T = 4v I
(theory/ab_assemble, 347 layout).  Restart-safe JSONL, one line per configuration.
usage: python3 theory/ab_qaut_run.py OUT.jsonl [--cs 4,8] [--vs ...] [--maxfib 65536] [--maxblock 2000000]
                                     [--part i/W] [--secs 1800] [--memgb 3.0]"""
import sys, os, json, time, signal, resource, itertools, hashlib
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_gauss as AGs, ab_mono_border as MB, ab_assemble as AS
from ab_qaut import build
from ab_lib import factorint, AG
from ab_sweep import abelian_types, unit_subgroup_gens
from ab_small_exact import solve_y

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN = [571, 599, 643, 659, 719, 739, 751, 823, 839, 859, 863, 883, 907, 919, 947, 955, 971, 991]


def block_vectors(sl, cap):
    """all y on the G'-orbits with block energy <= T1(0) - 3 * min energy, as an array (n x L1)"""
    w = sl.bsz.astype(np.int64); L = [np.array(x, np.int64) for x in sl.letters]; nb = len(L)
    T0 = int(sl.T1[0]); emin = [int(w[j] * (Lj ** 2).min()) for j, Lj in enumerate(L)]
    Emax = T0 - 3 * sum(emin)
    if Emax < sum(emin): return np.zeros((0, nb), np.int64)
    rest = np.cumsum(emin[::-1])[::-1].tolist() + [0]
    cur = np.zeros((1, 0), np.int64); E = np.zeros(1, np.int64)
    for j in range(nb):
        opts = L[j]; e = w[j] * opts ** 2
        newE = (E[:, None] + e[None, :]).ravel(); ok = newE + rest[j + 1] <= Emax
        cur = np.concatenate([np.repeat(cur, len(opts), 0), np.tile(opts, len(cur))[:, None]], 1)[ok]
        E = newE[ok]
        if len(cur) > cap: raise RuntimeError(f"block list > {cap}")
    return cur


def stage1(sl, cap, Qmods, b):
    Y0 = block_vectors(sl, cap)
    if len(Y0) == 0: return [], dict(blocks=0)
    full = Y0[:, sl.borb]
    # Q-level bound: |z^(chi)|^2 <= 4b for chi != 1 (z = Q-compression of y)
    nQ = len(Qmods); Qg = AG(Qmods)
    qidx = Qg.idx(sl.Gp.co[:, :nQ])
    Z = np.zeros((len(full), Qg.n), np.int64)
    for qq in range(Qg.n): Z[:, qq] = full[:, qidx == qq].sum(1)
    Zh = np.abs(np.fft.fftn(Z.reshape((-1,) + tuple(Qmods)).astype(float), axes=tuple(range(1, nQ + 1)))) ** 2
    keep = (Zh.reshape(len(Z), -1)[:, 1:] <= 4 * b + 1e-6).all(1)
    full = full[keep]; Y0 = Y0[keep]
    if len(full) == 0: return [], dict(blocks=0)
    # canonical sign: first nonzero entry positive (blocks may be negated)
    P = np.rint(sl.Gp.paf(full)).astype(np.int64)
    T1 = sl.T1
    rng = np.random.default_rng(7); wv = rng.integers(1, 2 ** 61, size=P.shape[1], dtype=np.int64)
    h = P @ wv
    n = len(h)
    if n * (n + 1) // 2 > 6 * 10 ** 7: raise RuntimeError(f"pair list too large ({n} block vectors)")
    I, J = np.triu_indices(n)
    ph = h[I] + h[J]
    o = np.argsort(ph); ps = ph[o]
    need = np.int64(T1 @ wv) - ph
    lo = np.searchsorted(ps, need, 'left'); hi = np.searchsorted(ps, need, 'right')
    sols = set()
    for p in np.nonzero(hi > lo)[0]:
        for qi in o[lo[p]:hi[p]]:
            if qi < p: continue
            if not (P[I[p]] + P[J[p]] + P[I[qi]] + P[J[qi]] == T1).all(): continue
            sols.add(tuple(sorted((int(I[p]), int(J[p]), int(I[qi]), int(J[qi])))))
    out = [np.array([full[i] for i in s_]) for s_ in sorted(sols)]
    return out, dict(blocks=int(n))


def configs(vs, cs):
    for v in vs:
        for c in cs:
            if (v + 1) % c: continue
            h = (v + 1) // c - 1
            qks = {4: ['z2sq3'], 8: ['singer7', 'gl3_3']}[c]
            for q, e in factorint(h).items():
                if e != 1 or q < 5: continue
                for Ht in abelian_types(h // q):
                    Hopts = [()]
                    for ai, n_ in enumerate(Ht):
                        for (o, u) in unit_subgroup_gens(n_): Hopts.append(((ai, u),))
                    for s in sorted(d for d in range(2, q) if (q - 1) % d == 0):
                        for qk in qks:
                            for Hm in Hopts:
                                yield dict(v=v, c=c, H1=list(Ht), q=q, s=s, f=(q - 1) // s, qaut=qk, Hmult=[list(x) for x in Hm])


class TO(Exception): pass


def run_one(cfg, maxfib, maxblock, coredir, log):
    sl = build(cfg['v'], cfg['c'], tuple(cfg['H1']), cfg['q'], cfg['s'], cfg['qaut'], [tuple(x) for x in cfg['Hmult']])
    rec = dict(orbits=int(sl.r), L1=len(sl.breps))
    Qmods = {4: (2, 2), 8: (2, 2, 2)}[cfg['c']]; b = cfg['c'] - 1
    sols, info = stage1(sl, maxblock, Qmods, b)
    rec.update(stage1_blocks=info['blocks'], stage1_solutions=len(sols))
    reps, _ = AGs.lag_reps(sl.G, sl.perms)
    gw = np.random.default_rng(5).integers(1, 2 ** 62, size=len(reps), dtype=np.int64)
    nQ = len(Qmods); Qg = AG(Qmods); qidx = Qg.idx(sl.Gp.co[:, :nQ])
    nb = 0; nund = 0; ncore = 0; nhad = 0
    for Y in sols:
        z = np.zeros((4, Qg.n), np.int64)
        for qq in range(Qg.n): z[:, qq] = Y[:, qidx == qq].sum(1)
        br = MB.border_for(z, Qmods)
        if br is None: continue
        nb += 1
        hits, m = solve_y(sl, Y, reps, gw, maxfib)
        if hits is None: nund += 1; continue
        for Xc in hits:
            ncore += 1
            ok = bool((sl.G.paf(Xc).sum(0) == sl.T).all())
            Zc = AS.compressions(Xc, sl.G.mods, nQ)
            Cf = AS.forced_corner(Zc, Qmods, MB.LAYOUT, br['K'], br['Kc'])
            H = AS.assemble(Xc, sl.G.mods, nQ, MB.LAYOUT, br['K'], br['Kc'], Cf) if Cf is not None else None
            had = H is not None and AS.is_hadamard(H)
            nhad += had
            if ncore <= 20:
                tag = f"v{cfg['v']}_c{cfg['c']}_q{cfg['q']}_s{cfg['s']}_{cfg['qaut']}_{ncore}"
                np.save(os.path.join(coredir, f'core_{tag}.npy'), Xc.astype(np.int8))
                if had: np.savez_compressed(os.path.join(coredir, f"H{4 * cfg['v']}_{tag}.npz"), H=np.packbits(H > 0, axis=1))
                log(dict(CORE=tag, core_ok=ok, hadamard=bool(had), order=4 * cfg['v']))
    rec.update(bordered=nb, stage2_undecided=nund, cores=ncore, hadamard=nhad)
    return rec


def main():
    a = sys.argv[1:]; outp = a[0]
    opt = lambda k, d: a[a.index(k) + 1] if k in a else d
    vs = [int(x) for x in opt('--vs', ','.join(map(str, OPEN))).split(',')]
    cs = [int(x) for x in opt('--cs', '4,8').split(',')]
    maxfib = int(opt('--maxfib', '65536')); maxblock = int(opt('--maxblock', '2000000'))
    part, nW = (int(x) for x in opt('--part', '0/1').split('/'))
    secs = int(opt('--secs', '1800')); memgb = float(opt('--memgb', '3.0'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30),) * 2)
    coredir = os.path.join(ROOT, 'logs/ab/qaut/cores'); os.makedirs(coredir, exist_ok=True)
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.dumps(json.loads(l)['cfg']))
            except Exception: pass
    fo = open(outp, 'a')
    def log(rec):
        fo.write(json.dumps(rec) + "\n"); fo.flush(); print(json.dumps(rec), flush=True)
    def handler(signum, frame): raise TO()
    signal.signal(signal.SIGALRM, handler)
    for i, cfg in enumerate(configs(vs, cs)):
        if i % nW != part or json.dumps(cfg) in done: continue
        rec = dict(cfg=cfg); t0 = time.time(); signal.alarm(secs)
        try: rec.update(run_one(cfg, maxfib, maxblock, coredir, log))
        except TO: rec['timeout'] = secs
        except (Exception, MemoryError) as e: rec['error'] = repr(e)[:160]
        signal.alarm(0); rec['secs'] = round(time.time() - t0, 1)
        log(rec)


if __name__ == "__main__":
    main()
