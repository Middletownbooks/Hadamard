r"""Paley-type c = 4 (b = 3) and c = 6 (b = 5) slices of every open order, decided exactly [M].

Configuration: v with 4 | v + 1, h = (v + 1)/4 - 1, Q = Z_2^2, a prime q || h with the GP index s on Z_q (s = 2: the
Paley classes; s = 3, 4 also tried when f >= 2), H' = Z_{h/q} split into its prime-power cyclic factors, and U = one
cyclic multiplier group on one cyclic factor of H' (all of them), or none.  Only configurations with at most
`--max-orbits` G'-orbits are taken (the block fibres of stage 2 are then at most 2^(G'-orbits) for s = 2).
Pipeline per configuration (restart-safe, one JSONL line each, logs/ab/paley_c4.jsonl):
  screen A (exact Q-level) -> exact stage 1 (theory/ab_s1_direct.py; skipped above `--max-vectors` block vectors)
  -> exact stage 2 for every stage-1 solution (theory/ab_small_exact.py, fibre cap 2^16)
  -> for every core: the c = 4 border (theory/ab_border_c4.py) and the full H(4v), checked; saved under logs/ab/paley_c4/.
c = 6 (--c 6): Q = Z_2 x Z_3, stage-1 solutions are kept only if their Q-compression is bordered by the 347 template
(theory/ab_borderfit.fits); cores are completed with that frame and the forced corner (theory/ab_assemble.py).
usage: python3 theory/ab_paley_c4.py [--c 4|6] [--vs 751,823,...] [--max-orbits 16] [--max-vectors 20000000]"""
import sys, os, json, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD
import ab_s1_direct as D, ab_small_exact as SE
from ab_lib import factorint

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN = [571, 599, 643, 659, 719, 739, 751, 823, 839, 859, 863, 883, 907, 917, 919, 947, 955, 971, 991]


def configs(v, max_orbits, c=4):
    if (v + 1) % c: return []
    h = (v + 1) // c - 1; out = []
    Qm = [2, 2] if c == 4 else [2, 3]
    for q, e in factorint(h).items():
        if e != 1 or q < 5: continue
        hp = h // q
        factors = [p ** k for p, k in sorted(factorint(hp).items())] if hp > 1 else []
        for s in (2, 3, 4):
            if (q - 1) % s or (q - 1) // s < 2: continue
            f = (q - 1) // s
            Uopts = [[]]
            for ai, n in enumerate(factors):
                for (o, u) in W.unit_subgroup_gens(n): Uopts.append([[2 + ai, o, u]])
            for U in Uopts:
                cfg = dict(v=v, c=c, b=c - 1, h=h, m=c * h, Q=Qm, H1=factors, q=q, s=s, f=f,
                           mods=Qm + factors + [q], nQ=2, U=U)
                try: sl = W.build(SD.cfg_from_json(json.dumps(cfg)))
                except Exception: continue
                if len(sl.breps) <= max_orbits: out.append((cfg, sl))
    return out


def main():
    a = sys.argv[1:]
    vs = [int(x) for x in a[a.index('--vs') + 1].split(',')] if '--vs' in a else OPEN
    mo = int(a[a.index('--max-orbits') + 1]) if '--max-orbits' in a else 16
    mv = int(a[a.index('--max-vectors') + 1]) if '--max-vectors' in a else 2 * 10 ** 7
    cc = int(a[a.index('--c') + 1]) if '--c' in a else 4
    logf = os.path.join(ROOT, f'logs/ab/paley_c{cc}.jsonl'); outdir = os.path.join(ROOT, f'logs/ab/paley_c{cc}')
    os.makedirs(outdir, exist_ok=True)
    done = set()
    if os.path.exists(logf):
        for l in open(logf):
            try: done.add(json.dumps(json.loads(l)['cfg'], sort_keys=True))
            except Exception: pass
    def log(rec):
        with open(logf, 'a') as fh: fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
    for v in vs:
        for cfg, sl in configs(v, mo, cc):
            key = json.dumps(cfg, sort_keys=True)
            if key in done: continue
            t0 = time.time(); rec = dict(cfg=cfg, orbits=int(sl.r), L1=len(sl.breps))
            c2 = SD.cfg_from_json(json.dumps(cfg))
            try:
                A = W.screen_A(sl, c2); rec['A'] = bool(A.get('A'))
                if not A.get('A'): log(dict(rec, secs=round(time.time() - t0, 1))); continue
                T1 = sl.T1; Gp = sl.Gp
                reps, _ = AGs.lag_reps(Gp, [Gp.act(u) for u in sl.Sp])
                w = sl.bsz.astype(np.int64)
                Emin = int(sum(w[j] * min(y * y for y in sl.letters[j]) for j in range(len(sl.breps))))
                Y, E = D.block_vectors(sl, int(T1[0]) - 3 * Emin, cap=mv, nQ=2, b=cfg['b'])
                if Y is None: log(dict(rec, stage1='skipped (block vectors > cap)', secs=round(time.time() - t0, 1))); continue
                full = Y[:, sl.borb]
                P = np.zeros((len(Y), len(reps)), np.int64)
                for c0 in range(0, len(Y), 1 << 16): P[c0:c0 + (1 << 16)] = Gp.paf(full[c0:c0 + (1 << 16)])[:, reps]
                keys, inv = np.unique(P, axis=0, return_inverse=True); inv = inv.ravel()
                groups = {}
                for i, g in enumerate(inv.tolist()): groups.setdefault(g, []).append(i)
                quads = D.ejoin4_sym(keys, T1[reps])
                sols = set()
                for qd in quads:
                    for a0 in groups[qd[0]]:
                        for a1 in groups[qd[1]]:
                            for a2 in groups[qd[2]]:
                                for a3 in groups[qd[3]]: sols.add(tuple(sorted((a0, a1, a2, a3))))
                ys = [full[list(s_)] for s_ in sorted(sols)]
                ys = [y for y in ys if (Gp.paf(y).sum(0) == T1).all()]
                rec.update(block_vectors=len(Y), stage1_solutions=len(ys))
                fitmap = {}
                if cc == 6:
                    import ab_borderfit as BF
                    qidx = Gp.co[:, 0] * 3 + Gp.co[:, 1]
                    keep = []
                    for i_, y in enumerate(ys):
                        z = np.zeros((4, 6), np.int64)
                        for k in range(4): np.add.at(z[k], qidx, y[k])
                        ft = BF.fits(z, 1388)
                        if ft: keep.append(y); fitmap[len(keep) - 1] = ft[0]
                    ys = keep; rec['bordered_stage1'] = len(ys)
                # stage 2
                greps, _ = AGs.lag_reps(sl.G, sl.perms)
                gw = np.random.default_rng(5).integers(1, 2 ** 62, size=len(greps), dtype=np.int64)
                ncore = 0; nskip = 0; first = None; firstfit = None
                for yi, y in enumerate(ys):
                    hits, m = SE.solve_y(sl, y, greps, gw)
                    if hits is None: nskip += 1; continue
                    for Xc in hits:
                        ncore += 1
                        if ncore <= 20: np.save(os.path.join(outdir, f"core_v{v}_q{cfg['q']}_s{cfg['s']}_{len(done)}_{ncore}.npy"), Xc.astype(np.int8))
                        if first is None: first = Xc; firstfit = fitmap.get(yi)
                rec.update(stage2_cores=ncore, stage2_skipped_fibre_cap=nskip)
                if first is not None and cc == 6:
                    import ab_assemble as AS, ab_borderfit as BF
                    T = BF.template(1388); perm, sg = firstfit
                    X6 = np.array([sg[i] * first[perm[i]] for i in range(4)], np.int64)
                    Z = AS.compressions(X6, tuple(cfg['mods']), 2)
                    Cf = AS.forced_corner(Z, (2, 3), T['layout'], T['K'], T['Kc'])
                    H = AS.assemble(X6, tuple(cfg['mods']), 2, T['layout'], T['K'], T['Kc'], Cf) if Cf is not None else None
                    rec['hadamard'] = bool(H is not None and AS.is_hadamard(H))
                    if rec['hadamard']:
                        tag = f"v{v}_q{cfg['q']}_s{cfg['s']}"
                        np.savez_compressed(os.path.join(outdir, f'H{4 * v}_{tag}.npz'), H=np.packbits(H > 0, axis=1),
                                            K=T['K'], Kc=T['Kc'], C=Cf, core=X6.astype(np.int8), mods=np.array(cfg['mods']))
                        rec['saved'] = f'logs/ab/paley_c6/H{4 * v}_{tag}.npz'
                elif first is not None:
                    import ab_border_c4 as B4
                    r = B4.border_c4(first, tuple(cfg['mods']), 2)
                    rec['hadamard'] = r is not None
                    if r is not None:
                        H, K, Kc, C = r
                        tag = f"v{v}_q{cfg['q']}_s{cfg['s']}"
                        np.savez_compressed(os.path.join(outdir, f'H{4 * v}_{tag}.npz'), H=np.packbits(H > 0, axis=1),
                                            K=K, Kc=Kc, C=C, core=first.astype(np.int8), mods=np.array(cfg['mods']))
                        rec['saved'] = f'logs/ab/paley_c4/H{4 * v}_{tag}.npz'
            except (MemoryError, RuntimeError) as e:
                rec['error'] = repr(e)[:200]
            rec['secs'] = round(time.time() - t0, 1)
            log(rec); done.add(key)


if __name__ == "__main__":
    main()
