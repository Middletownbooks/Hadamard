r"""Exact stage 2 over a stream of stage-1 solutions (JSON lines with key "y": 4 x L1 orbit values, as written by
theory/ab_s1_cpsat.py), for small b = 3 (c = 4) Paley-type slices [M].

Per stage-1 y (deduplicated up to block order and block sign, which leave stage 2 invariant): the exact one-shot
2+2 join of theory/ab_small_exact.solve_y (fibre cap 2^16 per block).  Every core is checked exactly
(sum_k PAF_G = T), saved, bordered by theory/ab_border_c4.border_c4 and the full H(4v) is checked (H H^T = 4v I).
Restart-safe: processed y's are logged to OUT (JSON lines); re-running skips them.
usage: python3 theory/ab_s2_stream.py CFG.json S1.jsonl OUT.jsonl [--cap 65536] [--watch]"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs, ab_s1sd as SD, ab_small_exact as SE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def canon(Yo):
    rows = []
    for r in Yo:
        r = list(r)
        rows.append(tuple(r) if r[0] > 0 else tuple(-x for x in r))
    return tuple(sorted(rows))


def main():
    a = sys.argv[1:]
    cfgf, s1f, outf = a[0], a[1], a[2]
    cap = int(a[a.index('--cap') + 1]) if '--cap' in a else 1 << 16
    cfg = SD.cfg_from_json(open(cfgf).read()); sl = W.build(cfg)
    reps, _ = AGs.lag_reps(sl.G, sl.perms)
    gw = np.random.default_rng(5).integers(1, 2 ** 62, size=len(reps), dtype=np.int64)
    tag = os.path.basename(cfgf).replace('cfg_', '').replace('.json', '')
    outdir = os.path.join(ROOT, 'logs/ab/paley6/cores'); os.makedirs(outdir, exist_ok=True)
    done = set()
    if os.path.exists(outf):
        for l in open(outf):
            try: done.add(tuple(tuple(r) for r in json.loads(l)['canon']))
            except Exception: pass
    ncore = 0; t0 = time.time(); pos = 0
    while True:
        new = 0
        if os.path.exists(s1f):
            with open(s1f) as fh:
                fh.seek(pos)
                for l in fh:
                    pos += len(l.encode())
                    try: z = json.loads(l)
                    except Exception: continue
                    if not z.get('verified', True): continue
                    c = canon(z['y'])
                    if c in done: continue
                    done.add(c); new += 1
                    Yo = np.array(c, np.int64)
                    Y = Yo[:, np.asarray(sl.borb)]                      # 4 x |G'|
                    hits, m = SE.solve_y(sl, Y, reps, gw, cap)
                    rec = dict(canon=[list(r) for r in c], fibres=m, cores=None if hits is None else len(hits),
                               secs=round(time.time() - t0, 1))
                    if hits:
                        for Xc in hits:
                            ok = bool((sl.G.paf(Xc).sum(0) == sl.T).all())
                            ncore += 1
                            np.save(os.path.join(outdir, f'core_{tag}_{ncore}.npy'), Xc.astype(np.int8))
                            rec.setdefault('core_checks', []).append(ok)
                            if ncore == 1:
                                import ab_border_c4 as B4
                                r = B4.border_c4(Xc, tuple(cfg['mods']), 2)
                                rec['hadamard'] = r is not None
                                if r is not None:
                                    H, K, Kc, C = r
                                    np.savez_compressed(os.path.join(outdir, f'H{4 * cfg["v"]}_{tag}.npz'),
                                                        H=np.packbits(H > 0, axis=1), K=K, Kc=Kc, C=C,
                                                        core=Xc.astype(np.int8), mods=np.array(cfg['mods']))
                    with open(outf, 'a') as f: f.write(json.dumps(rec) + "\n")
                    if hits: print(json.dumps(dict(stage='CORE', tag=tag, n=ncore, hadamard=rec.get('hadamard'))), flush=True)
        if '--watch' not in a: break
        if new == 0: time.sleep(30)
    print(json.dumps(dict(stage='end', tag=tag, ys=len(done), cores=ncore, secs=round(time.time() - t0, 1))), flush=True)


if __name__ == "__main__":
    main()
