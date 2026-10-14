r"""Rank GP-type slices by their staged Gaussian figures [E]: stage 1 (exact lag-0/row count x lattice Gaussian) and
stage 2 for typical y (exact H2, log-det constraint), plus the plain whole-slice figure.
usage: python3 theory/ab_rank.py CFGFILE.jsonl OUT.jsonl   (one config per line, as in theory/ab_sweep.py)"""
import sys, os, json, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_sweep as W, ab_gauss as AGs
from ab_s1sd import cfg_from_json

def rank_one(cfg, N1=3000, N2=5000, ny=3):
    sl = W.build(cfg); t0 = time.time()
    rec = dict(v=cfg['v'], c=cfg['c'], Q=cfg['Q'], H1=cfg['H1'], q=cfg['q'], s=cfg['s'], f=cfg['f'],
               U=[list(u) for u in cfg['U']], orbits=int(sl.r), Gp_orbits=len(sl.breps))
    wh = AGs.whole(sl, N=4000); rec['whole_log2E'] = wh.get('log2E')
    s1 = AGs.stage1(sl, N=N1)
    rec.update(stage1_log2Ny=s1.get('log2Ny'), stage1_bits=float(s1['bits']) if s1.get('bits') is not None else None,
               stage1_mahal=s1.get('mahal'), stage1_log2E=float(s1['stage1_log2E']) if s1.get('stage1_log2E') is not None else None)
    draw, info = AGs.stage1_sampler(sl, seed=3)
    if draw is not None:
        res = [AGs.stage2_typical(sl, draw(), N=N2, seed=i) for i in range(ny)]
        rec['stage2_H2'] = [r['H2'] for r in res]; rec['stage2_C2'] = [r['C2'] for r in res]
        rec['stage2_log2E'] = [r['stage2_log2E'] for r in res]
    rec['secs'] = round(time.time() - t0, 1)
    return rec

if __name__ == "__main__":
    import resource
    memgb = float(os.environ.get('AB_MEMGB', '3'))
    resource.setrlimit(resource.RLIMIT_AS, (int(memgb * 2 ** 30), int(memgb * 2 ** 30)))
    done = set()
    if os.path.exists(sys.argv[2]):
        for l in open(sys.argv[2]):
            try:
                d = json.loads(l)
                done.add(json.dumps([d.get('v'), d.get('c'), d.get('Q'), d.get('H1'), d.get('q'), d.get('s'), d.get('U')]))
            except Exception: pass
    out = open(sys.argv[2], 'a')
    for line in open(sys.argv[1]):
        line = line.strip()
        if not line: continue
        cfg = cfg_from_json(line)
        key = json.dumps([cfg['v'], cfg['c'], list(cfg['Q']), list(cfg['H1']), cfg['q'], cfg['s'], [list(u) for u in cfg['U']]])
        if key in done: continue
        try: rec = rank_one(cfg)
        except (Exception, MemoryError) as e: rec = dict(v=cfg['v'], c=cfg['c'], Q=list(cfg['Q']), H1=list(cfg['H1']), q=cfg['q'], s=cfg['s'], U=[list(u) for u in cfg['U']], error=repr(e)[:200])
        s = json.dumps(rec); print(s, flush=True); out.write(s + "\n"); out.flush()
