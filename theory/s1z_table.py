r"""Table of the 112 distinct skipped b >= 3 configurations (117 log entries; theory/s1z_budget_batch.skipped_configs):
stage-1 z-screen (logs/s1z/zscreen_skipped.jsonl), symmetry order (theory/s1z_lib.log2_symmetry), stage-2 typical-y
budget (logs/s1z/budget_skipped.jsonl), exact compression screens of stage 2 (logs/s1z/agg/, logs/s1z/c719/), verdict.
Writes logs/s1z/table.json and prints a markdown table.
usage: python3 theory/s1z_table.py [--md FILE]"""
import sys, os, json, math, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import log2_symmetry
from s1z_budget_batch import skipped_configs, cfgkey

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUDGET_MAX = 47.0          # prompt: drop configurations whose stage 2 costs more than 2^47 per core
SYM_DEAD = -4.0            # calibration: every exactly empty rung was estimated at <= 2^-4 solutions per symmetry class


def jl(path):
    out = {}
    if os.path.exists(path):
        for l in open(path):
            try:
                r = json.loads(l); out[json.dumps(r['cfgkey'])] = r
            except Exception: pass
    return out


def agg_results():
    """per configuration tag: multisets, dead (with level), alive, errors (cheap levels) and the per-type full level"""
    res = {}
    for f in glob.glob(os.path.join(ROOT, 'logs/s1z/agg/agg_*.jsonl')):
        tag = os.path.basename(f)[4:-6]
        recs = {}
        nms = None
        for l in open(f):
            r = json.loads(l)
            if r.get('stage') == 'setup': nms = r['multisets']
            if r.get('stage') == 'ms': recs[r['i']] = r
        res[tag] = dict(multisets=nms, done=len(recs), dead=sum(1 for r in recs.values() if r.get('dead_at') is not None),
                        alive=sum(1 for r in recs.values() if r.get('dead_at') is None and 'error' not in r),
                        errors=sum(1 for r in recs.values() if 'error' in r))
    return res


def tag_of(cfg):
    v, c, q, s, H1, U = cfgkey(cfg)
    return f"{v}_c{c}_q{q}_s{s}_H{'x'.join(map(str, H1))}_U{json.dumps(U).replace(' ', '').replace(chr(34), '')}"


def main():
    a = sys.argv[1:]
    B = jl(os.path.join(ROOT, 'logs/s1z/budget_skipped.jsonl'))
    Zs = jl(os.path.join(ROOT, 'logs/s1z/zscreen_skipped.jsonl'))
    AG = agg_results()
    rows = []
    for cfg in skipped_configs():
        key = json.dumps(cfgkey(cfg)); b = B.get(key, {}); z = Zs.get(key, {})
        c2 = SD.cfg_from_json(json.dumps(cfg))
        sy = log2_symmetry(c2)
        r = dict(key=cfgkey(cfg), v=cfg['v'], c=cfg['c'], q=cfg['q'], s=cfg['s'], f=cfg['f'], H1=cfg['H1'], U=cfg['U'],
                 L1=cfg['L1'], orbits=cfg['orbits'], log2_sym=sy['log2_sym'])
        r['zscreen'] = {k: z.get(k) for k in ('multisets', 'ordered_zquads', 'log2_N0_total', 'log2_stage1_total', 'mode',
                                               'stage1_empty', 'exact_multisets', 'exact_zero', 'error', 'timeout') if k in z}
        S1 = z.get('log2_stage1_total')
        r['S1'] = S1; r['S1_sym'] = round(S1 - sy['log2_sym'], 1) if S1 is not None else None
        r['budget'] = {k: b.get(k) for k in ('H2', 'C2', 'stage2_log2E', 'D_def', 'C_S', 'h_D', 'log2_cost_per_core', 'energy_ok') if k in b}
        per_y = b.get('stage2_log2E')
        if S1 is not None and per_y is not None:
            r['log2_cores_per_sym'] = round(S1 + per_y - sy['log2_sym'] - math.log2(cfg['s']), 1)
        t = tag_of(cfg)
        if t in AG: r['exact_stage2'] = AG[t]
        if cfg['v'] == 719 and cfg['q'] == 17 and cfg['s'] == 2 and cfg['U'] == [['neg', [1]]]:
            r['exact_stage2'] = dict(multisets=39, done=39, dead=39, alive=0, errors=0, note='per-type screen, all 242,977 type quadruples')
        # verdict
        es = r.get('exact_stage2', {})
        cost = b.get('log2_cost_per_core')
        if b.get('energy_ok') is False or z.get('stage1_empty'):
            v_ = 'stage 1 empty [M]'
        elif es and es.get('multisets') and es.get('dead') == es.get('multisets'):
            v_ = 'no core [M] (exact stage-2 compression screens)'
        elif r['S1_sym'] is not None and r['S1_sym'] <= SYM_DEAD and S1 is not None:
            v_ = 'stage 1 < 2^-4 per symmetry class [E]'
        elif S1 is None and z.get('multisets'):
            v_ = 'stage 1 empty [E] (target outside the lattice of every multiset)'
        elif cost is not None and cost > BUDGET_MAX:
            v_ = f'stage 2 over budget [E] (2^{cost:.0f} per core)'
        else:
            v_ = 'open'
        r['verdict'] = v_
        rows.append(r)
    json.dump(rows, open(os.path.join(ROOT, 'logs/s1z/table.json'), 'w'), indent=1)
    lines = ['| v, c | q, s (f) | H\' | U | L1 | stage 1 [M/E] (log2 ordered; per sym. class) | stage 2 per y, D_def, cost/core [E] | exact stage-2 screens | verdict |',
             '|---|---|---|---|---|---|---|---|---|']
    def fmt(x, d=1):
        return '–' if x is None else (f'{x:.{d}f}' if isinstance(x, float) else str(x))
    for r in sorted(rows, key=lambda r: (r['v'], r['c'], r['q'], r['s'], r['L1'])):
        z = r['zscreen']
        if r['verdict'].startswith('stage 1 empty [M]'): s1 = 'empty [M]'
        elif r['S1'] is None: s1 = '−∞ [E]' if z.get('multisets') else '–'
        else: s1 = f"{r['S1']:.1f}; {r['S1_sym']:+.1f}"
        b = r['budget']
        s2 = f"{fmt(b.get('stage2_log2E'))}, {fmt(b.get('D_def'))}, 2^{fmt(b.get('log2_cost_per_core'))}" if b.get('D_def') is not None else '–'
        es = r.get('exact_stage2')
        ex = f"{es['dead']}/{es['multisets']} multisets dead" if es and es.get('multisets') else '–'
        U = json.dumps(r['U']).replace('"', '')
        lines.append(f"| {r['v']}, {r['c']} | {r['q']}, {r['s']} ({r['f']}) | {'×'.join(map(str, r['H1']))} | {U} | {r['L1']} | {s1} | {s2} | {ex} | {r['verdict']} |")
    md = '\n'.join(lines)
    if '--md' in a: open(a[a.index('--md') + 1], 'w').write(md + '\n')
    print(md)


if __name__ == '__main__':
    main()
