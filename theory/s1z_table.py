r"""Table of the 112 distinct skipped b >= 3 configurations (117 log entries; theory/s1z_budget_batch.skipped_configs):
stage-1 z-screen (logs/s1z/zscreen_skipped.jsonl), symmetry order (theory/s1z_lib.log2_symmetry), stage-2 typical-y
budget (logs/s1z/budget_skipped.jsonl), exact compression screens of stage 2 (logs/s1z/agg/, logs/s1z/c719/), verdict.
Writes logs/s1z/table.json and prints a markdown table.
usage: python3 theory/s1z_table.py [--md FILE]"""
import sys, os, json, math, glob
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_s1sd as SD
from s1z_lib import log2_symmetry
from s1z_budget_batch import skipped_configs, cfgkey

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUDGET_MAX = 47.0          # prompt: drop configurations whose stage 2 costs more than 2^47 per core
SYM_DEAD = -4.0            # calibration: every exactly empty rung was estimated at <= 2^-4 solutions per symmetry class


def jl(path):
    """last record per configuration, a record without an error preferred"""
    out = {}
    if os.path.exists(path):
        for l in open(path):
            try:
                r = json.loads(l); k = json.dumps(r['cfgkey'])
                if k not in out or 'error' not in r or 'error' in out[k]: out[k] = r
            except Exception: pass
    return out


def agg_results():
    """per configuration tag: multisets, dead (with level), alive, errors (cheap levels) and the per-type full level"""
    res = {}
    for f in glob.glob(os.path.join(ROOT, 'logs/s1z/agg/agg_*.jsonl')):
        tag = os.path.basename(f)[4:-6]
        recs = {}
        nms = None; levels = None
        for l in open(f):
            r = json.loads(l)
            if r.get('stage') == 'setup': nms = r['multisets']; levels = r.get('levels')
            if r.get('stage') == 'ms': recs[r['i']] = r
        und = [r for r in recs.values() if r.get('dead_at') is None]
        res[tag] = dict(multisets=nms, done=len(recs), dead=sum(1 for r in recs.values() if r.get('dead_at') is not None),
                        alive=sum(1 for r in und if 'error' not in r), errors=sum(1 for r in und if 'error' in r),
                        dead_at=dict(Counter(json.dumps(r['dead_at']) for r in recs.values() if r.get('dead_at') is not None)),
                        passed=min((len(r.get('hits') or []) for r in und), default=None), levels=levels)
    return res


def char_results():
    out = {}
    p = os.path.join(ROOT, 'logs/s1z/agg/char_screens.jsonl')
    if os.path.exists(p):
        for l in open(p):
            r = json.loads(l)
            if r.get('multisets') is not None: out[r['tag']] = r
    return out


def u_text(cfg):
    """readable multiplier group: axes of G' = Q x H' (Q's axes first)"""
    nQ = len(cfg['Q']); mods = list(cfg['Q']) + list(cfg['H1'])
    def name(a):
        return (f"Q's Z{mods[a]}" if a < nQ else f"Z{mods[a]}") + ("" if mods.count(mods[a]) == 1 or a < nQ else f"(axis {a - nQ + 1})")
    if not cfg['U']: return '1'
    out = []
    for g in cfg['U']:
        if g[0] == 'neg': out.append('−1 on ' + ', '.join(name(a) for a in g[1] if mods[a] > 2))
        else:
            a, o, u = g
            out.append(('−1' if u == mods[a] - 1 else f'×{u}') + ' on ' + name(a))
    return '; '.join(out)


def level_name(cfg, keep):
    """the quotient Qbar x Z_q of a compression level (keep = Q's axes kept)"""
    return '×'.join([f"Z{cfg['Q'][a]}" for a in keep] + [f"Z{cfg['q']}"]) if len(keep) < len(cfg['Q']) else f"Q×Z{cfg['q']}"


def tag_of(cfg):
    v, c, q, s, H1, U = cfgkey(cfg)
    return f"{v}_c{c}_q{q}_s{s}_H{'x'.join(map(str, H1))}_U{json.dumps(U).replace(' ', '').replace(chr(34), '')}"


def main():
    a = sys.argv[1:]
    B = jl(os.path.join(ROOT, 'logs/s1z/budget_skipped.jsonl'))
    Zs = jl(os.path.join(ROOT, 'logs/s1z/zscreen_skipped.jsonl'))
    AG = agg_results(); CH = char_results()
    rows = []
    for cfg in skipped_configs():
        key = json.dumps(cfgkey(cfg)); b = B.get(key, {}); z = Zs.get(key, {})
        c2 = SD.cfg_from_json(json.dumps(cfg))
        sy = log2_symmetry(c2)
        r = dict(key=cfgkey(cfg), v=cfg['v'], c=cfg['c'], q=cfg['q'], s=cfg['s'], f=cfg['f'], Q=cfg['Q'], H1=cfg['H1'], U=cfg['U'],
                 L1=cfg['L1'], orbits=cfg['orbits'], log2_sym=sy['log2_sym'])
        r['zscreen'] = {k: z.get(k) for k in ('multisets', 'ordered_zquads', 'log2_N0_total', 'log2_stage1_total', 'mode',
                                               'evaluated', 'rel_se', 'stage1_empty', 'exact_multisets', 'exact_zero',
                                               'error', 'timeout') if k in z}
        S1 = z.get('log2_stage1_total')
        r['S1'] = S1; r['S1_sym'] = round(S1 - sy['log2_sym'], 1) if S1 is not None else None
        r['budget'] = {k: b.get(k) for k in ('H2', 'C2', 'stage2_log2E', 'D_def', 'C_S', 'h_D', 'log2_cost_per_core', 'energy_ok') if k in b}
        per_y = b.get('stage2_log2E')
        if S1 is not None and per_y is not None:
            r['log2_cores_per_sym'] = round(S1 + per_y - sy['log2_sym'] - math.log2(cfg['s']), 1)
        t = tag_of(cfg)
        if t in AG and AG[t].get('done'): r['exact_stage2'] = AG[t]
        if t in CH: r['char_screens'] = dict(multisets=CH[t]['multisets'], dead=CH[t]['dead'])
        if cfg['v'] == 719 and cfg['q'] == 17 and cfg['s'] == 2 and cfg['U'] == [['neg', [1]]]:
            r['exact_stage2'] = dict(multisets=39, done=39, dead=39, alive=0, errors=0,
                                     note='all 242,977 type quadruples; §5')
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
    lines = ['| v, c | q, s (f) | H\' | U | L1 | stage 1: log2 ordered; per sym. class | stage 2 per y, D_def, cost/core [E] | cores per class [E] | exact stage-2 screens | verdict |',
             '|---|---|---|---|---|---|---|---|---|---|']
    def fmt(x, d=1):
        return '–' if x is None else (f'{x:.{d}f}' if isinstance(x, float) else str(x))
    for r in sorted(rows, key=lambda r: (r['v'], r['c'], r['q'], r['s'], r['L1'])):
        z = r['zscreen']
        if r['verdict'].startswith('stage 1 empty [M]'): s1 = 'empty [M]'
        elif r['S1'] is None: s1 = '−∞ [E]' if z.get('multisets') else '–'
        else: s1 = f"{r['S1']:.1f}; {r['S1_sym']:+.1f}"
        b = r['budget']
        s2 = f"{fmt(b.get('stage2_log2E'))}, {fmt(b.get('D_def'))}, 2^{fmt(b.get('log2_cost_per_core'))}" if b.get('D_def') is not None else '–'
        es = r.get('exact_stage2'); ch = r.get('char_screens')
        if es and es.get('multisets'):
            ex = f"{es['dead']}/{es['multisets']} multisets dead"
            if es.get('note'): ex += f" ({es['note']})"
            elif es.get('levels') and es['dead'] < es['multisets'] and es.get('passed') is not None:
                lv = [level_name(r, k) for k in es['levels']]
                ex += f"; {es['multisets'] - es['dead']} pass {', '.join(lv[:es['passed']])}" + \
                      (f" ({lv[es['passed']]} too large)" if es['passed'] < len(lv) and es.get('errors') else '')
        elif ch and ch.get('multisets'):
            ex = f"per character: {ch['dead']}/{ch['multisets']} dead"
        else: ex = '–'
        U = u_text(r)
        cp = r.get('log2_cores_per_sym')
        cps = '–' if cp is None else f"{cp:+.0f}"
        lines.append(f"| {r['v']}, {r['c']} | {r['q']}, {r['s']} ({r['f']}) | {'×'.join(map(str, r['H1'])) or '1'} | {U} | {r['L1']} | {s1} | {s2} | {cps} | {ex} | {r['verdict']} |")
    md = '\n'.join(lines)
    if '--md' in a: open(a[a.index('--md') + 1], 'w').write(md + '\n')
    print(md)


if __name__ == '__main__':
    main()
