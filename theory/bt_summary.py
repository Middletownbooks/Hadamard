r"""Census tables for report/b3_push_answer.md from the logs (no computation; reproducible from logs/bt/).

Reads the census logs (logs/bt/census/k4_P4_w*.jsonl, logs/bt/census2/*.jsonl), the Theorem alpha screen
(logs/bt/census/thm_alpha.jsonl), the gate ladders (logs/bt/ladder_run*.jsonl), the level-A gates
(logs/bt/level1_*.jsonl), the exhaustive chains (logs/bt/decide_run1.log), the border census
(logs/bt/border_census.jsonl) and the next-gate costs (logs/bt/cost_next_gate.jsonl), and prints Markdown tables.
usage: python3 theory/bt_summary.py [--top 5]
"""
import sys, os, json, glob, argparse, collections

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def jl(path):
    out = []
    for f in sorted(glob.glob(os.path.join(ROOT, path))):
        for L in open(f):
            L = L.strip()
            if L.startswith('{'):
                try:
                    out.append(json.loads(L))
                except Exception:
                    pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--top', type=int, default=5)
    a = ap.parse_args()
    census = {}
    for d in jl('logs/bt/census/k4_P4_w*.jsonl') + jl('logs/bt/census2/*.jsonl'):
        d.setdefault('fam', 'k4')
        census[d['key']] = d
    alpha = {d['key']: d['dead'] for d in jl('logs/bt/census/thm_alpha.jsonl')}
    ladder = {}
    for d in jl('logs/bt/ladder_run*.jsonl'):
        ladder[d['key']] = d
    level1 = collections.defaultdict(list)
    for d in jl('logs/bt/level1_*.jsonl'):
        level1[d['key']].append(d)
    decide = {d['key']: d for d in jl('logs/bt/decide_run1.log') if 'key' in d}
    border = {d['key']: d for d in jl('logs/bt/border_census.jsonl')}
    cost = {d['key']: d for d in jl('logs/bt/cost_next_gate.jsonl')}

    def status(d):
        v = d.get('verdict', '')
        k = d['key']
        if v.startswith('no core [P] (coset'):
            return 'P:cp'
        if 'Theorem alpha' in v or 'order-3 element' in v or alpha.get(k):
            return 'P:alpha'
        if v.startswith('skipped'):
            return 'pre'
        if 'Q-level empty' in v:
            return 'M:Q'
        if not v.startswith('passes'):
            return 'other'
        if k in decide and decide[k]['verdict'].startswith('no core'):
            return 'M:chain'
        if k in ladder and ladder[k].get('verdict', '').startswith('no core'):
            return 'M:level'
        if any(x.get('levelA_solutions_weighted') == 0 for x in level1.get(k, [])):
            return 'M:level'
        return 'alive'

    groups = collections.defaultdict(list)
    for d in census.values():
        groups[(d['n'], d['fam'])].append(d)
    print('| n | family | structures | [P] coset-Parseval | [P] Theorem alpha | prefilter (4r - 6.5e < -45) | '
          '[M] Q-level empty | [M] deeper level / chain empty | alive | best alive: S, r, e, log2 E_Q |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    alive_all = []
    for (n, fam) in sorted(groups):
        c = collections.Counter(status(d) for d in groups[(n, fam)])
        alive = [d for d in groups[(n, fam)] if status(d) == 'alive']
        alive_all += alive
        best = max(alive, key=lambda d: (d.get('log2E_Q') if d.get('log2E_Q') is not None else -1e9), default=None)
        bs = '—' if best is None else \
            f"S{best['order']}, {best['r']}, {best['e']}, {best.get('log2E_Q') if best.get('log2E_Q') is not None else '−∞'}"
        print(f"| {n} | {fam} | {len(groups[(n, fam)])} | {c['P:cp']} | {c['P:alpha']} | {c['pre']} | {c['M:Q']} | "
              f"{c['M:level'] + c['M:chain']} | {len(alive)} | {bs} |")
    print()
    print(f'Top {a.top} alive structures per order and family (log2 E_Q = Q-level stratified count at the true target; '
          'next gate = cheapest exact level beyond Q, log2 of the 2+2 join pairs; border = Q-level classes with a border):')
    print()
    print('| n | family | G = P x O | S order | r | e | log2 E_Q | exact levels passed | next gate (log2 pairs) | border |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    by = collections.defaultdict(list)
    for d in alive_all:
        by[(d['n'], d['fam'])].append(d)
    for k in sorted(by):
        for d in sorted(by[k], key=lambda d: -(d.get('log2E_Q') if d.get('log2E_Q') is not None else -1e9))[:a.top]:
            lv = ladder.get(d['key'], {}).get('levels', [])
            passed = ', '.join(f"N{l['N']}: {l['solutions']}" for l in lv if l.get('solutions'))
            capped = [l['N'] for l in lv if l.get('capped')]
            if capped:
                passed += (' ' if passed else '') + f"(capped: N{', N'.join(map(str, capped))})"
            co = cost.get(d['key'], {}).get('log2_cost_next')
            b = border.get(d['key'])
            bs = '—' if b is None else (f"{b.get('with_border')}/{b.get('classes')}" if 'classes' in b else b['status'])
            print(f"| {d['n']} | {d['fam']} | {d['Pmods']} x {d['Omods']} | {d['order']} | {d['r']} | {d['e']} | "
                  f"{d.get('log2E_Q')} | {passed or '—'} | {co if co is not None else '—'} | {bs} |")


if __name__ == '__main__':
    main()
