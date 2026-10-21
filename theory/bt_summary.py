r"""Census tables for report/b3_push_answer.md from the logs (no computation; reproducible from logs/bt/).

Reads the census logs (logs/bt/census/k4_P4_w*.jsonl, split_k4.jsonl, logs/bt/census2/*.jsonl), the screens
(logs/bt/census/thm_alpha.jsonl, thm_rigid.jsonl), the gate ladders (logs/bt/ladder_run*.jsonl), the level-A gates
(logs/bt/level1_*.jsonl), the exhaustive chains (logs/bt/decide_run1.log), the border census
(logs/bt/border_census.jsonl), the next-gate costs (logs/bt/cost_next_gate.jsonl) and the cheapest-first gates
(logs/bt/gates_run*.jsonl).  `python3 theory/bt_summary.py` prints the tables; `--fill` rewrites the generated
blocks of report/b3_push_answer.md (between <!-- NAME --> and <!-- /NAME --> markers).
usage: python3 theory/bt_summary.py [--top 3] [--fill]
"""
import sys, os, json, glob, argparse, collections, math, subprocess

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


def rstr(r):
    return '/'.join(map(str, r)) if isinstance(r, list) else str(r)


def rmax(r):
    return max(r) if isinstance(r, list) else r


class Data:
    def __init__(self):
        self.census = {}
        for d in jl('logs/bt/census/k4_P4_w*.jsonl') + jl('logs/bt/census/split_k4.jsonl') + \
                [d for d in jl('logs/bt/census2/*.jsonl') if not d.get('verdict', '').startswith('skipped')]:
            d.setdefault('fam', 'k4')
            if 'key' in d:
                self.census[d['key']] = d
        for d in jl('logs/bt/census2/*.jsonl'):              # prefilter-skipped, unless re-run by bt_close.py
            if d.get('verdict', '').startswith('skipped') and d['key'] not in self.census:
                self.census[d['key']] = d
        self.alpha = {d['key']: d['dead'] for d in jl('logs/bt/census/thm_alpha.jsonl')}
        self.rigid = {d['key']: d['dead'] for d in jl('logs/bt/census/thm_rigid.jsonl')}
        self.ladder = {d['key']: d for d in jl('logs/bt/ladder_run*.jsonl')}
        self.level1 = collections.defaultdict(list)
        for d in jl('logs/bt/level1_*.jsonl'):
            self.level1[d['key']].append(d)
        self.decide = {d['key']: d for d in jl('logs/bt/decide_run1.log') if 'key' in d}
        self.border = {d['key']: d for d in jl('logs/bt/border_census.jsonl')}
        self.gates = collections.defaultdict(list)
        for d in jl('logs/bt/gates_run*.jsonl'):
            self.gates[d['key']].append(d)
        self.cost = {d['key']: d for d in jl('logs/bt/cost_next_gate.jsonl')}

    def status(self, d):
        v = d.get('verdict', '')
        k = d['key']
        if v.startswith('no core [P] (coset'):
            return 'P:cp'
        if 'Theorem alpha' in v or 'order-3 element' in v or self.alpha.get(k):
            return 'P:alpha'
        if 'Lemma R' in v or self.rigid.get(k):
            return 'P:rigid'
        if v.startswith('skipped'):
            return 'pre'
        if 'Q-level empty' in v:
            return 'M:Q'
        if not v.startswith('passes'):
            return 'other'
        if k in self.decide and self.decide[k]['verdict'].startswith('no core'):
            return 'M:chain'
        if k in self.ladder and self.ladder[k].get('verdict', '').startswith('no core'):
            return 'M:level'
        if any(x.get('levelA_solutions_weighted') == 0 for x in self.level1.get(k, [])):
            return 'M:level'
        if any(x.get('verdict', '').startswith('no core') for x in self.gates.get(k, [])):
            return 'M:level'
        return 'alive'

    def groups(self):
        g = collections.defaultdict(list)
        for d in self.census.values():
            g[(d['n'], d['fam'])].append(d)
        return g

    def levels_passed(self, d):
        lv = self.ladder.get(d['key'], {}).get('levels', [])
        sols = {l['N']: l['solutions'] for l in lv if l.get('solutions')}
        for x in self.gates.get(d['key'], []):
            if x.get('solutions_weighted'):
                sols[x['N']] = x['solutions_weighted']
        for x in self.level1.get(d['key'], []):
            if x.get('levelA_solutions_weighted'):
                sols[x['NA']] = x['levelA_solutions_weighted']
        capped = sorted({l['N'] for l in lv if l.get('capped')} - set(sols), reverse=True)
        return sols, capped

    def next_cost(self, d):
        """cheapest level (log2 pairs) among the levels not yet passed"""
        sols, _ = self.levels_passed(d)
        c = self.cost.get(d['key'])
        if not c:
            return None, None
        best = None
        for l in c.get('levels', []):
            if l['N'] in sols or l.get('log2_pairs') is None:
                continue
            if best is None or l['log2_pairs'] < best[0]:
                best = (l['log2_pairs'], l['N'])
        return best if best else (None, None)


def e_of(d):
    return d.get('log2E_Q') if d.get('log2E_Q') is not None else -1e9


def fmt_e(x):
    return '−∞' if x <= -1e8 else f'{x:.1f}'


def table_status(D):
    out = ['| n | family | structures | [P] coset-Parseval | [P] Theorem α | [P] Lemma R | prefilter only (4r − 6.5e < −45) | '
           '[M] Q-level empty | [M] deeper level / chain empty | alive | alive with log2 E_Q ≥ 0 | best alive: S, r, e, log2 E_Q |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for (n, fam), L in sorted(D.groups().items()):
        c = collections.Counter(D.status(d) for d in L)
        alive = [d for d in L if D.status(d) == 'alive']
        pos = [d for d in alive if e_of(d) >= 0]
        best = max(alive, key=e_of, default=None)
        bs = '—' if best is None else \
            f"S{best['order']}, {rstr(best['r'])}, {best['e']}, {best['log2E_Q'] if best.get('log2E_Q') is not None else '−∞'}"
        out.append(f"| {n} | {fam} | {len(L)} | {c['P:cp']} | {c['P:alpha']} | {c['P:rigid']} | {c['pre']} | {c['M:Q']} | "
                   f"{c['M:level'] + c['M:chain']} | {len(alive)} | {len(pos)} | {bs} |")
    return '\n'.join(out)


def table_top(D, top=3):
    out = ['| n | family | G = P × O | S order | r | e | log2 E_Q | exact levels passed (N: weighted solutions) | '
           'cheapest unrun level (log2 pairs) | border (classes with a border) |',
           '|---|---|---|---|---|---|---|---|---|---|']
    for (n, fam), L in sorted(D.groups().items()):
        alive = sorted([d for d in L if D.status(d) == 'alive'], key=lambda d: -e_of(d))[:top]
        for d in alive:
            sols, capped = D.levels_passed(d)
            passed = ', '.join(f"N{N}: {v}" for N, v in sorted(sols.items(), key=lambda t: -t[0])) or '—'
            nc, nN = D.next_cost(d)
            b = D.border.get(d['key'])
            bs = '—' if b is None else (f"{b.get('with_border')}/{b.get('classes')}" if 'classes' in b else 'not checked (k = 8)')
            out.append(f"| {n} | {fam} | {d['Pmods']} × {d['Omods']} | {d['order']} | {rstr(d['r'])} | {d['e']} | "
                       f"{d['log2E_Q'] if d.get('log2E_Q') is not None else '−∞'} | {passed} | "
                       f"{'—' if nc is None else f'N{nN}: {nc}'} | {bs} |")
    return '\n'.join(out)


def table_window(D):
    out = ['| n | family | best alive with r ≤ 60: log2 E_Q (S, r) | best alive with r ≤ 80 | best alive overall |',
           '|---|---|---|---|---|']
    for (n, fam), L in sorted(D.groups().items()):
        alive = [d for d in L if D.status(d) == 'alive' and d.get('log2E_Q') is not None]
        cells = []
        for lim in (60, 80, 10 ** 6):
            c = [d for d in alive if rmax(d['r']) <= lim]
            b = max(c, key=e_of, default=None)
            cells.append('—' if b is None else f"{b['log2E_Q']} (S{b['order']}, {rstr(b['r'])})")
        out.append(f"| {n} | {fam} | " + ' | '.join(cells) + ' |')
    return '\n'.join(out)


def table_rich(D, thr=40.0, per=3):
    out = ['| n | family | S, r, e | log2 E_Q | exact levels passed | capped levels | cheapest unrun level (log2 pairs) |',
           '|---|---|---|---|---|---|---|']
    for (n, fam), L in sorted(D.groups().items()):
        rich = sorted([d for d in L if D.status(d) == 'alive' and e_of(d) >= thr], key=lambda d: -e_of(d))[:per]
        for d in rich:
            sols, capped = D.levels_passed(d)
            passed = ', '.join(f"N{N}: {v}" for N, v in sorted(sols.items(), key=lambda t: -t[0])) or '—'
            nc, nN = D.next_cost(d)
            out.append(f"| {n} | {fam} | S{d['order']}, {rstr(d['r'])}, {d['e']} | {d['log2E_Q']} | {passed} | "
                       f"{', '.join('N%d' % x for x in capped) or '—'} | {'—' if nc is None else f'N{nN}: {nc}'} |")
    return '\n'.join(out)


def closed_families(D):
    lines = []
    for (n, fam), L in sorted(D.groups().items()):
        c = collections.Counter(D.status(d) for d in L)
        alive = [d for d in L if D.status(d) == 'alive']
        pre = c['pre']
        if not alive and not pre:
            lines.append(f"- **{n}, {fam}: closed [P, M].** All {len(L)} structures are dead: {c['P:cp']} coset-Parseval, "
                         f"{c['P:alpha']} Theorem α, {c['P:rigid']} Lemma R, {c['M:Q']} Q-level, "
                         f"{c['M:level'] + c['M:chain']} deeper level.")
        elif not alive:
            lines.append(f"- **{n}, {fam}:** no alive structure; {pre} structures are excluded only by the prefilter "
                         f"(not yet re-run exactly).")
        elif all(e_of(d) < -20 for d in alive):
            lines.append(f"- **{n}, {fam}: effectively closed.** {len(alive)} of {len(L)} structures pass the gates run, "
                         f"all with log2 E_Q < −20 (best {fmt_e(max(e_of(d) for d in alive))})"
                         + (f"; {pre} prefiltered" if pre else '') + '.')
    return '\n'.join(lines)


def screen_counts(D):
    tot = collections.Counter(D.status(d) for d in D.census.values())
    return (f"Over {len(D.census)} structures: coset-Parseval {tot['P:cp']}, Theorem α {tot['P:alpha']}, Lemma R "
            f"{tot['P:rigid']} (the screens [P]); Q-level empty {tot['M:Q']}, deeper level or chain empty "
            f"{tot['M:level'] + tot['M:chain']} [M]; prefilter only {tot['pre']}; alive {tot['alive']}.")


def border_counts(D):
    k = [b for b in D.border.values() if 'classes' in b]
    per = collections.Counter((b['n'], b['fam']) for b in k)
    return (f"{len(k)} structures ({', '.join(f'{n} {fam}: {c}' for (n, fam), c in sorted(per.items()))}), "
            f"{sum(b['classes'] for b in k)} Q-level classes, {sum(b['with_border'] for b in k)} with a border")


def fill(D, top):
    p = os.path.join(ROOT, 'report', 'b3_push_answer.md')
    s = open(p).read()
    calib = subprocess.run([sys.executable, os.path.join(ROOT, 'theory', 'bt_calib_table.py')], capture_output=True,
                           text=True).stdout.strip()
    alpha_n = sum(1 for v in D.alpha.values() if v)
    rigid_n = sum(1 for v in D.rigid.values() if v)
    blocks = {'CENSUS-TABLE': table_status(D), 'TOP-TABLE': table_top(D, top), 'WINDOW-TABLE': table_window(D),
              'RICH-TABLE': table_rich(D), 'CLOSED-FAMILIES': closed_families(D), 'SCREEN-COUNTS': screen_counts(D),
              'BORDER-COUNTS': border_counts(D), 'CALIB-TABLE': calib, 'ALPHA-COUNTS': str(alpha_n),
              'RIGID-COUNTS': f'{rigid_n} structures killed'}
    inline = {'ALPHA-COUNTS', 'RIGID-COUNTS', 'BORDER-COUNTS'}
    for name, text in blocks.items():
        a, b = f'<!-- {name} -->', f'<!-- /{name} -->'
        if a not in s:
            continue
        sep = '' if name in inline else '\n'
        i = s.index(a)
        if b in s:
            j = s.index(b) + len(b)
            s = s[:i] + a + sep + text + sep + b + s[j:]
        else:
            s = s[:i] + a + sep + text + sep + b + s[i + len(a):]
    open(p, 'w').write(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--top', type=int, default=3); ap.add_argument('--fill', action='store_true')
    a = ap.parse_args()
    D = Data()
    if a.fill:
        fill(D, a.top)
        return
    print(table_status(D)); print()
    print(table_top(D, a.top)); print()
    print(table_window(D)); print()
    print(table_rich(D)); print()
    print(closed_families(D)); print()
    print(screen_counts(D)); print(border_counts(D))


if __name__ == '__main__':
    main()
