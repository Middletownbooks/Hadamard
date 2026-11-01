r"""Tables of the Z_q census from the logs (no computation) [M gates, E counts].
Reads logs/zq/census/{fam}_{n}.jsonl (theory/zq_census.py) and logs/zq/gates/{fam}_{n}.jsonl
(theory/bt_gates_levels.py).  A structure is dead if a screen [P] or the Q-level gate [M] killed it, or if any exact
level G/N it was run at has no solution [M]; otherwise alive.
usage: python3 theory/zq_summary.py [--keys fam_n OUT]   (alive keys of one order, for bt_gates_levels --keys @OUT)
"""
import sys, os, json, glob, math
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def jl(p):
    out = []
    if os.path.exists(p):
        for L in open(p):
            L = L.strip()
            if L.startswith('{'):
                try: out.append(json.loads(L))
                except Exception: pass
    return out


def load(tag):
    C = {d['key']: d for d in jl(f'{ROOT}/logs/zq/census/{tag}.jsonl')}
    G = {}
    for d in jl(f'{ROOT}/logs/zq/gates/{tag}.jsonl'):
        G.setdefault(d['key'], []).append(d)
    return C, G


def status(d, g):
    v = d.get('verdict', '')
    if 'coset-Parseval' in v: return 'cp'
    if 'Lemma R' in v: return 'rigid'
    if 'Q-level empty' in v: return 'qempty'
    if v.startswith('error') or 'cap' in v: return 'undecided'
    if any(x.get('verdict', '').startswith('no core') for x in g): return 'deeper'
    return 'alive'


def table(tags):
    rows = []
    for tag in tags:
        C, G = load(tag)
        if not C:
            continue
        fam, n = tag.rsplit('_', 1)
        cnt = {k: 0 for k in ('cp', 'rigid', 'qempty', 'deeper', 'undecided', 'alive')}
        alive = []
        for k, d in C.items():
            s = status(d, G.get(k, [])); cnt[s] += 1
            if s == 'alive':
                alive.append(d)
        pos = [d for d in alive if d.get('log2E_Q') is not None and d['log2E_Q'] >= 0]
        best = max(alive, key=lambda d: d['log2E_Q'] if d.get('log2E_Q') is not None else -1e9) if alive else None
        # next gate of the best alive: the levels run on it
        nxt = None
        if best is not None:
            lv = G.get(best['key'], [])
            pas = [f"N{x['N']}:{x.get('solutions_weighted')}" for x in lv if x.get('verdict') == 'passes']
            cap = [f"N{x['N']}:{x['verdict']}" for x in lv if not x.get('verdict', '').startswith(('passes', 'no core'))]
            nxt = dict(passed=pas, open=cap)
        rows.append(dict(n=int(n), fam=fam, structures=len(C), **cnt, alive_pos=len(pos),
                         best=None if best is None else
                         f"S{best['order']}, r{best['r']}, e{best['e']}, {best.get('log2E_Q')}", next=nxt))
    return rows


def md(rows):
    out = ['| 4n | n | family | structures | [P] coset-Parseval | [P] Lemma R | [M] Q-level empty | '
           '[M] finer level empty | undecided | alive | alive, log2 E_Q ≥ 0 | best alive: S, r, e, log2 E_Q | '
           'best alive: levels passed / open |',
           '|' + '---|' * 13]
    for r in sorted(rows, key=lambda r: (r['n'], r['fam'])):
        nx = '—' if not r['next'] else (', '.join(r['next']['passed']) or '—') + ' / ' + \
            (', '.join(r['next']['open']) or '—')
        out.append(f"| {4 * r['n']} | {r['n']} | {r['fam']} | {r['structures']} | {r['cp']} | {r['rigid']} | "
                   f"{r['qempty']} | {r['deeper']} | {r['undecided']} | {r['alive']} | {r['alive_pos']} | "
                   f"{r['best'] or '—'} | {nx} |")
    return '\n'.join(out)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--keys':
        C, G = load(sys.argv[2])
        ks = [k for k, d in C.items() if status(d, G.get(k, [])) == 'alive']
        ks.sort(key=lambda k: -(C[k].get('log2E_Q') or -1e9))
        open(sys.argv[3], 'w').write(','.join(ks))
        print(len(ks))
    else:
        tags = sorted(os.path.basename(p)[:-6] for p in glob.glob(f'{ROOT}/logs/zq/census/*.jsonl')
                      if 'regress' not in p)
        print(md(table(tags)))
