"""Build the part-2 table of report/h4892_answer.md section 4 from the logs in logs/h4892/part2/ (prints markdown).
Per key: the tool's A and B results (A,B queue logs NAME.jsonl, B-only reruns B_NAME.jsonl; a missing final record
with rc 124 in done.jsonl = 3600 s timeout), and the independent second method (any_NAME.out)."""
import json, os, glob
D = 'logs/h4892/part2'
keys = [l.rstrip('\n').split('\t') for l in open(f'{D}/keys.txt')]
done = {}
for l in open(f'{D}/done.jsonl'):
    r = json.loads(l); done[r['name']] = r
def tool(path):
    if not os.path.exists(path): return None
    rec = None
    for l in open(path):
        r = json.loads(l)
        if r.get('phase') in ('decision', 'final'): rec = r
    return rec
def fmt_m(r, m):
    if r is None or m not in r: return None
    x = r[m]
    if 'solutions' in x:
        pl = [p.get(m) for p in r.get('planted') or []] if r.get('phase') == 'final' else None
        ps = '' if pl is None else f", planted {sum(1 for p in pl if p)}/{len(pl)}"
        return f"{x['solutions']} ({x.get('blocks')} blocks, {x.get('secs')} s{ps})"
    if 'undecided' in x:
        u = x['undecided']
        return None if 'not run' in u else 'undecided: ' + u.split('(')[0]
    return None
def anyres(name):
    p = f'{D}/any_{name}.out'
    if not os.path.exists(p): return None
    L = [json.loads(l) for l in open(p) if l.startswith('{')]
    if not L: return None
    r = L[0]; pl = L[1]['planted'] if len(L) > 1 else []
    if 'solutions' in r:
        ok = sum(1 for x in pl if x.get('recovered')); return f"{r['solutions']} ({r['blocks']} blocks, planted {ok}/{len(pl)})"
    return 'capped (' + r.get('undecided', '?').split('(')[0] + ')'
print('| 4v | key | A | B | independent (h4892_stage1_any) |')
print('|---|---|---|---|---|')
for name, key in keys:
    r1 = tool(f'{D}/{name}.jsonl'); r2 = tool(f'{D}/B_{name}.jsonl') or tool(f'{D}/{name}_toolB.jsonl')
    A = fmt_m(r1, 'A') or ('timeout 3600 s' if done.get(name, {}).get('rc') == 124 else '–')
    B = fmt_m(r1, 'B') or fmt_m(r2, 'B') or ('timeout 3600 s' if done.get('B_' + name, {}).get('rc') == 124 else
                                              ('not reached' if 'B_' + name not in done and r2 is None else '–'))
    print(f"| {name.split('_')[0]} | `{key}` | {A} | {B} | {anyres(name) or '–'} |")
