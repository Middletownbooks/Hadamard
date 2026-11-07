"""Step 1: coverage map of Kumar's skew work for every open census value v (prompt §1, §4 step 1).

Sources in K = vatskishan03/hadamard (read-only clone, commit recorded):
  search/NOTES.md (every paragraph mentioning v as a whole number, in a section on skew matrices or the census),
  skew/NEGACYCLIC.md, skew/PALEY_SKEW.md, search/skewplug/**/README.md, search/negaskew/README.md,
  search/skewplug/scan38/data/shapes38.json (554 Paley-product shapes with his decided/undecided status),
  search/skewplug/scan38/data/census_shapes15.json.
Output: logs/skew/coverage.jsonl, one record per open v:
  {v, prime, notes:[{sec, title, kind, snippet}], plug_shapes:{n, decided, undecided}, lead, methods}
kind is a keyword class of the paragraph (exhaustive / estimate / heuristic-run / solved-other / mention).
Lead = values his notes name as running, queued or a lead (spin primes of §14, §115 budgets, 'queued', 'pending').
This is a text map: a paragraph match is evidence of contact, not a decision; 'kind' is a heuristic label.
  python3 theory/skew_coverage.py K
"""
import json, re, sys, os, glob, subprocess
K = sys.argv[1] if len(sys.argv) > 1 else '/home/user/vatskishan03/hadamard'
here = os.path.dirname(os.path.abspath(__file__))
openv = json.load(open(f'{here}/../data/skew_open_census.json'))['open_v']
isprime = lambda n: n > 1 and all(n % p for p in range(2, int(n ** .5) + 1))

SKEWSECS = {10, 12, 14, 15, 16, 20, 21, 22, 23, 26, 34, 35, 36, 39, 48, 64, 66, 73, 87, 88, 89, 90, 91, 93, 94,
            107, 110, 111, 112, 115, 92}
def paragraphs(path, secfilter=None):
    sec, title, buf = None, '', []
    out = []
    for line in open(path, encoding='utf8'):
        m = re.match(r'^## (\d+)\. (.*)', line)
        if m:
            if buf: out.append((sec, title, ' '.join(buf)))
            sec, title, buf = int(m.group(1)), m.group(2).strip(), []
            continue
        if line.strip() == '' or re.match(r'^\s*[-*] ', line):
            if buf: out.append((sec, title, ' '.join(buf)))
            buf = [line.strip()] if line.strip() else []
        else:
            buf.append(line.strip())
    if buf: out.append((sec, title, ' '.join(buf)))
    return out

def kind(t):
    tl = t.lower()
    if re.search(r'\bfound\b|skew order \d+|settled|is solved|\bhit\b', tl) and not re.search(r'no hit|without a hit', tl): k = 'solved-or-hit-context'
    else: k = None
    if re.search(r'exhaustiv|ruled out|exactly empty|exactly dead|infeasible|no (admissible|feasible)|stage 1 is empty|dies|die by|empty', tl): return 'exact-negative' if k is None else k + '+exact-negative'
    if re.search(r'no solution|nothing|budget ran out|sampled|unknown|undecided|timeout|not finished|stopped', tl): return 'heuristic-negative' if k is None else k + '+heuristic-negative'
    if re.search(r'estimate|10\^|predicted|log', tl): return 'estimate' if k is None else k
    return k or 'mention'

srcs = [(f'{K}/search/NOTES.md', True)] + [(p, False) for p in
        [f'{K}/skew/NEGACYCLIC.md', f'{K}/skew/PALEY_SKEW.md'] + glob.glob(f'{K}/search/skewplug/**/README.md', recursive=True)
        + glob.glob(f'{K}/search/negaskew/**/README.md', recursive=True) + glob.glob(f'{K}/search/listgen/README.md')]
paras = []
for p, isnotes in srcs:
    for sec, title, t in paragraphs(p):
        if isnotes and sec not in SKEWSECS and 'skew' not in t.lower(): continue
        paras.append((os.path.relpath(p, K) + (f' §{sec}' if isnotes else ''), title, t))

shapes = json.load(open(f'{K}/search/skewplug/scan38/data/shapes38.json'))
LEADS = {1951: '§14 spin prime (scheduler queue)', 2269: '§14 spin prime', 2437: '§14 spin prime',
         661: '§115 dd_ktree 30 h budget, no solution (being run)', 883: '§115 dd_ktree 44 h budget, no solution (being run)'}
commit = subprocess.run(['git', '-C', K, 'log', '-1', '--format=%h'], capture_output=True, text=True).stdout.strip()
out = open(f'{here}/../logs/skew/coverage.jsonl', 'w')
summ = {'any': 0, 'exact-negative': 0, 'none': 0}
for v in openv:
    pat = re.compile(r'(?<![\d.,·⁰-⁹])' + str(v) + r'(?![\d⁰-⁹]|\s?·\s?10|,\d)')
    hits = []
    for src, title, t in paras:
        for m in pat.finditer(t):
            s = t[max(0, m.start() - 160): m.end() + 200]
            hits.append(dict(src=src, title=title[:60], kind=kind(s), snippet=s))
            break
    sh = [s for s in shapes if s['v'] == v]
    lead = LEADS.get(v)
    for h in hits:
        if re.search(r'queued|pending|scheduled', h['snippet'].lower()): lead = (lead or '') + ' queued/pending: ' + h['src']
    rec = dict(v=v, prime=isprime(v), n_paragraphs=len(hits), kinds=sorted({h['kind'] for h in hits}),
               plug_shapes=dict(n=len(sh), decided=sum(1 for s in sh if s.get('decided')),
                                undecided=[(s['m'], s['q'], s['order'], s['nu']) for s in sh if not s.get('decided')][:8]),
               lead=lead, kumar_commit=commit, notes=hits)
    out.write(json.dumps(rec, ensure_ascii=False) + '\n')
    summ['any'] += bool(hits or sh); summ['none'] += not (hits or sh)
    summ['exact-negative'] += any('exact-negative' in h['kind'] for h in hits)
print(summ, 'open', len(openv))
